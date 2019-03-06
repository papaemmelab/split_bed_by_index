let doc = """
Split bed intervals by estimated read counts from bam index

Usage:
  split_bed_by_index [options] <BED> <BAM> <OUT>

Arguments:
  <BED>     Input bed file
  <BAM>     Input bam file
  <OUT>     Output bed file

Options:
  -h --help                 Show this screen.
  -v --version              Show version.
  -c --count <count>        Number of reads to split  [default: 1000000].
  -m --merge                Merge overlapping regions before split.
"""

import algorithm
import docopt
import hts
import math
import os
import sequtils
import streams
import strutils
import tables
import times

const chunk_size = 16_384

proc `/`(s1: openArray[int|float], s2: openArray[int|float]): seq[float] =
    let len = s1.len
    assert len == s2.len
    newSeq(result, len)
    for i in 0 ..< len:
        if s2[i] != 0:
            result[i] = s1[i] / s2[i]
        else:
            result[i] = 0

proc `/`(s: openArray[int|float], v: int|float): seq[float] =
    assert v != 0
    let d = float(v)
    newSeq(result, s.len)
    for i in 0 ..< s.len:
        result[i] = float(s[i]) / d

type
    region_t = ref object
        chrom: string
        start: int
        stop: int
        name: string
        count: int

# from https://github.com/brentp/hts-nim-tools
proc bed_line_to_region(line: string): region_t =
    var
        cse = line.strip().split('\t', 5)

    if len(cse) < 3:
        stderr.write_line("skipping bad bed line:", line.strip())
        return nil
    var
        s = parse_int(cse[1])
        e = parse_int(cse[2])
        reg = region_t(chrom: cse[0], start: s, stop: e, count:0)
    if len(cse) > 3:
        reg.name = cse[3]
    return reg

# from https://github.com/brentp/hts-nim-tools
proc bed_to_table(bed: string): TableRef[string, seq[region_t]] =
    var bed_regions = newTable[string, seq[region_t]]()
    var hf = hts.hts_open(cstring(bed), "r")
    var kstr: hts.kstring_t
    kstr.l = 0
    kstr.m = 0
    kstr.s = nil
    while hts_getline(hf, cint(10), addr kstr) > 0:
        if ($kstr.s).startswith("track "):
            continue
        if $kstr.s[0] == "#":
            continue
        var v = bed_line_to_region($kstr.s)
        if v == nil: continue
        discard bed_regions.hasKeyOrPut(v.chrom, new_seq[region_t]())
        bed_regions[v.chrom].add(v)

    # since it is read into mem, can also well sort.
    for chrom, ivs in bed_regions.mpairs:
        sort(ivs, proc (a, b: region_t): int = a.start - b.start)

    hts.free(kstr.s)
    return bed_regions

proc readIndex(fs: FileStream, sizes: var seq[seq[int]], total_bytes: var seq[int], total_reads: var seq[int]) =
    # read index and set sizes, total_bytes, and total_reads
    discard fs.readStr(4)
    var n_ref = fs.readInt32()
    for chr in 0..<n_ref:
        sizes.add(newSeq[int]())
        var n_bin = fs.readInt32()
        for bins in 0..<n_bin:
            var bin = fs.readUint32()
            var n_chunk = fs.readInt32()
            for chunk in 0..<n_chunk:
                var chunk_beg = fs.readUint64()
                var chunk_end = fs.readUint64()
                if bin == 37450 and chunk == 1:
                    total_reads.add(int(chunk_beg) + int(chunk_end))
        if total_reads.len < chr+1:
            total_reads.add(0)
        var n_intv = fs.readInt32()
        var offset = 0
        for intv in 0..<n_intv:
            var ioffset = fs.readUint64()
            if offset != 0:
                sizes[chr].add(int(ioffset) - int(offset))
            offset = int(ioffset)
        total_bytes.add(sum(sizes[chr]))
    var n_no_coor = fs.readUint64()
    fs.close()

proc splitRegion(chrom: string, contig: int, start_pos: int, end_pos: int, sizes: seq[seq[int]], count: int, bytes_per_read: float, out_bed: File) =
    # split a region and write to file
    var pos = start_pos
    var counter = 0.0
    let start_chunk = int(start_pos / chunk_size)
    for i,chunk in sizes[contig][int(start_pos / chunk_size)..int(end_pos / chunk_size)]:
        let chunk_ind = int(start_pos / chunk_size) + i
        let est_reads = float(chunk)/bytes_per_read
        counter += est_reads
        if counter > float(count):
            echo chrom & ":" & $pos & "-" & $((chunk_ind+1)*chunk_size - 1) & "\t" & $counter
            out_bed.writeLine(chrom & "\t" & $pos & "\t" & $((chunk_ind+1)*chunk_size - 1) & "\t" & $counter)
            pos = (chunk_ind+1)*chunk_size
            counter = 0.0
    echo chrom & ":" & $pos & "-" & $end_pos & "\t" & $counter
    out_bed.writeLine(chrom & "\t" & $pos & "\t" & $end_pos & "\t" & $counter)

proc mergeRegions(regions: TableRef[string, seq[region_t]]): TableRef[string, seq[region_t]] =
    # adapted from https://www.geeksforgeeks.org/merging-intervals/
    # we already know the structure is sorted
    result = newTable[string, seq[region_t]]()
    for chrom, intervals in regions.pairs():
        # put first interval in stack
        result[chrom] = @[intervals[0]]
        if intervals.len > 1:
            # start from next interval and merge if needed
            for interval in intervals[1..intervals.len-1]:
                # get last element
                let top = result[chrom][result[chrom].len-1]
                # if not overlapping add to stack
                if top.stop < interval.start:
                    result[chrom].add(interval)
                # else update stop if current greater than top
                elif top.stop < interval.stop:
                    top.stop = interval.stop
                    discard result[chrom].pop()
                    result[chrom].add(top)


proc main() =
    let args = docopt(doc, version = "0.1.0")
    echo "BED: ", $args["<BED>"]
    echo "BAM: ", $args["<BAM>"]
    echo "OUT: ", $args["<OUT>"]
    echo "Count: ", $args["--count"], " reads"
    echo ""

    var
        bam:Bam
        sizes = newSeq[seq[int]]() # matrix of bytes per chunk per contig
        total_bytes = newSeq[int]() # bytes per contig
        total_reads = newSeq[int]() # reads per contig
        regions = bed_to_table($args["<BED>"])
        idx_path: string

    let
        time = now()
        split = parseInt($args["--count"])
        out_bed = open($args["<OUT>"], fmWrite)
        bam_path = $args["<BAM>"]

    open(bam, bampath, index=true)

    # find index
    if os.fileExists(bam_path & ".bai"):
        idx_path = bam_path & ".bai"
    elif os.fileExists(bam_path.replace(".bam",".bai")):
        idx_path = bam_path.replace(".bam",".bai")
    else:
        raise newException(ValueError, "Bam index could not be found")

    # resolve overlaps
    if args["--merge"]:
        regions = mergeRegions(regions)

    # read index
    let idx = newFileStream(idx_path, fmRead)
    readIndex(idx, sizes, total_bytes, total_reads)

    # calculate average bytes per read
    let bytes_per_read = sum(total_bytes/total_reads)/float((total_bytes/total_reads).len)

    # loop through regions and split
    for target in targets(bam.hdr):
        var chrom = target.name
        if not regions.contains(chrom) or regions[chrom].len == 0:
            continue
        for region in regions[chrom]:
            splitRegion(chrom, target.tid, region.start, region.stop, sizes, split, bytes_per_read, out_bed)
    echo ""
    echo "Elapsed time: " & $(now() - time)
    close(bam)
    close(out_bed)

when isMainModule:
  main()
