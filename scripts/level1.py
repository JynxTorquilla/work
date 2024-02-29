#!/usr/bin/env python

import sys
import time
import os
import tempfile

import pysam
from numpy import *


__version__ = "2015.04.06"


def find_assembly(handle):
    assembly = None
    species = None
    for row in handle.header['SQ']:
        if 'AS' in row:
            if assembly:
                if assembly != row['AS']:
                    raise ValueError("Different genome assembly found (found %s, expected %s)" % (row['AS'], assembly))
            else:
                assembly = row['AS']
        if 'SP' in row:
            if species:
                if species != row['SP']:
                    raise ValueError("Different species found (found %s, expected %s)" % (row['SP'], species))
            else:
                species = row['SP']
    return assembly

def find_condition(filename, handle):
    # Use the information in the SAM/BAM header if available, otherwise use
    # the file name directly
    for row in handle.header.get('RG', []):
        library = row.get('LB')
        sample = row.get('SM')
        condition = "%s.%s" % (library, sample)
        break
    else:
        basename = os.path.basename(filename)
        root, extension = os.path.splitext(basename)
        if extension.lower() in ('.bam', '.sam'):
            condition = root
        else:
            condition = basename
    return condition

def open_input_file(input_filename):
    if input_filename.endswith(".bam"):
        mode = 'rb'
    else:
        mode = 'r'
    handle = pysam.Samfile(input_filename, mode)
    return handle

def open_output_file(output_filename):
    if output_filename.endswith(".gz"):
        import gzip
        output = gzip.open(output_filename, 'w')
    elif output_filename.endswith(".bz2"):
        import bz2
        output = bz2.BZ2File(output_filename, 'w')
    else:
        output = open(output_filename, 'w')
    return output

def write_meta_data(output, input_filenames, conditions, assembly, chromosomes,
                    threshold, flag_select, flag_skip,
                    quality, identity, process_fingerprint):
    # General meta-data
    output.write("##FileFormat = OSCtable\n")
    output.write("##ProtocolREF = OP-CAGE-Tag-Clustering-v1.1\n")
    output.write("##ParameterValue[level1_script_version] = %s\n" % __version__)
    if assembly:
        output.write("##ParameterValue[genome_assembly] = %s\n" % assembly)
    output.write("##ParameterValue[threshold] = %f\n" % threshold)
    output.write("##ParameterValue[flag_select] = 0x%x\n" % flag_select)
    output.write("##ParameterValue[flag_skip] = 0x%x\n" % flag_skip)
    output.write("##ParameterValue[quality] = %f\n" % quality)
    output.write("##ParameterValue[identity] = %f\n" % (identity*100))
    if process_fingerprint:
        output.write("##ParameterValue[fingerprint] = True\n")
    output.write("##Date = %s\n" % time.strftime("%Y-%m-%d"))
    output.write("##ContactName = LSA support staff\n")
    output.write("##ContactEmail = lsa-support@gsc.riken.jp\n")
    for input_filename in input_filenames:
        input_filename = os.path.basename(input_filename)
        output.write("##InputFile = %s\n" % input_filename)
    output.write("##ChromosomeNameOrder = %s\n" % ",".join(chromosomes))
    output.write("##ColumnVariable[id] = identifier of the level-1 promoter\n")
    output.write("##ColumnVariable[chrom] = name of the chromosome, identical to what is written in genome assembly\n")
    output.write('##ColumnVariable[start] = start genomic coordinate of the promoter, starting from 0\n')
    output.write("##ColumnVariable[end] = end genomic coordinate of the promoter, always equal to start+1\n")
    output.write("##ColumnVariable[strand] = genomic strand on which the promoter is located\n")
    for condition in conditions:
        output.write("##ColumnVariable[raw.%s] = rescued tag count of sample %s\n" % (condition, condition))
    for condition in conditions:
        output.write("##ColumnVariable[norm.%s] = normalized expression for sample %s\n" % (condition, condition))
    output.write("id\tchrom\tstart.0base\tend\tstrand")
    for condition in conditions:
        output.write("\traw.%s" % condition)
    for condition in conditions:
        output.write("\tnorm.%s" % condition)
    output.write("\n")

def merge(output, tempfiles, chromosomes, name, totals):
    for c, chromosome in enumerate(chromosomes):
        for strand in '+-':
            handle = tempfiles[c][strand]
            handle.flush()
            handle.seek(0)
            for line in handle:
                words = line.split()
                position = int(words[0])
                values = map(float, words[1:])
                promoter = "L1_%s_%s_%s_%d" % (name, chromosome, strand, position)
                output.write("%s\t%s\t%d\t%d\t%s" % (promoter, chromosome, position, position + 1, strand))
                for value in values:
                    output.write("\t%f" % value)
                for value, total in zip(values, totals):
                    tpm = 1.e6 * value / total
                    output.write("\t%f" % tpm)
                output.write("\n")
            handle.close()

def write_tempfile_line(handle, position, values):
    handle.write("%d" % position)
    for value in values:
        handle.write("\t%f" % value)
    handle.write("\n")

def find_tag_count(read):
    count = 1.0
    weight = 1.0
    for key, value in read.tags:
        if key == 'XC':
            count = value
        elif key == 'XW':
            weight = value
    count *= weight
    return count

def get_next_mapped_read(handle, flag_select=0, flag_skip=4, quality=0,
                         identity=0):
    for read in handle:
        if read.flag & flag_select != flag_select: # skip if not selected
            continue
        if read.flag & flag_skip: # skip unwanted tags (default: unmapped).
            continue
        if read.mapq < quality:
            # skip tags with a mapping quality below the cutoff
            continue
        if identity > 0:
            length = sum([c for i, c in read.cigar if "MIDNSHP"[i] in "MDN"])
            for key, value in read.tags:
                if key == 'NM':
                    errors = value
                    break
            else:
                raise Exception(
                    "Filtering on percentage identity but the edit distance (NM tag) was not found in the input file")
            if length - errors < length * identity:
                continue
        return read
    return None

def write_temporary_files(handles, chromosomes, threshold,
                          flag_select, flag_skip, quality, identity,
                          process_fingerprint=None):
    n = len(handles)
    totals = zeros(n)
    tempfiles = [{'+': tempfile.SpooledTemporaryFile("r+w"),
                  '-': tempfile.SpooledTemporaryFile("r+w")}
                 for chromosome in chromosomes]
    minusstrand_expression_cache = {}
    minusstrand_fingerprint = [dict() for i in range(n)]
    reads = [get_next_mapped_read(handle, flag_select, flag_skip, quality,
                                  identity) for handle in handles]
    while True:
        loci = [(read.rname, read.pos) for read in reads if read]
        if not loci:
            break
        chromosome, position = min(loci)
        plusstrand_expression = zeros(n)
        for i in range(n):
            read = reads[i]
            plusstrand_fingerprint = set()
            while True:
                if not read:
                    break
                elif read.rname > chromosome:
                    break
                elif read.pos > position:
                    break
                if process_fingerprint:
                    # remove PCR duplicates using fingerprints
                    fingerprint = find_fingerprint(read)
                count = find_tag_count(read)
                if read.flag & 16:  # minus strand
                    length = 0
                    for operation, value in read.cigar:
                        if operation in (0, 2, 3):
                            length += value
                    tss = position + length - 1
                    key = (chromosome, tss)
                    if not key in minusstrand_expression_cache:
                        minusstrand_expression_cache[key] = zeros(n)
                    if process_fingerprint:
                        if not key in minusstrand_fingerprint[i]:
                            minusstrand_fingerprint[i][key] = set()
                        if fingerprint in minusstrand_fingerprint[i][key]:
                            read = get_next_mapped_read(handles[i], flag_select, flag_skip, quality, identity)
                            continue
                        else:
                            minusstrand_fingerprint[i][key].add(fingerprint)

                    minusstrand_expression_cache[key][i] += count
                else:  # plus strand
                    if process_fingerprint:
                        if fingerprint in plusstrand_fingerprint:
                            read = get_next_mapped_read(handles[i], flag_select, flag_skip, quality, identity)
                            continue
                        else:
                            plusstrand_fingerprint.add(fingerprint)
                    plusstrand_expression[i] += count
                read = get_next_mapped_read(handles[i], flag_select, flag_skip,
                                            quality, identity)
            reads[i] = read
        # First write out previous positions on the minus strand
        loci = minusstrand_expression_cache.keys()
        loci.sort()
        for locus in loci:
            if locus > (chromosome, position):
                break
            values = minusstrand_expression_cache[locus]
            del minusstrand_expression_cache[locus]
            totals += values
            if max(values) < threshold:
                continue
            write_tempfile_line(tempfiles[locus[0]]['-'], locus[1], values)
        # Then write out the current position on the plus strand,
        # followed by the minus strand
        for strand in '+-':
            if strand == '+':
                values = plusstrand_expression
            elif strand == '-':
                if not position in minusstrand_expression_cache:
                    break
                values = minusstrand_expression_cache[position]
                del minusstrand_expression_cache[position]
            totals += values
            if max(values) < threshold:
                continue
            write_tempfile_line(tempfiles[chromosome][strand], position, values)
    # Finally write any remaining positions on the minus strand
    loci = minusstrand_expression_cache.keys()
    loci.sort()
    for locus in loci:
        values = minusstrand_expression_cache[locus]
        del minusstrand_expression_cache[locus]
        totals += values
        if max(values) < threshold:
            continue
        chromosome, position = locus
        write_tempfile_line(tempfiles[chromosome]['-'], position, values)
    return tempfiles, totals


def run(input_filenames, output_filename, name, threshold, flag_select,
        flag_skip, quality, identity, process_fingerprint):
    assembly, chromosomes = None, None
    handles = []
    for filename in input_filenames:
        handle = open_input_file(filename)
        handles.append(handle)
        try:
            current_assembly = find_assembly(handle)
        except ValueError, message:
            print "Inconsistent genome assembly or species found in %s: %s" % (filename, message)
            sys.exit(2)
        if assembly is None:
            assembly = current_assembly
        elif assembly != current_assembly:
            print "Input files have inconsistent genome assemblies"
            sys.exit(2)
        if chromosomes is None:
            chromosomes = [row['SN'] for row in handle.header['SQ']]
        elif chromosomes != [row['SN'] for row in handle.header['SQ']]:
            print "Input files have inconsistent chromosome sets"
            sys.exit(2)
    if name is None:
        name = assembly
    conditions = [find_condition(input_filename, handle)
                  for input_filename, handle in zip(input_filenames, handles)]
    tempfiles, totals = write_temporary_files(handles, chromosomes, threshold,
                                              flag_select, flag_skip, quality,
                                              identity, process_fingerprint)
    for handle in handles:
        handle.close()
    # Write the final output
    output = open_output_file(output_filename)
    write_meta_data(output, input_filenames, conditions, assembly, chromosomes,
                    threshold, flag_select, flag_skip, quality, identity,
                    process_fingerprint)
    merge(output, tempfiles, chromosomes, name, totals)
    output.close()

def find_fingerprint(read):
    """find and return the fingerprint id from the read's name"""
    return read.qname.split('FP:')[1].split(';')[0]

def version():
    print """\
Level-1 promoter clustering version %s.
Copyright (C) 2010-2015 Michiel JL de Hoon (mdehoon@gsc.riken.jp),
                        Mickael Mendez (mickael.mendez@riken.jp),
RIKEN Center for Life Science Technologies, Yokohama, Japan.
    """ % __version__


def usage():
    version()
    print """\
Usage:

python level1.py -o <level1.osc[.gz|.bz2]>
       [-t <threshold>] [-n <name>] [-f <flag>] [-F <flag>] [-q <quality>] [-i <identity>] [<fingerprint>]
       samfile1.[sam|bam] samfile2.[sam|bam] ...

     Options:
        -o <level1.osc[.gz|.bz2]>
                           name of the output file with the level-1
                           promoters in OSC Table format; if the file name
                           ends in .gz or .bz2, a gzipped or bzipped file
                           is written.
        -t <threshold>     the minimum (absolute) tag count required in
                           at least one of the experimental conditions
                           (default: 1)
        -f <flag_select>   SAM/BAM flag to filter for. Use either a decimal
                           integer, or a hexadecimal value starting with 0x.
                           (default: 0, keep everything).
        -F <flag_skip>     SAM/BAM flag to filter against. Use either a decimal
                           integer, or a hexadecimal value starting with 0x.
                           (default: 4, unmapped)
        -q <quality>       quality cutoff; all reads with a mapping quality
                           less than this value are ignored (default: 0)
        -i <identity>      identity cutoff; all reads with a percentage
                           identical nucleotides to the genome less than this
                           value are ignored (default: 0)
        -n <name>          the identifier to use in promoter names, which
                           are formatted as L1_name_chromosome_strand_position
                           (default: the assembly name as found in the input
                           files)
        <fingerprint>      Flag to remove PCR duplicate using fingerprints.
                           If this flag is on, the program will search the
                           fingerprint id in the read's name with the following
                           pattern: "FP:XXX;" where "XXX" is the fingerprint id
        <samfile.sam>      names of the input files containing the CAGE tag
                           mapping in SAM or BAM format.


This script creates level-1 promoters from the SAM/BAM mapping files.
"""


def main(argv):
    import getopt
    try:
        opts, args = getopt.getopt(argv, "vho:t:f:F:q:i:n:",
                                   ["version", "help", "output", "threshold",
                                    "flag_select", "flag_skip", "quality",
                                    "identity", "name", "fingerprint"])
    except getopt.GetoptError, error:
        print str(error)
        usage()
        sys.exit(2)
    output_filename = None
    flag_select = 0
    flag_skip = 4
    quality = 0.0
    identity = 0.0
    threshold = 1.0
    name = None
    process_fingerprint = None
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-v", "--version"):
            version()
            sys.exit()
        elif opt in ("-o", "--output"):
            output_filename = arg
        elif opt in ("-t", "--threshold"):
            threshold = float(arg)
        elif opt in ("-f", "--flag_select"):
            flag_select = int(arg, 0)
        elif opt in ("-F", "--flag_skip"):
            flag_skip = int(arg, 0)
        elif opt in ("-q", "--quality"):
            quality = float(arg)
        elif opt in ("-i", "--identity"):
            identity = float(arg) / 100.0
        elif opt in ("-n", "--name"):
            name = arg
        elif opt == "--fingerprint":
            process_fingerprint = True
    if not output_filename:
        usage()
        sys.exit(2)
    input_filenames = args
    run(input_filenames, output_filename, name, threshold, flag_select,
        flag_skip, quality, identity, process_fingerprint)

if __name__ == "__main__":
    main(sys.argv[1:])
