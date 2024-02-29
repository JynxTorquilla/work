#!/usr/bin/env python

import sys
import time
import os
from numpy import *


__version__ = "2013.06.27"


def analyze_meta_data(handles):
    assembly = None
    chromosomes = None
    samples = {'raw': [], 'norm': []}
    for handle in handles:
        if assembly==None:
            assembly = handle.assembly
        else:
            if assembly != handle.assembly:
                message = "Different genome assembly found (found %s, expected %s)\n" % (handle.assembly, assembly)
                sys.stderr.write(message)
                sys.exit(2)
        if chromosomes==None:
            chromosomes = handle.chromosomes
        else:
            if chromosomes != handle.chromosomes:
                message = "Different chromosome sets found (found %s, expected %s)\n" % (",".join(handle.chromosomes), ",".join(chromosomes))
                sys.stderr.write(message)
                sys.exit(2)
        # current line is the header
        raw_samples = take(handle.samples, handle.raw_indices)
        tpm_samples = take(handle.samples, handle.tpm_indices)
        samples['raw'].extend(raw_samples)
        samples['norm'].extend(tpm_samples)
    return assembly, chromosomes, samples

def write_meta_data(output, assembly, chromosomes, input_filenames, samples):
    # General meta-data
    output.write("##FileFormat = OSCtable\n")
    output.write("##ProtocolREF = OP-CAGE-Tag-Clustering-v1.1\n")
    output.write("##ParameterValue[level2_script_version] = %s\n" % __version__)
    if assembly:
        output.write("##ParameterValue[genome_assembly] = %s\n" % assembly)
    output.write("##ChromosomeNameOrder = %s\n" % ",".join(chromosomes))
    output.write("##Date = %s\n" % time.strftime("%Y-%m-%d"))
    output.write('##ContactName = LSA support staff\n')
    output.write('##ContactEmail = lsa-support@gsc.riken.jp\n')
    # Input file names
    for input_filename in input_filenames:
        basename = os.path.basename(input_filename)
        output.write("##InputFile = %s\n" % basename)
    # Column variables
    output.write('##ColumnVariable[id] = identifier of this tss cluster\n')
    output.write('##ColumnVariable[chrom] = name of the chromosome, identical to what is written in genome assembly\n')
    output.write('##ColumnVariable[start.0base] = start genomic coordinate of a cluster, using zero-based coordinates, always same or smaller than "end"\n')
    output.write('##ColumnVariable[end] = end genomic coordinate of this cluster\n')
    output.write('##ColumnVariable[strand] = genomic strand of this cluster\n')
    output.write('##ColumnVariable[pos] = genomic coordinate of the most highly expressed position of this cluster\n')
    for sample in samples['raw']:
        output.write('##ColumnVariable[raw.%s] = rescued tag count of sample %s\n' % (sample, sample))
    for sample in samples['norm']:
        output.write('##ColumnVariable[norm.%s] = normalized expression by total tag count for sample %s\n' % (sample, sample))
    # Write header line
    output.write("id\tchrom\tstart.0base\tend\tstrand\tpos")
    for word in samples['raw']:
        output.write("\traw.%s" % word)
    for word in samples['norm']:
        output.write("\tnorm.%s" % word)
    output.write("\n")


def make_line(library, chromosome, strand, tss, start, end, raw, tpm):
    line = "L2_%s_%s_%s_%d\t%s\t%d\t%d\t%s\t%d" % (library, chromosome, strand, tss, chromosome, start, end, strand, tss)
    for value in raw:
        line += "\t" + str(value)
    for value in tpm:
        line += "\t" + str(value)
    line += "\n"
    return line

class TSS(tuple):
    """tuple of chromosome_index, strand, position, with an attribute values
    containing the expression values"""
    def get_chromosome(self):
        return self[0]
    chromosome = property(get_chromosome)
    def get_strand(self):
        return self[1]
    strand = property(get_strand)
    def get_position(self):
        return self[2]
    position = property(get_position)

class Level1File:
    def __init__(self, filename):
        self.assembly = None
        sys.stderr.write("Opening %s\n" % filename)
        if filename.endswith(".gz"):
            import gzip
            self.handle = gzip.open(filename)
        elif filename.endswith(".bz2"):
            import bz2
            self.handle = bz2.BZ2File(filename)
        else:
            self.handle = open(filename)
        for line in self.handle:
            if not line.startswith("#"):
                break
            if line.startswith("##ParameterValue[genome_assembly]"):
                key, value = line.strip().split("=", 1)
                self.assembly = value.strip()
            elif line.startswith("##ChromosomeNameOrder ="):
                key, value = line.strip().split("=", 1)
                self.chromosomes = value.strip().split(",")
        # Current line is the header file
        words = line.strip().split("\t")
        if words[0]!="id":
            raise ValueError("Expected 'id' as first column")
        if words[1]!="chrom":
            raise ValueError("Expected 'chrom' as second column")
        if words[2] not in ("start", "start.0base"):
            raise ValueError("Expected 'start' or 'start.0base' as third column")
        if words[3]!="end":
            raise ValueError("Expected 'end' as fourth column")
        if words[4]!="strand":
            raise ValueError("Expected 'strand' as fifth column")
        self.samples = []
        self.raw_indices = []
        self.tpm_indices = []
        words = words[5:]
        n = len(words)
        for i in range(n):
            word = words[i]
            if word.startswith("raw."):
                sample = word[4:]
                self.samples.append(sample)
                self.raw_indices.append(i)
            elif word.startswith("norm."):
                sample = word[5:]
                self.samples.append(sample)
                self.tpm_indices.append(i)
            else:
                raise ValueError("Unknown column %s" % word)
    def __iter__(self):
        return self
    def next(self):
        line = self.handle.next()
        words = line.strip().split('\t')
        position = int(words[2])
        assert int(words[3])==position+1
        chromosome = self.chromosomes.index(words[1])
        strand = words[4]
        tss = TSS([chromosome, strand, position])
        tss.values = map(float, words[5:])
        return tss
    def close(self):
        self.handle.close()

def create_clusters(handles, chromosomes, separation, tpm_indices):
    chromosome, strand, representative, start, end = None, None, None, None, None
    nhandles = len(handles)
    tsss = [handle.next() for handle in handles]
    values = zeros(sum([len(tss.values) for tss in tsss]))
    current_values = zeros(sum([len(tss.values) for tss in tsss]))
    while True:
        candidates = [tss for tss in tsss if tss]
        if not candidates:
            break
        closest = min(candidates)
        if closest.chromosome!=chromosome or closest.strand!=strand or closest.position >= end + separation-1:
            yield (chromosome, strand, representative, start, end, values)
            chromosome, strand, start = closest
            end = start+1
            values[:] = 0.0
            maximum = 0
        current_values[:] = 0.0
        i = 0
        for j in range(nhandles):
            tss = tsss[j]
            handle = handles[j]
            n = len(handle.samples)
            if tss==closest:
                # This checks the locus only, not the expression values
                current_values[i:i+n] += tss.values
                try:
                    tss = handle.next()
                except StopIteration:
                    tss = None
                tsss[j] = tss
            i += n
        sum_tpm = sum(take(current_values, tpm_indices))
        if sum_tpm > maximum:
            maximum = sum_tpm
            representative = closest.position
        end = closest.position+1
        values += current_values
    yield (chromosome, strand, representative, start, end, values)

def gather_tpm_indices(handles):
    # Find the indices at which the tags-per-millions are stored
    offset = 0
    tpm_indices = []
    for handle in handles:
        for index in handle.tpm_indices:
            tpm_indices.append(index+offset)
        offset += len(handle.raw_indices) + len(handle.tpm_indices)
    tpm_indices = array(tpm_indices, int)
    return tpm_indices

def gather_raw_indices(handles):
    # Find the indices at which the raw counts are stored
    offset = 0
    raw_indices = []
    for handle in handles:
        for index in handle.raw_indices:
            raw_indices.append(index+offset)
        offset += len(handle.raw_indices) + len(handle.raw_indices)
    raw_indices = array(raw_indices, int)
    return raw_indices

def open_input_files(input_filenames):
    input_files = []
    for input_filename in input_filenames:
        try:
            input_file = Level1File(input_filename)
        except IOError:
            message = "Failed to open the input file %s\n" % input_filename
            sys.stderr.write(message)
            sys.exit(2)
        input_files.append(input_file)
    return input_files

def run(input_filenames, output_filename, removed_filename, separation, threshold, identifier):
    input_files = open_input_files(input_filenames)
    assembly, chromosomes, samples = analyze_meta_data(input_files)
    if identifier==None:
        if assembly:
            identifier = assembly
        else:
            identifier = ""
    output_file = open(output_filename, 'w')
    write_meta_data(output_file, assembly, chromosomes, input_filenames, samples)
    if removed_filename:
        removed_file = open(removed_filename, 'w')
    else:
        removed_file = None
    tpm_indices = gather_tpm_indices(input_files)
    raw_indices = gather_raw_indices(input_files)
    clusters = create_clusters(input_files, chromosomes, separation, tpm_indices)
    for cluster in clusters:
        (chromosome_index, strand, tss, start, end, values) = cluster
        if chromosome_index==None or strand==None:
            continue
        chromosome = chromosomes[chromosome_index]
        raw = take(values, raw_indices)
        tpm = take(values, tpm_indices)
        if max(tpm) >= threshold:
            line = make_line(identifier, chromosome, strand, tss, start, end, raw, tpm)
            output_file.write(line)
        elif removed_file:
            line = make_line(identifier, chromosome, strand, tss, start, end, raw, tpm)
            removed_file.write(line)
    for input_file in input_files:
        input_file.close()
    output_file.close()
    if removed_file:
        removed_file.close()


def version():
    print """\
Level-2 promoter clustering version %s.
Copyright (C) 2010-2012 Michiel JL de Hoon (mdehoon@gsc.riken.jp),
RIKEN Omics Science Center, Yokohama, Japan.
    """ % __version__

def usage():
    version()
    print """\
Usage:

python level2.py -o <level2.osc> [-r <removed.osc>] [-t <threshold>]
                 [-s <separation>] [-n <name>]
                 <level1file1.osc[.gz|.bz2]> <level1file2.osc[.gz|.bz2]> ...

     Options:
        -o <level2.osc>    name of the output file with the level-2
                           promoters in OSC Table format.
        -r <removed.osc>   name of the output file to which the level-2
                           promoters are written that did not pass the
                           threshold.
        -t <threshold>     the minimum tags-per-million required in at least
                           one of the experimental conditions
                           (default: 10 tpm)
        -s <separation>    the minimum distance in base pairs between
                           adjacent clusters on the genome.
                           (default: 20)
        -n <name>          the identifier to use in promoter names, which
                           are formatted as L2_name_chromosome_strand_position
                           (default: the assembly name as found in the input
                           files)
        <level1file.osc[.gz|.bz2]>
                           names of the input files containing the level-1
                           promoters in OSC Table format. Optionally, the
                           input file can be gzipped or bzipped.

This script performs level-1 to level-2 promoter clustering.
"""

def main(argv):
    import getopt
    try:
        opts, args = getopt.getopt(argv, "vho:r:t:s:n:", ["version", "help", "name"])
    except getopt.GetoptError, error:
        sys.stderr.write("%s\n" % error)
        usage()
        sys.exit(2)
    input_filenames = None
    output_filename = None
    removed_filename = None
    separation = 20
    threshold = 10.0
    name = None
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-v", "--version"):
            version()
            sys.exit()
        elif opt=="-o":
            output_filename = arg
        elif opt=="-r":
            removed_filename = arg
        elif opt=="-t":
            threshold = float(arg)
        elif opt=="-s":
            separation = int(arg)
        elif opt in ("-n", "--name"):
            name = arg
    input_filenames = args
    if not input_filenames or not output_filename:
        usage()
        sys.exit(2)
    run(input_filenames, output_filename, removed_filename, separation, threshold, name)


if __name__ == "__main__":
    main(sys.argv[1:])
