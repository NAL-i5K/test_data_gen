#!/usr/env python
import argparse
import re
import os
import errno

parser = argparse.ArgumentParser(description='Generate random datasets from real data.')
parser.add_argument('-a', type=str, nargs=1,
                   help='Genome assembly file')
parser.add_argument('-g', type=str, nargs=1,
                    help='OGS Gff3')
parser.add_argument('-c', type=str, nargs=1,
                    help='OGS CDS fasta')
parser.add_argument('-d', type=str, nargs=1,
                    help='OGS cDNA fasta')
parser.add_argument('-p', type=str, nargs=1,
                    help='OGS peptide fasta')
parser.add_argument('-o', type=str, nargs=1,
                    help='Output directory')
#args = parser.parse_args()
args = parser.parse_args()


# Part 1 generate unique ID for new project
import subprocess
id = subprocess.check_output(["date", "+%s%N"])
id.rstrip('\r\n')
id = id.replace("\n", "")
print id
# Set up output path
out_dir = genome = args.o[0]   # Output directory
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

mkdir_p(out_dir)

genome_out_file = out_dir + '/' + id + '_assembly.fa'
gff_out_file = out_dir + '/' + id + '_ogs.gff'
cds_out_file = out_dir + '/' + id + '_ogs_cds.fa'
cdna_out_file = out_dir + '/' + id + '_ogs_cdna.fa'
pep_out_file = out_dir + '/' + id + '_ogs_pep.fa'
# Set up lookup table files

genome_fa_lookup_file = out_dir + '/' + id + '_assembly_lookup.txt'
ogs_gff_lookup_file = out_dir + '/' + id + '_ogs_gff_lookup.txt'
cds_fa_lookup_file = out_dir + '/' + id + '_ogs_cds_fa_lookup.txt'
cdna_fa_lookup_file = out_dir + '/' + id + '_ogs_cdna_fa_lookup.txt'
pep_fa_lookup_file = out_dir + '/' + id + '_ogs_pep_fa_lookup.txt'
project_lookup_file = out_dir + '/' + id + '_project_lookup.txt'


# Part 2 
# organize incoming datafiles, make sure that everything is present.
# Required files:
genome = args.a[0]   # Genome assembly
ogs_gff = args.g[0]  # OGS.gff3
ogs_cds = args.c[0]  # OGS_CDS.fa
ogs_cdna = args.d[0] # OGS_cDNA.fa
ogs_pep = args.p[0]  # OGS_peptide.fa
seq_lookup = {}

# Part 3
# open files, and start proccessing. As we process, we need to keep track of the IDs to maintain fidelity.At the end, write two files out:
#
# Part 3a. subroutine for processing fasta files

def process_fasta( input_fa, output_fa, output_table, type):
    seq_count = 0
    print input_fa, output_fa, output_table
    fa_file = open(input_fa, "r")
    fa_out = open(output_fa, "w")
    fa_lookup = open(output_table, "w")
    for x in fa_file:
        m = re.search('>' '(.*)', x)
        if m:
            seq_count += 1
            new_seq_id = id + '_' + type + '_' + str(seq_count)
            fa_out.write('>' + new_seq_id + '\n')
            seq_lookup[m.group(1)] = new_seq_id
            fa_lookup.write( m.group(1) + '\t' + new_seq_id + '\n')
        else:
            fa_out.write(x.rstrip() + '\n')
    fa_file.close()
    fa_out.close()
    fa_lookup.close()
    return [seq_lookup]

# Part 3b Process the fastas and build our table for the gff file
genome_lookup = process_fasta(genome, genome_out_file, genome_fa_lookup_file, 'scaffold')
cds_lookup =  process_fasta(ogs_cds, cds_out_file, cds_fa_lookup_file, 'cds')
cdna_lookup =  process_fasta(ogs_cdna, cdna_out_file, cdna_fa_lookup_file, 'cdna')
pep_lookup =  process_fasta(ogs_pep, pep_out_file, pep_fa_lookup_file, 'pep')

# Part 4 process the gff
def process_gff(in_gff, out_gff, output_table):
    print in_gff, out_gff, output_table
    gff_file = open(in_gff, "r")
    gff_out = open(out_gff, "w")
    gene_lookup = open(output_table, "w")
    seq_count = {}
    gff_lookup = {}
    for x in gff_file:
        m = re.search('ID=([^;]*)', x)
        if m is not None:
            line = x.split()
            seq = line[0]
            x = x.replace(seq, seq_lookup[seq])
            type = line[2]
            if type in seq_count:
                seq_count[type] += 1
            else:
                seq_count[type] = 1
            orig_id = m.group(1).rstrip()
            if orig_id in seq_lookup:
                new_seq_id = seq_lookup[orig_id]
            else:
                new_seq_id = id + '_' + type + '_' + str(seq_count[type])
            x = x.replace(orig_id, new_seq_id)
            seq_lookup[orig_id] = new_seq_id
            m2 = re.search('Parent=([^;]*)?', x)
            if m2 is not None:
                parent = ''
                parent = m2.group(1).rstrip()
                x = x.replace(parent, seq_lookup[parent])
            if orig_id not in gff_lookup:
                gff_lookup[orig_id] = new_seq_id
                gene_lookup.write( orig_id + '\t' + new_seq_id + '\n')
            gff_out.write(x)
    gff_file.close()
    gff_out.close()
    gene_lookup.close()
    return [seq_lookup]

gene_lookup = process_gff(ogs_gff, gff_out_file, ogs_gff_lookup_file)

# Part 6 build the global lookup table 
project_lookup = open(project_lookup_file , "w")
for item in sorted(seq_lookup):
    project_lookup.write(item + '\t' + seq_lookup[item] + '\n')

project_lookup.close()
