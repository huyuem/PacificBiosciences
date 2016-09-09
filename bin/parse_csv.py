"""
this script parse the sample.csv file for pacbio isoseq pipeline
and generate bash source files
"""
import sys
import os
import csv
from collections import defaultdict

if len(sys.argv) != 2:
    print("usage: {} sample.csv".format(sys.argv[1]))
    sys.exit(1)
sample_info = defaultdict(dict)
with open(sys.argv[1], 'rb') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        sample_name = '.'.join(['-'.join([row['year'], row['month'], row['day']]),
                                row['genome'], row['tissue'], row['treatment'], row['size_bin']])
        if not os.path.exists(row['cell_path']):
            sys.stderr.write("cannot access {}, please double check\n".format(row['cell_path']))
            continue
        if row['barcode_path'] != 'NA' and not os.path.exists(row['barcode_path']):
            sys.stderr.write("cannot access {}, please double check\n".format(row['barcode_path']))
            continue
        sample_info[sample_name] = row
# write to a bash source file
print('declare -xi SampleSize={}'.format(len(sample_info)))

print('declare -xa SampleNames=({})'.format(
    ' '.join(["'" + sample + "'" for sample, _ in sample_info.items()])))

print('declare -xa Tissues=({})'.format(
    ' '.join(["'" + row['tissue'] + "'" for sample, row in sample_info.items()])))

print('declare -xa SampleCellPaths=({})'.format(
    ' '.join(["'" + row['cell_path'] + "'" for sample, row in sample_info.items()])))

print('declare -xa Treatments=({})'.format(
    ' '.join(["'" + row['treatment'] + "'" for sample, row in sample_info.items()])))

print('declare -xa SizeBins=({})'.format(
    ' '.join(["'" + row['size_bin'] + "'" for sample, row in sample_info.items()])))

print('declare -xa BarcodeFiles=({})'.format(
    ' '.join(["'" + row['barcode_path'] + "'" for sample, row in sample_info.items()])))

print('declare -xa BarcodeNumbers=({})'.format(
    ' '.join(["'" + row['barcode'] + "'" for sample, row in sample_info.items()])))

print('declare -xa Genomes=({})'.format(
    ' '.join(["'" + row['genome'] + "'" for sample, row in sample_info.items()])))
