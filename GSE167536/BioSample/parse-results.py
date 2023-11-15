#!/usr/bin/env python3
import re
import csv
import argparse
import sys

pat_sn = r'Sample name: (\w+);'
pat_attr = r'/(.+)="(.+)"'

def parse(lines):
    output = []
    cursample = None
    for line in lines:
        if m := re.search(pat_sn, line):
            if cursample is not None:
                output.append(cursample)
            cursample = dict(SampleName=m.group(1))

        if cursample is not None:
            if m := re.search(pat_attr, line):
                cursample.update({m.group(1): m.group(2)})

    
    if cursample is not None:
        output.append(cursample)
    return output

def output_csv(parsed, fd):
    header = {k for r in parsed for k in r.keys()} - {'SampleName'}
    header = ['SampleName'] + list(sorted(header))

    wtr = csv.DictWriter(fd, header)
    wtr.writeheader()
    wtr.writerows(parsed)

def cli():
    p = argparse.ArgumentParser('Parse BioSample text file')
    p.add_argument('input', help='text file from BioSample')
    p.add_argument('output', nargs='?', default=None, help='Output file, or stdout if not present')

    return p.parse_args()

if __name__ == '__main__':
    args = cli()

    with open(args.input) as f:
        parsed = parse(f.readlines())
    
    if args.output:
        with open(args.output, 'wt') as f:
            output_csv(parsed, f)
    else:
        output_csv(parsed, sys.stdout)