#! /usr/bin/env python

import argparse
import sys

parser = argparse.ArgumentParser(description='Simple tool to extract UMIs from the 5\' or 3\' end of fastq reads and '
                                             'add them to the end of the first word of the sequence name, preceded by '
                                             'an underscore (_)')
parser.add_argument('-i', '--infile', type=argparse.FileType('rU'), default=sys.stdin,
                    help='Input fastq file (default: STDIN)')
parser.add_argument('-o', '--outfile', type=argparse.FileType('w'), default=sys.stdout,
                    help='Output fastq file (default: STDOUT)')
parser.add_argument('-u', '--umilen', type=int, default=5, help='Length of UMI (default: 5)')
parser.add_argument('-t', '--trim', type=int, default=4, help='Extra bases to trim after/before the UMI (default: 4)')
parser.add_argument('-m', '--minlen', type=int, default=10,
                    help='Minimum remaining read length after trimming (default: 10)')
parser.add_argument('-5', '--fiveprime', action='store_true', help='UMI is at 5\' end of read (default is 3\')')
parser.add_argument('-v', '--verbose', action='store_true', help='Report statistics to STDERR')
# parser.add_argument('-p', '--numproc', type=int, default=1, help='Number of processes')
opts = parser.parse_args()

if opts.umilen <= 0:
    raise ValueError('UMILEN must be at least 1 nt (currently %d)' % opts.umilen)

if opts.trim < 0:
    raise ValueError('TRIM must be non-negative (currently %d)' % opts.trim)

if opts.minlen <= 0:
    raise ValueError('MINLEN must be at least 1 nt (currently %d)' % opts.minlen)

numreported = 0
numtooshort = 0
while True:
    try:
        name = next(opts.infile).strip()
        seq = next(opts.infile).strip()
        plus = next(opts.infile).strip()
        qual = next(opts.infile).strip()
    except StopIteration:
        if opts.verbose:
            sys.stderr.write('umiextract reported %d reads; %d rejected for being shorter than %d nt\n'
                             % (numreported, numtooshort, opts.minlen))
        break
    if len(seq) >= opts.minlen + opts.trim + opts.umilen:
        name_part = name.partition(' ')
        if opts.fiveprime:
            name_new = name_part[0] + '_' + seq[opts.trim:opts.trim+opts.umilen] + name_part[1] + name_part[2]
            seq_new = seq[opts.trim + opts.umilen:]
            qual_new = qual[opts.trim + opts.umilen:]
        else:
            umistart = len(seq)-opts.trim-opts.umilen
            name_new = name_part[0] + '_' + seq[umistart:umistart+opts.umilen] + name_part[1] + name_part[2]
            seq_new = seq[:umistart]
            qual_new = qual[:umistart]
        opts.outfile.write('\n'.join([name_new, seq_new, plus, qual_new])+'\n')
        numreported += 1
    else:
        numtooshort += 1
