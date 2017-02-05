"""
   lasig.py - extract signature alleles
   Copyright (C) 2013 Giulio Genovese

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   Written by Giulio Genovese <giulio.genovese@gmail.com>
"""

import argparse, sys, resource, gzip, csv, vcf, os, numpy as np, warnings

def lasig(argv):
    parser = argparse.ArgumentParser(description='latools.py sig: extract signature alleles', add_help=False, usage='latools.py sig [options] <sig.csv> [in.vcf]')
    parser.add_argument('sig', metavar='sig.csv', type=str, help='ancestry signatures file (can be gzipped)')
    parser.add_argument('i', metavar='in.vcf', nargs='?', type=str, default='/dev/stdin', help='VCF file with ancestry signatures [stdin]')
    parser.add_argument('-d', metavar='DIR', type=str, default='.', help='directory for output files')
    parser.add_argument('-z', action='store_true', help='whether to output compressed CSV files')

    try:
        parser.error = parser.exit
        args = parser.parse_args(argv)
    except SystemExit:
        parser.print_help()
        exit(2)

    if args.sig == '-':
        sigfile = sys.stdin
        csv_reader = csv.reader(sigfile)
        row = next(csv_reader)
    else:
        try:
            sigfile = gzip.open(args.sig,'rb')
            csv_reader = csv.reader(sigfile)
            row = next(csv_reader)
        except IOError:
            sigfile.close()
            sigfile = open(args.sig,'rt')
            csv_reader = csv.reader(sigfile)
            row = next(csv_reader)
    
    # build ancestry signatures dictionary
    sig = {row[0] + ':' + row[1] : row[3:]}
    for row in csv_reader:
        sig[row[0] + ':' + row[1]] = row[3:]
    
    numpop = len(row)-3
    
    # open the input VCF file
    vcf_reader = vcf.Reader(filename=args.i)
    numsmpl = len(vcf_reader.samples)
    
    if resource.getrlimit(resource.RLIMIT_NOFILE)[0] < numsmpl:
        resource.setrlimit(resource.RLIMIT_NOFILE, (numsmpl+5,numsmpl+5))
    
    # open output files with ancestry signatures
    outfiles = []
    for i in range(numsmpl):
        if args.z:
            outfiles.append(csv.writer(gzip.open(args.d + os.sep + vcf_reader.samples[i] + '.sig.gz','wb'),delimiter='\t'))
        else:
            outfiles.append(csv.writer(open(args.d + os.sep + vcf_reader.samples[i] + '.sig','wt'),delimiter='\t'))
    
    phased = np.zeros(numsmpl,dtype=int)
    mat = np.zeros((numsmpl,numpop),dtype=int)
    pat = np.zeros((numsmpl,numpop),dtype=int)
    
    # flag for detection of phased data
    flag = True
    
    # record = next(vcf_reader)
    # iterate over the lines in the VCF file
    for record in vcf_reader:
        key = record.CHROM + ':' + str(record.POS)
        if key in sig:
            value = sig[key]
            for i in range(numpop):
                if value[i] in record.alleles:
                    for j in range(len(record.alleles)):
                        if value[i]==record.alleles[j]:
                            allele = j
                    for j in range(numsmpl):
                        if record.samples[j].called:
                            if record.samples[j].phased:
                                phased[j] = True
                            elif flag:
                                warnings.warn('Unphased data detected',Warning)
                                flag = False
                            if int(record.samples[j].gt_alleles[0]) == allele:
                                mat[j][i] = 1
                            if int(record.samples[j].gt_alleles[1]) == allele:
                                pat[j][i] = 1
        for j in range(numsmpl):
            if np.any(mat[j]) or np.any(pat[j]):
                if record.ID:
                    ID = record.ID
                else:
                    ID = '.'
                outfiles[j].writerow([record.CHROM, record.POS, ID, phased[j]] + list(mat[j]) + list(pat[j]))
                phased[j]=0
                mat[j]*=0
                pat[j]*=0
