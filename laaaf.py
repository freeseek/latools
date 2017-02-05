"""
   laaaf.py - estimate ancestral allele frequencies
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

import os, argparse, sys, ctypes, gzip, numpy as np, vcf, collections, re, functools, scipy.optimize, math
import libla

def laaaf(argv):
    libpath = os.path.abspath(os.path.dirname(sys.argv[0]) + os.sep + 'latools.so')
    parser = argparse.ArgumentParser(description='latools.py map: admixture-map SNPs', add_help=False, usage='latools.py map [options] <la.bg> [in.vcf]')
    parser.add_argument('la', metavar='la.bg', type=str, help='local ancestry file in bedgraph format (can be gzipped)')
    parser.add_argument('i', metavar='in.vcf', nargs='?', type=str, default='/dev/stdin', help='input VCF file (can be gzipped) [stdin]')
    parser.add_argument('-o', metavar='FILE', type=argparse.FileType('w'), default=sys.stdout, help='output VCF file [stdout]')
    parser.add_argument('-e', metavar='FILE', type=argparse.FileType('r'), help='list of training samples of European ancestry')
    parser.add_argument('-a', metavar='FILE', type=argparse.FileType('r'), help='list of training samples of African ancestry')
    parser.add_argument('-n', metavar='FILE', type=argparse.FileType('r'), help='list of training samples of Native American ancestry')
    parser.add_argument('-s', action='store_true', help='output VCF with no samples genotypes')
    parser.add_argument('-g', metavar='INT', type=int, default=20, help='minimum number of genotypes [20]')
    parser.add_argument('-p', metavar='FILE', type=str, default=libpath, help='library for likelihood computations [' + libpath + ']')

    try:
        parser.error = parser.exit
        args = parser.parse_args(argv)
    except SystemExit:
        parser.print_help()
        exit(2)

    clib = ctypes.CDLL(args.p)

    # open the BedGraph file and stores it in a numpy array
    la_samples, chrom, start, end, la = libla.getla(args.la,clib)

    # compute ancestry proportions
    la_ancXX = np.sum(la>0,axis=0).astype(ctypes.c_double)
    la_ancEE = np.sum(la==1,axis=0).astype(ctypes.c_double) / la_ancXX
    la_ancEA = np.sum(la==2,axis=0).astype(ctypes.c_double) / la_ancXX
    la_ancAA = np.sum(la==3,axis=0).astype(ctypes.c_double) / la_ancXX
    la_ancEN = np.sum(la==4,axis=0).astype(ctypes.c_double) / la_ancXX
    la_ancAN = np.sum(la==5,axis=0).astype(ctypes.c_double) / la_ancXX
    la_ancNN = np.sum(la==6,axis=0).astype(ctypes.c_double) / la_ancXX

    # add training samples from ancestral populations
    tr_samples = list(la_samples)
    tr_ancEE, tr_ancEA, tr_ancAA, tr_ancEN, tr_ancAN, tr_ancNN = la_ancEE, la_ancEA, la_ancAA, la_ancEN, la_ancAN, la_ancNN
    # tr_ancE, tr_ancA, tr_ancN = la_ancE, la_ancA, la_ancN
    trfile = [args.e, args.a, args.n]
    for i in range(3):
        if trfile[i] != None:
            tr_samples += [sample.rstrip('\r\n') for sample in trfile[i].readlines()]
            trfile[i].close()
            tr_ancEE = np.concatenate((tr_ancEE,np.tile(1.0 if i==0 else 0.0,len(tr_samples)-len(tr_ancEE))))
            tr_ancEA = np.concatenate((tr_ancEA,np.tile(0.0,len(tr_samples)-len(tr_ancEA))))
            tr_ancAA = np.concatenate((tr_ancAA,np.tile(1.0 if i==1 else 0.0,len(tr_samples)-len(tr_ancAA))))
            tr_ancEN = np.concatenate((tr_ancEN,np.tile(0.0,len(tr_samples)-len(tr_ancEN))))
            tr_ancAN = np.concatenate((tr_ancAN,np.tile(0.0,len(tr_samples)-len(tr_ancAN))))
            tr_ancNN = np.concatenate((tr_ancNN,np.tile(1.0 if i==2 else 0.0,len(tr_samples)-len(tr_ancNN))))

    # open the input VCF file
    vcf_reader = vcf.Reader(filename=args.i)
    vcf_samples = vcf_reader.samples

    # add INFO fields to the VCF structure
    _Info = collections.namedtuple('Info', ['id', 'num', 'type', 'desc'])
    vcf_reader.infos['EUR_AAF'] = _Info('EUR_AAF','1','Float','Estimated Ancestral Allele Frequency for EUR Samples')
    vcf_reader.infos['AFR_AAF'] = _Info('AFR_AAF','1','Float','Estimated Ancestral Allele Frequency for AFR Samples')
    vcf_reader.infos['NAT_AAF'] = _Info('NAT_AAF','1','Float','Estimated Ancestral Allele Frequency for AMR Samples')

    # this code is extremely slow, but does the job
    tr2vcf = [(x, y) for x in range(len(tr_samples)) for y in range(len(vcf_samples)) if tr_samples[x] == vcf_samples[y]]
    tr2vcf_loc = [x[0] for x in tr2vcf]
    vcf2tr_loc = [x[1] for x in tr2vcf]
    la2vcf = [(x, y) for x in range(len(la_samples)) for y in range(len(vcf_samples)) if la_samples[x] == vcf_samples[y]]
    la2vcf_loc = [x[0] for x in la2vcf]
    vcf2la_loc = [x[1] for x in la2vcf]

    # iterate over the lines in the VCF file
    for record in vcf_reader:
        gl0, gl1, gl2, ind = libla.getGL(record,vcf2tr_loc,clib)
        n = len(ind)

        if n < args.g: vcf_writer.write_record(record); continue
        frqE = frqA = frqN = (sum(gl1)/2+sum(gl2))/len(gl0)

        tr_ind = [tr2vcf_loc[i] for i in ind]
        ancEE = (ctypes.c_double * n)(); ancEE[:] = tr_ancEE[tr_ind]
        ancEA = (ctypes.c_double * n)(); ancEA[:] = tr_ancEA[tr_ind]
        ancAA = (ctypes.c_double * n)(); ancAA[:] = tr_ancAA[tr_ind]
        ancEN = (ctypes.c_double * n)(); ancEN[:] = tr_ancEN[tr_ind]
        ancAN = (ctypes.c_double * n)(); ancAN[:] = tr_ancAN[tr_ind]
        ancNN = (ctypes.c_double * n)(); ancNN[:] = tr_ancNN[tr_ind]

        def f(x,y,z): return -libla.sumlog(libla.lkl3(n,x,y,z,ancEE,ancEA,ancAA,ancEN,ancAN,ancNN,gl0,gl1,gl2,clib),clib)

        # estimate ancestral allele frequencies
        for i in range(3):
            frqE = scipy.optimize.fminbound(functools.partial(f,y=frqA,z=frqN),args.m,1-args.m)
            frqA = scipy.optimize.fminbound(functools.partial(f,frqE,z=frqN),args.m,1-args.m)
            frqN = scipy.optimize.fminbound(functools.partial(f,frqE,frqA),args.m,1-args.m)

        record.add_info('EUR_AAF',round(1e4*frqE)/1e4)
        record.add_info('AFR_AAF',round(1e4*frqA)/1e4)
        record.add_info('NAT_AAF',round(1e4*frqN)/1e4)
