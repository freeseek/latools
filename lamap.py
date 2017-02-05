"""
   lamap.py - admixture-map SNPs
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

import os, argparse, sys, ctypes, numpy as np, vcf, collections, functools, scipy.optimize
import libla

# sys.argv = '/home/genovese/Documents/python/latools/latools.py map /home/genovese/Documents/python/latools/la.bg.gz /tmp/in.vcf'.split(' ')
# (args, argv) = parser.parse_known_args(sys.argv[1:])

def lamap(argv):
    libpath = os.path.abspath(os.path.dirname(sys.argv[0]) + os.sep + 'latools.so')
    parser = argparse.ArgumentParser(description='latools.py map: admixture-map SNPs', add_help=False, usage='latools.py map [options] <la.bg> [in.vcf]')
    parser.add_argument('la', metavar='la.bg', type=str, help='local ancestry file in bedgraph format (can be gzipped)')
    parser.add_argument('i', metavar='in.vcf', nargs='?', type=str, default='/dev/stdin', help='input VCF file (can be gzipped) [stdin]')
    parser.add_argument('-o', metavar='FILE', type=argparse.FileType('w'), default=sys.stdout, help='output VCF file [stdout]')
    parser.add_argument('-e', metavar='FILE', type=argparse.FileType('r'), help='list of training samples of European ancestry')
    parser.add_argument('-a', metavar='FILE', type=argparse.FileType('r'), help='list of training samples of African ancestry')
    parser.add_argument('-n', metavar='FILE', type=argparse.FileType('r'), help='list of training samples of Native American ancestry')
    parser.add_argument('-f', action='store_true', help='compute ancestral allele frequencies only')
    parser.add_argument('-r', metavar='FLOAT', type=float, default=0, help='LOD prior for fixed polymorphism')
    parser.add_argument('-s', action='store_true', help='output VCF with no samples genotypes')
    parser.add_argument('-g', metavar='INT', type=int, default=20, help='minimum number of genotypes [20]')
    parser.add_argument('-m', metavar='FLOAT', type=float, default=.002, help='minimum estimated ancestral allele frequency [0.002]')
    parser.add_argument('-t', metavar='FLOAT', type=float, default=.05, help='minimum minor allele frequency to attempt mapping [0.05]')
    parser.add_argument('-l', metavar='INT', type=int, default=10000, help='number of LOD scores computed for each SNP [10000]')
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
    if not args.f:
        vcf_reader.infos['MAP1'] = _Info('MAP1','1','String','Coordinates of first best admixture-mapping')
        vcf_reader.infos['MAP2'] = _Info('MAP2','1','String','Coordinates of second best admixture-mapping')

    vcf_writer = vcf.Writer(args.o,vcf_reader)

    # this code is slow, but does the job
    tr2vcf = [(x, y) for x in range(len(tr_samples)) for y in range(len(vcf_samples)) if tr_samples[x] == vcf_samples[y]]
    tr2vcf_loc = [x[0] for x in tr2vcf]
    vcf2tr_loc = [x[1] for x in tr2vcf]
    la2vcf = [(x, y) for x in range(len(la_samples)) for y in range(len(vcf_samples)) if la_samples[x] == vcf_samples[y]]
    la2vcf_loc = [x[0] for x in la2vcf]

    # iterate over the lines in the VCF file
    for record in vcf_reader:
        # extracts genotype from the VCF file
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

        def f(x,y,z): return -libla.sumlog(libla.lkl3(n,x,y,z,ancEE,ancEA,ancAA,ancEN,ancAN,ancNN,gl0,gl1,gl2,clib),clib) \
            -( (x==0 or x==1) + (y==0 or y==1) + (z==0 or z==1) ) * args.r

        # estimate ancestral allele frequencies
        for i in range(3):
            frqE = scipy.optimize.fminbound(functools.partial(f,y=frqA,z=frqN),args.m,1-args.m)
            frqA = scipy.optimize.fminbound(functools.partial(f,frqE,z=frqN),args.m,1-args.m)
            frqN = scipy.optimize.fminbound(functools.partial(f,frqE,frqA),args.m,1-args.m)

        record.add_info('EUR_AAF',round(1e4*frqE)/1e4)
        record.add_info('AFR_AAF',round(1e4*frqA)/1e4)
        record.add_info('NAT_AAF',round(1e4*frqN)/1e4)

        # compute ancestral allele frequencies only
        if args.f:
            vcf_writer.write_record(record)
            continue

        # skip rare SNPs
        if frqE < args.t and frqA < args.t and frqN < args.t or \
            frqE > 1 - args.t and frqA > 1 - args.t and frqN > 1 - args.t:
            vcf_writer.write_record(record)
            continue


        la_ind = [la2vcf_loc[i] for i in ind if i<len(la_samples)]
        n = len(la_ind)

        # skip SNPs with insufficient genotype
        if n < args.g: vcf_writer.write_record(record); continue

        ancEE = (ctypes.c_double * n)(); ancEE[:] = tr_ancEE[la_ind]
        ancEA = (ctypes.c_double * n)(); ancEA[:] = tr_ancEA[la_ind]
        ancAA = (ctypes.c_double * n)(); ancAA[:] = tr_ancAA[la_ind]
        ancEN = (ctypes.c_double * n)(); ancEN[:] = tr_ancEN[la_ind]
        ancAN = (ctypes.c_double * n)(); ancAN[:] = tr_ancAN[la_ind]
        ancNN = (ctypes.c_double * n)(); ancNN[:] = tr_ancNN[la_ind]

        lods = [0, 0]
        inds = [-1, -1]
        c_la = (ctypes.c_int8 * n)()
        step = max(1,int(round(len(la)/args.l)))
        lkl = libla.lkl3(n,frqE,frqA,frqN,ancEE,ancEA,ancAA,ancEN,ancAN,ancNN,gl0,gl1,gl2,clib)

        # identify the two best mappings
        for i in range(1,len(la),step): # first locus is skipped
            c_la[:] = la[i][la_ind]
            lod = libla.lod3(frqE,frqA,frqN,c_la,gl0,gl1,gl2,lkl,clib)
            if lod > 0:
                if lod > lods[0]:
                    if inds[0] == -1 or chrom[i] != chrom[inds[0]]:
                        lods[1], inds[1] = lods[0], inds[0]
                    lods[0], inds[0] = lod, i
                elif lod > lods[1] and chrom[i] != chrom[inds[0]]:
                    lods[1], inds[1] = lod, i

        # refine the two best mappings
        for k in [x for x in range(2) if inds[x] != -1]:
            for i in range(max(0,inds[k]-5*step),min(inds[k]+5*step,len(la))):
                if chrom[i] == chrom[inds[k]]:
                    c_la[:] = la[i][la_ind]
                    lod = libla.lod3(frqE,frqA,frqN,c_la,gl0,gl1,gl2,lkl,clib)
                    if lod > lods[k]:
                        lods[k], inds[k] = lod, i
            record.add_info('MAP' + str(k+1),str(chrom[inds[k]]) + ':' + str(start[inds[k]]) + '-' + str(end[inds[k]]) + ',' + str(np.round(100.0*lods[k])/100.0))

        if args.s:
            record.samples = []
        vcf_writer.write_record(record)
