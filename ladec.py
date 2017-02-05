"""
   ladec.py - deconvolve local ancestry
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

import os, argparse, sys, ctypes, numpy as np
import libla

def ladec(argv):
    # setup information about chromsome starts and ends
    chromName = [ '1', '2', '3', '4', '5', '6',
              '7', '8', '9','10','11','12',
             '13','14','15','16','17','18',
             '19','20','21','22', 'X', 'Y']
    chromStart = [10000,    10000,    60000,    10000, 10000, 60000,
                  10000,    10000,    10000,    60000, 60000, 60000,
               19020000, 19000000, 20000000,    60000,     0, 10000,
                  60000,    60000,  9411193, 16050000, 10000, 10000]
    chromEnd = [249251621, 243199373, 198022430, 191154276, 180915260, 171115067,
                159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
                115169878, 107349540, 102531392,  90354753,  81195210,  78077248,
                 59128983,  63026520,  48129895,  51305566, 155270560,  59373566]
#    tel5 = {chromName[i]: chromStart[i] for i in range(len(chromName))}
#    tel3 = {chromName[i]: chromEnd[i] for i in range(len(chromName))}

    libpath = os.path.abspath(os.path.dirname(sys.argv[0]) + os.sep + 'latools.so')
    parser = argparse.ArgumentParser(description='latools.py dec: deconvolve local ancestry', add_help=False, usage='latools.py dec [options] [in.sig]')
    parser.add_argument('i', metavar='in.sig', nargs='?', type=str, default='/dev/stdin', help='ancestry signatures file (can be gzipped) [stdin]')
    parser.add_argument('-o', metavar='FILE', type=argparse.FileType('w'), default=sys.stdout, help='output BED file [stdout]')
    parser.add_argument('-a', metavar='FLOAT', type=float, default=15, help='ancestry switch cost [15]')
    parser.add_argument('-b', metavar='FLOAT', type=float, default=3, help='phase switch cost [3]')
    parser.add_argument('-c', metavar='FLOAT', type=float, default=2, help='signature error cost [2]')
    parser.add_argument('-p', metavar='FILE', type=str, default=libpath, help='library for Viterbi computations [' + libpath + ']')

    try:
        parser.error = parser.exit
        args = parser.parse_args(argv)
    except SystemExit:
        parser.print_help()
        exit(2)

    clib = ctypes.CDLL(args.p)

    if args.a < args.b:
        raise Exception('Ancestry switch cost should be greater than phase switch cost')

    # open the signatures file and stores it in a numpy array
    numpop, chrom, pos, obs = libla.getsig(args.i,clib)

    # build transition matrix
    trans = args.a / args.c * (np.kron(np.ones((numpop,numpop)),1-np.eye(numpop)) \
                            + np.kron(1-np.eye(numpop),np.ones((numpop,numpop))))
    for i in range(numpop):
        for j in range(i+1,numpop):
            trans[numpop*i+j,numpop*j+i] = min(args.b,2*args.a) / args.c
            trans[numpop*j+i,numpop*i+j] = min(args.b,2*args.a) / args.c
            for k in range(j+1,numpop):
                trans[numpop*i+j,numpop*j+k] = min(args.a + args.b,2*args.a) / args.c
                trans[numpop*j+k,numpop*i+j] = min(args.a + args.b,2*args.a) / args.c
                trans[numpop*i+j,numpop*k+i] = min(args.a + args.b,2*args.a) / args.c
                trans[numpop*k+i,numpop*i+j] = min(args.a + args.b,2*args.a) / args.c
                trans[numpop*i+k,numpop*j+i] = min(args.a + args.b,2*args.a) / args.c
                trans[numpop*j+i,numpop*i+k] = min(args.a + args.b,2*args.a) / args.c
                trans[numpop*i+k,numpop*k+j] = min(args.a + args.b,2*args.a) / args.c
                trans[numpop*k+j,numpop*i+k] = min(args.a + args.b,2*args.a) / args.c
                trans[numpop*j+i,numpop*k+j] = min(args.a + args.b,2*args.a) / args.c
                trans[numpop*k+j,numpop*j+i] = min(args.a + args.b,2*args.a) / args.c
                trans[numpop*j+k,numpop*k+i] = min(args.a + args.b,2*args.a) / args.c
                trans[numpop*k+i,numpop*j+k] = min(args.a + args.b,2*args.a) / args.c

    path2anc = [int(1+min(i,j)+max(i,j)*(max(i,j)+1)/2) for i in range(numpop) for j in range(numpop)]

    for c in range(len(chromName)):
        # identifies indices related to the chromosome of interest
        ind = [i for i in range(len(chrom)) if chrom[i]==c+1]
        if not ind:
            continue

        # run Viterbi algorithm
        T = len(ind)

        path = (ctypes.c_byte * T)()
        phased = np.array([obs[ind[i]][0] for i in range(T)],dtype=ctypes.c_int8)
        matobs = np.array([obs[ind[i]][1:numpop+1] for i in range(T)],dtype=ctypes.c_int8)
        patobs = np.array([obs[ind[i]][numpop+1:2*numpop+1] for i in range(T)],dtype=ctypes.c_int8)

        clib.getpath.restype = None;
        clib.getpath(path,ctypes.c_long(T),ctypes.c_long(numpop),
                     ctypes.c_void_p(trans.ctypes.data),
                     ctypes.c_void_p(phased.ctypes.data),
                     ctypes.c_void_p(matobs.ctypes.data),
                     ctypes.c_void_p(patobs.ctypes.data))

        anc = [path2anc[path[i]] for i in range(T)]

        # when a switch happens, the mid-point should be taken
        start = chromStart[c]
        code = anc[0]
        for t in range(1,T):
            if anc[t] != code:
                args.o.write("%s\t%d\t%d\t%d\n" % (chromName[c],start,int((pos[ind[t-1]]+pos[ind[t]])/2),code))
                start = int((pos[ind[t-1]]+pos[ind[t]])/2)
                code = anc[t]
        args.o.write("%s\t%d\t%d\t%d\n" % (chromName[c],start,chromEnd[c],code))
