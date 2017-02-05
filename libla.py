"""
   liblatools.py - library for local ancestry computations
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

import sys, gzip, ctypes, numpy as np, re, math

# retrieve the local ancestry matrix
def getla(filename,clib):
    if filename == '-':
        infile = sys.stdin
        line = infile.readline()
    else:
        infile = gzip.open(filename,'rb')
        try:
            line = infile.readline()
        except IOError:
            infile.close()
            infile = open(filename,'rt')
            line = infile.readline()
    lines = infile.readlines()
    infile.close()
    samples = bytes.decode(line).split()[3:]
    n = len(lines)
    m = len(samples)
    chrom = (ctypes.c_long * n)()
    start = (ctypes.c_long * n)()
    end = (ctypes.c_long * n)()
    la = np.zeros((n,m),dtype=ctypes.c_int8)
    lines_p = (ctypes.c_char_p * n)()
    lines_p[:] = lines

    clib.getla.restype = None;
    clib.getla(chrom,start,end,la.ctypes.data_as(ctypes.c_void_p),
               ctypes.c_long(n),ctypes.c_long(m),lines_p)
    return samples, chrom, start, end, la


# retrieve the ancestry signatures matrix
def getsig(filename,clib):
    infile = gzip.open(filename,'rb')
    try:
        lines = infile.readlines()
    except IOError:
        infile.close()
        infile = open(filename,'rt')
        lines = infile.readlines()
    infile.close()
    T = len(lines)
    n = len(bytes.decode(lines[0]).split()[3:])
    chrom = (ctypes.c_long * T)()
    pos = (ctypes.c_long * T)()
    obs = np.zeros((T,n),dtype=ctypes.c_int8)
    lines_p = (ctypes.c_char_p * T)()
    lines_p[:] = lines
    clib.getsig.restype = None;
    clib.getsig(chrom,pos,obs.ctypes.data_as(ctypes.c_void_p),
                ctypes.c_long(T),ctypes.c_long(n),lines_p)
    return n/2, chrom, pos, obs

# retrieve the genotype likelihoods
def getGL(record,loc,clib):
    if 'GT' not in record.FORMAT.split(':'):
        raise Exception("getGL: Missing genotype")
    else:
        ind = np.array([record.samples[i].gt_nums!=None for i in loc]).nonzero()[0]
    if 'GL' in record.FORMAT.split(':'):
        gl0 = np.power(10.0,np.array([record.samples[loc[i]].data.GL[0] for i in ind],dtype=ctypes.c_double))
        gl1 = np.power(10.0,np.array([record.samples[loc[i]].data.GL[1] for i in ind],dtype=ctypes.c_double))
        gl2 = np.power(10.0,np.array([record.samples[loc[i]].data.GL[2] for i in ind],dtype=ctypes.c_double))
    elif 'PL' in record.FORMAT.split(':'):
        gl0 = np.power(10.0,-np.array([record.samples[loc[i]].data.PL[0] for i in ind],dtype=ctypes.c_double)/10.0)
        gl1 = np.power(10.0,-np.array([record.samples[loc[i]].data.PL[1] for i in ind],dtype=ctypes.c_double)/10.0)
        gl2 = np.power(10.0,-np.array([record.samples[loc[i]].data.PL[2] for i in ind],dtype=ctypes.c_double)/10.0)
    else:
        ref = re.compile('0[/|]0')
        het = re.compile('0[/|][1-9]|[1-9][/|]0')
        hom = re.compile('[1-9][/|][1-9]')
        gl0 = np.array([not not ref.match(record.samples[loc[i]].gt_nums) for i in ind],dtype=ctypes.c_double)
        gl1 = np.array([not not het.match(record.samples[loc[i]].gt_nums) for i in ind],dtype=ctypes.c_double)
        gl2 = np.array([not not hom.match(record.samples[loc[i]].gt_nums) for i in ind],dtype=ctypes.c_double)
    glsum = gl0 + gl1 + gl2
    n = len(ind)
    c_gl0 = (ctypes.c_double * n)(); c_gl0[:] = gl0 / glsum
    c_gl1 = (ctypes.c_double * n)(); c_gl1[:] = gl1 / glsum
    c_gl2 = (ctypes.c_double * n)(); c_gl2[:] = gl2 / glsum
    return c_gl0, c_gl1, c_gl2, ind

# compute genotype likelihoods from ancestry proportions
def lkl3(n,frqA,frqB,frqC,ancAA,ancAB,ancBB,ancAC,ancBC,ancCC,gl0,gl1,gl2,clib):
    if frqA<0 or frqA>1 or frqB<0 or frqB>1 or frqC<0 or frqC>1:
        raise Exception("lkl3: Frequencies out of bounds")
    if len(ancAA)<n or len(ancAB)<n or len(ancBB)<n or len(ancAC)<n or len(ancBC)<n or len(ancCC)<n or len(gl0)<n or len(gl1)<n or len(gl2)<n:
        raise Exception("lkl3: Vector lengths insufficient")
    lkl = (ctypes.c_double * n)()

    clib.lkl3.restype = None;
    clib.lkl3(lkl,ctypes.c_long(n),ctypes.c_double(frqA),
              ctypes.c_double(frqB),ctypes.c_double(frqC),
              ancAA,ancAB,ancBB,ancAC,ancBC,ancCC,gl0,gl1,gl2)
    return lkl

# compute sum(log(.)) in an efficient way
def sumlog(x,clib):
    answer = ctypes.c_double()
    clib.sumlog.restype = ctypes.c_int;
    if clib.sumlog(x,ctypes.c_int(len(x)),
                   ctypes.byref(answer))!=0:
        raise Exception("sumlog: Error computing sumlog")
    return answer.value;

# compute genotype LOD score for local ancestry vs. ancestry proportions
def lod3(frqA,frqB,frqC,la,gl0,gl1,gl2,lkl,clib):
    if frqA<0 or frqA>1 or frqB<0 or frqB>1 or frqC<0 or frqC>1:
        raise Exception("lod3: Allele frequency out of bounds")
    n = len(lkl)
    if len(la)!=n or len(gl0)<n or len(gl1)<n or len(gl2)<n:
        raise Exception("lod3: Vector length insufficient")
    clib.odd3.restype = ctypes.c_double;
    odd = clib.odd3(ctypes.c_long(n),ctypes.c_double(frqA),
                    ctypes.c_double(frqB),ctypes.c_double(frqC),
                    la,gl0,gl1,gl2,lkl)
    if odd == 0:
        return float("-inf")
    else:
        return math.log10(odd)
