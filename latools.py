#!/usr/bin/env python3
"""
   lamap.py - map ancestry informative markers from genotype and local ancestry inference
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

import os, argparse, sys
sys.path.append(os.path.dirname(sys.argv[0]))
import lamap, lasig, ladec, laaaf

parser = argparse.ArgumentParser(description='latools.py: local ancestry tools', add_help=False, usage='latools.py <subcommand> [options]')
subparsers = parser.add_subparsers(title='subcommands',help='valid subcommands')

parser_map = subparsers.add_parser('map', help='Admixture-map SNPs', add_help=False)
parser_map.set_defaults(func=lamap.lamap)

parser_map = subparsers.add_parser('aaf', help='Estimate ancestral allele frequencies', add_help=False)
parser_map.set_defaults(func=laaaf.laaaf)

parser_sig = subparsers.add_parser('sig', help='Extract signature alleles', add_help=False)
parser_sig.set_defaults(func=lasig.lasig)

parser_dec = subparsers.add_parser('dec', help='Deconvolve local ancestry', add_help=False)
parser_dec.set_defaults(func=ladec.ladec)

try:
    parser.error = parser.exit
    (args, argv) = parser.parse_known_args(sys.argv[1:])
except SystemExit:
    parser.print_help()
    exit(2)

args.func(argv)
