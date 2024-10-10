#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import gzip
import textwrap
from pathlib import Path
from typing import Iterator, List
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"



def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=Path,
                        default=Path("OTU.fasta"), help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file: Path, minseqlen: int) -> Iterator[str]:
    """Read a compressed fasta and extract all fasta sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :return: A generator object that provides the Fasta sequences (str).
    """
    with gzip.open(amplicon_file, "rt") as fasta:
        amplicon = ""
        first_line = True
        for line in fasta:
            if first_line:
                amplicon = ""
                first_line = False
                continue
            if line.startswith(">"):
                if len(amplicon) >= minseqlen:
                    yield amplicon
                    amplicon = ""
                else:
                    amplicon = ""
            else:
                amplicon += line.strip()
        yield amplicon



def dereplication_fulllength(amplicon_file: Path, minseqlen: int, mincount: int) -> Iterator[List]:
    """Dereplicate the set of sequence

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :param mincount: (int) Minimum amplicon count
    :return: A generator object that provides a (list)[sequences, count] of sequence 
            with a count >= mincount and a length >= minseqlen.
    """
    sequences_dict = {}
    for sequence in read_fasta(amplicon_file, minseqlen):
        if sequence in sequences_dict:
            sequences_dict[sequence] += 1
        else:
            sequences_dict[sequence] = 1
    sequences_dict = sorted(sequences_dict.items(), key=lambda item: item[1], reverse=True)
    for sequence, occurence in sequences_dict:
        if occurence >= mincount:
            yield [sequence, occurence]


def get_identity(alignment_list: List[str]) -> float:
    """Compute the identity rate between two sequences

    :param alignment_list: (list) A list of aligned sequences 
            in the format ["SE-QUENCE1", "SE-QUENCE2"]
    :return: (float) The rate of identity between the two sequences.
    """
    seq_i = alignment_list[0]
    seq_j = alignment_list[1]
    same_nucl = 0
    for i in range(len(seq_i)):
        if seq_i[i] == seq_j[i]:
            same_nucl += 1
    return (same_nucl/len(seq_i))*100


def abundance_greedy_clustering(amplicon_file: Path, minseqlen: int, mincount: int, chunk_size: int, kmer_size: int) -> List:
    """Compute an abundance greedy clustering regarding sequence count and identity.
    Identify OTU sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length.
    :param mincount: (int) Minimum amplicon count.
    :param chunk_size: (int) A fournir mais non utilise cette annee
    :param kmer_size: (int) A fournir mais non utilise cette annee
    :return: (list) A list of all the [OTU (str), count (int)] .
    """
    otu_dict = {}
    for seq_occ in dereplication_fulllength(amplicon_file, minseqlen, mincount):
        seq_i = seq_occ[0]
        occ_i = seq_occ[1]
        if len(otu_dict) == 0:
            otu_dict[seq_i] = occ_i
        else:
            new_otu = True
            for seq_j in otu_dict.keys():
                if seq_j != seq_i:
                    alignment = nw.global_align(seq_i, seq_j, gap_open=-1, gap_extend=-1,
                    matrix=str(Path(__file__).parent / "MATCH"))
                    identity = get_identity(list(alignment))
                    if identity >= 97:
                        new_otu = False
                        break
            if new_otu:
                otu_dict[seq_i] = occ_i
    list_otu = [[seq, occ] for seq, occ in otu_dict.items()]
    return list_otu



def write_OTU(OTU_list: List, output_file: Path) -> None:
    """Write the OTU sequence in fasta format.

    :param OTU_list: (list) A list of OTU sequences
    :param output_file: (Path) Path to the output file
    """
    with open(output_file, "w+", encoding="utf-8") as out:
        count = 1
        for otu in OTU_list:
            seq = otu[0]
            occ = otu[1]
            out.write(f">OTU_{count} occurrence:{occ}\n")
            out.write(textwrap.fill(seq, width=80))
            out.write("\n")
            count += 1


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    amplicon = args.amplicon_file
    output = args.output_file
    minlen = args.minseqlen
    mincount = args.mincount

    chunk=100
    kmer=8

    outus = abundance_greedy_clustering(amplicon, minlen, mincount, chunk, kmer)

    write_OTU(outus, output)


if __name__ == '__main__':
    main()
