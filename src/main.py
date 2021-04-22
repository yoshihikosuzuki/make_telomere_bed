from typing import Dict
import argparse
from os.path import splitext, isfile
from logzero import logger
import bits.seq as bs
import bits.util as bu


def split_fasta(fasta_fname: str, size: int = 1000000) -> str:
    out_fname = f"{splitext(fasta_fname)[0]}.split.fasta"
    if isfile(out_fname):
        logger.info(f"Output file ({out_fname}) already exists. Skip splitting contigs.")
        return out_fname
    logger.info(f"Splitting sequences in {fasta_fname} into {size} bp substrings")
    contigs = bs.load_fasta(fasta_fname)
    split_contigs = [bs.FastaRecord(name=f"{contig.name}/{i*size}_{min((i+1)*size,contig.length)}",
                                    seq=seq)
                     for contig in contigs
                     for i, seq in enumerate(bs.split_seq(contig.seq, width=size))]
    bs.save_fasta(split_contigs, out_fname)
    return out_fname


def run_trf(fasta_fname: str, trf_path: str) -> str:
    out_fname = f"{splitext(fasta_fname)[0]}.trf"
    if isfile(out_fname):
        logger.info(f"Output file ({out_fname}) already exists. Skip running TRF.")
        return out_fname
    command = f"{trf_path} {fasta_fname} 2 7 7 80 10 50 50 -h -ngs > {out_fname}"
    logger.info(f"Running TRF: $ {command}")
    logger.info(f"NOTE: If TRF takes forever, stop and run this program again with `-s` option.")
    bu.run_command(command)
    return out_fname


def parse_trf(trf_fname: str, unit_seq: str, split_contigs: bool) -> None:
    out_fname = f"{splitext(trf_fname)[0]}.telomere.bed"
    logger.info(f"Wrting results to {out_fname}")
    unit_len = len(unit_seq)
    unit_seq = unit_seq.lower()
    r = bs.EdlibRunner(mode="global", revcomp=True, cyclic=True)
    with open(out_fname, 'w') as g:
        with open(trf_fname, 'r') as f:
            for line in f:
                data = line.strip().split()
                if line.startswith('@'):
                    cname = data[0][1:]
                    if split_contigs:
                        cname, _offset = cname.split('/')
                        offset = _offset.split('_')[0]
                    else:
                        offset = 0
                else:
                    b, e, _, ncopy = data[:4]
                    useq = data[13].lower()
                    if len(useq) == unit_len and r.align(unit_seq, useq).diff == 0.:
                        b, e = str(offset + int(b)), str(offset + int(e))
                        g.write('\t'.join([cname, b, e, ncopy]) + '\n')


def main():
    args = parse_args()
    if args.split_contigs:
        args.contig_fasta = split_fasta(args.contig_fasta)
    trf_fname = run_trf(args.contig_fasta, args.trf_path)
    parse_trf(trf_fname, args.unit_sequence, args.split_contigs)
    logger.info("Finished")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "contig_fasta",
        type=str,
        help="A fasta file name.")
    parser.add_argument(
        "unit_sequence",
        type=str,
        help="Unit sqeuence like `TTAGGG`. Both capital and small are OK.")
    parser.add_argument(
        "-t",
        "--trf_path",
        type=str,
        default="trf",
        help="Path to TRF's executable. [trf]")
    parser.add_argument(
        "-s",
        "--split_contigs",
        action="store_true",
        help="Split contigs in `contig_fasta` into 1 Mbp subsequences. TRF sometimes freezes for long contigs. Use this option if TRF takes forever. [False]")
    parser.add_argument(
        "-I",
        "--ignore_exist",
        action="store_true",
        help="Ignore existing files and generate them again (except the final .bed file). [False]")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    main()
