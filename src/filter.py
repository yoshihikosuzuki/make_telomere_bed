from dataclasses import dataclass
from typing import NamedTuple, Optional, List
import argparse
import sys
from os.path import splitext, isfile
from logzero import logger
import bits.seq as bs
import bits.util as bu


class BedRecord(NamedTuple):
    chrom: str
    begin: int
    end: int
    score: float

    def to_str(self) -> str:
        return '\t'.join(list(map(str, self)))


def load_bed(in_bed_fname: str) -> List[BedRecord]:
    with open(in_bed_fname, 'r') as f:
        return [BedRecord(c, int(b), int(e), float(s))
                for c, b, e, s in list(map(lambda x: x.strip().split('\t'),
                                           f.readlines()))]


def save_bed(bed_records: List[BedRecord],
             out_bed_fname: Optional[str] = None) -> None:
    f = open(out_bed_fname, 'w') if out_bed_fname is not None else sys.stdout
    for r in bed_records:
        f.write(f"{r.to_str()}\n")
    if out_bed_fname is not None:
        f.close()


def filter_bed(bed_records: List[BedRecord],
               max_gap_len: int,
               min_score: float) -> None:
    grouped_bed_records = []
    prev_r = BedRecord("", None, None, None)
    for r in bed_records:
        if (prev_r.chrom == r.chrom
                and r.begin - prev_r.end < max_gap_len):
            grouped_bed_records[-1].append(r)
        else:
            grouped_bed_records.append([r])
        prev_r = r
    filtered_bed_records = []
    for rs in grouped_bed_records:
        tot_score = sum([r.score for r in rs])
        if tot_score >= min_score:
            filtered_bed_records += rs
    logger.info(
        f"# of records: {len(bed_records)} -> {len(filtered_bed_records)}")
    return filtered_bed_records


def main():
    args = parse_args()
    bed_records = load_bed(args.in_bed)
    filtered_bed_records = filter_bed(bed_records,
                                      args.max_gap_len,
                                      args.min_score)
    save_bed(filtered_bed_records)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "in_bed",
        type=str,
        help="Input .bed file (of telomere motifs).")
    parser.add_argument(
        "-l",
        "--max_gap_len",
        type=int,
        default=3000,
        help="Max distance between two adjacent positions to be merged. [3000]")
    parser.add_argument(
        "-m",
        "--min_score",
        type=float,
        default=100,
        help="Threshold of the score. Every adjacent positions within <`max_gap_len` bp are merged and the total score among them is used. [100]")
    args = parser.parse_args()
    assert isfile(args.in_bed), \
        f"{args.in_bed} does not exist"
    return args


if __name__ == "__main__":
    main()
