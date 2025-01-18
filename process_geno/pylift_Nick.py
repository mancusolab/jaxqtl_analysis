#! /usr/bin/env python
import argparse as ap
import gzip
import logging
import os
import sys

from collections import defaultdict

import sqlite3


def get_logger(name, path=None):
  logger = logging.getLogger(name)
  if not logger.handlers:
    # Prevent logging from propagating to the root logger
    logger.propagate = 0
    console = logging.StreamHandler()
    logger.addHandler(console)

    log_format = "[%(asctime)s - %(levelname)s] %(message)s"
    date_format = "%Y-%m-%d %H:%M:%S"
    formatter = logging.Formatter(fmt=log_format, datefmt=date_format)
    console.setFormatter(formatter)

    if path is not sys.stdout and path is not None:
      disk_log_stream = open("{}.log".format(path), "w")
      disk_handler = logging.StreamHandler(disk_log_stream)
      logger.addHandler(disk_handler)
      disk_handler.setFormatter(formatter)

  return logger



def build_rsidlist(batch):
  # CHR    POS   REF   ALT   SNP   BETA  SE   P    NCASE  NCONTROL    N
  for row in batch:
    rsid = row[4]
    if rsid.startswith('rs'):
      yield rsid[2:]
    else:
      continue

  return

def process_batch(batch, dbconn, drop_miss=False):

  # push to iterator for faster/lower mem req
  rsids = ", ".join(build_rsidlist(batch))

  c = dbconn.cursor()
  query = "SELECT DISTINCT chrom, coord, rsid FROM rsid_to_coord WHERE rsid IN ({:s})".format(rsids)

  lookup = dict()
  for result in c.execute(query):
    #('22', 16875195, 25274)
    chrom, pos, rsid_int = result
    lookup[f"rs{rsid_int}"] = (chrom, str(pos))

  # CHR    POS   REF   ALT   SNP   BETA  SE   P    NCASE  NCONTROL    N
  for row in batch:
    chrom = row[0]
    pos = row[1]
    rsid = row[4]
    if rsid in lookup:
      chrom, pos = lookup[rsid]
    elif drop_miss:
      # if drop_miss is set, and we've reached this point
      # we know this variant doesn't have db lookup hit
      # and skip outputting it
      continue

    row[0] = chrom
    row[1] = pos

    yield row

  return

def main(args):
  argp = ap.ArgumentParser(description="")
  argp.add_argument("gwas")
  argp.add_argument("chrom")
  argp.add_argument("rsidx")
  argp.add_argument("--chunk-size", type=int, default=int(1e6))
  argp.add_argument("--drop-miss", action="store_true", default=False)
  argp.add_argument("-v", "--verbose", action="store_true", default=False)
  argp.add_argument("-o", "--output", type=ap.FileType("w"),
           default=sys.stdout)

  args = argp.parse_args(args)

  log = get_logger(__name__)
  if args.verbose:
    log.setLevel(logging.DEBUG)
  else:
    log.setLevel(logging.INFO)

  dbconn = sqlite3.connect(args.rsidx)

  with gzip.open(args.gwas, "rt") as gwas:
    # CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P EXTRA
    header = gwas.readline()
    args.output.write(header)

    row = gwas.readline().strip().split()
    start = int(row[1])
    batch = [row]
    for line in gwas:
      row = line.strip().split()
      pos = int(row[1])
      batch.append(row)
      if pos - start >= args.chunk_size:
        query = f"{args.chrom}:{start}-{pos}"
        log.info(f"Processing batch {query}")
        for entry in process_batch(batch, dbconn, args.drop_miss):
          args.output.write("\t".join(entry) + os.linesep)
        batch = []
        start = pos

  # any remaining
  if batch:
    query = f"{args.chrom}:{start}-{pos}"
    log.info(f"Processing batch {query}")
    for entry in process_batch(batch, dbconn, args.drop_miss):
      args.output.write("\t".join(entry) + os.linesep)

  return 0


if __name__ == "__main__":
  sys.exit(main(sys.argv[1:]))
