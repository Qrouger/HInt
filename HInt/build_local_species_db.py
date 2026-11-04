#!/usr/bin/env python3
"""Build the local accession -> taxid database read instead of UniProt.

ColabFold ships, next to uniref30_2302, the taxonomy it uses to pair MSAs:
  <db>_seq_h (+ .index)  internal key -> header "UniRef100_<accession> ..."
  <db>_mapping           internal key -> NCBI taxid (binary MMseqs2 "TAXM")
Joining both gives accession -> taxid for the same UniRef release as the MSAs,
including accessions deleted from UniProt since then. The mapping taxid is
used; the header "TaxID=" is the fallback.

Run once (20 to 90 min, ~14 GB, plus as much temporary space next to --out):
  build_local_species_db.py --colabfold-db-dir /path/to/database \
      --out /path/to/database/uniref30_2302_species.sqlite
then point AlphaPulldown to it:
  export ALPHAPULLDOWN_LOCAL_SPECIES_DB=/path/to/database/uniref30_2302_species.sqlite

Nothing in the ColabFold database directory is modified.
"""

from __future__ import annotations

import argparse
import mmap
import os
import sqlite3
import sys
import time

import numpy as np

# Must match alphapulldown/utils/mmseqs_species_identifiers.py.
LOCAL_SPECIES_DB_TABLE = "local_species_ids"
TAXM_MAGIC = 0x0C170013


def log(message: str) -> None:
    print(time.strftime("%H:%M:%S"), message, flush=True)


def open_taxonomy_mapping(path: str) -> tuple[np.ndarray, np.ndarray]:
    """Memory-map an MMseqs2 binary taxonomy mapping as (keys, taxids)."""
    size = os.path.getsize(path)
    with open(path, "rb") as handle:
        head = handle.read(64)
    if int.from_bytes(head[:4], "little") != TAXM_MAGIC:
        raise ValueError(f"{path}: not a binary MMseqs2 taxonomy mapping")
    # The (key, taxid) pairs follow a short header; keys read 0, 1, 2, ...
    for offset in range(4, 17):
        chunk = head[offset:offset + 24]  # first three (key, taxid) pairs
        if len(chunk) < 24 or (size - offset) % 8:
            continue
        pairs = np.frombuffer(chunk, dtype="<u4")
        if list(pairs[0::2]) == [0, 1, 2]:
            break
    else:
        raise ValueError(f"{path}: cannot locate the (key, taxid) pairs")
    data = np.memmap(path, dtype="<u4", mode="r", offset=offset).reshape(-1, 2)
    keys = data[:, 0]
    sample = keys[:: max(1, len(keys) // 100000)].astype(np.int64)
    if not np.all(np.diff(sample) > 0):
        raise ValueError(f"{path}: keys are not sorted")
    return keys, data[:, 1]


def header_taxid(header: bytes) -> int:
    start = header.find(b"TaxID=")
    if start < 0:
        return 0
    start += 6
    end = start
    while end < len(header) and 48 <= header[end] <= 57:
        end += 1
    return int(header[start:end]) if end > start else 0


def build_local_species_db(
    colabfold_db_dir: str,
    out_path: str,
    *,
    db_name: str = "uniref30_2302_db",
    progress_every: int = 20_000_000,
) -> int:
    """Build the accession -> taxid SQLite database; return its row count."""
    prefix = os.path.join(colabfold_db_dir, db_name)
    for suffix in ("_seq_h", "_seq_h.index", "_mapping"):
        if not os.path.exists(prefix + suffix):
            raise FileNotFoundError(prefix + suffix)
    keys, taxids = open_taxonomy_mapping(prefix + "_mapping")
    log(f"taxonomy mapping: {len(keys):,} entries")
    direct = len(keys) > 0 and int(keys[-1]) == len(keys) - 1

    def mapped_taxid(key: int) -> int:
        if direct:
            return int(taxids[key]) if key < len(keys) else 0
        row = int(np.searchsorted(keys, key))
        return int(taxids[row]) if row < len(keys) and keys[row] == key else 0

    out_path = os.path.abspath(os.path.expanduser(out_path))
    out_dir = os.path.dirname(out_path)
    os.makedirs(out_dir, exist_ok=True)
    temp_path = out_path + ".building"
    for leftover in (temp_path, temp_path + "-journal"):
        if os.path.exists(leftover):
            os.remove(leftover)
    # Sorting for the index needs temporary space: keep it next to the output
    # rather than in /tmp, which is often small on clusters.
    os.environ.setdefault("SQLITE_TMPDIR", out_dir)
    connection = sqlite3.connect(temp_path)
    connection.execute("PRAGMA journal_mode=OFF")
    connection.execute("PRAGMA synchronous=OFF")
    connection.execute(
        f"CREATE TABLE {LOCAL_SPECIES_DB_TABLE} "
        "(accession TEXT NOT NULL, taxid INTEGER NOT NULL)"
    )
    insert = f"INSERT INTO {LOCAL_SPECIES_DB_TABLE} VALUES (?, ?)"
    rows: list[tuple[str, int]] = []
    count = scanned = 0
    started = time.time()
    with open(prefix + "_seq_h", "rb") as data_file, \
            open(prefix + "_seq_h.index", "rb") as index_file:
        data = mmap.mmap(data_file.fileno(), 0, access=mmap.ACCESS_READ)
        try:
            for line in index_file:
                scanned += 1
                key, offset, length = line.split(b"\t")
                start = int(offset)
                header = data[start:start + int(length)]
                space = header.find(b" ")
                name = header[: space if space > 0 else len(header)].rstrip(b"\x00\n")
                if name.startswith(b"UniRef100_"):
                    name = name[10:]
                taxid = mapped_taxid(int(key)) or header_taxid(header)
                if taxid:
                    rows.append((name.decode(), taxid))
                    count += 1
                if len(rows) >= 200_000:
                    connection.executemany(insert, rows)
                    rows.clear()
                if progress_every and scanned % progress_every == 0:
                    log(f"{scanned:,} headers read, {count:,} rows "
                        f"({(time.time() - started) / 60:.0f} min)")
        finally:
            data.close()
    if rows:
        connection.executemany(insert, rows)
    connection.commit()
    log(f"indexing {count:,} rows ...")
    connection.execute(
        f"CREATE INDEX {LOCAL_SPECIES_DB_TABLE}_accession "
        f"ON {LOCAL_SPECIES_DB_TABLE}(accession)"
    )
    connection.commit()
    connection.close()
    os.replace(temp_path, out_path)
    log(f"local species database written: {out_path} "
        f"({count:,} rows, {(time.time() - started) / 60:.0f} min)")
    return count


def check(out_path: str) -> bool:
    """Compare a few entries with the taxids UniProt gives."""
    expected = {"Q96UW3": "5270", "A0A2K6Q9X5": "61622"}
    connection = sqlite3.connect(f"file:{out_path}?mode=ro", uri=True)
    ok = True
    for accession, taxid in expected.items():
        row = connection.execute(
            f"SELECT taxid FROM {LOCAL_SPECIES_DB_TABLE} WHERE accession = ?",
            (accession,),
        ).fetchone()
        found = str(row[0]) if row else None
        ok &= found == taxid
        log(f"check {accession}: {found} "
            + ("OK" if found == taxid else f"unexpected (UniProt says {taxid})"))
    connection.close()
    return ok


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("\n\n", 1)[1],
    )
    parser.add_argument("--colabfold-db-dir", required=True,
                        help="folder holding uniref30_2302_db_seq_h and uniref30_2302_db_mapping")
    parser.add_argument("--out", required=True, help="SQLite file to create")
    parser.add_argument("--db-name", default="uniref30_2302_db")
    parser.add_argument("--force", action="store_true",
                        help="rebuild even if --out already exists")
    args = parser.parse_args(argv)

    out_path = os.path.abspath(os.path.expanduser(args.out))
    if os.path.exists(out_path) and not args.force:
        log(f"{out_path} already exists (use --force to rebuild)")
        return 0
    try:
        build_local_species_db(args.colabfold_db_dir, out_path, db_name=args.db_name)
    except (FileNotFoundError, ValueError) as exc:
        log(f"error: {exc}")
        return 1
    check(out_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
