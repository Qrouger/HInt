#!/usr/bin/env python3

import json
import zipfile
import hashlib
import io
import os
import time
from collections import defaultdict, Counter
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import requests
import logging
import pickle
import gzip

# Configure global logger
logging.basicConfig(
    filename="./log_file/HInt.log",  # Log file name
    level=logging.INFO,  # Log level
    format="%(asctime)s - %(levelname)s - %(message)s"  # Log format
)

# Use the logger configuration from HInt.py
logger = logging.getLogger()


def HInt_to_QSProteome (ppi, AF_version) :
    """
    Create a Zip file to send to the database.

    Parameters :
    ----------
    ppi : str
    AF_version : str

    """
    folder = Path(ppi)
    name = folder.name
    zip_path = folder / f"{name}_qsproteome_upload.zip"

    if AF_version == "3" :
        AF_version = "3.0.2"
        num_recycles = 10

        RENAME = {
            "_model.cif": "_model_0.cif",
            "_summary_confidences.json": "_summary_confidences_0.json",
            "_confidences.json": "_full_data_0.json",
        }

        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf :
            for old, new in RENAME.items() :
                zf.write(folder / f"{name}{old}", f"{name}/{name}{new}")

            data = json.loads((folder / f"{name}_data.json").read_text())
            chains = [
                {"proteinChain": {"sequence": s["protein"]["sequence"], "count": len(s["protein"]["id"]) if isinstance(s["protein"]["id"], list) else 1}}
                for s in data["sequences"] if "protein" in s
            ]

            job = [{"name": name, "modelSeeds": [str(x) for x in data["modelSeeds"]], "sequences": chains, "dialect": "alphafold3", "version": AF_version, "num_recycles": num_recycles}]
            zf.writestr(f"{name}/{name}_job_request.json", json.dumps(job))

    if AF_version == "2" :
        AF_version = "2.3.2"
        num_recycles = 3
        AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
       "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
       "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
        # AF2 has no summary/confidences/data json: build them from ranked_0.pdb and the best result_*.pkl
        best = json.loads((folder / "ranking_debug.json").read_text())["order"][0]

        pkl = folder / f"result_{best}.pkl"
        with (open(pkl, "rb") if pkl.exists() else gzip.open(f"{pkl}.gz", "rb")) as f :
            result = pickle.load(f)
        num_recycles = int(result.get("num_recycles", num_recycles))

        pdb = (folder / "ranked_0.pdb").read_text()
        atoms = [l for l in pdb.splitlines() if l.startswith("ATOM")]
        residues = list(dict.fromkeys((l[21], int(l[22:26]), l[17:20]) for l in atoms))
        seqs = {}
        for chain, _, res in residues :
            seqs[chain] = seqs.get(chain, "") + AA3.get(res, "X")

        summary = {"iptm": float(result["iptm"]) if "iptm" in result else None, "ptm": float(result["ptm"]),
                   "ranking_score": float(result["ranking_confidence"]), "num_recycles": num_recycles}
        full_data = {"atom_chain_ids": [l[21] for l in atoms], "atom_plddts": [float(l[60:66]) for l in atoms],
                     "pae": result["predicted_aligned_error"].round(2).tolist(),
                     "token_chain_ids": [c for c, _, _ in residues], "token_res_ids": [r for _, r, _ in residues]}

        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf :
            zf.writestr(f"{name}/{name}_model_0.pdb", pdb)
            zf.writestr(f"{name}/{name}_summary_confidences_0.json", json.dumps(summary))
            zf.writestr(f"{name}/{name}_full_data_0.json", json.dumps(full_data))

            chains = [{"proteinChain": {"sequence": s, "count": n}} for s, n in Counter(seqs.values()).items()]

            job = [{"name": name, "modelSeeds": [], "sequences": chains, "dialect": "alphafold2", "version": AF_version, "num_recycles": num_recycles}]
            zf.writestr(f"{name}/{name}_job_request.json", json.dumps(job))
    return zip_path


def screenshot(model_bytes) :
    """
    Simple PNG of the Calpha trace, coloured by chain (the database requires an image)

    Parameters :
    ----------
    model_bytes : bytes
        CIF or PDB file data

    Returns :
    -------
    bytes : PNG image data

    """
    cols, chains = [], defaultdict(list)
    for line in model_bytes.decode().splitlines() :
        if line.startswith("_atom_site.") :
            cols.append(line.split(".")[1].strip())
        elif line.startswith("ATOM") and cols :  # mmCIF (AF3)
            f = dict(zip(cols, line.split()))
            if f["label_atom_id"] == "CA" :
                chains[f["label_asym_id"]].append((float(f["Cartn_x"]), float(f["Cartn_y"])))
        elif line.startswith("ATOM") and line[12:16].strip() == "CA" :  # PDB (AF2)
            chains[line[21]].append((float(line[30:38]), float(line[38:46])))
    fig, ax = plt.subplots(figsize=(4, 4))
    for xy in chains.values() :
        ax.plot(*zip(*xy), lw=1.5)
    ax.set_aspect("equal")
    ax.axis("off")
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=100)
    plt.close(fig)
    return buf.getvalue()


def content_type (name) :
    """
    Determine the content type based on the file extension.

    Parameters :
    ----------
    name : str

    Returns :
    -------
    str : Content type

    """
    if name.endswith(".json"):
        return "application/json"
    if name.endswith(".cif"):
        return "application/cif"
    if name.endswith(".pdb"):
        return "chemical/x-pdb"
    return "application/octet-stream"


def upload_run (s, run_id, files, credits) :
    """
    Upload a single run to the QSProteome database.

    Parameters :
    ----------
    s : requests.Session
    run_id : str
    files : dict
    credits : dict

    """
    TIMEOUT = 40
    BASE = os.environ.get("QSPROTEOME_URL", "https://qsproteome.org")
    model_name = next(n for n in files if "model_0" in n)
    model = files[model_name]

    digest = hashlib.sha256(model).hexdigest()
    if s.post(f"{BASE}/api/rds/check-structure-hash", json={"hash": digest}, timeout=TIMEOUT).json().get("exists") :
        logger.info(f"[{run_id}] already in database, skipping upload")
        return

    job_id = f"af3ws_{int(time.time() * 1000)}_{run_id}"
    urls = {}
    for name, blob, ctype in [("screenshot.png", screenshot(model), "image/png")] + \
                             [(n, b, content_type(n)) for n, b in files.items()]:
        key = f"{job_id}/{name}"
        r = s.get(f"{BASE}/api/s3/get-presigned-url", params={"fileName": key, "fileType": ctype}, timeout=TIMEOUT)
        r.raise_for_status()
        requests.put(r.json()["url"], data=blob, headers={"Content-Type": ctype}, timeout=TIMEOUT).raise_for_status()
        urls[name] = f"https://qsproteome.s3.us-west-1.amazonaws.com/uploads/{key}"

    def url(k) :
        return next((u for n, u in urls.items() if k in n), None)

    s.post(f"{BASE}/api/rds/complex-data/upload", json={
        "complex_id": job_id,
        "structure_file": url("model_0"),
        "config_file": url("job_request"),
        "metrics_file": url("full_data_0"),
        "summary_file": url("summary_confidences_0"),
        "image_file": url("screenshot.png"),
        "credits": credits,
    }, timeout=TIMEOUT).raise_for_status()
    logger.info(f"[{run_id}] send -> {BASE}/structure/{job_id}")


def send_to_QSproteome (zip_file) :
    """
    Upload a zip file containing AlphaFold predictions to the QSProteome database, then remove it.

    Parameters :
    ----------
    zip_file : str

    """
    BASE = os.environ.get("QSPROTEOME_URL", "https://qsproteome.org")
    KEYS = ["model_0", "summary_confidences_0", "full_data_0", "job_request"]
    credits = {"firstName": "HInt", "lastName": "PPI", "affiliation": "University of Rennes", "publicationLink": ""}

    try :
        requests.get(BASE, timeout=10).raise_for_status()
    except requests.RequestException as e :
        logger.error(f"QSProteome unreachable, upload skipped ({e}). Zip kept : {zip_file}")
        return

    runs = defaultdict(dict)  # one run per top-level folder, as on the site
    with zipfile.ZipFile(zip_file) as z :
        for name in z.namelist() :
            if any(k in name for k in KEYS) :
                runs[name.split("/")[0]][name] = z.read(name)

    with requests.Session() as s :
        for run_id, files in runs.items() :
            upload_run(s, run_id, files, credits)

    Path(zip_file).unlink(missing_ok=True)

