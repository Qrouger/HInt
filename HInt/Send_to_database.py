#!/usr/bin/env python3

import json
import zipfile
import hashlib
import io
import os
import time
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import requests
import logging

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
    if AF_version == "3" :
        AF_version = "3.0.2"
        num_recycles = 10

        RENAME = {
            "_model.cif": "_model_0.cif",
            "_summary_confidences.json": "_summary_confidences_0.json",
            "_confidences.json": "_full_data_0.json",
        }
    if AF_version == "2" :
        AF_version = "2.3.2"
        num_recycles = 3


    
    folder = Path(ppi)
    name = folder.name
    zip_path = folder / f"{name}_qsproteome_upload.zip"

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

    return zip_path

def screenshot(cif_bytes) :
    """
    Simple PNG of the Calpha trace, coloured by chain (the database requires an image)

    Parameters :
    ----------
    cif_bytes : bytes
        CIF file data

    Returns :
    -------
    bytes : PNG image data

    """
    cols, chains = [], defaultdict(list)
    for line in cif_bytes.decode().splitlines() :
        if line.startswith("_atom_site.") :
            cols.append(line.split(".")[1].strip())
        elif line.startswith("ATOM") and cols :
            f = dict(zip(cols, line.split()))
            if f["label_atom_id"] == "CA" :
                chains[f["label_asym_id"]].append((float(f["Cartn_x"]), float(f["Cartn_y"])))
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
    BASE = os.environ.get("QSPROTEOME_URL", "https://qsproteome.org")
    model_name = next(n for n in files if "model_0" in n)
    model = files[model_name]

    digest = hashlib.sha256(model).hexdigest()
    if s.post(f"{BASE}/api/rds/check-structure-hash", json={"hash": digest}).json().get("exists") :
        logger.info(f"[{run_id}] already in database, skipping upload")
        return

    job_id = f"af3ws_{int(time.time() * 1000)}_{run_id}"
    urls = {}
    for name, blob, ctype in [("screenshot.png", screenshot(model), "image/png")] + \
                             [(n, b, content_type(n)) for n, b in files.items()]:
        key = f"{job_id}/{name}"
        r = s.get(f"{BASE}/api/s3/get-presigned-url", params={"fileName": key, "fileType": ctype})
        r.raise_for_status()
        requests.put(r.json()["url"], data=blob, headers={"Content-Type": ctype}).raise_for_status()
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
    }).raise_for_status()

    logger.info(f"[{run_id}] send -> {BASE}/structure/{job_id}")


def send_to_QSproteome (zip_file) :
    """
    Upload a zip file containing AlphaFold predictions to the QSProteome database.

    Parameters :
    ----------
    zip_file : str

    """
    KEYS = ["model_0", "summary_confidences_0", "full_data_0", "job_request"]
    credits = {"firstName": "HInt", "lastName": "PPI", "affiliation": "University of Rennes"}

    runs = defaultdict(dict)  # one run per top-level folder, as on the site
    with zipfile.ZipFile(zip_file) as z :
        for name in z.namelist() :
            if not name.endswith("/") and any(k in name for k in KEYS) :
                run_id = name.split("/")[0] if "/" in name else "single"
                runs[run_id][name] = z.read(name)

    with requests.Session() as s :
        for run_id, files in runs.items() :
            if not any("model_0" in n for n in files) :
                logger.info(f"[{run_id}] no file model_0, ignored")
                continue
            try :
                upload_run(s, run_id, files, credits)
            except requests.RequestException as e :
                logger.info(f"[{run_id}] Error : {e}")


