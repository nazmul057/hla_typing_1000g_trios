# Stage 1: Download reference files
from pathlib import Path
from datetime import datetime, timezone
import urllib.request
import urllib.error
import json
import logging
import requests

DATA_DIR = Path("families/trios")
DATA_DIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(DATA_DIR / "collect_trios_py.log", mode="w"),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger(__name__)

# https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt
# https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index
# https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_698_related_high_coverage.sequence.index

BASE_URL = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage"
SOURCES = {
    "pedigree":   f"{BASE_URL}/20130606_g1k_3202_samples_ped_population.txt",
    "index_2504": f"{BASE_URL}/1000G_2504_high_coverage.sequence.index",
    "index_698":  f"{BASE_URL}/1000G_698_related_high_coverage.sequence.index",
}

# Reference only — superpopulation → SpecHLA mapping.
# Not used in code; per-member mapping uses POPULATION_TO_SPECHLA below.
# SPECHLA_POP_MAP = {
#     "EAS": "Asian",
#     "SAS": "Asian",
#     "AFR": "Black",
#     "EUR": "Caucasian",
#     "AMR": "Unknown",
# }

POPULATION_TO_SPECHLA = {
    # East Asian → Asian
    "CHB": "Asian", "JPT": "Asian", "CHS": "Asian", "CDX": "Asian",
    "KHV": "Asian", "CHD": "Asian",
    # South Asian → Asian
    "GIH": "Asian", "PJL": "Asian", "BEB": "Asian", "STU": "Asian", "ITU": "Asian",
    # European → Caucasian
    "CEU": "Caucasian", "TSI": "Caucasian", "GBR": "Caucasian",
    "FIN": "Caucasian", "IBS": "Caucasian",
    # African → Black
    "YRI": "Black", "LWK": "Black", "GWD": "Black", "MSL": "Black", "ESN": "Black",
    # African-diaspora in the Americas → Black
    "ASW": "Black", "ACB": "Black",
    # Admixed American → Unknown
    "MXL": "Unknown", "PUR": "Unknown", "CLM": "Unknown", "PEL": "Unknown",
}

def _get_file_size_over_http(url: str) -> int:
    try:
        r = requests.head(url, allow_redirects=True, timeout=10)
    except requests.RequestException as e:
        raise ValueError(f'Network error reaching {url}: {e}') from e
    if r.status_code != 200:
        raise ValueError(f'HTTP {r.status_code} for {url}')
    if 'Content-Length' not in r.headers:
        raise ValueError(f'No Content-Length header for {url}')
    return int(r.headers['Content-Length'])

def ebi_to_aws_url(ebi_url: str, is_related: bool = False) -> tuple[str, int, int]:
    """
    Convert an EBI CRAM URL to its AWS S3 equivalent and verify it by
    fetching the sizes of the .cram and .cram.crai files.

    Parameters
    ----------
    ebi_url : str
        EBI URL like 'ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram'
    is_related : bool
        False for the original 2,504 unrelated samples.
        True for the additional 698 trio-completing samples.

    Returns
    -------
    (aws_url, cram_size_bytes, crai_size_bytes)

    Raises
    ------
    ValueError
        If the EBI URL doesn't match the expected layout, or if either the
        .cram or .cram.crai cannot be reached on AWS.
    """
    AWS_HOST = 'https://1000genomes.s3.us-east-1.amazonaws.com'
    BASE = '1000G_2504_high_coverage'
    SUBDIR = 'additional_698_related/' if is_related else ''

    path = ebi_url.split('://', 1)[-1]
    parts = path.strip('/').split('/')
    if len(parts) < 6 or parts[0] != 'ftp.sra.ebi.ac.uk' or parts[1:3] != ['vol1', 'run']:
        raise ValueError(f'Not a recognized EBI run-archive URL: {ebi_url!r}')

    err_id = parts[-2]
    filename = parts[-1]
    if not err_id.startswith('ERR') or not filename.endswith('.final.cram'):
        raise ValueError(f'Unexpected ERR id or filename in: {ebi_url!r}')

    aws_url = f'{AWS_HOST}/{BASE}/{SUBDIR}data/{err_id}/{filename}'

    cram_size = _get_file_size_over_http(aws_url)
    crai_size = _get_file_size_over_http(aws_url + '.crai')

    return aws_url, cram_size, crai_size

def main_run():
    paths = {}
    for name, url in SOURCES.items():
        dest = DATA_DIR / Path(url).name
        if dest.exists():
            log.info(f"[cached]   {dest.name} ({dest.stat().st_size:,} bytes)")
        else:
            log.info(f"[download] {url}")
            urllib.request.urlretrieve(url, dest)
            log.info(f"           -> {dest} ({dest.stat().st_size:,} bytes)")
        paths[name] = dest

    # ===============================================================================================
    # Stage 2: Build trios

    with open(paths["pedigree"]) as f:
        header = f.readline().strip().split()
        header_cols_stripped = [h.strip() for h in header]
        col = {name: idx for idx, name in enumerate(header_cols_stripped)}
        rows = [line.strip().split() for line in f if line.strip()]
        for r in rows:
            if len(r) != len(header):
                raise ValueError(f"Malformed row (got {len(r)} fields, expected {len(header)}): {r}")

    sample_ids = {row[col["SampleID"]].strip() for row in rows}

    trios = {}
    n_duos = 0

    for row in rows:
        fid = row[col["FamilyID"]].strip()
        sid = row[col["SampleID"]].strip()
        pid = row[col["FatherID"]].strip()
        mid = row[col["MotherID"]].strip()
        sex = row[col["Sex"]].strip()
        pop = row[col["Population"]].strip()
        spop = row[col["Superpopulation"]].strip()

        if pid == "0" and mid == "0":
            continue

        if (pid == "0") or (mid == "0") or (pid not in sample_ids) or (mid not in sample_ids):
            n_duos += 1
            continue

        key = f"{fid}_{pid}_{mid}_{sid}"
        if key in trios:
            raise ValueError(f"Duplicate trio key: {key}")

        trios[key] = {
            "family_id": fid,
            "population": pop,
            "superpopulation": spop,
            "members": {
                "child":  {"sample_id": sid, "sex": "M" if sex == "1" else "F"},
                "father": {"sample_id": pid, "sex": "M"},
                "mother": {"sample_id": mid, "sex": "F"},
            },
        }

    log.info(f"Trios: {len(trios)}")
    log.info(f"Duos:  {n_duos}")

    # ===============================================================================================
    # Building cram_index

    cram_index = {}

    samples_2504 = []
    samples_698 = []

    for index_file_path in [paths["index_2504"], paths["index_698"]]:
        
        is_related = False
        if (index_file_path == paths["index_698"]):
            is_related = True
        
        with open(index_file_path) as f:
            header = None
            for line in f:
                if line.startswith("##"):
                    continue
                if line.startswith("#"):
                    header = line.lstrip("#").strip().split("\t")
                    break
                raise ValueError(f"{index_file_path.name}: hit data line before finding header")

            if header is None:
                raise ValueError(f"{index_file_path.name}: no header line found")

            header_cols_stripped = [h.strip() for h in header]
            col = {name: idx for idx, name in enumerate(header_cols_stripped)}
            rows = [ln.rstrip("\n").split("\t") for ln in f if ln.strip()]

        log.info(f"{index_file_path.name}: {len(rows)} rows, creating cram file indexes...")

        cmplt_c = 0
        for row in rows:
            cmplt_c += 1
            sid               = row[col["SAMPLE_NAME"]].strip()
            raw_url           = row[col["ENA_FILE_PATH"]].strip()
            if not raw_url.startswith("ftp://"):
                raise ValueError(f"Expected ftp:// URL for {sid}, got: {raw_url}")
            cram_url          = raw_url.replace("ftp://", "https://", 1)
            cram_md5          = row[col["MD5SUM"]].strip()
            sample_population = row[col["POPULATION"]].strip()

            cram_aws_url, cram_size, crai_size = ebi_to_aws_url(raw_url, is_related=is_related)
            
            if cmplt_c % 25 == 0:
                print(f"  {cmplt_c}/{len(rows)} checked")

            if sid not in cram_aws_url:
                raise ValueError(f"Sample ID {sid} not found in CRAM URL: {cram_aws_url}")
            
            s_hla_pop = POPULATION_TO_SPECHLA[sample_population]

            new_entry = {
                "sample_id": sid,
                "cram_url": cram_url,
                "cram_aws_url": cram_aws_url,
                "cram_size": cram_size,
                "crai_size": crai_size,
                "cram_md5": cram_md5,
                "sample_population": sample_population,
                "sample_spechla_population": s_hla_pop
            }

            if sid in cram_index:
                raise ValueError(
                    f"Conflicting entries for {sid}:\n  existing: {cram_index[sid]}\n  new:      {new_entry}"
                )

            cram_index[sid] = new_entry

            tup = (sid, cram_aws_url, s_hla_pop)
            if is_related:
                samples_698.append(tup)
            else:
                samples_2504.append(tup)

        log.info(f"  Done. Running total: {len(cram_index)}")
    
    expected = 3202
    combined = len(samples_2504) + len(samples_698)
    if combined != len(cram_index) or combined != expected:
        raise ValueError(
            f"Sample count mismatch: 2504={len(samples_2504)}, "
            f"698={len(samples_698)}, sum={combined}, "
            f"cram_index={len(cram_index)}, expected={expected}"
        )
    
    log.info(f"Combined: {len(cram_index)} unique samples (expected 3202)")
    log.info(f"One Sample dict check: {cram_index.get('HG00443')}")

    # ===============================================================================================
    # Enrich Family data

    unique_samples = set()

    for key, trio in trios.items():
        for role, member in trio["members"].items():
            sid = member["sample_id"]
            if sid not in cram_index:
                raise ValueError(f"Sample {sid} ({role} in {key}) not found in cram_index")
            sample_info = cram_index[sid]

            member["cram_url"]     = sample_info["cram_url"]
            member["cram_aws_url"] = sample_info["cram_aws_url"]
            member["cram_size"]    = sample_info["cram_size"]
            member["crai_size"]    = sample_info["crai_size"]
            member["cram_md5"]     = sample_info["cram_md5"]
            
            member["member_spechla_population"] = sample_info["sample_spechla_population"]
            member["member_sample_population"] = sample_info["sample_population"]
            unique_samples.add(sid)

    log.info(f"Enriched {len(trios)} trios with CRAM URLs")
    log.info(f"Unique samples: {len(unique_samples)}")

    # ===============================================================================================
    # Write Output file

    output = {
        "metadata": {
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "n_trios": len(trios),
            "n_unique_samples": len(unique_samples),
            "n_duos_skipped": n_duos,
            "source_files": {
                "pedigree":   paths["pedigree"].name,
                "index_2504": paths["index_2504"].name,
                "index_698":  paths["index_698"].name,
            },
        },
        "family_trios": trios,
    }

    out_path = DATA_DIR / "families_1000G_trios.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    
    samples_2504_path = DATA_DIR / "samples_2504.tsv"
    with open(samples_2504_path, "w") as f:
        for sid, url, pop in samples_2504:
            f.write(f"{sid}\t{url}\t{pop}\n")
    log.info(f"Wrote {samples_2504_path} ({len(samples_2504)} samples)")

    samples_698_path = DATA_DIR / "samples_698.tsv"
    with open(samples_698_path, "w") as f:
        for sid, url, pop in samples_698:
            f.write(f"{sid}\t{url}\t{pop}\n")
    log.info(f"Wrote {samples_698_path} ({len(samples_698)} samples)")

    log.info(f"Wrote {out_path} ({out_path.stat().st_size:,} bytes)")
    log.info(f"  Trios:        {output['metadata']['n_trios']}")
    log.info(f"  Unique Samples:  {output['metadata']['n_unique_samples']}")
    log.info(f"  Duos skipped: {output['metadata']['n_duos_skipped']}")


if __name__ == "__main__":
    try:
        main_run()
    except Exception:
        log.exception("Pipeline failed")
        raise


