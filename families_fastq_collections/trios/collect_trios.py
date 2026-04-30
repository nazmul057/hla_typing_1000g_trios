# Stage 1: Download reference files
from pathlib import Path
from datetime import datetime, timezone
import urllib.request
import urllib.error
import json
import logging

DATA_DIR = Path("families_fastq_collections/trios")
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

    for index_file_path in [paths["index_2504"], paths["index_698"]]:
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

        i = 0
        for row in rows:
            i += 1
            sid               = row[col["SAMPLE_NAME"]].strip()
            raw_url           = row[col["ENA_FILE_PATH"]].strip()
            if not raw_url.startswith("ftp://"):
                raise ValueError(f"Expected ftp:// URL for {sid}, got: {raw_url}")
            cram_url          = raw_url.replace("ftp://", "https://", 1)
            crai_url          = cram_url + ".crai"
            cram_md5          = row[col["MD5SUM"]].strip()
            sample_population = row[col["POPULATION"]].strip()

            
            # Verified all .crai URLs exist on 2026-04-30; re-enable for new cohorts.
            # req = urllib.request.Request(crai_url, method="HEAD")
            # try:
            #     with urllib.request.urlopen(req, timeout=30) as resp:
            #         if resp.status != 200:
            #             raise ValueError(f"CRAI not OK for {sid}: HTTP {resp.status} at {crai_url}")
            # except urllib.error.HTTPError as e:
            #     raise ValueError(f"CRAI missing for {sid}: HTTP {e.code} at {crai_url}") from e
            # except urllib.error.URLError as e:
            #     raise ValueError(f"CRAI unreachable for {sid}: {e.reason} at {crai_url}") from e
            
            # if i % 25 == 0:
            #     print(f"  {i}/{len(rows)} checked")
            

            if sid not in cram_url:
                raise ValueError(f"Sample ID {sid} not found in CRAM URL: {cram_url}")
            
            new_entry = {
                "sample_id": sid,
                "cram_url":  cram_url,
                "crai_url":  crai_url,
                "cram_md5":  cram_md5,
                "sample_population": sample_population,
            }

            if sid in cram_index:
                raise ValueError(
                    f"Conflicting entries for {sid}:\n  existing: {cram_index[sid]}\n  new:      {new_entry}"
                )

            cram_index[sid] = new_entry

        log.info(f"  Done. Running total: {len(cram_index)}")

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
            member["cram_url"] = sample_info["cram_url"]
            member["crai_url"] = sample_info["crai_url"]
            member["cram_md5"] = sample_info["cram_md5"]
            member["member_spechla_pop"] = POPULATION_TO_SPECHLA[sample_info["sample_population"]]
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


