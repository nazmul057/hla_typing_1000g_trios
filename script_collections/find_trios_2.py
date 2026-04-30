import csv

input_path = "raw_downloads/20130606_g1k_3202_samples_ped_population.txt"
output_path = "processed_data/1kGP_trios.csv"

# Step 1: read every row into a list of dictionaries, keyed by column name
people = []
with open(input_path) as f:
    header = f.readline().split()
    for line in f:
        fields = line.split()
        if len(fields) == len(header):
            people.append(dict(zip(header, fields)))

print(f"Read {len(people)} people from the pedigree file")

# Step 2: build a lookup so we can check whether a given ID exists in the cohort
sample_ids = {p["SampleID"] for p in people}

# Step 3: find every person who has BOTH parents listed AND both parents are in the cohort
trios = []
for p in people:
    father_id = p["FatherID"]
    mother_id = p["MotherID"]
    if father_id != "0" and mother_id != "0":
        if father_id in sample_ids and mother_id in sample_ids:
            trios.append({
                "FamilyID": p["FamilyID"],
                "ChildID": p["SampleID"],
                "FatherID": father_id,
                "MotherID": mother_id,
                "Sex": "M" if p["Sex"] == "1" else "F",
                "Population": p["Population"],
                "Superpopulation": p["Superpopulation"],
            })

print(f"Found {len(trios)} complete trios")

# Step 4: write them to CSV
with open(output_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=[
        "FamilyID", "ChildID", "FatherID", "MotherID",
        "Sex", "Population", "Superpopulation"
    ])
    writer.writeheader()
    writer.writerows(trios)

print(f"Saved to: {output_path}")

# Step 5: quick summary so you can sanity-check
from collections import Counter
by_superpop = Counter(t["Superpopulation"] for t in trios)
print("\nTrios by super-population:")
for pop, count in sorted(by_superpop.items()):
    print(f"  {pop}: {count}")