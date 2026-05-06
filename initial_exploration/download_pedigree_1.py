import requests

url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt"
output_path = "initial_exploration/20130606_g1k_3202_samples_ped_population.txt"

response = requests.get(url)
response.raise_for_status()  # stops here if the download fails

with open(output_path, "wb") as f:
    f.write(response.content)

print(f"Saved to: {output_path}")