import pandas as pd

# SNVs
snvs = pd.read_csv(
    "../../data/SNVs/allSNVs_Masterfile_REACH_SSC_VariantsOnly.tsv",
    sep = '\t', header = None
)
snvs.columns = ["Location", "Mutant"]
snvs[["Chromosome", "Location"]] = (snvs["Location"]
    .str
    .split(":", expand = True))
snvs["Location"] = snvs["Location"].apply(lambda x: x.split("-")[0])
snvs["Mutant"] = snvs["Mutant"].apply(lambda x: x[0])
snvs["Location"] = snvs["Location"].astype(int)
print(snvs.head())

# GNOMAD results
gnomad = pd.read_csv(
    "../../data/SNVs/allSNVs_gnomad.txt",
    sep = '\t', header = None
)
gnomad.columns = [
    "Chromosome", "Location", "dbSNP", "Reference", "Mutant", 
    "Allele_Frequency", "Filter", "Info"
]
gnomad.insert(len(gnomad.columns), "GNOMAD", True)
print(gnomad.loc[:, ["Chromosome", "Location", "Mutant", "GNOMAD"]].head())

snvs_gnomad = pd.merge(
    snvs,
    gnomad.loc[:, ["Chromosome", "Location", "Mutant", "GNOMAD"]],
    on = ["Chromosome", "Location", "Mutant"],
    how = "left"
)

print(snvs.shape)
print(snvs_gnomad.shape)

snvs_gnomad.to_csv(
    "../../data/SNVs/allSNVs_Masterfile_REACH_SSC_VariantsOnly_gnomadAnnot.tsv",
    sep = '\t'
)
