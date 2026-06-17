import glob
import numpy as np
import pandas as pd
import random
import os

dataDir = "/datapool/life-zhanghk/niulk/my_project/18_mammalian/03_upset/output"
os.makedirs(f"{dataDir}/datasets", exist_ok=True)

#function to convert a value to int or float
def getIntorStrchr(chr):
    try:
        return int(chr)
    except:
        return str(chr)

#for generating a unique ID
generated_numbers = set()
def generate_unique_random_number():
    while True:
        # Generate a 10-digit random number
        new_number = random.randint(1000000000, 9999999999)
        if new_number not in generated_numbers:
            generated_numbers.add(new_number) # Add the number to the set
            return new_number

def joinDuplicationsfromTranlocations(dfDiag):
    '''
    This function is used to join the
    duplicated G4s which might happen due to
    tranlocations such as 4,19 in gorilla. It takes
    duplicated G4s from two IDs and combines them.
    '''
    dfDiagDup = dfDiag[dfDiag["is_duplicated"] == True]

    dfDiag.set_index("ID", inplace=True)
    extractAll = dfDiag.loc[dfDiagDup["ID"]]
    extractAll.reset_index(inplace=True)

    dupGroups = []
    for idx, group in extractAll.groupby(by="ID"):
        if len(group) > 1:
            dupGroups.append(idx)

    realExtract = dfDiag.loc[dupGroups]
    realExtract.reset_index(inplace=True)

    tobeWorkedOn = realExtract.loc[realExtract["is_duplicated"] == True].copy()
    tobeWorkedOn.sort_values(by=["SPECIES", "CHR", "START", "STRAND"], inplace=True)

    tobeWorkedOn.drop(columns=["is_duplicated"], inplace=True)
    groupDups = tobeWorkedOn.groupby(by=["SPECIES", "CHR", "START", "STRAND","LENGTH","SCORE"], as_index=False).agg(list)

    tobeMergdeddf = pd.DataFrame(columns=dfDiag.columns[:-1])

    for combination in groupDups["ID"]:
        newcombined = dfDiag.loc[combination]
        newcombined.drop_duplicates(inplace=True)
        newcombined.reset_index(inplace=True)
        newcombined.drop(columns=["is_duplicated"], inplace=True)
        pick = newcombined["ID"].unique()[0]
        newcombined["ID"] = pick
        tobeMergdeddf = pd.concat([tobeMergdeddf, newcombined], ignore_index=True)

    tobeMergdeddf = tobeMergdeddf[tobeMergdeddf.columns[-1:].tolist() + tobeMergdeddf.columns[:-1].tolist()]
    tobeMergdeddf.drop_duplicates(inplace=True)

    removedDups = dfDiag.drop(tobeWorkedOn["ID"])
    removedDups.reset_index(inplace=True)
    removedDups.drop(columns=["is_duplicated"], inplace=True)
    finalmergeddfDiag = pd.concat([removedDups, tobeMergdeddf], ignore_index=True)
    finalmergeddfDiag.drop_duplicates(subset=finalmergeddfDiag.columns[1:], inplace=True)

    return finalmergeddfDiag

print("=== Starting makeSharednUniqueDatabase ===")

# Step 1: make all the graphs from each chromosome into a set
print("1. Collecting G4s from all graph files...")
graphs = glob.glob(f"{dataDir}/datasets/hsa*/hsa*.graph")
graphValues = set()

for graph in graphs:
    try:
        df = pd.read_csv(graph, header=None, sep="\t", low_memory=False, na_values="NA", dtype=str)

        # Extract all non-null values
        for row in df.values:
            for item in row:
                if item is not None and str(item) != 'nan' and str(item) != '' and str(item) != 'NA':
                    graphValues.add(str(item))
    except Exception as e:
        print(f"Warning: Error processing file {graph}: {e}")
        continue

print(f"  Collected graphValues: {len(graphValues):,} unique G4s")

# Step 2: get the g4s which are annotated as unique but also aligned
print("2. Collecting alignedUnique nodes...")
nodefile = glob.glob(f"{dataDir}/datasets/hsa*/hsa*_alignedUnique.nds")
alignednodes = set()
for nodes in nodefile:
    df = pd.read_csv(nodes, header=None, low_memory=False, dtype=str)
    nodeList = list(df.values.flatten())
    alignednodes.update([x for x in nodeList if str(x) != 'nan'])

print(f"  Collected alignednodes: {len(alignednodes):,} nodes")

# Step 3: get the g4s which are aligned & unique but not in shared
print("3. Computing alignedUnique but not shared G4s...")
uniqueAligned = alignednodes - graphValues
print(f"  Found {len(uniqueAligned):,} alignedUnique but not shared G4s")

with open(f"{dataDir}/datasets/alignedUniquehsaG.egs", "w") as f:
    for uniq in uniqueAligned:
        f.write(f"{uniq}\n")

# Step 4: collect all the beds of all the species and make a set
print("4. Collecting bed files from all species...")
speciesSName = ["Homo_sapiens","Pan_troglodytes", "Pan_paniscus",  "Gorilla_gorilla", "Pongo_pygmaeus", "Pongo_abelii" ]
bedValues = set()

for nos, species in enumerate(speciesSName):
    print(f"  Processing species {species}...")
    paths = glob.glob(f"/datapool/life-zhanghk/niulk/my_project/18_mammalian/03_upset/GreatApeT2T-G4s-main/datasets/pG4s/{species}/chr*.pqsfinder.filtered.bed")
    speciesBED = []

    for file in paths:
        if 'Unknown' not in file:  # Exclude chrUnknown
            try:
                bed = pd.read_csv(file, header=None, sep="\t")
                speciesBED.append(bed)
            except Exception as e:
                print(f"    Warning: Cannot read {file}: {e}")

    if speciesBED:
        speciesBED = pd.concat(speciesBED, axis=0, ignore_index=True)
        speciesBED["VAL"] = f'{nos+1}|' + speciesBED[0].str[3:] + '|' + speciesBED[1].astype(str) + '|' + speciesBED[5] + '|' + speciesBED[4].astype(str) + '|' + speciesBED[3].astype(str)
        bedValues.update(list(speciesBED["VAL"].values.flatten()))

print(f"  Collected bedValues: {len(bedValues):,} unique G4s")

# Step 5: get the g4s which are not in alignment graphs or are not shared (i.e. unique)
print("5. Computing unique G4s...")
uniqueValues = bedValues - graphValues
print(f"  Found {len(uniqueValues):,} unique G4s")

with open(f"{dataDir}/datasets/uniquehsaG.egs", "w") as f:
    for uniq in uniqueValues:
        f.write(f"{uniq}\n")

print("Generated set of Unique G4s")

# Step 6: Continue with subsequent processing
speciesmap = dict({"Homo_sapiens": "hs1", "Pan_troglodytes": "pan", "Pan_paniscus": "pan",
                   "Gorilla_gorilla": "gor", "Pongo_pygmaeus": "pon", "Pongo_abelii": "pon"})

hsamap = pd.read_csv("/datapool/life-zhanghk/niulk/my_project/18_mammalian/03_upset/GreatApeT2T-G4s-main/datasets/dataForAnalysis/maps/align.hsa.map", sep="\t", header=0, index_col=0)

uniqueValues = pd.read_csv(f"{dataDir}/datasets/uniquehsaG.egs", header=None, names=["MARKER"])
uniqueValues[["SPECIES","CHR","START","STRAND","LENGTH","SCORE"]] = uniqueValues["MARKER"].str.split("|", expand=True)

print("6. Processing each chromosome...")
excludechr = [2, 5, 17]
for hsa in list(range(1, 23)) + ['X', 'Y']:

    print(f"  Processing chromosome HSA {hsa}")

    markers = []
    for nos, species in enumerate(speciesSName):
        genus = speciesmap[species]
        if hsa not in excludechr:
            chr = getIntorStrchr(hsamap[f"{hsa}"][f"{genus}"])
            locate = uniqueValues[(uniqueValues["SPECIES"] == f"{nos+1}") & (uniqueValues["CHR"] == f"{chr}")]["MARKER"]
            markers.extend(list(locate))

        elif hsa == 2:
            if genus == "hs1":
                chr = 2
                locate = uniqueValues[(uniqueValues["SPECIES"] == f"{nos+1}") & (uniqueValues["CHR"] == f"{chr}")]["MARKER"]
                markers.extend(list(locate))
            elif genus == "pan":
                for chr in [12, 13]:
                    locate = uniqueValues[(uniqueValues["SPECIES"] == f"{nos+1}") & (uniqueValues["CHR"] == f"{chr}")]["MARKER"]
                    markers.extend(list(locate))
            else:
                for chr in [11, 12]:
                    locate = uniqueValues[(uniqueValues["SPECIES"] == f"{nos+1}") & (uniqueValues["CHR"] == f"{chr}")]["MARKER"]
                    markers.extend(list(locate))

        elif hsa == 5 or hsa == 17:
            if genus != "gor":
                chr = getIntorStrchr(hsamap[f"{hsa}"][f"{genus}"])
                locate = uniqueValues[(uniqueValues["SPECIES"] == f"{nos+1}") & (uniqueValues["CHR"] == f"{chr}")]["MARKER"]
                markers.extend(list(locate))
            else:
                for chr in [4, 19]:
                    locate = uniqueValues[(uniqueValues["SPECIES"] == f"{nos+1}") & (uniqueValues["CHR"] == f"{chr}")]["MARKER"]
                    markers.extend(list(locate))

    with open(f"{dataDir}/datasets/hsa{hsa}/uniquehsa{hsa}.egs", "w") as f:
        for uniq in markers:
            f.write(f"{uniq}\n")

    # Read the shared graph for each chromosome
    graph = pd.read_csv(f"{dataDir}/datasets/hsa{hsa}/hsa{hsa}.graph", header=None, sep="\t", low_memory=False, na_values="NA", dtype=str)

    # Make the unique G4s as the same shape as the graph of the shared g4s in that chromosome
    uniquegraph = pd.read_csv(f"{dataDir}/datasets/hsa{hsa}/uniquehsa{hsa}.egs", header=None, low_memory=False, dtype=str)
    nancolumns = pd.DataFrame(np.nan, index=uniquegraph.index, columns=range(1, len(graph.columns)))
    uniquegraph = pd.concat([uniquegraph, nancolumns], axis = 1)

    # Append the unique G4s to the shared graph
    graph = pd.concat([graph, uniquegraph], axis=0, ignore_index=True)
    graph.to_csv(f"{dataDir}/datasets/hsa{hsa}/allhsa{hsa}.graph", header=False, index=False, sep="\t", na_rep="NA")

    print(f"   Generated dataframe for HSA {hsa}")

    # Make the dataframe from the graph
    allsharedset = []
    for idx in range(len(graph)):
        sharedset = pd.DataFrame()
        # Filter out NaN values
        non_nan_values = [x for x in graph.iloc[idx].dropna().tolist() if str(x) != 'nan']
        if non_nan_values:
            sharedset[["SPECIES","CHR","START","STRAND","LENGTH","SCORE"]] = pd.DataFrame(non_nan_values)[0].str.split('|', expand=True)
            sharedset[["SPECIES","START","LENGTH","SCORE"]] = sharedset[["SPECIES","START","LENGTH","SCORE"]].apply(pd.to_numeric)
            sharedset.sort_values(by=["SPECIES","CHR","START"], inplace=True)
            sharedset['ID'] = f'#{generate_unique_random_number()}'
            sharedset.reset_index(drop=True, inplace=True)
            allsharedset.append(sharedset)

    # Concatenate all the dataframes
    if allsharedset:
        allsharedset = pd.concat(allsharedset, axis = 0, ignore_index=True)

        # Reorder the columns, to bring the ID columns to the first
        cols = allsharedset.columns.tolist()
        new_cols = [cols[-1]] + cols[:-1]
        allsharedset = allsharedset[new_cols]

        # Save the dataframe
        allsharedset.to_csv(f"{dataDir}/datasets/hsa{hsa}/allhsa{hsa}.graph.df", header=True, index=False, sep="\t")
    else:
        print(f"    Warning: No data for HSA {hsa}")

print("7. Merging data from all chromosomes...")
ids = glob.glob(f"{dataDir}/datasets/hsa*/allhsa*.graph.df")

# Read all the graphs from all the chromosomes
genomesharedset = []
for graphdffile in ids:
    print(f"  Reading {graphdffile}")
    try:
        graphdf = pd.read_csv(graphdffile, header=0, sep="\t")
        genomesharedset.append(graphdf)
    except Exception as e:
        print(f"    Warning: Cannot read {graphdffile}: {e}")

# Concatenate all the dataframes to create the final database
if genomesharedset:
    genomesharedset = pd.concat(genomesharedset, axis=0, ignore_index=True)
    genomesharedset["is_duplicated"] = genomesharedset.duplicated(subset=genomesharedset.columns[1:], keep=False)

    # Join the duplicated G4s from translocations and save the final database
    genomesharedset = joinDuplicationsfromTranlocations(genomesharedset)
    genomesharedset.to_csv(f"{dataDir}/datasets/allhsaG.graph.df", header=True, index=False, sep="\t")
    print(f"Complete! Final database saved to: {dataDir}/datasets/allhsaG.graph.df")
else:
    print("Error: No data files found!")
