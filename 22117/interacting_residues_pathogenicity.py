# Read data from the dataset.txt file
data_file = "BRCA1.txt"

# Initialize variables
interacting_residues = set(["S1655", "G1656", "K1702", "E1698", "R1699", "K1702", "Q1811", "D1813"])
pathogenic_mutations = {}

with open(data_file, "r") as f:
    next(f)  # Skip the header line
    for line in f:
        _, mutation, pathogenicity_score, _ = line.strip().split("\t")  # extract data for each mutation
        aminoacid = mutation[:-1]

        if aminoacid in interacting_residues:
            pathogenicity_score = float(pathogenicity_score)
            if aminoacid not in pathogenic_mutations:
                pathogenic_mutations[aminoacid] = (mutation, pathogenicity_score)
            else:
                _, current_score = pathogenic_mutations[aminoacid]
                if pathogenicity_score > current_score:
                    pathogenic_mutations[aminoacid] = (mutation, pathogenicity_score)

# Sort the amino acids by highest pathogenicity score
sorted_aminoacids = sorted(pathogenic_mutations.keys(), key=lambda x: pathogenic_mutations[x][1], reverse=True)

# Print the top 3 amino acids with highest pathogenicity scores
print("Top 3 amino acids with highest pathogenicity scores:")
for aminoacid in sorted_aminoacids[:3]:
    best_mutation, score = pathogenic_mutations[aminoacid]
    print(f"{aminoacid}: Mutation: {best_mutation}, Pathogenicity Score: {score}")


