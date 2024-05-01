# Read data from the dataset.txt file
data_file = "BRCA1.txt"

# Initialize a dictionary to store the counts of pathogenic mutations per amino acid
amino_acid_counts = {}
linecounter = 0

# Read each line from the data file
with open(data_file, "r") as f:
    for line in f:
        linecounter += 1
        _, amino_acid, _, mutation_type = line.strip().split("\t")
        if mutation_type == "pathogenic":
            # Extract the amino acid position (number) from the second column
            amino_acid_position = amino_acid[1:-1]
            # Increment the count for the corresponding amino acid position
            amino_acid_counts[amino_acid_position] = amino_acid_counts.get(amino_acid_position, 0) + 1
        else:
            # Extract the amino acid position (number) from the second column
            amino_acid_position = amino_acid[1:-1]
            # Increment the count for the corresponding amino acid position
            amino_acid_counts[amino_acid_position] = amino_acid_counts.get(amino_acid_position, 0)

# Print the results for each amino acid
print("Amino Acid Position\tPathogenic Mutations")
total_pathogenic_mutations = 0
for position, count in sorted(amino_acid_counts.items()):
    print(f"{position}\t\t\t{count}")
    total_pathogenic_mutations += count

# Calculate the total number of amino acids analyzed
total_amino_acids = len(amino_acid_counts)
print("\nTotal amino acids analyzed:", total_amino_acids)
print("\n",total_pathogenic_mutations, "out of ",linecounter, "(",round(total_pathogenic_mutations/linecounter*100,2),"%) of the aminoacids in this sequence are pathogenic")
