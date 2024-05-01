# Read data from the dataset.txt file
data_file = "alphamissense_BRCA1.txt"

# cd ../../mnt/c/Users/cecil/22117/

# Initialize variables
linecounter = 0

# C termianl
c_start = 1756
c_end = 1855
c_pathogenic_mutation_counts = 0
c_total_mutations = 0
c_aminoacids_met = []
c_missense_mutations = {}

# Rest of sequence
rest_pathogenic_mutation_counts = 0
rest_total_mutations = 0
rest_aminoacids_met = []


# Read each line from the data file
with open(data_file, "r") as f:
    for line in f:
        linecounter += 1
        _, mutation, pathogenicity_score, mutation_type = line.strip().split("\t")  # extract data for each mutation
        aminoacid = str(mutation[1:-1])

        if linecounter == 1:    # pass the header
            pass

        # Extract data if the mutations are occuring in the C termianl
        elif int(aminoacid) in range(c_start, c_end):
            c_total_mutations += 1
            if mutation_type == "pathogenic":
                c_missense_mutations[mutation] = pathogenicity_score
                c_pathogenic_mutation_counts += 1
            if aminoacid not in c_aminoacids_met:
                c_aminoacids_met.append(aminoacid)

        # Extract data for rest of the sequence
        else:
            rest_total_mutations += 1
            if aminoacid not in rest_aminoacids_met:
                rest_aminoacids_met.append(aminoacid)
            if mutation_type == "pathogenic":
                rest_pathogenic_mutation_counts += 1
print("\n")
print("C terminal data: ")
print("Legnth of sequence: ", len(c_aminoacids_met))
print("Total amount of mutations: ", c_total_mutations)
print("Total amount of pathogenic mutations: ", c_pathogenic_mutation_counts)
print("Percentage of pathogenic mutations: ", round(c_pathogenic_mutation_counts / c_total_mutations * 100, 2), "%")
print("\n")
print("Protein without C terminal data: ")
print("Legnth of sequence: ", len(rest_aminoacids_met))
print("Total amount of mutations:  ", rest_total_mutations)
print("Total amount of pathogenic mutations: ", rest_pathogenic_mutation_counts)
print("Percentage of pathogenic mutations: ", round(rest_pathogenic_mutation_counts / rest_total_mutations * 100, 2), "%")
print("\n")

# Sort the dictionary by values in descending order
sorted_dict = dict(sorted(c_missense_mutations.items(), key=lambda item: float(item[1]), reverse=True))
# Print the top five key-value pairs
print("The five most pathogenic missense mutations in the C terminal are:")
for key, value in list(sorted_dict.items())[:5]:
    print(f"{key}: {value}")
print("\n")

# python3 PathogenicMissenseMutations_Analysis.py
