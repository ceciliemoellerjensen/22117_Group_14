# Read data from the dataset.txt file
data_file = "alphamissense_BRCA1.txt"

# cd ../../mnt/c/Users/cecil/22117/

# Initialize variables
linecounter = 0

# C terminal variables
c_terminus_start = 1756
c_terminus_end = 1855
c_terminus_pathogenic_mutation_counts = 0
c_terminus_total_mutations = 0
c_aminoacids_met = []
c_terminal_missense_mutations = {}

# N terminal variables
n_terminus_start = 224
n_terminus_end = 500
n_terminus_pathogenic_mutation_counts = 0
n_terminus_total_mutations = 0
n_aminoacids_met = []
n_terminal_missense_mutations = {}

# rest of amino acids 
rest_pathogenic_mutation_counts = 0
rest_total_mutations = 0
rest_aminoacids_met = []


# Read each line from the data file
with open(data_file, "r") as f:
    for line in f:
        linecounter += 1
        _, mutation, pathogenicity_score, mutation_type = line.strip().split("\t")
        aminoacid = str(mutation[1:-1])

        if linecounter == 1:    # header 
            pass

        elif int(aminoacid) in range(n_terminus_start, n_terminus_end):
            n_terminus_total_mutations += 1
            if mutation_type == "pathogenic":
                n_terminal_missense_mutations[mutation] = pathogenicity_score
                n_terminus_pathogenic_mutation_counts += 1
            if aminoacid not in n_aminoacids_met:
                n_aminoacids_met.append(aminoacid)

        elif int(aminoacid) in range(c_terminus_start, c_terminus_end):
            c_terminus_total_mutations += 1
            if mutation_type == "pathogenic":
                c_terminal_missense_mutations[mutation] = pathogenicity_score
                c_terminus_pathogenic_mutation_counts += 1
            if aminoacid not in c_aminoacids_met:
                c_aminoacids_met.append(aminoacid)
        else:
            rest_total_mutations += 1
            if aminoacid not in rest_aminoacids_met:
                rest_aminoacids_met.append(aminoacid)
            if mutation_type == "pathogenic":
                rest_pathogenic_mutation_counts += 1
print("\n")
print("C terminal data: ")
print("Amount of amino acids in C terminal: ", len(c_aminoacids_met))
print("Total amount of mutations in C terminal: ", c_terminus_total_mutations)
print("Total amount of pathogenic mutations: ", c_terminus_pathogenic_mutation_counts)
print("Percentage of pathogenic mutations: ", round(c_terminus_pathogenic_mutation_counts / c_terminus_total_mutations * 100, 2), "%")
print("\n")
print("N terminal data: ")
print("Amount of amino acids in N terminal: ", len(n_aminoacids_met))
print("Total amount of mutations in N terminal: ", n_terminus_total_mutations)
print("Total amount of pathogenic mutations: ", n_terminus_pathogenic_mutation_counts)
print("Percentage of pathogenic mutations: ", round(n_terminus_pathogenic_mutation_counts / n_terminus_total_mutations * 100, 2), "%")
print("\n")
print("Protein without C or N terminal data: ")
print("Amount of amino acids not in C or N terminal: ", len(rest_aminoacids_met))
print("Total amount of mutations not in C or N terminal: ", rest_total_mutations)
print("Total amount of pathogenic mutations: ", rest_pathogenic_mutation_counts)
print("Percentage of pathogenic mutations: ", round(rest_pathogenic_mutation_counts / rest_total_mutations * 100, 2), "%")
print("\n")

# Sort the dictionary by values in descending order
sorted_dict = dict(sorted(c_terminal_missense_mutations.items(), key=lambda item: float(item[1]), reverse=True))
# Print the top five key-value pairs
print("The five most pathogenic missense mutations in the C terminal are:")
for key, value in list(sorted_dict.items())[:5]:
    print(f"{key}: {value}")
print("\n")

# python3 PathogenicMissenseMutations_Analysis_Optimized.py
