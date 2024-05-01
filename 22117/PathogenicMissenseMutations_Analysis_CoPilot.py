def process_mutation_data(data_file):
    c_start, c_end = 1756, 1855
    c_pathogenic_mutation_counts, c_total_mutations = 0, 0
    c_aminoacids_met, c_missense_mutations = [], {}

    rest_pathogenic_mutation_counts, rest_total_mutations = 0, 0
    rest_aminoacids_met = []

    with open(data_file, "r") as f:
        for linecounter, line in enumerate(f, start=1):
            if linecounter == 1:  # Skip the header
                continue

            _, mutation, pathogenicity_score, mutation_type = line.strip().split("\t")
            aminoacid = str(mutation[1:-1])

            if int(aminoacid) in range(c_start, c_end):
                c_total_mutations += 1
                if mutation_type == "pathogenic":
                    c_missense_mutations[mutation] = pathogenicity_score
                    c_pathogenic_mutation_counts += 1
                if aminoacid not in c_aminoacids_met:
                    c_aminoacids_met.append(aminoacid)
            else:
                rest_total_mutations += 1
                if aminoacid not in rest_aminoacids_met:
                    rest_aminoacids_met.append(aminoacid)
                if mutation_type == "pathogenic":
                    rest_pathogenic_mutation_counts += 1

    return c_missense_mutations, c_aminoacids_met, c_total_mutations, c_pathogenic_mutation_counts, rest_aminoacids_met, rest_total_mutations, rest_pathogenic_mutation_counts


data_file = "BRCA1.txt"
c_missense_mutations, c_aminoacids_met, c_total_mutations, c_pathogenic_mutation_counts, rest_aminoacids_met, rest_total_mutations, rest_pathogenic_mutation_counts = process_mutation_data(data_file)

print("\nC terminal data:")
print(f"Length of sequence: {len(c_aminoacids_met)}")
print(f"Total amount of mutations: {c_total_mutations}")
print(f"Total amount of pathogenic mutations: {c_pathogenic_mutation_counts}")
print(f"Percentage of pathogenic mutations: {round(c_pathogenic_mutation_counts / c_total_mutations * 100, 2)}%")

print("\nProtein without C terminal data:")
print(f"Length of sequence: {len(rest_aminoacids_met)}")
print(f"Total amount of mutations: {rest_total_mutations}")
print(f"Total amount of pathogenic mutations: {rest_pathogenic_mutation_counts}")
print(f"Percentage of pathogenic mutations: {round(rest_pathogenic_mutation_counts / rest_total_mutations * 100, 2)}%")
print("\n")
# Sort the dictionary by values in descending order
sorted_dict = dict(sorted(c_missense_mutations.items(), key=lambda item: float(item[1]), reverse=True))
# Print the top five key-value pairs
print("The five most pathogenic missense mutations in the C terminal are:")
for key, value in list(sorted_dict.items())[:5]:
    print(f"{key}: {value}")
