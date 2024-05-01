# Read data from the dataset.txt file
data_file = "BRCA1.txt"

# Initialize a dictionary to store the counts of pathogenic mutations per amino acid
c_terminus_start = 1755
c_terminus_end = 1856

c_terminus_pathogenic_mutation_counts = 0
rest_pathogenic_mutation_counts = 0

c_terminus_total_mutations = 0
rest_total_mutations = 0

rest_aminoacids_met = []
c_rest_aminoacids_met = []

# Read each line from the data file
with open(data_file, "r") as f:
    for line in f:
        total_mutations += 1
        _, mutation, _, mutation_type = line.strip().split("\t")
        aminoacid = str(mutation[1:-1])

        if int(aminoacid) in range(c_terminus_start,c_terminus_end):
        	c_terminus_total_mutations += 1
	        if mutation_type == "pathogenic":
	        	c_terminus_pathogenic_mutation_counts += 1
	        if aminoacid not in c_aminoacids_met:
        		c_aminoacids_met.append("aminoacid")
		else:
	    	rest_total_mutations += 1
	        if aminoacid not in rest_aminoacids_met:
	        	rest_aminoacids_met.append("aminoacid")
	        if mutation_type == "pathogenic":
	        	rest_pathogenic_mutation_counts += 1

print("Amount of aminoacids in C terminal: ", len(c_aminoacids_met))
print("Total amount of mutations in C terminal: ", c_terminus_total_mutations)
print("Total amount of pathogenic mutations: ", c_terminus_pathogenic_mutation_counts)
print("Percentage of pathogenic mutations: ", round(c_terminus_pathogenic_mutation_counts/c_terminus_total_mutations*100,2), "%")



#print("Amount of aminoacids in in protein without C terminal: ", len(rest_aminoacids_met))


