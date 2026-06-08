input_file = open("data/sequences.fasta", "r")
output_file = open("results/gc_results.txt", "w")

for line in input_file:
    line = line.strip()

    if line == "":
        continue
    elif line.startswith(">"):
        output_file.write("sequence_ID: " + line + "\n")
    else:
        length = len(line)
        gc = (line.count("G") + line.count("C")) / length * 100
        output_file.write("Length: " + str(length) + "\n")
        output_file.write("GC Content: " + str(round(gc, 2)) + "%\n")
        output_file.write("---\n")

input_file.close()
output_file.close()

print("Done. Results saved to results/gc_results.txt")