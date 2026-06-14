input_file = open("data/tp53.fasta", "r")
output_file = open("results/tp53_parsed.txt", "w")

current_id = ""
current_seq = ""

for line in input_file:
    line = line.strip()
    
    if line == "":
        continue
    elif line.startswith(">"):
        if current_id != "":
            length = len(current_seq)
            gc = (current_seq.count("G") + current_seq.count("C")) / length * 100
            output_file.write("ID:" + current_id + "\n")
            output_file.write("Length:" + str(length) + "\n")
            output_file.write("GC Content:" + str(round(gc, 2)) + "%\n")
            output_file.write("---\n")
        current_id = line
        current_seq = ""
    else:
        current_seq += line

# Write the last sequence if it exists
if current_id != "":
    length = len(current_seq)
    gc = (current_seq.count("G") + current_seq.count("C")) / length * 100 if length > 0 else 0
    output_file.write("ID:" + current_id + "\n")
    output_file.write("Length:" + str(length) + "\n")
    output_file.write("GC Content:" + str(round(gc, 2)) + "%\n")
    output_file.write("---\n")

input_file.close()
output_file.close()

print("Done. Results saved to results/parsed_results.txt")