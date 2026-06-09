def calculate_gc(sequence):
    length = len(sequence)
    if length == 0:
        return 0
    gc_count = (sequence.count("G") + sequence.count("C")) / length * 100
    return round(gc_count, 2)

def parse_fasta(input_path, output_path):
    input_file = open(input_path, "r")
    output_file = open(output_path, "w")

    current_id = ""
    current_seq = ""

    for line in input_file:
        line = line.strip()
        
        if line == "":
            continue
        elif line.startswith(">"):
            if current_id != "":
                gc = calculate_gc(current_seq)
                output_file.write("ID: " + current_id + "\n")
                output_file.write("Length: " + str(len(current_seq)) + "\n")
                output_file.write("GC Content: " + str(gc) + "%\n")
                output_file.write("---\n")
            current_id = line
            current_seq = ""
        else:
            current_seq += line

    if current_id != "":
        gc = calculate_gc(current_seq)
        output_file.write("ID: " + current_id + "\n")
        output_file.write("Length: " + str(len(current_seq)) + "\n")
        output_file.write("GC Content: " + str(gc) + "%\n")
        output_file.write("---\n")

    input_file.close()
    output_file.close()
    print("Done. Results saved to " + output_path)

def main():
    parse_fasta("data/tp53.fasta", "results/tp53_analyzed.txt")

main()
          
                