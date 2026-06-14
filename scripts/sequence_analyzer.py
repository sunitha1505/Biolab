def calculate_gc(sequence):
    length = len(sequence)
    if length == 0:
        return 0
    gc = (sequence.count("G") + sequence.count("C")) / length * 100
    return round(gc, 2)

def calculate_at(sequence):
    length = len(sequence)
    if length == 0:
        return 0
    at = (sequence.count("A") + sequence.count("T")) / length * 100
    return round(at, 2)

def analyze_fasta(input_path, output_path):
    input_file = open(input_path, "r")

    sequences = []
    current_id = ""
    current_seq = ""

    for line in input_file:
        line = line.strip()
        if line == "":
            continue
        elif line.startswith(">"):
            if current_id != "":
                sequences.append({
                    "id": current_id,
                    "length": len(current_seq),
                    "gc": calculate_gc(current_seq),
                    "at": calculate_at(current_seq)
                })
            current_id = line
            current_seq = ""
        else:
            current_seq += line

    if current_id != "":
        sequences.append({
            "id": current_id,
            "length": len(current_seq),
            "gc": calculate_gc(current_seq),
            "at": calculate_at(current_seq)
        })

    input_file.close()

    highest_gc = max(sequences, key=lambda x: x["gc"])
    lowest_gc = min(sequences, key=lambda x: x["gc"])

    output_file = open(output_path, "w")
    output_file.write("=" * 60 + "\n")
    output_file.write("SEQUENCE ANALYSIS REPORT\n")
    output_file.write("=" * 60 + "\n\n")

    for seq in sequences:
        output_file.write("ID: " + seq["id"] + "\n")
        output_file.write("Length: " + str(seq["length"]) + " bp\n")
        output_file.write("GC Content: " + str(seq["gc"]) + "%\n")
        output_file.write("AT Content: " + str(seq["at"]) + "%\n")
        output_file.write("-" * 40 + "\n")

    output_file.write("\nSUMMARY\n")
    output_file.write("=" * 60 + "\n")
    output_file.write("Total sequences analyzed: " + str(len(sequences)) + "\n")
    output_file.write("Highest GC: " + highest_gc["id"] + " (" + str(highest_gc["gc"]) + "%)\n")
    output_file.write("Lowest GC:  " + lowest_gc["id"] + " (" + str(lowest_gc["gc"]) + "%)\n")

    output_file.close()
    print("Done. Report saved to " + output_path)

import sys
def main():
    if len(sys.argv) != 3:
        print("Usage: python sequence_analyzer.py <input.fasta> <output.txt>")
    else:
        analyze_fasta(sys.argv[1], sys.argv[2])

main()        