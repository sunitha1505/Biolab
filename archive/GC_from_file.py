file = open("data/sequences.fasta", "r")

for line in file:
    line = line.strip()

    if line == "":
        continue
    elif line.startswith(">"):
        print("Sequence ID:", line)
    else:
        length = len(line)
        gc = (line.count("G") + line.count("C")) / length * 100
        print("Length:", length)
        print("GC Content:", round(gc, 2), "%")
        print("---")

file.close()        