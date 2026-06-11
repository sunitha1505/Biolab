from Bio import SeqIO

for record in SeqIO.parse("data/tp53.fasta", "fasta"):
    print("ID:", record.id)
    print("Length:", len(record.seq))
    print("GC Content:", round((record.seq.count("G") + record.seq.count("C")) / len(record.seq) * 100, 2), "%")
    print("First 20 bases:", record.seq[:20])
    print("Complement:", record.seq[:20].complement())
    print("Reverse Complement:", record.seq[:20].reverse_complement())    