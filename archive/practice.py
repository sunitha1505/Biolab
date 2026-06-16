sequence = "ATGCGTAACCGT"

length = len(sequence)
gc_content = ((sequence.count("G") + sequence.count("C")) / length) * 100

print("sequence:", sequence)
print("length:", length)
print("gc_content:", gc_content)

if gc_content >= 50:
    print("GC content is high.")
else:
    print("GC content is low.")    
