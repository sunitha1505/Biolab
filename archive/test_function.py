def calculate_gc(sequence):
    length = len(sequence)
    gc = (sequence.count("G") + sequence.count("C")) / length * 100 if length > 0 else 0
    return round(gc, 2)

print(calculate_gc("ATGCGCGC"))
print(calculate_gc("AAAATTTT"))