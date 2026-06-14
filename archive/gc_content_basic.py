gene = "ATGCGTAACCGT"

length = len(gene)
g_count = gene.count("G")
c_count = gene.count("C")

gc_content = (g_count + c_count) / length * 100

print("sequece:", gene)
print("length:", length)
print("GC Count:", gc_content, "%")
