f = open("h1n1-ha.fasta")
g = open("rlinput", "w")

lines = f.readlines()
s = ""
for i in range(1, len(lines)):
    s += lines[i].strip()

g.write(lines[0])
g.write(s)
