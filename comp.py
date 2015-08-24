f = open("comp.txt")
flines = f.readlines()
count = [0, 0]
for i in range(len(flines[0])):
    if flines[0][i] != flines[1][i]:
        count[0] += 1
print count
