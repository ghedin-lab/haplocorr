import time

def consensus(corRef):
    f = open(corRef)
    lines = f.readlines()
    final = []
    for i in range(1, len(lines), 2):
        final.append(lines[i].strip())
    return final

def createProbs(consensus):
    fired = time.clock()
    f = open("preprocess", "w")
    insert = []
    for i in range(len(consensus[0])):
        insert.append([])
        for j in range(i):
            insert[i].append([[0, 0, 0, 0], [0, 0, 0, 0],
                              [0, 0, 0, 0], [0, 0, 0, 0]])
    print len(consensus)
    for i in range(len(consensus)):
        if i%100 == 0:
            print str(i) + " time = " + str(time.clock() - fired)
        for j in range(len(consensus[i])):
            chars = ["A", "G", "C", "T"]
            if consensus[i][j] in chars:
                ind = chars.index(consensus[i][j])
            for l in range(j):
                if consensus[i][l] in chars:
                    ind2 = chars.index(consensus[i][l])
                    try:
                        insert[j][l][ind][ind2] += 1
                    except:
                        print j
                        print len(insert[j])
                        print l
                        print ind
                        print ind2
                        raise ValueError
    for i in range(len(insert)):
        for j in range(len(insert[i])):
            for l in range(len(insert[i][j])):
                for m in range(len(insert[i][j][l])):
                    f.write(str(insert[i][j][l][m]) + ",")
                f.write(" ")
            f.write("\t")
        f.write("\n")
    f.close()
    print "total = " + str(time.clock() - fired)

createProbs(consensus("correctionReference"))
