import HMM

g = open("rlinput")
h = open("comp.txt", "w")
lines = g.readlines()
hmm = HMM.HMM(lines[1].strip(), "correctionReference", 10)
h.write(str(hmm.walk()) + "\n")
g.close()
h.close()
