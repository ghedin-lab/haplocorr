import linecache
import math
import time
import copy



def levdist(str1, str2):
    """
    requires: two strings of equal length
    insures: the number of distinct indices between the two strings
    """
    count = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            count += 1
    return count

def normMemo(memo):
    """
    requires: a memo mapping strigs to (probability, string)
    ensures: normalizes the remaining probability to 1
    """
    total = 0.0
    for i in memo:
        total += 10**memo[i][0]
    for i in memo:
        memo[i] = (memo[i][0] - math.log10(total), memo[i][1])
    return memo

class HMM:
    def __init__(self, reads, corRef, k):
        """
        requires: reads is a dna sequence read by minION
                  corRef is a reference database with reads of the same types
                      as reads in fasta format
                  k is the length of the reference k-mers
        ensures: returns a hidden markov model using the reference data and
                     the input data to form transition and emission 
                     probabilites
        """
        self.reads = reads
        self.k = k
        self.codes = {"A":"AMRWVHDN", "G":"CMSYVHBN", 
                      "C":"GRSKVDBN", "T":"TYKHDBN"}
        self.consensus = self.buildConsensus(corRef)
        self.gradient = self.buildGradient()
        print self.gradient[0:10]
        self.links = self.setLinks("preprocess")
        self.YTab = self.buildYTab()
        self.XTab = self.buildXTab()

    def buildConsensus(self, corRef):
        """
        requires: corRef is a reference database with reads of the same types
            as reads in fasta format
        ensures: returns a list of all of the reads in the correction
            reference
        """
        f = open(corRef)
        lines = f.readlines()
        consensus = []
        for i in range(1, len(lines), 2):
            read = lines[i]
            if "N" not in read:
                consensus.append(read.strip())
        return consensus

    def buildGradient(self):
        """
            ensures: returns the probability of each base occuring at a given
                         position
        """
        fired = time.clock()
        #keeping position constant
        chars = ["A", "G", "C", "T"]
        gradient = []
        for i in range(len(self.consensus[1])):
            #currently takes about 2 seconds per 100 reads
            if i%100 == 0:
                print str(i) + ": " + str(time.clock() - fired)
            total = 0
            #setting the list
            gradient.append([0, 0, 0, 0])
            insert = []
            #going through all reads at a position
            for j in range(len(self.consensus)):
                parts = 0
                newinsert = []
                for k in self.codes:
                    #adding all possible configurations given ambiguity codes
                    if self.consensus[j][i] in self.codes[k]:
                        newinsert.append(chars.index(k))
                        parts += 1
                #pairing what index to add to and how many combinations to
                #add to the total
                for k in range(len(newinsert)):
                    newinsert[k] = (newinsert[k], parts)
                insert += newinsert
                total += 1.0
            #adds the number of instances of each base at a given position
            #with respect to the number of total reads in the reference
            for index in range(len(insert)):
                gradient[i][insert[index][0]] += 1.0/insert[index][1]/total
        print "this took " + str(time.clock() - fired) + " seconds"
        return gradient

    def checkYs(self, tab, read):
        """
        requires: tab is a table with strings of length k as keys
                      read is a string of length k
        ensures: returns a list of all keys in tab that could have a transition
                     edge into read
        """
        final = []
        for i in "ACGT":
            #if a possibilities are found, they are returned. 
            #If not, nothing happens.
            try:
                tab[i + read[:-1]]
                final.append(i + read[:-1])
            except KeyError:
                pass
        return final

    def checkPossible(self, kmer):
        """
        requires: kmer is a string of length k comprised of bases and IUPAC
                      ambiguity codes
        ensures: returns all possible versions of the string with ambiguity
                     codes converted to bases
        """
        final = [""]
        for i in range(len(kmer)):
            add = ""
            for j in self.codes:
                #creating a string of possibilities
                if kmer[i] in self.codes[j]:
                    add += j
            temp = []
            #adding each of those possibilities to all current possibilities
            for j in range(len(add)):
                for l in range(len(final)):
                    temp.append(final[l] + add[j])
            final = copy.copy(temp)
        return final

    def buildYTab(self):
        """
        builds the transition edges between reference kmers
        """
        f = open("lens.txt", "a")
        k = self.k
        kmers = {}
        fired = time.clock()
        for i in range(len(self.consensus)):
            #yay timing. Currently takes ~6s per 100 reads
            if i%100 == 0:
                print "y.i = " + str(i) + " time = " + str(time.clock() - fired)
            for j in range(len(self.consensus[i]) - self.k + 1):
                #finding all possible iterations of a kmer from reference
                check = self.consensus[i][j:j + k]
                insert = self.checkPossible(check)
                for possible in insert:
                        #putting reads in table
                        kmers[possible] = []
        for i in kmers:
            #finding potential inedges
            kmers[i] = (self.checkYs(kmers, i))
            f.write(str(len(kmers[i])) + "\n")
        return kmers

    def emissionProb(self, inp, emit, start):
        """
        requires: inp is a kmer from the input read
                  emit is a string that can be emitted
                  start is the position of the first base in inp
        ensures: returns the emission probability from inp to emit
        """
        prob = 1.0
        tsprob = .000025
        tvprob = .0000012
        chars = ["A", "G", "C", "T"]
        for i in range(len(inp)):
            #probability that read was accurate (the number currently here is
            #nonsense and will need to be updated as better data comes out or
            #someone with a better biological intuition than me figures it out
            #This does not seem to be a uniform probability but treating it as
            #such here.
            if inp[i] == emit[i]:
                prob *= 90
                prob /= 100
            else:
                #takes probability of a base showing up at this position
                possible = copy.deepcopy(self.gradient[start + i])
                adjust = []
                for j in range(len(possible)):
                    #adjusting for transition and tranversion probabilites
                    #setting how much to adjust by. Does not adjust totally
                    #conserved regions
                    if possible[j] != 1:
                        ts = tsprob*possible[j]
                        tv = tvprob*possible[j]
                        possible[j] -= ts
                        possible[j] -= tv
                        adjust.append((ts, tv))
                    else:
                        adjst.append((0, 0))
                for j in range(len(possible)):
                    for k in range(len(possible)):
                        if k == len(possible) - j:
                            possible[k] += adjust[j][0]
                        elif k != j:
                            possible[k] += adjust[j][1]
                prob *= (possible[chars.index(emit[i])]/len(self.consensus))
                #normalizing probabilities based on chance of inaccurate read
                prob *= 10
                prob /= 100
        prob = math.log10(prob)
        return prob

    def buildXTab(self):
        """
        ensures: builds a table with emission probabilites from the input to
                     the reference kmers and with transition probabilities
                     between input kmers
        """
        k = self.k
        fired = time.clock()
        table = []
        for i in range(len(self.reads) - k + 1):
            if i%100 == 0:
                print "at i = " + str(i) + " the time = " + str(time.clock() - fired)
            pos = []
            #taking kmer out of read and adding it to its position
            read = self.reads[i:i+k]
            pos.append(read)
            emission = {}
            for j in self.YTab:
                #creating table of emission -> emission prob
                if levdist(read, j) <= 3:
                    eprob = self.emissionProb(read, j, i)
                    emission[j] = eprob
            #adding the self edge
            emission[read] = self.emissionProb(read, read, i)
            #putting emission table in position
            pos.append(emission)
            #adding inedges
            if i > 0:
                pos.append(self.checkYs(table[-1][1], read))
            else:
                pos.append([])
            table.append(pos)
        return table

    def unList(self, list):
        """
        requires: list is a list of length 1
        ensures: returns the value inside of the list outside of a list
        """
        return int(list[0])

    def setLinks(self, fil):
        """
        requires: fil is a preprocessed file detailing linkages
        ensures: returns positional links between bases
        """
        fired = time.clock()
        f = open(fil)
        #split of linkages to a position is along lines
        lines = f.readlines()
        chars = ["A", "G", "C", "T"]
        final = []
        for i in range(len(lines)):
            if i%100 == 0:
                print "i = " + str(i) + " and t = " + str(time.clock() - fired)
            #split between links to position i between positions <i are along
            #tabs
            pos = lines[i].split("\t")
            ch = []
            for ind in range(len(pos)):
                spot = []
                #split between linkages to a character at position i are along
                #spaces
                a = pos[ind].split()
                for j in range(len(a)):
                    #split between linkages to character at position i given a
                    #character at a position less than i are split along commas
                    b = a[j].split(",")
                    #the newline character creates an empty list at the end
                    #which causes bugs. since there are only four possibilities
                    #This harmlessly solves that problem
                    b = b[:4]
                    char = {}
                    for k in range(len(b) - 1):
                        #checking for significant linkage. This can probably
                        #be specified to control for a little bit of noise, but
                        #I don't know how to do that yet
                        if int(b[k]) > (self.gradient[ind][j] * 
                           self.gradient[i][j]) * len(self.consensus):
                            char[chars[k]] = int(b[k])
                    spot.append(char)
                ch.append(spot)
            final.append(ch[:-1])
        return final
                    
    def TProbs(self, str1, str2, i):
        """
        requires: str1 is the first i - 1 bases determined to be in the read
                  str2 is the kmer ending at position i
                  i is the position of the final base in str2
        ensures: returns the transition probability from str1 to str2
        """
        prob = 0
        chars = ["A", "G", "C", "T"]
        link = self.links[i]
        for ind in range(len(link)):
            #The dictionary structure worked several thousand times (not an
            #exaggeration) faster, which forced this try and except structure
            try:
                #In order, getting linkage, frequency of both chars, and
                #expectation to determine linkage, then if link is there,
                #adjusting probability to control for likelihood
                num = link[ind][chars.index(str2[-1])][str1[ind]]
                expfreq = (self.gradient[ind][chars.index(str1[ind])]
                          *self.gradient[i][chars.index(str2[-1])])
                if expfreq != 0:
                    exptot = len(self.consensus) * expfreq
                    if num >exptot:
                        prob += (math.log10(num) - math.log10(exptot))
            except KeyError:
                pass
        if prob < 0:
            print "oops"
        return prob

    def walk(self):
        """
        ensures: returns the highest probability traversal through the read
        """
        fired = time.clock()
        #Much easier to start the memo with a starting point
        memo = self.XTab[0][1]
        for i in memo:
            memo[i] = (memo[i], i)
        k = self.k
        for i in range(1, len(self.XTab)):
            #taking current emission table
            temp = self.XTab[i][1]
            newmemo = {}
            for j in temp:
                #getting inedges for possibility from table
                if j in self.YTab:
                    possible = [pry for pry in self.YTab[j] if pry in memo]
                else:
                    possible = [prx for prx in self.XTab[i][2] if prx in memo]
                if i == 25:
                    print possible
                #going through inedges and finding most likely path in. picks
                #first option as base and then updates when it finds a more
                #likely path
                if len(possible) > 0:
                    if i == 25:
                        print "hi"
                    check = (self.TProbs(memo[possible[0]][1], j, i + k - 1)
                            + temp[j]
                            + memo[possible[0]][0])
                    checkind = 0
                    if i == 25:
                        print "(" + str(check) + ", " + str(checkind) + ")"
                    for poss in range(1, len(possible)):
                        newprob = (self.TProbs(memo[possible[poss]][1], 
                                              j, i + k - 1)
                                  + temp[j]
                                  + memo[possible[poss]][0])
                        if newprob > check:
                            check = newprob
                            checkind = poss
                    if i == 25 :
                        print "adding"
                        print "(" + str(check) + ", " + str(checkind) + ")"
                    #adds to the new memo
                    newmemo[j] = (check,
                                  memo[possible[checkind]][1] + j[-1])
                if i == 25:
                    print newmemo
            #overwrites previous memo
            if i == 25:
                print "newmemo = " + str(newmemo)
            memo = copy.deepcopy(normMemo(newmemo))
            if i == 25:
                print "memo = " + str(memo)
            if len(memo) > 0:
                print self.XTab[i][0] in memo
                print "at i = " + str(i) + ", len(memo) = " + str(len(memo))
        return memo
