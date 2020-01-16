import sys

class Seeding:
    def __init__(self,str1,str2,searchlen):
        super().__init__()
        self.initseed(str1,str2,searchlen)

    def initseed(self,str1,str2,searchlen):
        self.searchlen = searchlen
        self.str1hashs = []
        for i in range(0,len(str1)+1-searchlen):
            tempsubs = str1[i:i+searchlen]
            temphash = hash(tempsubs)
            self.str1hashs.append(temphash)

        #print(self.str1hashs)
        #print('Seeding...')

        #str1veri = []
        #for i in range(0,len(str1)+1-searchlen):
        #    tempsubs = str1[i:i+searchlen]
        #    verihash = hash(tempsubs)
        #    str1veri.append(verihash == self.str1hashs[i])

        #print(str1veri)

        self.str2result = []
        for i in range(0,len(str2)+1-searchlen):
            tempsubs = str2[i:i+searchlen]
            temphash = hash(tempsubs)
            try:
                tempindex = self.str1hashs.index(temphash)
            except ValueError:
                tempindex = -1

            self.str2result.append(tempindex)

        #print(self.str2result)
        #print(self.str2result[241])
        return self.str2result

    def seedindexs(self):
        self.indexpairs = []
        for i in range(len(self.str2result)):
            if self.str2result[i] != -1:
                self.indexpairs.append([self.str2result[i],i])

        if len(self.indexpairs) == 0:
            return self.indexpairs

        #print(self.indexpairs)

        self.seedsresults = []
        self.seedsresults.append([self.indexpairs[0][0],self.indexpairs[0][1],self.searchlen])
        for i in range(1,len(self.indexpairs)):
            this = self.indexpairs[i]
            last = self.indexpairs[i-1]
            seedlen = self.searchlen
            if this[0] == last[0]+1 and this[1] == last[1]+1:
                self.seedsresults[-1][2] += 1
            else:
                self.seedsresults.append([this[0],this[1],seedlen])

        #print('Seeding Complete. Total seed number: %d'%len(self.seedsresults))
        #print(self.seedsresults)
        sorted(self.seedsresults,key = lambda item: item[2])
        if len(self.seedsresults) > 50:
            self.seedsresults = self.seedsresults[0:50]

        return self.seedsresults
