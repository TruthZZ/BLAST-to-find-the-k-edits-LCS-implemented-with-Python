import sys
import numpy as np
from Seed import Seeding
from BandedDP import BandDP
import threading

class MyBlast:
    def __init__(self,str1,str2,minlength,errtol):
        super().__init__()
        self.initBlast(str1,str2,minlength,errtol)

    def initBlast(self,str1,str2,minlength,errtol):
        self.str1 = str1
        self.str2 = str2
        self.minlength = minlength
        self.errtol = errtol
        
        self.searchseedlen = minlength//(errtol + 1)
        self.doublechecked = False

    def startseeding(self):
        self.seedins = Seeding(self.str1,self.str2,self.searchseedlen)
        self.seedsresults = self.seedins.seedindexs()
        while(len(self.seedsresults) == 0):
            self.doublechecked = True
            if self.searchseedlen < 20:
                self.searchseedlen = self.searchseedlen - 2
            else:
                self.searchseedlen = self.searchseedlen//2
            self.seedins = Seeding(self.str1,self.str2,self.searchseedlen)
            self.seedsresults = self.seedins.seedindexs()

        return self.seedsresults

    def printseeds(self):
        print('total seeds number: %d'%(len(self.seedsresults)))
        for i in range(0,len(self.seedsresults)):
            print('Seed %d:'%(i+1))
            print('Seed length: %d'%(self.seedsresults[i][2]))
            print('String 1 index: %d - %d'%(self.seedsresults[i][0],self.seedsresults[i][0]+self.seedsresults[i][2]))
            print('String 2 index: %d - %d'%(self.seedsresults[i][1],self.seedsresults[i][1]+self.seedsresults[i][2]))
            print('Seed: '+self.str1[self.seedsresults[i][0]:self.seedsresults[i][0]+self.seedsresults[i][2]])
            #print('Seed: '+self.str2[self.seedsresults[i][1]:self.seedsresults[i][1]+self.seedsresults[i][2]])

    def initDPrange(self,seedinfo):
        seedstart1 = seedinfo[0]
        seedstart2 = seedinfo[1]
        seedlength = seedinfo[2]

        start1 = max(0,seedstart1-self.errtol*seedlength)
        end1 = min(len(self.str1),seedstart1+(self.errtol+1)*seedlength)
        sub1 = self.str1[start1:end1]

        start2 = max(0,seedstart2-self.errtol*seedlength)
        end2 = min(len(self.str2),seedstart2+(self.errtol+1)*seedlength)
        sub2 = self.str2[start2:end2]

        if seedstart1-self.errtol*seedlength >= 0:
            lost1 = 0
        else:
            lost1 = 0 - (seedstart1-self.errtol*seedlength)

        if seedstart2-self.errtol*seedlength >= 0:
            lost2 = 0
        else:
            lost2 = 0 - (seedstart2-self.errtol*seedlength)

        basediff = lost1 - lost2
        rightward = 0
        downward = 0
        if basediff > 0:
            rightward = basediff
        elif basediff < 0:
            downward = -basediff

        return [sub1,sub2,rightward,downward,start1,start2]

    def startDP(self):
        dp_threads = []
        self.matchresults = []
        seednumber = 1
        for item in self.seedsresults:
            #print('Extending seed: %d / %d'%(seednumber,len(self.seedsresults)),end = '\r')
            #[sub1,sub2,rw,dw,sstart1,sstart2] = self.initDPrange(item)
            #DPinstance = BandDP(sub1,sub2,self.errtol,rw,dw)
            #DPinstance.Analyze()
            #thisthread = threading.Thread(target = DPinstance.Analyze)
            #[str1start,str2start,maxlen1,maxlen2] = DPinstance.searchMax()
            #tstart1 = sstart1 + str1start
            #tstart2 = sstart2 + str2start

            #self.matchresults.append([tstart1,tstart2,maxlen1,maxlen2])
            #self.singleDP(item)
            thisthread = threading.Thread(target = self.singleDP,args = (item,))
            dp_threads.append(thisthread)
            thisthread.start()
            seednumber = seednumber + 1

        for thread in dp_threads:
            thread.join()

        self.firresult = []
        self.firresult.append(self.matchresults[0])
        for i in range(len(self.matchresults)-1):
            this = self.matchresults[i]
            next = self.matchresults[i+1]
            if this[0] != next[0] and this[1] != next[1] and this[2] != next[2] and this[3] != next[3]:
                self.firresult.append(next)

        return self.firresult

    def singleDP(self,item):
        [sub1,sub2,rw,dw,sstart1,sstart2] = self.initDPrange(item)
        DPinstance = BandDP(sub1,sub2,self.errtol,rw,dw)
        DPinstance.Analyze()
        [str1start,str2start,maxlen1,maxlen2] = DPinstance.searchMax()
        tstart1 = sstart1 + str1start
        tstart2 = sstart2 + str2start

        #print('extending done','  ',tstart1,'  ',maxlen1,'  ',tstart2,'  ',maxlen2)

        self.matchresults.append([tstart1,tstart2,maxlen1,maxlen2])


    def printfirstresult(self):
        maxlennum = 0
        firmaxlen = 0
        firmaxlen1 = 0
        firmaxlen2 = 0
        i = 0
        maxindex = 0
        #print('\n'+'Total LCS number: %d'%len(self.firresult))
        for item in self.firresult:
            #print('\n'+'LCS NO. %d'%(i+1))
            tstart1 = item[0]
            tstart2 = item[1]
            maxlen1 = item[2]
            maxlen2 = item[3]
            
            if min(maxlen1,maxlen2) >= self.minlength:
                maxlennum += 1

            commonindex1 = str(tstart1)+' - '+str(tstart1+maxlen1-1)
            commonindex2 = str(tstart2)+' - '+str(tstart2+maxlen2-1)
            commonstr1 = self.str1[tstart1:tstart1+maxlen1]
            commonstr2 = self.str2[tstart2:tstart2+maxlen2]

            if min(maxlen1,maxlen2) >= firmaxlen:
                firmaxlen = max(min(maxlen1,maxlen2),firmaxlen)
                maxindex = i
                firmaxlen1 = maxlen1
                firmaxlen2 = maxlen2
                maxsubindex1 = commonindex1
                maxsubindex2 = commonindex2
                maxLCS1 = commonstr1
                maxLCS2 = commonstr2

            #print('Max length in string1: %d'%maxlen1)
            #print('Index in string1: '+commonindex1)
            #print('common substring1: '+commonstr1)
            #print('Max length in string2: %d'%maxlen2)
            #print('Index in string2: '+commonindex2)
            #print('common substring2: '+commonstr2)

            i += 1

        #print('\n'+'First Report:')
        #if maxlennum == 0:
            #print('Sorry. No LCS longer than %d'%self.minlength)
        #else:
        #    print('Total num of LCS longer than %d: %d'%(self.minlength,maxlennum))
        #print('The longest LCS length in string1 is: %d'%firmaxlen1)
        #print('Longest LCS index in string1: '+maxsubindex1)
        #print('Longest LCS in string1: '+maxLCS1)
        #print('The longest LCS length in string2 is: %d'%firmaxlen2)
        #print('Longest LCS index in string2: '+maxsubindex2)
        #print('Longest LCS in string2: '+maxLCS2)

        self.firsstmaxlen = [firmaxlen1,firmaxlen2]

        return [firmaxlen1,firmaxlen2]

    def seconditer(self):
        firmaxlen = self.firstmaxlength()
        #print(firmaxlen)
        if firmaxlen < self.minlength and self.doublechecked == False:
            #print('second check!')
            self.doublechecked = True
            self.secondblast = MyBlast(self.str1,self.str2,min(self.firsstmaxlen),self.errtol)
            self.secondblast.startseeding()
            #self.secondblast.printseeds()
            self.finalresult = self.secondblast.startDP()
            secondmaxlen = self.secondmaxlength()
            finalmaxlen = self.printfinalresult()
            self.finalmaxlength = secondmaxlen
            print('\n'+'Search Complete!')

            return self.finalresult

        else:
            print('\n'+'Search Complete!')
            self.finalmaxlength = firmaxlen
            return self.firresult

    def printfinalresult(self):
        maxlennum = 0
        finalmaxlen = 0
        finalmaxlen1 = 0
        finalmaxlen2 = 0
        i = 0
        maxindex = 0
        #print('\n'+'Total LCS number: %d'%len(self.finalresult))
        for item in self.finalresult:
            #print('\n'+'LCS NO. %d'%(i+1))
            tstart1 = item[0]
            tstart2 = item[1]
            maxlen1 = item[2]
            maxlen2 = item[3]
            
            if min(maxlen1,maxlen2) >= self.minlength:
                maxlennum += 1

            commonindex1 = str(tstart1)+' - '+str(tstart1+maxlen1-1)
            commonindex2 = str(tstart2)+' - '+str(tstart2+maxlen2-1)
            commonstr1 = self.str1[tstart1:tstart1+maxlen1]
            commonstr2 = self.str2[tstart2:tstart2+maxlen2]

            if min(maxlen1,maxlen2) >= finalmaxlen:
                finalmaxlen = max(min(maxlen1,maxlen2),finalmaxlen)
                maxindex = i
                finalmaxlen1 = maxlen1
                finalmaxlen2 = maxlen2
                maxsubindex1 = commonindex1
                maxsubindex2 = commonindex2
                maxLCS1 = commonstr1
                maxLCS2 = commonstr2

            #print('Max length in string1: %d'%maxlen1)
            #print('Index in string1: '+commonindex1)
            #print('common substring1: '+commonstr1)
            #print('Max length in string2: %d'%maxlen2)
            #print('Index in string2: '+commonindex2)
            #print('common substring2: '+commonstr2)

            i += 1

        #print('\n'+'Final Report:')
        #if maxlennum == 0:
        #    print('Sorry. No LCS longer than %d'%self.minlength)
        #else:
        #    print('Total num of LCS longer than %d: %d'%(self.minlength,maxlennum))
        #print('The longest LCS length in string1 is: %d'%finalmaxlen1)
        #print('Longest LCS index in string1: '+maxsubindex1)
        #print('Longest LCS in string1: '+maxLCS1)
        #print('The longest LCS length in string2 is: %d'%finalmaxlen2)
        #print('Longest LCS index in string2: '+maxsubindex2)
        #print('Longest LCS in string2: '+maxLCS2)

        self.finalmaxlength = [finalmaxlen1,finalmaxlen2]

        return [finalmaxlen1,finalmaxlen2]

    def firstmaxlength(self):
        firstmaxlength1 = 0
        firstmaxlength2 = 0
        for item in self.firresult:
            firstmaxlength1 = max(firstmaxlength1,item[2])
            firstmaxlength2 = max(firstmaxlength2,item[3])

        return min(firstmaxlength1,firstmaxlength2)

    def secondmaxlength(self):
        secondmaxlength1 = 0
        secondmaxlength2 = 0
        for item in self.finalresult:
            secondmaxlength1 = max(secondmaxlength1,item[2])
            secondmaxlength2 = max(secondmaxlength2,item[3])

        return min(secondmaxlength1,secondmaxlength2)




