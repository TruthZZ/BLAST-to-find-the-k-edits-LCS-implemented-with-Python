import sys
import numpy as np

class BandDP:
    def __init__(self,str1,str2,errtol,rw,dw):
        super().__init__()
        self.initBandDP(str1,str2,errtol,rw,dw)

    def initBandDP(self,str1,str2,errtol,rw,dw):
        self.str1 = str1
        self.str2 = str2
        self.errtol = errtol
        self.rw = rw
        self.dw = dw
        self.len1 = len(str1)
        self.len2 = len(str2)

        self.lenMatrix = np.zeros((self.len1,self.len2,self.errtol+1,2))
        #print(self.lenMatrix.shape)

    def Analyze(self):
        #print('dw:'+str(self.dw))
        for enum in range(0,self.errtol+1):
            for i in range(self.dw-self.errtol,self.len1):
                left = max(0,self.rw+i-self.errtol-self.dw)
                right = min(self.len2,self.rw+i+self.errtol+1-self.dw)
                for j in range(left,right):
                    [self.lenMatrix[i,j,enum,0],self.lenMatrix[i,j,enum,1]] = self.callen(i,j,enum)
                    #if enum >= 2 and self.lenMatrix[i,j,enum] != self.lenMatrix[i,j,enum-2]:
                    #    self.lenMatrix[i,j,enum] -= 1
                    #print(self.lenMatrix[i,j,enum])
                    #self.lenMatrix[i,j,enum] = 1

            #print('Complete error number: %d'%enum,end = '\r')

        #print(self.str1[241:241+125])
        #print(self.str2[243:243+125])
        
        #print(self.lenMatrix[0:12,0:12,0])
        #print(self.lenMatrix[0:12,0:12,1])
        #print(self.lenMatrix[0:12,0:12,2])
        #print(self.lenMatrix[12:24,12:24,0])
        #print(self.lenMatrix[241:251,243:253,0])
        #print(self.lenMatrix[241:251,243:253,1])
        #print(self.lenMatrix[241:251,243:253,2])

        return self.lenMatrix
        #print(self.len1)
        #print(self.len2)
        #print(self.dw)
        #print(self.rw)
        #print(self.errtol)
        #print(self.lenMatrix.shape)

    def callen(self,i,j,errnum):
        #print([i,j])
        k = 0
        length = 0
        length1 = 0
        length2 = 0
        if i>0 and j>0 and self.str1[i-1]==self.str2[j-1] and self.str1[i] == self.str2[j]:
            return [self.lenMatrix[i-1,j-1,errnum,0]-1,self.lenMatrix[i-1,j-1,errnum,1]-1]
            #return [0,0]
        while(k<1 and i+length < self.len1 and j+length < self.len2):
            #if i+length > 0 and j+length > 0 and self.str1[i+length-1] == self.str2[j+length-1] and self.str1[i+length] == self.str2[j+length]:
            #    return [self.lenMatrix[i+length-1,j+length-1,errnum,0],self.lenMatrix[i+length-1,j+length-1,errnum,1]]
            if self.str1[i+length] == self.str2[j+length]:
                length = length + 1
                length1 = length
                length2 = length
            else:
                k = k + 1
                errleft = errnum - k
                nextlength1 = 0
                nextlength2 = 0
                if k<=errnum:
                    if i+length1 < self.len1-1 and j+length2 < self.len2-1:
                        #if j == max(self.rw+i-self.errtol-self.dw,0):
                        #    nextlength = max(self.lenMatrix[i+length+1,j+length+1,errleft],self.lenMatrix[i+length,j+length+1,errleft])
                        #elif j == min(self.len2,self.rw+i+self.errtol+1-self.dw):
                        #    nextlength = max(self.lenMatrix[i+length+1,j+length+1,errleft],self.lenMatrix[i+length+1,j+length,errleft])
                        #else:
                        #    nextlength = max(self.lenMatrix[i+length+1,j+length+1,errleft],self.lenMatrix[i+length,j+length+1,errleft],self.lenMatrix[i+length+1,j+length,errleft])
                        nextlength1 = max(self.lenMatrix[i+length1+1,j+length2+1,errleft,0],self.lenMatrix[i+length1,j+length2+1,errleft,0],self.lenMatrix[i+length1+1,j+length2,errleft,0])
                        nextlength2 = max(self.lenMatrix[i+length1+1,j+length2+1,errleft,1],self.lenMatrix[i+length1,j+length2+1,errleft,1],self.lenMatrix[i+length1+1,j+length2,errleft,1])
                        if nextlength1 > self.lenMatrix[i+length1+1,j+length2+1,errleft,0] and nextlength1 > self.lenMatrix[i+length1+1,j+length2,errleft,0]:
                            nextlength1 = nextlength1 - 1
                        if nextlength2 > self.lenMatrix[i+length1+1,j+length2+1,errleft,1] and nextlength2 > self.lenMatrix[i+length1,j+length2+1,errleft,1]:
                            nextlength2 = nextlength2 - 1
                    elif i+length1 == self.len1-1 or j+length2 == self.len2-1:
                        nextlength1 = 1
                        nextlength2 = 1

                    length1 = length1 + nextlength1 + 1
                    length2 = length2 + nextlength2 + 1

        #print([length1,length2])
        return [length1,length2]

    def searchMax(self):
        maxindex1 = int(np.argmax(self.lenMatrix[:,:,self.errtol,0]))
        maxlength1 = int(np.max(self.lenMatrix[:,:,self.errtol,0]))
        maxindex2 = int(np.argmax(self.lenMatrix[:,:,self.errtol,1]))
        maxlength2 = int(np.max(self.lenMatrix[:,:,self.errtol,1]))
        rowindex = int(maxindex1//(self.len2))
        linindex = int(maxindex2%(self.len2))
        #commonindex1 = str(rowindex)+' - '+str(rowindex+maxlength-1)
        #commonindex2 = str(linindex)+' - '+str(linindex+maxlength-1)
        #commonstr1 = self.str1[rowindex:rowindex+maxlength]
        #commonstr2 = self.str2[linindex:linindex+maxlength]

        #print('Max length: %d'%maxlength)
        #print('Index in string1: '+commonindex1)
        #print('common substring1: '+commonstr1)
        #print('Index in string2: '+commonindex2)
        #print('common substring2: '+commonstr2)

        return [rowindex,linindex,maxlength1,maxlength2]
        


