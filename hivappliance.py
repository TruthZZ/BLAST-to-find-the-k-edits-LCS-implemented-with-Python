import sys
import os
from MyBlastPy import MyBlast
import numpy  as np
import threading
import time
from matplotlib import pyplot as plt
from scipy.interpolate import spline
import gc

def align(str1,str2,minlen,errtol,longestlist,exceedlist,longestinfo,exceedinfo,threshold,i,j):
    print('Start comparing: NO.%d and NO.%d'%(i+1,j+1))
    #print(len(str1))
    #print(len(str2))
    testblast = MyBlast(str1,str2,minlen,errtol)
    testblast.startseeding()
    #testblast.printseeds()
    firstresult = testblast.startDP()
    firstmaxlen = testblast.printfirstresult()
    secondresult = testblast.seconditer()
    thisexceedinfo = []
    del testblast
    gc.collect()

    longest = 0
    for item in secondresult:
        if min(item[2],item[3]) >= threshold:
            exceedlist.append(min(item[2],item[3]))
            thisexceedinfo.append(item)

        if min(item[2],item[3]) > longest:
            longestitem = item
        longest = max(longest,min(item[2],item[3]))

    longestlist.append(longest)

    longestinfo[i,j] = longestitem
    longestinfo[j,i] = longestitem
    exceedinfo[i,j] = thisexceedinfo
    exceedinfo[j,i] = thisexceedinfo

    print('Complete comparing: NO.%d and NO.%d'%(i+1,j+1))

def alignonetomulti(str1,str2list,minlen,errtol,longestlist,exceedlist,longestinfo,exceedinfo,threshold,i):
    print('Start comparing: No.%d in dataset1'%i)
    for str2 in str2list:
        testblast = MyBlast(str1,str2,minlen,errtol)
        testblast.startseeding()
        #testblast.printseeds()
        firstresult = testblast.startDP()
        firstmaxlen = testblast.printfirstresult()
        secondresult = testblast.seconditer()

        longest = 0
        for item in secondresult:
            if min(item[2],item[3]) >= threshold:
                exceedlist.append(min(item[2],item[3]))

            longest = max(longest,min(item[2],item[3]))

        longestlist.append(longest)

    print('Complete comparing: No.%d in dataset1'%i)

def writereport(longestinfo,exceedinfo,names,totaldatanum,timeused,minimumlength,errtol):
    resultfile = r'./resultreport.txt'
    with open(resultfile,'w+') as rf:
        rf.write('Minimum length required: %d\n'%minimumlength)
        rf.write('Edits tolerence: %d\n'%errtol)
        rf.write('Time consumption: %.2f seconds\n'%timeused)
        rf.write('\n')
        for i in range(0,totaldatanum-1):
            for j in range(i+1,totaldatanum):
                longestitem = longestinfo[i,j]
                longest = min(longestitem[2],longestitem[3])
                thisexceedinfo = exceedinfo[i,j]
                
                rf.write('NO.%d and NO.%d\n'%(i+1,j+1))
                rf.write('Name '+str(i+1)+': '+names[i])
                rf.write('Name '+str(j+1)+': '+names[j])
                rf.write('Longest LCS:\n')
                rf.write('LCS length: %d\n'%longest)
                rf.write('Index in genome %d: %d - %d\n'%(i+1,longestitem[0],longestitem[0]+longestitem[2]))
                rf.write('Index in genome %d: %d - %d\n'%(j+1,longestitem[1],longestitem[1]+longestitem[3]))
                rf.write('\n')

                if len(thisexceedinfo) == 0:
                    rf.write('No LCS exceed the minimum length!\n')

                for ipie in range(0,len(thisexceedinfo)):
                    rf.write('Exceed LCS NO.%d:\n'%(ipie+1))
                    rf.write('Length in genome %d: %d\n'%(i+1,thisexceedinfo[ipie][2]))
                    rf.write('Index in genome %d: %d - %d\n'%(i+1,thisexceedinfo[ipie][0],thisexceedinfo[ipie][0]+thisexceedinfo[ipie][2]))
                    rf.write('Length in genome %d: %d\n'%(j+1,thisexceedinfo[ipie][3]))
                    rf.write('Index in genome %d: %d - %d\n'%(j+1,thisexceedinfo[ipie][1],thisexceedinfo[ipie][1]+thisexceedinfo[ipie][3]))
                    rf.write('\n')

                rf.write('\n')

def findcore(exceedinfo,totaldatanum):
    exceednums = np.zeros((1,totaldatanum))
    for i in range(0,totaldatanum):
        exceednums[0,i] = sum(len(exceedinfo[i,j]) for j in range(0,totaldatanum) if j != i)

    coreindex = np.argmax(exceednums)
    coreedgenum = np.max(exceednums)
    return [coreindex,coreedgenum]

def findmarker(exceedinfo,coreindex,coreglength,matchedgid,totaldatanum):
    markercounter = np.zeros((1,coreglength))
    for i in range(0,totaldatanum):
        if i < coreindex and i not in matchedgid:
            for item in exceedinfo[coreindex,i]:
                markerstart = item[1]
                markerend = item[1]+item[3]

                for j in range(markerstart,min(markerend,coreglength)):
                    markercounter[0,j] += 1

        elif i > coreindex and i not in matchedgid:
            for item in exceedinfo[coreindex,i]:
                markerstart = item[0]
                markerend = item[0]+item[2]

                for j in range(markerstart,min(markerend,coreglength)):
                    markercounter[0,j] += 1

    markerstart = np.argmax(markercounter)
    markertimes = np.max(markercounter)
    mlength = 0
    while markerstart+mlength < coreglength and markercounter[0,markerstart+mlength] == markertimes:
        mlength += 1
    markerend = markerstart+mlength

    return [markerstart,mlength,markertimes]

def findallmarkers(exceedinfo,coreindex,coreglength,totaldatanum):
    matchedgid = []
    allmarkers = []
    matchedgid.append(coreindex)
    added = [True]
    while len(matchedgid) < totaldatanum and True in added:
        added = []
        [markerstart,markerlength,markertimes] = findmarker(exceedinfo,coreindex,len(dataset[coreindex]),matchedgid,totaldatanum)
        #print('find one time')
        if markertimes > 0:
            allmarkers.append([markerstart,markerlength,markertimes])
            #print('append one time')
        for i in range(0,totaldatanum):
            if i < coreindex and i not in matchedgid:
                #if len(exceedinfo[coreindex,i]) == 0:
                #    print('Still Empty!')
                for j in range(0,len(exceedinfo[coreindex,i])):
                    thisitem = exceedinfo[coreindex,i][j]
                    thismarkerstart = thisitem[1]
                    thismarkerend = thisitem[1] + thisitem[3]
                    if thismarkerstart <= markerstart and thismarkerend >= markerstart+markerlength:
                        matchedgid.append(i)
                        added.append(True)
                        break
                    else:
                        added.append(False)
                        continue

            elif i > coreindex and i not in matchedgid:
                #if len(exceedinfo[coreindex,i]) == 0:
                #    print('Still Empty!')
                for j in range(0,len(exceedinfo[coreindex,i])):
                    thisitem = exceedinfo[coreindex,i][j]
                    thismarkerstart = thisitem[0]
                    thismarkerend = thisitem[0] + thisitem[2]
                    if thismarkerstart <= markerstart and thismarkerend >= markerstart+markerlength:
                        matchedgid.append(i)
                        added.append(True)
                        break
                    else:
                        added.append(False)
                        continue

            #print(matchedgid)

    return allmarkers

def compmarkers(coreindex,exceedinfo,longestinfo,totaldatanum):
    newexceedinfo = exceedinfo
    for i in range(0,totaldatanum):
        if i != coreindex:
            thisitem = newexceedinfo[coreindex,i]
            if len(thisitem) == 0:
                #print('Empty!')
                thisitem.append(longestinfo[coreindex,i])
                newexceedinfo[coreindex,i] = thisitem

    return newexceedinfo

def writecorereport(coreindex,coreedgenum,exceedinfo,corename,names,totaldatanum,allmarkers):
    corefile = r'./corereport.txt'
    with open(corefile,'w+') as cf:
        cf.write('Core genome is NO.%d\n'%coreindex)
        cf.write('Core genome name: '+corename)
        cf.write('Core genome edge number: %d\n'%coreedgenum)
        cf.write('\n')
        cf.write('Markers:\n')
        for i in range(0,len(allmarkers)):
            cf.write('Maker NO.%d\n'%(i+1))
            item = allmarkers[i]
            markerlength = item[1]
            markerstart = item[0]
            markertimes = item[2]
            cf.write('Marker length: %d\n'%markerlength)
            cf.write('Marker index: %d - %d\n'%(markerstart,markerstart+markerlength))
            cf.write('Aligned times: %d\n'%markertimes)
            cf.write('\n')

        for i in range(0,totaldatanum):
            if i < coreindex:
                cf.write('with other genomes NO.%d\n'%(i+1))
                cf.write('Name '+str(i+1)+': '+names[i])
                for j in range(0,len(exceedinfo[coreindex,i])):
                    thisitem = exceedinfo[coreindex,i][j]
                    thismarkerstart = thisitem[1]
                    thismarkerend = thisitem[1]+thisitem[3]
                    for k in range(0,len(allmarkers)):
                        thismarker = allmarkers[k]
                        markerstart = thismarker[0]
                        markerlength = thismarker[1]
                        markertimes = thismarker[2]
                        if thismarkerstart <= markerstart and thismarkerend >= markerstart+markerlength:
                            cf.write('LCS NO.%d\n'%(j+1))
                            cf.write('Aligned with marker NO.%d\n'%(k+1))
                            cf.write('Length in core genome: %d\n'%(thisitem[3]))
                            cf.write('Index in core genome: %d - %d\n'%(thismarkerstart,thismarkerend))
                            cf.write('Length in aligned genome: %d\n'%(thisitem[2]))
                            cf.write('Index in aligned genome: %d - %d\n'%(thisitem[0],thisitem[0]+thisitem[2]))
                            cf.write('\n')

            elif i > coreindex:
                cf.write('with other genomes NO.%d\n'%(i+1))
                cf.write('Name '+str(i+1)+': '+names[i])
                for j in range(0,len(exceedinfo[coreindex,i])):
                    thisitem = exceedinfo[coreindex,i][j]
                    thismarkerstart = thisitem[0]
                    thismarkerend = thisitem[0]+thisitem[2]
                    for k in range(0,len(allmarkers)):
                        thismarker = allmarkers[k]
                        markerstart = thismarker[0]
                        markerlength = thismarker[1]
                        markertimes = thismarker[2]
                        if thismarkerstart <= markerstart and thismarkerend >= markerstart+markerlength:
                            cf.write('LCS NO.%d\n'%(j+1))
                            cf.write('Aligned with marker NO.%d\n'%(k+1))
                            cf.write('Length in core genome: %d\n'%(thisitem[2]))
                            cf.write('Index in core genome: %d - %d\n'%(thismarkerstart,thismarkerend))
                            cf.write('Length in aligned genome: %d\n'%(thisitem[3]))
                            cf.write('Index in aligned genome: %d - %d\n'%(thisitem[1],thisitem[1]+thisitem[3]))
                            cf.write('\n')

            cf.write('\n')


def plotresult(longestlist,exceedlist):
    longestresult = []
    longestrefer = []
    for item in longestlist:
        if item in longestrefer:
            index = longestrefer.index(item)
            longestresult[index][1] += 1
        else:
            longestrefer.append(item)
            longestresult.append([item,1])

    longestresult = sorted(longestresult)
    print(longestresult)

    exceedresult = []
    exceedrefer = []
    for item in exceedlist:
        if item in exceedrefer:
            index = exceedrefer.index(item)
            exceedresult[index][1] += 1
        else:
            exceedrefer.append(item)
            exceedresult.append([item,1])

    exceedresult = sorted(exceedresult)
    print(exceedresult)

    longx = []
    longy = []
    excex = []
    excey = []

    for item in longestresult:
        longx.append(item[0])
        longy.append(item[1])

    for item in exceedresult:
        excex.append(item[0])
        excey.append(item[1])

    longxbar = ['0-250','251-500','501-1000','1001-2000','2001-4000','>4000']
    longybar = [0,0,0,0,0,0]
    for i in range(0,len(longx)):
        if longx[i] < 251:
            longybar[0] += longy[i]
        elif longx[i] < 501:
            longybar[1] += longy[i]
        elif longx[i] < 1001:
            longybar[2] += longy[i]
        elif longx[i] < 2001:
            longybar[3] += longy[i]
        elif longx[i] < 4001:
            longybar[4] += longy[i]
        else:
            longybar[5] += longy[i]

    excexbar = ['250-500','501-750','751-1000','1001-1500','1501-2000','2001-3000','3001-4000','>4000']
    exceybar = [0,0,0,0,0,0,0,0]
    for i in range(0,len(excex)):
        if excex[i] < 501:
            exceybar[0] += excey[i]
        elif excex[i] < 751:
            exceybar[1] += excey[i]
        elif excex[i] < 1001:
            exceybar[2] += excey[i]
        elif excex[i] < 1501:
            exceybar[3] += excey[i]
        elif excex[i] < 2001:
            exceybar[4] += excey[i]
        elif excex[i] < 3001:
            exceybar[5] += excey[i]
        elif excex[i] < 4001:
            exceybar[6] += excey[i]
        else:
            exceybar[7] += excey[i]

    plt.figure(1)
    plt.bar(range(len(longybar)),longybar)
    plt.xticks(range(len(longybar)),longxbar)
    for x,y in enumerate(longybar):
        plt.text(x,y,'%s'%y,ha = 'center',va = 'bottom')
    plt.title('Longest Length Plot')
    plt.xlabel('Length')
    plt.ylabel('Number')
    plt.savefig('./longestplot_bar.jpg')

    plt.figure(2)
    plt.bar(range(len(exceybar)),exceybar)
    plt.xticks(range(len(exceybar)),excexbar)
    for x,y in enumerate(exceybar):
        plt.text(x,y,'%s'%y,ha = 'center',va = 'bottom')
    plt.title('Exceeded Length Plot')
    plt.xlabel('Length')
    plt.ylabel('Number')
    plt.savefig('./exceededplot_bar.jpg')


    #print('Pass 1!!!')
    #print(longx)

    #longxarray = np.array(longx)
    ##print(longxarray)
    #longyarray = np.array(longy)
    #longxnew = np.linspace(longxarray.min(),longxarray.max(),50)
    ##print(longxnew)
    #longynew = spline(longxarray,longyarray,longxnew)
    ##print(longynew)

    #excexarray = np.array(excex)
    #exceyarray = np.array(excey)
    #excexnew = np.linspace(excexarray.min(),excexarray.max(),50)
    #exceynew = spline(excexarray,exceyarray,excexnew)

    #print('Pass 2!!!')

    #plt.figure(1)
    #plt.plot(longxnew,longynew)
    #plt.title('Longest Length Plot')
    #plt.xlabel('Length')
    #plt.ylabel('Number')
    #plt.savefig('./longestplot.jpg')
    ##plt.show()

    #plt.figure(2)
    #plt.plot(excexnew,exceynew)
    #plt.title('Exceeded Length Plot')
    #plt.xlabel('Length')
    #plt.ylabel('Number')
    #plt.savefig('./exceededplot.jpg')

    plt.figure(3)
    plt.plot(longx,longy)
    plt.title('Longest Length Plot')
    plt.xlabel('Length')
    plt.ylabel('Number')
    plt.savefig('./longestplot_ori.jpg')

    plt.figure(4)
    plt.plot(excex,excey)
    plt.title('Exceeded Length Plot')
    plt.xlabel('Length')
    plt.ylabel('Number')
    plt.savefig('./exceededplot_ori.jpg')
    plt.show()


if __name__ == '__main__':
    path = './hivgenomes'
    
    files = os.listdir(path)
    names = []
    dataset = []
    gnum = 0

    for file in files:
        f = open(path+'/'+file)
        ls = []
        oristr = ''
        for line in f:
            if not line.startswith('>'):
                ls.append(line.replace('\n',''))
            else:
                names.append(line)
                gnum += 1
                if len(ls) != 0:
                    for subs in ls:
                        oristr += subs

                    dataset.append(oristr)

                    ls = []
                    oristr = ''

        for subs in ls:
            oristr += subs
        dataset.append(oristr)

        f.close()

    print('Genomes number in dataset: ',len(dataset),'  ',gnum,'  ',len(names))
    #print(dataset1[0])
    #print('')
    #print(dataset2[0])

    testnumber = int(input('Please enter the number of genomes you want to test: '))
    totaldatanum = testnumber
    minlength = int(input('Please enter the minimum length you wish to find: '))
    errtol = int(input('Please enter the edits number you wish to include: '))
    plotthreshold = int(input('Please enter the plot threshold: '))
    print('Please choose the paralllel method: (enter 1 or 2)')
    parallelchoice = int(input('1. one to one compare  2. one to multiple compare: '))
    #print(len(dataset3))

    #print(dataset1)
    #print(len(dataset1[0]))
    #print(len(dataset2[0]))

    longestlist = []
    exceedlist = []
    longestinfo = np.ndarray((totaldatanum,totaldatanum),dtype = list)
    exceedinfo = np.ndarray((totaldatanum,totaldatanum),dtype = list)
    align_threads = []

    #align(dataset1[0],dataset2[0],250,2,longestlist,exceedlist)

    start = time.clock()

    if parallelchoice == 1:
        for i in range(0,totaldatanum-1):
            for j in range(i+1,totaldatanum):
                data1 = dataset[i]
                data2 = dataset[j]
                athread = threading.Thread(target = align,args = (data1,data2,minlength,errtol,longestlist,exceedlist,longestinfo,exceedinfo,plotthreshold,i,j))
                align_threads.append(athread)
                athread.start()
                time.sleep(8)
        #i = 1
        #for data1 in dataset[0:10]:
        #    j = 1
        #    for data2 in dataset[i+1:60]:
        #        athread = threading.Thread(target = align,args = (data1,data2,minlength,errtol,longestlist,exceedlist,plotthreshold,i,j))
        #        align_threads.append(athread)
        #        athread.start()
        #        j += 1
        #    i += 1

        for athread in align_threads:
            athread.join()

    elif parallelchoice == 2:
        i = 1
        for data1 in dataset1[0:2]:
            aomthread = threading.Thread(target = alignonetomulti,args = (data1,dataset2[0:2],minlength,errtol,longestlist,exceedlist,longestinfo,exceedinfo,plotthreshold,i))
            align_threads.append(aomthread)
            aomthread.start()
            i += 1

        for aomthread in align_threads:
            aomthread.join()

    timeused = time.clock()-start

    #print(longestlist)
    #print(exceedlist)
    print('Time used: ',timeused)

    plotresult(longestlist,exceedlist)
    [coreindex,coreedgenum] = findcore(exceedinfo,totaldatanum)
    #print('Core genome: NO.%d'%(coreindex+1))
    #print('Core genome name: ',names[coreindex])
    #print('Core edge number: %d'%coreedgenum)
    #[markerstart,markerlength,markertimes] = findmarker(exceedinfo,coreindex,len(dataset[coreindex]),totaldatanum)
    #writecorereport(coreindex,coreedgenum,exceedinfo,names[coreindex],totaldatanum,markerstart,markerlength,markertimes)
    newexceedinfo = compmarkers(coreindex,exceedinfo,longestinfo,totaldatanum)
    allmarkers = findallmarkers(newexceedinfo,coreindex,len(dataset[coreindex]),totaldatanum)
    print(allmarkers)
    writereport(longestinfo,exceedinfo,names,totaldatanum,timeused,plotthreshold,errtol)
    writecorereport(coreindex,coreedgenum,newexceedinfo,names[coreindex],names,totaldatanum,allmarkers)

