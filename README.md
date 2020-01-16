# BLAST-to-find-the-k-edits-LCS-implemented-with-Python
BLAST with Pigeon Hole Theory to find the k-edits LCS implemented with Python

The code contains a API and an script to analyze the hiv dataset.

The 3 codes "MyBlastPy.py", "Seed.py","BandedDP.py" are the API I built.
It can compair two genome strings, with the following code:

    testblast = MyBlast(str1,str2,minlen,errtol)
    testblast.startseeding()
    firstresult = testblast.startDP()
    firstmaxlen = testblast.printfirstresult()
    secondresult = testblast.seconditer()

The input parameters are:
str1: the first genome string to compair
str2: the second genome string to compair
minlen: the minimum reqiured length of the LCS. It is used to search the seeds.
errtol: the edit tolerence of the LCS including insertion, deletion and mutation.

The output parameter is a python list that contains the information of all the LCS of which the length is longer than the minimum required length. And if there is not LCS longer than that, it would record the longest.
For each common string. The data structure is as follow:
[startindex1,startindex2,length1,length2]
startindex1 is the start point index in genome1.
the rest is as the first one.
So, the whole output list is like this:
[[startindex1_1,startindex1_2,length1_1,length1_2],[startindex2_1,startindex2_2,length2_1,length2_2],.......,[startindexn_1,startindexn_2,lengthn_1,lengthn_2]]
that represent the information of n common strings.

The script is to run all against all comparison in the dataset. It will find the distribution of common strings according to their length and show it in a bar chart. It can find the core genome in the dataset that is aligned longer than the minimum length for the most times. It will also find the markers in the core genome. The markers can be found in all the genomes in the dataset.

When running the script, the datas should be include in the fold "./hivgenomes".
When you run the script, there will be parameters that youcan choose.
The first line shows the total number of genomes in your dataset.
The second is to choose the number of genomes you need to run all against all comparison.
The third is the minimum length you require.
The forth is the edit tolerence.
The fifth is the plot threshold. Usually, it's the same as the third parameter.
The six is parallel method of the comparison. This you must choose "1".

The script will generate 4 plots and two txt files to record the result.
