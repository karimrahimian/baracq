from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
class FastQ():
    def __init__(self,filename):
        self.filename = filename
        self.monitoring_step = 10_000
        self.reading_limit = 500_000
        pass
    def Initialize(self):
        self.read_length = self.Get_Read_Length()
        self.read_lengths = []
        self.count_read = None
        empty_array = np.zeros((1,self.read_length))
        self.scores = np.array(empty_array)
    def Get_Read_Length(self):
        for i,record in enumerate(SeqIO.parse(self.filename, "fastq")):
            read = str(record.seq)
            return len(read)
    def ReadFastLengths(self):
        temp_array = []
        for i,record in enumerate(SeqIO.parse(self.filename, "fastq")):
            read_length = len(str(record.seq))
            temp_array.append(read_length)

            if (i%10000==0):
                print(i)
                self.read_lengths = np.append(self.read_lengths,temp_array)
                temp_array = []
            if (i > self.reading_limit==0):
                break
        print("Unique")
        print (np.unique(self.read_lengths))
    def GetFirstRead(self):
        for i, record in enumerate(SeqIO.parse(self.filename, "fastq")):
            id = record.id
            sequence = str(record.seq)
            score_number = record.letter_annotations["phred_quality"]
            #score_letter = recor
            break
        print(f"{id}\n{sequence}\n{score_number}")
    def GetReadCount(self):
        for i,record in enumerate(SeqIO.parse(self.filename, "fastq")):
            pass
        print (f"Read Count = {i}")

    def ReadFastQScores(self):
        print ("Reading Fastq file started")
        temp_array = []
        for i,record in enumerate(SeqIO.parse(self.filename, "fastq")):
            #seq = str(record.seq)

            score = record.letter_annotations["phred_quality"]
            temp_array.append(score)
            if (i % 5000 == 0):
                self.scores = np.append(self.scores,np.array(temp_array),axis=0)
                temp_array = []
            if (i % self.monitoring_step ==0):
                print(i)
            if (i>self.reading_limit):
                break
        print("Reading Fastq file ended")
    def FindMatch(self,sub_sequence):
        # Open the FASTQ file and loop through the records
        count = 0
        with open(self.filename, "r") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                # Count the number of occurrences of the subsequence in the record
                count += record.seq.count(sub_sequence)
        return count
    def DrawBoxPlot(self):
        colnames = []
        for i in range(self.read_length):
            colnames.append("Base "+str(i+1))
        fig = plt.figure(figsize=(15, 5), dpi=300)
        df = pd.DataFrame(self.scores,columns=colnames)
        boxplot = df.boxplot(column=colnames,showfliers=False,grid=False)
        boxplot.set_xticklabels(boxplot.get_xticklabels(), rotation=90)


        plt.title("Quality of each position in read")
        plt.xlabel("Bases")
        plt.ylabel("Quality Score")
        plt.savefig("box.png",bbox_inches="tight")
    def DrawDensityPlot(self):
        fig = plt.figure(figsize=(3, 3), dpi=100)
        sns.distplot(self.read_lengths,bins=1)
        #sns.kdeplot(self.read_lengths)
        #sns.despine()
        plt.xticks([100],["100"])
        plt.tight_layout()
        plt.savefig("density.png", bbox_inches="tight")
    def PrintFirstThosandLines(self):
        # Open the FASTQ file and read in the first 1000 lines
        for i,record in enumerate(SeqIO.parse(self.filename, "fastq")):
            if i >= 4000:  # stop reading after the 4th line of the 1000th record
                break
            read =  record.seq
            id = record.id
            print(f"{id}_{read}")


fastq = FastQ("../fastq/short_reads.fastq")
fastq.Initialize()
fastq.ReadFastLengths()
fastq.DrawDensityPlot()

fastq.ReadFastQScores()
fastq.DrawBoxPlot()
fastq.PrintFirstThosandLines()



fastq.GetReadCount()
fastq.GetFirstRead()
count = fastq.FindMatch("TTAAATGGAA")
print (f"Count of occur :{count}")


