import os
import re
import subprocess
import csv
import pandas as pd
import numpy as np
from statistics import stdev
from statistics import mean
path = "/home/paco/root_ex/TempdirSRR/"
#####
### Data download and dump
#####
with open("data_sources.sh", 'r') as file:
    lines = [line for line in file]
for i in range(168,170):
    TempdirSRR = os.path.exists("./TempdirSRR")
    if TempdirSRR == False:
        os.system("mkdir ./TempdirSRR")
    l = (lines[i])
    #print(l)
    SRA = re.search("SRR\\d+\\s", l)
    SRA = SRA.group()
    SRA = re.search("SRR\\d+", SRA)
    SRA = SRA.group()
    output = re.search("\\w+$", l)
    output = output.group()
    os.system("wget https://sra-pub-run-odp.s3.amazonaws.com/sra/"
         + SRA + "/" + SRA + " -O ./TempdirSRR/" + output)
    os.system("fastq-dump -v -F ./TempdirSRR/" + output + " --split-files -O ./TempdirSRR/")
#####
### Data download and dump finished!!!!
#####

#####
### Trim_galore
#####
    allfiles = sorted(os.listdir('/home/paco/root_ex/TempdirSRR/'))
    fq = []
    for file in allfiles:
        if re.search("(_[1-9]_[12]\.fastq)$", file):
            fq.append(file)
    #print(fq)
    #print(len(fq))
    if len(fq) == 1:
        p1 = (fq[0])
        p2 = (fq[0])
        print("Single-end reads detected: " + p1)
    elif len(fq) == 2:
        p1 = (fq[0])
        p2 = (fq[1])
        print("Paired-end reads detected: " + p1 + " & " + p2)
    if re.search("([1-9]_1\.fastq)$", p1) and re.search("([1-9]_2\.fastq)$", p2):
        os.system("trim_galore --paired ./TempdirSRR/" + p1 + " ./TempdirSRR/" + p2 +
        " -j 34 --fastqc -o " + path +
        " --clip_R1 5 --clip_R2 5 --three_prime_clip_R1 2 --three_prime_clip_R2 2")
    elif re.search("([1-9]_1\.fastq)$", p1) and re.search("([1-9]_1\.fastq)$", p2):
        os.system("trim_galore ./TempdirSRR/" + p1 + 
        " -j 34 --fastqc -o " + path + " --clip_R1 5 --three_prime_clip_R1 2")
    #####
    ### Trim_galore Finished!!!
    ##### 
    
    #####
    ### Kallisto
    #####
    #os.system(mkdir counts)
    path = "/home/paco/root_ex/TempdirSRR/"
    allfiles = sorted(os.listdir('/home/paco/root_ex/TempdirSRR/'))
    #print(allfiles)
    fq = []
    for file in allfiles:
        if re.search("(((val_[12])|(trimmed))\.((fastq)|(fq)))$", file):
            fq.append(file)
    #print(fq)
    if len(fq) == 1:
        p1 = (fq[0])
        p2 = (fq[0])
        print("Single-end reads detected: " + p1)
        print(p1)
        print(p2)
    elif len(fq) == 2:
        p1 = (fq[0])
        p2 = (fq[1])
        print(p1)
        print(p2)
    if re.search("(val_1\.((fastq)|(fq)))$", p1) and re.search("(val_2\.((fastq)|(fq)))$", p2):
        print("Paired-end reads detected:")
        os.system("kallisto quant -i /home/paco/root_ex/SRA/Arath.index -t 30 -o /home/paco/root_ex/counts/" + p1[:-11] + " ./TempdirSRR/" + p1 + " ./TempdirSRR/" + p2)
    elif re.search("(1_trimmed\.((fastq)|(fq)))$", p1) and (re.search("(1_trimmed\.((fastq)|(fq)))$", p2) or re.search("(val_1\.((fastq)|(fq)))$", p2)):
        print("Single-end reads detected:")
        os.system("unzip /home/paco/root_ex/TempdirSRR/\*.zip -d ./TempdirSRR/")
        l = subprocess.check_output("grep -zoP '(?=Length\sCount)(.|\\n)*?(?<=END_MODULE)' " + path + p1[:-3] + "_fastqc/fastqc_data.txt", shell=True)
        l = str(l)
        l = l.split('\\n')
        l = l[:-1]
        l.pop(0)
        l = pd.DataFrame(l)
        l = l.set_axis(['x'], axis=1)
        l[['length', 'counts']] = l['x'].str.split(pat = '\\\\t', expand = True)
        l = l.drop(['x'], axis=1)
        read_length = []
        read_count = []
        read_count_double = []
        read_count_single = []
        mean_l = []
        for i in range(0, l.shape[0]):
            if re.search("-", l.at[i, 'length']):
                read_length.append(l.at[i, 'length'])
                read_count.append(l.at[i, 'counts'])
                read_count_double.append(float(read_count[i])/2)
                read_count_double.append(float(read_count[i])/2)
            else:
                read_length.append(l.at[i, 'length'])
                read_count.append(l.at[i, 'counts'])
                read_count_single.append(float(read_count[i]))
        if len(read_count_double) == 0:
            print("calculating -l and -s from file " + p1)
            mean_l = np.repeat(read_length, read_count_single, axis = 0)
            mean_l = [int(n) for n in mean_l]
            os.system("kallisto quant -i /home/paco/root_ex/SRA/Arath.index -t 30 --single -o /home/paco/root_ex/counts/" +
                p1[:-11] +
                " -l " + str(mean(mean_l)) + 
                " -s " + str(stdev(mean_l)) + " ./TempdirSRR/" + p1
            )
        else:
            print("calculating -l and -s from file " + p1)
            read_length = [i.split('-', 1) for i in read_length]
            read_length = [item for sublist in read_length for item in sublist]
            mean_l = np.repeat(read_length, read_count_double, axis = 0)
            mean_l = [int(n) for n in mean_l]
            os.system("kallisto quant -i /home/paco/root_ex/SRA/Arath.index -t 30 --single -o /home/paco/root_ex/counts/" +
                p1[:-11] +
                " -l " + str(mean(mean_l)) + 
                " -s " + str(stdev(mean_l)) + " ./TempdirSRR/" + p1
            )
    print("deleting temporal files")
    os.system("rm -r ./TempdirSRR")