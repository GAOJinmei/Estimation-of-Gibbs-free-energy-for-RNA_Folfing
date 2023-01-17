import numpy as np
import pandas as pd
# for calculating the distance
import math
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
import argparse
parser = argparse.ArgumentParser(description="Generate 10 files with 20 scors in each file")
parser.add_argument("-d","--directory",type = str,help="The directory you put all the pdb file",default="pdb/")

args = parser.parse_args()
dir_path = args.directory


# create an empty dataframe
# N means the name of the nucleuotide
rna_df = pd.DataFrame(columns=['Atom', 'N', 'Chain', 'X','Y','Z'])
rna_df

###################Read the pdb files line by line ##################

# read the 10 pdb files and choose the C3' atoms and put them into the dataframe
i = 0
for path in os.listdir(dir_path):
    if 'pdb' in path:
        p = '{}{}'.format(dir_path,path)
        print(p)
        with open(p,'r') as pdbfile:
            for line in pdbfile:
                if line[:4] == 'ATOM' and line[13:16] == "C3'":
                    rna_data = [line[13:16],line[19:20],line[21:22]+str(i),line[32:38],line[40:46],line[48:54]]
                    # rna_df = rna_df.append(rna_data, ignore_index=True)
                    rna_df.loc[len(rna_df)] = rna_data
            i += 1
                    
# convert data type 
rna_df = rna_df.astype({'X':'float','Y':'float','Z':'float'})


# calculate the distance between each two C3 and put the results into an array
def distances(df):
    ns = ['A','C','G','U']
    dic_pair={}
    # create an empty dictionary with 10 keys
    for a in range(len(ns)):
        for b in range(a,len(ns)):
            key = '{}{}'.format(ns[a],ns[b])
            dic_pair[key]=[]
    # dic_pair = create_dic()
    # distance_arr = np.empty([r,r])
    r = df.shape[0]
    for i in range(r-4):
        for j in range(i+4,r):
            # check if the pair is belong to the same chain
            if df.loc[i,'Chain'] == df.loc[i,'Chain']:
                # calculate the distance between two C3
                cordinates_i = list(rna_df.iloc[i,3:6])
                cordinates_j = list(rna_df.iloc[j,3:6])
                distance = math.dist(cordinates_i,cordinates_j)
                pair = '{}{}'.format(rna_df.loc[i,'N'],rna_df.loc[j,'N'])
                if distance < 21:
                    if pair in dic_pair.keys():
                        dic_pair[pair].append(int(distance))
                    else:
                        pair = '{}{}'.format(rna_df.loc[j,'N'],rna_df.loc[i,'N'])
                        dic_pair[pair].append(int(distance))
    return dic_pair
print('--------------- Calculating the distances...... -----------')
d = distances(rna_df)

###################### Calculate observed and reference frequencies ######################
# calculate the observed probability (i.e. frequency)
# input is the distances of i and j
def observe_frequency(Dijs):
    dic_obs={}
    for i in range(21):
        # the numebr of distances in distance i/the total number of distance of pair [i,j]
        f = Dijs.count(i)/len(Dijs)
        dic_obs[i]=f
    # plt.bar(dic_obs.keys(),dic_obs.values(),width=0.9,align='edge')
    # plt.xlabel('Distance')
    # plt.ylabel('OBS')
    # plt.show()
    return dic_obs

# calculate the reference frequency
def reference_frequency(d):
    dic_ref = {}
    for i in range(21):
        n_xx_d = 0
        n_xx = 0
        for k,v in d.items():
            n_xx_d += v.count(i)
            n_xx += len(v)
        dic_ref[i] = n_xx_d/n_xx
    return dic_ref

dic_f_ref=reference_frequency(d)

# score
def score(dic_f_obs,dic_f_ref):
    dic_score= {}
    for k,v in dic_f_obs.items():
        s = 10
        if v !=0 and dic_f_ref[k] != 0:
            s = (-1) * math.log(v/dic_f_ref[k])
        dic_score[k]=s
    return dic_score

def write_txt(pair, dic):
    
    os.system('rm {}.txt'.format(pair))
    f = open('{}.txt'.format(pair),'x')
    for k,v in dic.items():
        f = open('{}.txt'.format(pair),'a')
        f.write('{}'.format(k) + '\t' + '{}'.format(v)+'\n')

        f.close()

def calculate_scores():
    for k,v in d.items():
        pair = '{}'.format(k)
        dic_f_obs = observe_frequency(v)
        # print('---------Paire {} --------'.format(k))
        dic_score = score(dic_f_obs,dic_f_ref)
        # print('---------Scores---------')
        tqdm(write_txt(pair, dic_score))
        # print(dic_score)

print('--------------- Calculating the distances...... -----------')
calculate_scores()
print('--------------- Program finished, you can check your directory -----------')





