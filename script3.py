
import numpy as np
import pandas as pd
import math
import argparse
import pathlib
parser = argparse.ArgumentParser(description="Test the  Gibbs free energy of one RNA")
parser.add_argument("-pdb","--pdb_file",type = pathlib.Path,help="your pdb file path",default="2lx1.pdb")

args = parser.parse_args()
pdb_path = args.pdb_file
################ Read pdb file #################

# create an empty dataframe
# N means the name of the nucleuotide
rna_df = pd.DataFrame(columns=['Atom', 'N', 'Chain', 'X','Y','Z'])
rna_df

# read the 10 pdb files and choose the C3' atoms and put them into the dataframe
with open(pdb_path,'r') as pdbfile:
    for line in pdbfile:
        if line[:4] == 'ATOM' and line[13:16] == "C3'":
            rna_data = [line[13:16],line[19:20],line[21:22],line[32:38],line[40:46],line[48:54]]
            rna_df.loc[len(rna_df)] = rna_data


# convert data type 
rna_df = rna_df.astype({'X':'float','Y':'float','Z':'float'})
rna_df.dtypes


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
                cordinates_i = rna_df.iloc[i,3:6]
                cordinates_j = rna_df.iloc[j,3:6]
                distance = math.dist(cordinates_i,cordinates_j)
                pair = '{}{}'.format(rna_df.loc[i,'N'],rna_df.loc[j,'N'])
                if distance < 20:
                    if pair in dic_pair.keys():
                        dic_pair[pair].append(distance)
                    else:
                        pair = '{}{}'.format(rna_df.loc[j,'N'],rna_df.loc[i,'N'])
                        dic_pair[pair].append(distance)
    return dic_pair

print('--------------- Calculating the distances...... -----------')
d = distances(rna_df)


################# Linear intepretion ##################

def linear_intepretion(dist_test,pair):
    df = pd.read_csv('{}.txt'.format(pair), header = None,sep='\t')
    if dist_test>=0 and dist_test<=0.5 or (dist_test - math.floor(dist_test)) > 0.5:
        x1, x2 = math.floor(dist_test), math.ceil(dist_test)
    else:
        x1, x2 = math.floor(dist_test)-1, math.ceil(dist_test)-1
    y1, y2 = df[1][x1],df[1][x2]
    E_score = y1 + (dist_test - x1) * ((y2 - y1) / (x2 - x1))
    return E_score
        

#iterate the distances
def calculate_energy():
    Gibbs_free_energy = 0
    for k,v in d.items():
        for dist in v:
            E_score = linear_intepretion(dist,k)
            Gibbs_free_energy += E_score
    return Gibbs_free_energy

print('--------------- Calculating the Free energy...... -----------')
Gibbs_free_energy = calculate_energy()

print('The predicted Gibbes free energy for the RNA is: {}'.format(Gibbs_free_energy))









