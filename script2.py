
import matplotlib.pyplot as plt
import pandas as pd


ns = ['A','C','G','U']
x = []
y = []
for a in range(len(ns)):
    for b in range(a,len(ns)):
        pair = '{}{}'.format(ns[a],ns[b])
        df = pd.read_csv('{}.txt'.format(pair), header = None, sep='\t')
        plt.plot(df[0], df[1],'--bo')
        plt.xlabel('Distance')
        plt.ylabel('Score')
        plt.savefig('{}.png'.format(pair),transparent = True,facecolor ="w")
        plt.show()





