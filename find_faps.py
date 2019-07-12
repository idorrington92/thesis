
# coding: utf-8

import numpy as np
import os
from collections import Counter
import pandas as pd

faps =pd.read_fwf('mva_pvalues.txt',index_col=0,header=None)
faps = faps.drop('170817')
faps = faps.drop('E264930')


sorted_faps = sorted(faps.values)
sorted_faps = [float(i) for i in sorted_faps]



cumsum = np.cumsum(sorted_faps)
frac = np.array(range(1,len(cumsum)+1))/float(len(cumsum))


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy.random as rnd
import scipy.special as spec


prob = rnd.uniform(low=0.,high=1.,size=(int(1e5),len(cumsum)))
sorted_prob = [sorted(row) for row in prob]

plow = np.percentile(sorted_prob,100*spec.erfc(2/np.sqrt(2))/2,axis=0)
phigh = np.percentile(sorted_prob,100*(1-spec.erfc(2/np.sqrt(2))/2),axis=0)


matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.serif'] = 'Computer Modern Roman'
matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']
matplotlib.rcParams.update({'font.size': 26})

fig, ax = plt.subplots(figsize=(10,10))
ax.loglog(sorted_faps,frac,'b-',label='Observed',marker='^',linewidth=0.5,markersize=5)
ax.loglog(plow,np.linspace(1/len(plow),1,len(plow)),'k:')
ax.loglog(phigh,np.linspace(1/len(phigh),1,len(phigh)),'k:')
ax.loglog(frac,frac,'k--',label='Expected')
ax.legend(loc=4)
plt.xlim([min(sorted_faps)*0.9,1])
#plt.xlim([0.8,1])
plt.ylim([min(frac)*0.9,1])
#plt.ylim([0.8,1])
plt.xlabel('$p$-value')
plt.ylabel('Fraction of GRBs')
plt.grid()

plt.tight_layout()
plt.savefig('mva_pvalue.png', ppi=100)

