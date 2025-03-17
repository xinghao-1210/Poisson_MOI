# -*- coding: utf-8 -*-
"""
Created on Sat Jul 21 01:05:07 2018

@author: Xinghao
"""

import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
import pandas as pd
moi=float(input('MOI:\n'))
n=np.arange(0,11)
p=st.poisson.pmf(n,moi)
p_df=pd.DataFrame((p*100).reshape(1,11),index=['%'],columns=n)

plt.figure()
plt.plot(n,p)
plt.title('Poisson:' + str(moi))
plt.xlabel('Particle/Cell')
plt.xticks(np.arange(0,11,1))
plt.ylabel('Probability of Particle/Cell')
plt.show()
sr=(1-p[0])*100
print('Survival Rate(%):\n',sr)
print('Probability(%):\n',p_df)
