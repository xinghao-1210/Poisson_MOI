# -*- coding: utf-8 -*-
"""
Created on Sat Jul 21 01:05:07 2018

@author: Xinghao
"""

import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import six
import pandas as pd
from scipy.optimize import curve_fit
from sympy import *
#Percentage transduction-MOI plot
sr=[]
moi=np.arange(0.01,8.1,0.01)
for i in moi:
    n=0
    p=st.poisson.pmf(n,i)
    sr.append((1-p)*100)

#MOI取x时，至少有一个转染的细胞百分比（数值曲线）
def func_sr(x):
    return (1-np.exp(-x))*100
#MOI取x时，至少有一个转染的细胞百分比（拟合曲线）
def func_fit(x,a,b,c,d):
    return a+b*np.exp(c*x**d)
popt,pcov=curve_fit(func_fit,moi,sr)
a=popt[0]
b=popt[1]
c=popt[2]
d=popt[3]
print('a=%d,b=%d,c=%d,d=%d'%(a,b,c,d))

#计算Percentage transduction
v_name=str.upper(input('Transduced virus:\n'))
y_select1=float(input('Transduced cell number with selection_1: '))
y_nselect1=float(input('Transduced cell number without selection_1: '))
y_select2=float(input('Transduced cell number with selection_2: '))
y_nselect2=float(input('Transduced cell number without selection_2: '))
y=np.mean([y_select1/y_nselect1*100,y_select2/y_nselect2*100])
print('\nPercentage transduction(%%) (%s): %f' %(v_name, y)) #Percentage-transduction （Cell number selection/NO selection）
#根据实验Percentage transduction计算实验MOI
x=symbols('x')
moi_x=float(solve(1-exp(-x)-y/100)[0])
print('MOI=%.4f'%moi_x)

#plt.xkcd() #xkcd style
plot_moi=plt.figure(figsize=(12,6))
plt.plot(moi,sr,markersize=5,marker='x',label='raw data')
plt.plot(moi,func_sr(moi),color='r',label='function curve')
plt.plot(moi,func_fit(moi,a,b,c,d),color='y',label='fit curve')
plt.plot(moi,y+moi*0,'--',color='g',label='Experimental Percentage')
plt.plot(moi_x,y,'o',markersize=7,color='r')
plt.plot([moi_x,moi_x],[0,y],color='r',ls='--',linewidth=2)
plt.plot([moi_x,0],[y,y],color='r',ls='--',linewidth=2)
plt.annotate(('(%.4f,%.4f)'%(moi_x,y)),
             xy=(moi_x,y), xycoords='data',
             xytext=(+10, -30), textcoords='offset points', fontsize=12,
             arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))

plt.title('Percentage transduction-MOI (%s)'%v_name)
plt.xlabel('MOI')
plt.xticks(np.arange(0,8.5,0.5),fontsize=10)
plt.ylabel('Percetage transduction(%)')
plt.legend(loc='lower right')
plt.show()

#根据实验MOI计算P(1)
n=np.arange(0,11)
p=st.poisson.pmf(n,moi_x)
p1=st.poisson.pmf(1,moi_x)
p_df=pd.DataFrame((p*100).reshape(1,11),index=['%'],columns=n)
pd.set_option('precision',4)

#plt.xkcd() #xkcd style
plot_poisson=plt.figure(figsize=(12,6))
plt.plot(n,p*100,'o--')
plt.plot(1,p1*100,'o',markersize=7,color='r')
plt.plot([1,1],[0,p1*100],color='r',ls='--',linewidth=2)
plt.plot([1,0],[p1*100,p1*100],color='r',ls='--',linewidth=2)
plt.annotate(('(%.4f,%.4f)'%(1,p1*100)),
            xy=(1,p1*100),xycoords='data',
            xytext=(+10, +20), textcoords='offset points', fontsize=12,
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))

plt.title('Poisson: MOI=%.4f (%s)'%(moi_x,v_name))
plt.xlabel('Particle/Cell')
plt.xticks(np.arange(0,11,1))
plt.ylabel('Probability of Particle/Cell(%)')
plt.show()
sr=(1-p[0])*100
print('Survival Rate(%%):\n %.4f'%sr)
print('Probability(%):\n',p_df)

plots_pdf=PdfPages('%s MOI_POISSON.pdf'%v_name)
plots_pdf.savefig(plot_moi)
plots_pdf.savefig(plot_poisson)
plots_pdf.close()
writer=pd.ExcelWriter('%s MOI_POISSON.xlsx'%v_name,engine='xlsxwriter')
p_df.to_excel(writer,sheet_name='%s MOI_POISSON'%v_name)
writer.save()
