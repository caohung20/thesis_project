import psutil
process = psutil.Process()
print("Memory usage:", process.memory_info().rss)

import time
time_start=time.time()
import numpy as np
import pandas as pd
import scipy as sp
from scipy.sparse import csc_matrix
import sys

def value(a,b):
    value1=Ybus.getrow(a).getcol(b).toarray()[0,0]
    return value1

j = 1j  # defining imaginary number as j 
linedata=pd.read_excel("G:\do_an_tot_nghiep\lineData.xlsx",sheet_name = 0)
nl = linedata['bus_nl'].values
nr = linedata['bus nr'].values
R = linedata['R'].values
X = linedata['X'].values
Bc = j*linedata['B'].values
a = linedata['Line code and tap setting'].values
nbr = len(linedata['bus_nl'])
nbus = max(max(nl), max(nr))
Z=np.ones(nbr,dtype=complex)
y = np.ones((nbr,1), dtype=complex)
subYbus=np.empty(2*nbr,dtype=complex)
for i in range(nbr):
    Z[i]= R[i] + j*X[i]
    # defining branch admittance as complex
    y[i]=y[i]/Z[i] 

    if a[i] <= 0:
        a[i] = 1
    else: 
        pass
    subYbus[i]=- y[i] / a[i]
    #Ybus matrix is symmetric, so double subYbus 
    # to calculate bottom half elements of diagonal
    subYbus[nbr+i]=subYbus[i]
b=0
for n in range(nbus):              
    #calculate diagonal elements
    for k in range(nbr):
        if nl[k] == n+1:
            b+= y[k] / (a[k]**2) + Bc[k]
        elif nr[k] == n+1:
            b += y[k] + Bc[k]
        else:
            pass
    subYbus=np.append(subYbus,b)
    b=0
for i in range(nbr):       
    #subtract to make nr and nl become indices
    nl[i]-=1
    nr[i]-=1
#Ybus matrix is symmetric, so double indices 
# to calculate bottom half elements of diagonal
snl= np.concatenate((nl, nr))
snr=np.concatenate((nr, nl))

for i in range(nbus):       
    #append diagonal element to nl and nr for indices
    snl=np.append(snl,i)
    snr=np.append(snr,i)
Ybus=csc_matrix((subYbus, (snl, snr)),shape = (nbus, nbus))#.toarray()
print(Ybus)
time_end=time.time()
print(time_end-time_start)
print(sys.getsizeof(Ybus))      #return required memory to run program 
print(np.angle(value(0,0),deg=False))