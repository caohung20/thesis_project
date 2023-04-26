import time
time_start=time.time()
import numpy as np
import pandas as pd
import sys
# defining imaginary number as j
j = 1j   
linedata=pd.read_excel("G:\do_an_tot_nghiep\Inputs12_2.xlsx",sheet_name = 2,skiprows=1)
nl = linedata['FROMBUS'].values
nr = linedata['TOBUS'].values
R = linedata['R(Ohm)'].values
X = linedata['X(Ohm)'].values
Bc = j*linedata['B(microSiemens)'].values
basemva=100
basevolt=12
zbase=basevolt**2/basemva
#a = linedata['Line code and tap setting'].values
nbr = len(linedata['FROMBUS'])
print(nbr)
nbus = max(max(nl), max(nr))
Z=np.ones(nbr,dtype=complex)
y = np.ones((nbr,1), dtype=complex)
for i in range(0,nbr):
    #convert to per unit system
    Z[i]= (R[i] + j*X[i])/zbase
    # defining branch admittance as complex
    y[i]=y[i]/Z[i]
    Bc[i]=Bc[i]/2*zbase*10**(-6)
#calculate with tap setting
'''
for n in range(nbr):
    if a[n] <= 0:
        a[n] = 1
    else: 
        pass  
'''
for n in range(nbr):
    # initialize Ybus to zero
    Ybus = np.zeros((int(nbus),int(nbus)), dtype=complex)  
    
    for k in range(nbr):
        #[nl[k]-1] => correct position in Ybus
        Ybus[nl[k]-1][nr[k]-1] = Ybus[nl[k]-1][nr[k]-1]  - y[k] #/ a[k] 
        Ybus[nr[k]-1][nl[k]-1] = Ybus[nl[k]-1][nr[k]-1]
#formation of the diagonal elements
for n in range(nbus):
    for k in range(nbr):
        if nl[k] == n+1:
            'Ybus[n,n] += y[k] / (a[k]**2) + Bc[k]'
            Ybus[n,n] += y[k] + Bc[k]
        elif nr[k] == n+1:
            Ybus[n,n] += y[k] + Bc[k]
        else:
            pass
print(Ybus)
time_end=time.time()
print(time_end-time_start)
print(sys.getsizeof(Ybus))
