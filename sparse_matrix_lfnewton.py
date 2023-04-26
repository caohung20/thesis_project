# Power flow solution by Newton-Raphson method
import time
time_start=time.time()
import psutil
process = psutil.Process()
print("Memory usage:", process.memory_info().rss)
import sys
import numpy as np
#import sparse_matrix_lfybus as spYbus
import pandas as pd
import math

def value(row,col):
    #return value of a specific cell in Ybus sparse matrix
    value=Ybus.getrow(row).getcol(col).toarray()[0,0]
    return value
def angle(a):
    #return angle of a complex
    phase_angle =np.angle(a,deg=False)
    return phase_angle
import sparse_matrix_lfybus as spYbus  

nbus=spYbus.nbus
Ybus=spYbus.Ybus
nbr=spYbus.nbr
nl=spYbus.nl
nr=spYbus.nr

basemva=100 #base power to convert to per unit
accuracy=0.001
maxiter=100 #maximum of iteration
ns=0
ng=0
Pgg=np.zeros(nbus,dtype=complex)
Qgg=np.zeros(nbus,dtype=complex)
ngs=np.zeros(nbus)
nss=np.zeros(nbus)
busdata=pd.read_excel("G:\do_an_tot_nghiep\myExample.xlsx",sheet_name = 0)
Vm=np.zeros(nbus,dtype=complex)
delta=np.zeros(nbus,dtype=complex)
yload=np.zeros(nbus,dtype=complex)
deltad=np.zeros(nbus,dtype=complex)
kb=np.zeros(nbus)
Pd=np.zeros(nbus)
Qd=np.zeros(nbus)
Pg=np.zeros(nbus,dtype=complex)
Qg=np.zeros(nbus,dtype=complex)
Qmin=np.zeros(nbus)
Qmax=np.zeros(nbus)
Qsh=np.zeros(nbus)
V=np.zeros(nbus,dtype=complex)
P=np.zeros(nbus,dtype=complex)
Q=np.zeros(nbus,dtype=complex)
S=np.zeros(nbus,dtype=complex)
for k in range(nbus):
    n=int(busdata.iloc[k][0])
    n-=1
    kb[n]=busdata.iloc[k][1]    #bus code
    Vm[n]=busdata.iloc[k][2]    #voltage magnitude
    delta[n]=busdata.iloc[k][3] #angle degree
    Pd[n]=busdata.iloc[k][4]    #P load
    Qd[n]=busdata.iloc[k][5]    #Q load
    Pg[n]=busdata.iloc[k][6]    #P generator
    Qg[n]=busdata.iloc[k][7]    #Q generator
    Qmin[n]=busdata.iloc[k][8]  #Qmin generator
    Qmax[n]=busdata.iloc[k][9]  #Qmax generator
    Qsh[n]=busdata.iloc[k][10]  #injected Q from shunt capacitor
    if Vm[n] <= 0:
        Vm[n] = 1.0             #estimate flat start
        V[n] = 1 + 1j*0
    else:
        delta[n] = np.pi/180*delta[n]   #convert angle from degree to radian
        V[n] = Vm[n]*(np.cos(delta[n]) + 1j*np.sin(delta[n]))
        P[n]=(Pg[n]-Pd[n])/basemva      #express P Bus in per units
        Q[n]=(Qg[n]-Qd[n]+ Qsh[n])/basemva  #express Q Bus in per units
        S[n] = P[n] + 1j*Q[n]           
for k in range(nbus):
    if kb[k] == 1:
        ns += 1         #count the number of slack bus 
    elif kb[k] == 2:
        ng += 1         #count the number of P-V bus
    ngs[k] = ng         
    nss[k] = ns
#Ym=np.abs(Ybus)         #return Y bus magnitude
#t = np.angle(Ybus,deg=False)    # return phase angle of Ybus
m=2*nbus-ng-2*ns        #ng bus voltage control => m equations involving delta Q and delta V => corresponding columns of the jacobian matrix are eliminated. => 2*nbus-2-ng components left
maxerror = 1
converge=1
iter = 0
# Start of iterations
DC = np.zeros(m,dtype=complex)      
J = np.zeros((m,m))
while maxerror >= accuracy and iter <= maxiter:
    # Test for max power mismatch
    A = np.zeros((m,m),dtype=complex)     #initializing Jacobian matrix
    iter+=1
    for n in range(1, nbus+1):
        nn = int(n - nss[n-1])            #when slack bus appears, it will eliminate slack bus components in jacobi matrix    
        lm = int(nbus + n - ngs[n-1] - nss[n-1] - ns)   #(nbus -ns)=> start J3; (n - ngs[n-1] - nss[n-1]) is similar to the above equation
        J11 = 0
        J22 = 0
        J33 = 0
        J44 = 0
        for i in range (1,nbr+1):
            if (nl[i-1]+1) ==n or (nr[i-1]+1) ==n:      #search in data if i match n => to eliminate cases which buses aren't connected 
                if  (nl[i-1]+1) ==n:
                    l = nr[i-1]+1
                if (nr[i-1]+1)==n:                  #already minus 1 in previous file
                    l = nl[i-1]+1
                J11 += Vm[n-1] * Vm[l-1] * abs(value(n-1,l-1)) * math.sin((angle(value(n-1,l-1)) - delta[n-1] + delta[l-1]).real)        #diagonal elements of J1
                J33 += Vm[n-1] * Vm[l-1] * abs(value(n-1,l-1)) * math.cos((angle(value(n-1,l-1)) - delta[n-1] + delta[l-1]).real)        #diagonal elements of J3
                if kb[n-1] != 1:            #slackbus doesn't exist in J2 and J4
                    J22 += Vm[l-1] * abs(value(n-1,l-1)) * math.cos((angle(value(n-1,l-1)) - delta[n-1] + delta[l-1]).real)              #diagonal elements of J2
                    J44 += Vm[l-1] * abs(value(n-1,l-1)) * math.sin((angle(value(n-1,l-1)) - delta[n-1] + delta[l-1]).real)              #diagonal elements of J4
                else:
                    pass
            
                if kb[n-1] != 1 and kb[l-1] != 1:        #check if both nodes are not slack bus
                    lk = int(nbus + l - ngs[l-1] - nss[l-1] - ns)
                    ll = int(l - nss[l-1])
                    # off diagonal elements of J1
                    A[nn-1][ll-1] = -Vm[n-1] * Vm[l-1] * abs(value(n-1,l-1)) * math.sin((angle(value(n-1,l-1)) - delta[n-1] + delta[l-1]).real)
                
                    if kb[l-1] == 0:
                        # off diagonal elements of J2
                         A[nn-1][lk-1] = Vm[n-1] * abs(value(n-1,l-1)) * math.cos((angle(value(n-1,l-1)) - delta[n-1] + delta[l-1]).real)
                
                    if kb[n-1] == 0:
                        # off diagonal elements of J3
                        A[lm-1][ll-1] = -Vm[n-1] * Vm[l-1] * abs(value(n-1,l-1)) * math.cos((angle(value(n-1,l-1)) - delta[n-1] + delta[l-1]).real)
                
                    if kb[n-1] == 0 and kb[l-1] == 0:
                        # off diagonal elements of J4
                        A[lm-1][lk-1] = -Vm[n-1] * abs(value(n-1,l-1)) * math.sin((angle(value(n-1,l-1)) - delta[n-1] + delta[l-1]).real)
                else:
                    pass
            else:
                pass
        Pk = Vm[n-1]**2 * abs(value(n-1,n-1)) * math.cos((angle(value(n-1,n-1))).real) + J33
        Qk = -Vm[n-1]**2 * abs(value(n-1,n-1)) * math.sin((angle(value(n-1,n-1))).real) - J11
        if kb[n-1] == 1: # Swing bus P
            P[n-1] = Pk
            Q[n-1] = Qk
        if kb[n-1] == 2:
            Q[n-1] = Qk 
            if Qmax[n-1] != 0:
                Qgc = Q[n-1] * basemva + Qd[n-1] - Qsh[n-1]
                if iter <= 7:               #Between the 2th & 6th iterations
                    if iter > 2:            #the Mvar of generator buses are
                        if Qgc < Qmin[n-1]:   #tested. If not within limits Vm(n)
                            Vm[n-1] += 0.01   #is changed in steps of 0.01 pu to
                        elif Qgc > Qmax[n-1]: #bring the generator Mvar within
                            Vm[n-1] -= 0.01   #the specified limits.
        if kb[n-1] != 1:
            A[nn-1][nn-1] = J11         #diagonal elements of J1
            DC[nn-1] = P[n-1] - Pk
        if kb[n-1] == 0:
            A[nn-1][lm-1] = 2 * Vm[n-1] * abs(value(n-1,n-1)) * math.cos(angle(value(n-1,n-1))) + J22    #diagonal elements of J2
            A[lm-1][nn-1] = J33         #diagonal elements of J3
            A[lm-1][lm-1] = -2 * Vm[n-1] * abs(value(n-1,n-1)) * math.sin(angle(value(n-1,n-1))) - J44   #diagonal elements of J4
            DC[lm-1] = Q[n-1] - Qk
    #matrix A left division Matrix DC
    DX = np.linalg.solve(A, DC.T) 
    for n in range(1,nbus+1):
        nn = int(n - nss[n-1])
        lm = int(nbus + n - ngs[n-1] - nss[n-1] - ns)
        if kb[n-1] != 1:
            delta[n-1] +=  DX[nn-1]
        if kb[n-1] == 0:
            Vm[n-1] = Vm[n-1] + DX[lm-1]
    maxerror = max(abs(DC))
    if iter == maxiter and maxerror > accuracy:
        print('\nWARNING: Iterative solution did not converge after', iter, 'iterations.\n\n')
        print('Press Enter to terminate the iterations and print the results \n')
        converge = 0
        input()
    else:
        pass
if converge != 1:
    tech = '                      ITERATIVE SOLUTION DID NOT CONVERGE'
else:
    tech = '                   Power Flow Solution by Newton-Raphson Method'

V = Vm * np.cos(delta) + 1j * Vm * np.sin(delta)
deltad = 180 / np.pi * delta
i = 1j
k = 0
for n in range(nbus):
    if kb[n] == 1:
        k += 1
        S[n] = P[n] + 1j * Q[n]
        Pg[n] = P[n] * basemva + Pd[n]
        Qg[n] = Q[n] * basemva + Qd[n] - Qsh[n]
        Pgg[k-1] = Pg[n]
        Qgg[k-1] = Qg[n]
    elif kb[n] == 2:
        k += 1
        S[n] = P[n] + 1j * Q[n]
        Qg[n] = Q[n] * basemva + Qd[n] - Qsh[n]
        Pgg[k-1] = Pg[n]
        Qgg[k-1] = Qg[n]
    yload[n] = (Pd[n] - 1j * Qd[n] + 1j * Qsh[n]) / (basemva * Vm[n]**2)
Pgt = sum(Pg)
Qgt = sum(Qg)
Pdt = sum(Pd)
Qdt = sum(Qd)
Qsht = sum(Qsh)
time_end=time.time()
print(time_end-time_start)

