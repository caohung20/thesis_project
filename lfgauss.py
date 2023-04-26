# Power flow solution by Newton-Raphson method
import psutil
process = psutil.Process()
print("Memory usage:", process.memory_info().rss)
import numpy as np
import lfybus
import pandas as pd
import cmath
nbus=lfybus.nbus
Ybus=lfybus.Ybus
nbr=lfybus.nbr
nl=lfybus.nl
nr=lfybus.nr
basemva=100 #base power to convert to per unit
accuracy=0.001
maxiter=100 #maximum of iteration
ns=0
ng=0
Pgg=np.zeros(nbus,dtype=complex)
Qgg=np.zeros(nbus,dtype=complex)
ngs=np.zeros(nbus)
nss=np.zeros(nbus)
busdata=pd.read_excel("G:\do an tot nghiep\myExample.xlsx",sheet_name = 0)
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
DV=np.zeros(nbus)
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
num = 0
AcurBus = 0
converge = 1
Vc = np.zeros((nbus, 1), dtype=complex)
Sc = np.zeros((nbus, 1), dtype=complex)   
#generate value if variable does not exist 
while 'accel' not in locals():
    accel = 1.3

while 'accuracy' not in locals():
    accuracy = 0.001

while 'basemva' not in locals():
    basemva = 100

while 'maxiter' not in locals():
    maxiter = 100
iter = 0
maxerror = 10
while maxerror >= accuracy and iter <= maxiter:
    iter += 1
    import numpy as np
    for n in range(nbus):
        YV = 0 + 1j * 0
        for L in range(nbr):
            if nl[L] == n:
                k = nr[L]
                YV += Ybus[n, k] * V[k]
            elif nr[L] == n:
                k = nl[L]
                YV += Ybus[n, k] * V[k]
        Sc = np.conj(V[n]) * (Ybus[n, n] * V[n] + YV)
        Sc = np.conj(Sc)
        DP[n] = P[n] - np.real(Sc)
        DQ[n] = Q[n] - np.imag(Sc)
        if kb[n] == 1:
            S[n] = Sc
            P[n] = np.real(Sc)
            Q[n] = np.imag(Sc)
            DP[n] = 0
            DQ[n] = 0
            Vc[n] = V[n]
        elif kb[n] == 2:
            Q[n] = np.imag(Sc)
            S[n] = P[n] + 1j * Q[n]
            if Qmax[n] != 0:
                Qgc = Q[n] * basemva + Qd[n] - Qsh[n]
                if abs(DQ[n]) <= 0.005 and iter >= 10:
                    if DV[n] <= 0.045:
                        if Qgc < Qmin[n]:
                            Vm[n] += 0.005
                            DV[n] += 0.005
                        elif Qgc > Qmax[n]:
                            Vm[n] -= 0.005
                            DV[n] += 0.005
                    else:
                        pass
                else:
                    pass
            else:
                pass
        if kb[n] != 1:
            Vc[n] = (np.conj(S[n]) / np.conj(V[n]) - YV) / Ybus[n, n]
        else:
            pass
        if kb[n] == 0:
            V[n] = V[n] + accel * (Vc[n] - V[n])
        elif kb[n] == 2:
            VcI = np.imag(Vc[n])
            VcR = np.sqrt(Vm[n] ** 2 - VcI ** 2)
            Vc[n] = VcR + 1j * VcI
            V[n] = V[n] + accel * (Vc[n] - V[n])
    maxerror = max(max(abs(np.real(DP))), max(abs(np.imag(DQ))))
    if iter == maxiter and maxerror > accuracy:
        print('\nWARNING: Iterative solution did not converge after', iter, 'iterations.\n\n')
        input('Press Enter to terminate the iterations and print the results \n')
        converge = 0
    else:
        pass
if converge != 1:
    tech = '                      ITERATIVE SOLUTION DID NOT CONVERGE'
else:
    tech = '                   Power Flow Solution by Gauss-Seidel Method'
k = 0
for n in range(nbus):
    Vm[n] = abs(V[n])
    deltad[n] = cmath.phase(V[n]) * 180 / cmath.pi
    if kb[n] == 1:
        S[n] = P[n] + 1j * Q[n]
        Pg[n] = P[n] * basemva + Pd[n]
        Qg[n] = Q[n] * basemva + Qd[n] - Qsh[n]
        k = k + 1
        Pgg[k] = Pg[n]
    elif kb[n] == 2:
        k = k + 1
        Pgg[k] = Pg[n]
        S[n] = P[n] + 1j * Q[n]
        Qg[n] = Q[n] * basemva + Qd[n] - Qsh[n]
    yload[n] = (Pd[n] - 1j * Qd[n] + 1j * Qsh[n]) / (basemva * Vm[n]**2)

Pgt = sum(Pg)
Qgt = sum(Qg)
Pdt = sum(Pd)
Qdt = sum(Qd)
Qsht = sum(Qsh)
busdata[:, 2] = Vm
busdata[:, 3] = deltad

    
