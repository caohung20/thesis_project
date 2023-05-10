# Power flow solution by Gauss-Seidel method
import psutil
process = psutil.Process()
print("Memory usage:", process.memory_info().rss)
import numpy as np
import pandas as pd
import math
import cmath
from scipy.sparse import csc_matrix
#base power to convert to per unit
basemva=100 
accuracy=0.001
#maximum of iteration
maxiter=100 
#base voltage to convert to per unit
basevolt=12
j = 1j   
linedata=pd.read_excel("G:\do_an_tot_nghiep\input_test.xlsx",sheet_name = 1,skiprows=1)
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
            Ybus[n,n] += y[k] + Bc[k]
        elif nr[k] == n+1:
            Ybus[n,n] += y[k] + Bc[k]
        else:
            pass

print(Ybus)
#start gauss seidel program
basemva=100 #base power to convert to per unit
accuracy=0.001
accel=1.8
maxiter=100 #maximum of iteration
busdata=pd.read_excel("G:\do_an_tot_nghiep\input_test.xlsx",sheet_name = 0,skiprows=1)
Pgg=np.zeros(nbus,dtype=complex)
Vm=np.zeros(nbus,dtype=complex)
delta=np.zeros(nbus,dtype=complex)
Qgg=np.zeros(nbus,dtype=complex)
yload=np.zeros(nbus,dtype=complex)
deltad=np.zeros(nbus,dtype=complex)
kb=busdata['CODE'].values
Pd=np.zeros(nbus)
Qd=np.zeros(nbus)
Pg=np.zeros(nbus,dtype=complex)
Qg=np.zeros(nbus,dtype=complex)
Qmin=busdata['QgenMin[kvar]'].values
Qmax=busdata['QgenMax[kvar]'].values
Qsh=np.zeros(nbus)
V=np.zeros(nbus,dtype=complex)
P=np.zeros(nbus,dtype=complex)
Q=np.zeros(nbus,dtype=complex)
S=np.zeros(nbus,dtype=complex)
DV=np.zeros(nbus)
DP=np.zeros(nbus,dtype=float)
DQ=np.zeros(nbus,dtype=float)
for k in range(nbus):
    n=int(busdata.iloc[k][0])
    n-=1
    Vm[n]=busdata.iloc[k][8]
    Pd[n]=busdata.iloc[k][4]/(10**3)    #P load
    Qd[n]=busdata.iloc[k][5]/(10**3)    #Q load
    Qsh[n]=busdata.iloc[k][6]/(10**3)  #injected Q from shunt capacitor
    if math.isnan(Qmax[n]):
        Qmax[n]=0
    if math.isnan(Qmin[n]):
        Qmin[n]=0
    #Vm is nan 
    if Vm[n] <= 0 or math.isnan(Vm[n]):
        #estimate flat start
        Vm[n] = 1.0             
        V[n] = 1 + 1j*0
    else:
        delta[n] = np.pi/180*delta[n]
        V[n] = Vm[n]*(np.cos(delta[n]) + 1j*np.sin(delta[n]))
        #express P Bus in per units
        P[n]=(Pg[n]-Pd[n])/basemva 
        #express Q Bus in per units     
        Q[n]=(Qg[n]-Qd[n]+ Qsh[n])/basemva  
        S[n] = P[n] + 1j*Q[n] 
num = 0
AcurBus = 0
converge = 1
Vc = np.zeros((nbus, 1), dtype=complex)   
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
    for n in range(nbus):
        YV = 0 + 1j * 0
        for L in range(nbr):
            if nl[L] == n+1:
                k = nr[L]-1
                YV += Ybus[n, k] * V[k]
            elif nr[L] == n+1:
                k = nl[L]-1
                YV += Ybus[n, k] * V[k]
        Sc = np.conj(V[n]) * (Ybus[n, n] * V[n] + YV)
        Sc = np.conj(Sc)
        DP[n] = P[n] - np.real(Sc)
        DQ[n] = Q[n] - np.imag(Sc)
        if kb[n] == 3:
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
        if kb[n] != 3:
            Vc[n] = (np.conj(S[n]) / np.conj(V[n]) - YV) / Ybus[n, n]
        else:
            pass
        if kb[n] == 1:
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
    if kb[n] == 3:
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
SLT = 0
#######################################################################
#print bus out
print(tech)
print('                      Maximum Power Mismatch = %g \n' % maxerror)
print('                             No. of Iterations = %g \n\n' % iter)

head = [    '    Bus  Voltage  Angle    ------Load------    ---Generation---   Injected',    '    No.  Mag.     Degree     MW       Mvar       MW       Mvar       Mvar ',    '                                                                          ']
print('\n'.join(head))

for n in range(nbus):
    print(' %5g' % n, end='')
    print(' %7.3f' % Vm[n].real, end=' ')
    print(' %8.3f' % deltad[n].real, end=' ')
    print(' %9.3f' % Pd[n], end=' ')
    print(' %9.3f' % Qd[n], end=' ')
    print(' %9.3f' % Pg[n].real, end=' ')
    print(' %9.3f ' % Qg[n].real, end=' ')
    print(' %8.3f' % Qsh[n])

print('      ')
print('    Total              ', end=' ')
print(' %9.3f' % Pdt, end=' ')
print(' %9.3f' % Qdt, end=' ')
print(' %9.3f' % Pgt.real, end=' ')
print(' %9.3f' % Qgt.real, end=' ')
print(' %9.3f\n\n' % Qsht)
#######################################################################################
print('\n')
print('                           Line Flow and Losses \n\n')
print('     --Line--  Power at bus & line flow    --Line loss--\n')
print('     from  to    MW      Mvar     MVA       MW      Mvar\n')

for n in range(nbus):
    busprt = 0
    for L in range(nbr):
        if busprt == 0:
            print('   \n'), 
            print('%6g' % (n+1), end='')
            print('      %9.3f' % (P[n]*basemva).real, end='')
            print('%9.3f' % (Q[n]*basemva).real, end='')
            print('%9.3f\n' % (abs(S[n]*basemva)), end='')

            busprt = 1
        else:
            pass

        if nl[L] == n+1:
            k = nr[L]
            In = (V[n] - V[k-1])*y[L] + Bc[L]*V[n]
            Ik = (V[k-1] - V[n])*y[L] + Bc[L]*V[k-1]
            Snk = V[n]*np.conj(In)*basemva
            Skn = V[k-1]*np.conj(Ik)*basemva
            SL  = Snk + Skn
            SLT = SLT + SL
        elif nr[L] == n+1:
            k = nl[L]
            In = (V[n] - V[k-1])*y[L] + Bc[L]*V[n]
            Ik = (V[k-1] - V[n])*y[L] + Bc[L]*V[k-1]
            Snk = V[n]*np.conj(In)*basemva
            Skn = V[k-1]*np.conj(Ik)*basemva
            SL  = Snk + Skn
            SLT = SLT + SL
        else:
            pass

        if nl[L] == n+1 or nr[L] == n+1:
            print('%12g' % k, end='')
            print('%9.3f' % Snk.real, end='')
            print('%9.3f' % Snk.imag, end='')
            print('%9.3f' % abs(Snk), end='')
            print('%9.3f' % SL.real, end='')

            if nl[L] == n+1:
                print('%9.3f\n' % SL.imag, end='')
            else:
                print('%9.3f\n' % SL.imag, end='')
        else:
            pass

SLT = SLT/2
print('   \n'), 
print('    Total loss                         ', end='')
print('%9.3f' % SLT.real, end='')
print('%9.3f\n' % SLT.imag, end='')
