# Power flow solution by Newton-Raphson method
"""
import psutil
process = psutil.Process()
print("Memory usage:", process.memory_info().rss)
"""
import numpy as np
import pandas as pd
import math
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
#start newton-raphson program
busdata=pd.read_excel("G:\do_an_tot_nghiep\input_test.xlsx",sheet_name = 0,skiprows=1)
ns=0
ng=0
Pgg=np.zeros(nbus,dtype=complex)
Qgg=np.zeros(nbus,dtype=complex)
ngs=np.zeros(nbus)
nss=np.zeros(nbus)
Vm=np.zeros(nbus,dtype=complex)
delta=np.zeros(nbus,dtype=complex)
yload=np.zeros(nbus,dtype=complex)
deltad=np.zeros(nbus,dtype=complex)
kb=busdata['CODE'].values
Pd=np.zeros(nbus,dtype=float)
Qd=np.zeros(nbus,dtype=float)
Pg=np.zeros(nbus)
Qg=np.zeros(nbus)
Qmin=busdata['QgenMin[kvar]'].values
Qmax=busdata['QgenMax[kvar]'].values
Qsh=np.zeros(nbus,dtype=float)
V=np.zeros(nbus,dtype=complex)
P=np.zeros(nbus,dtype=complex)
Q=np.zeros(nbus,dtype=complex)
S=np.zeros(nbus,dtype=complex)
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
for k in range(nbus):
    if kb[k] == 3:
        #count the number of slack bus 
        ns += 1         
    elif kb[k] == 2:
        #count the number of P-V bus
        ng += 1         
    ngs[k] = ng         
    nss[k] = ns
#return Y bus magnitude
Ym=np.abs(Ybus) 
# return phase angle of Ybus        
t = np.angle(Ybus,deg=False)  
#ng bus voltage control => m equations involving delta Q and delta V 
# => corresponding columns of the jacobian matrix are eliminated. 
# => 2*nbus-2-ng components left  
m=2*nbus-ng-2*ns        
maxerror = 1
converge=1
iter = 0
DC = np.zeros(m,dtype=complex)      
J = np.zeros((m,m))
# Start of iterations
while maxerror >= accuracy and iter <= maxiter:
    # Test for max power mismatch
    #initializing Jacobian matrix
    A = np.zeros((m,m),dtype=complex)     
    iter+=1
    for n in range(1, nbus+1):
        #when slack bus appears, it will eliminate slack bus components 
        # in jacobi matrix  
        nn = int(n - nss[n-1]) 
        #(nbus -ns)=> start J3; (n - ngs[n-1] - nss[n-1]) is similar 
        # to the above equation             
        lm = int(nbus + n - ngs[n-1] - nss[n-1] - ns)   
        J11 = 0
        J22 = 0
        J33 = 0
        J44 = 0
        for i in range(1, nbr+1):
            if nl[i-1] == n or nr[i-1] == n:
                if nl[i-1] == n:
                    l = nr[i-1]
                if nr[i-1]==n:
                    l = nl[i-1]
                J11 += Vm[n-1] * Vm[l-1] * Ym[n-1][l-1] * math.sin((t[n-1][l-1] - delta[n-1] + delta[l-1]).real)        #diagonal elements of J1
                J33 += Vm[n-1] * Vm[l-1] * Ym[n-1][l-1] * math.cos((t[n-1][l-1] - delta[n-1] + delta[l-1]).real)        #diagonal elements of J3
                if kb[n-1] != 3:
                    #slackbus doesn't exist in J2 and J4
                    J22 += Vm[l-1] * Ym[n-1][l-1] * math.cos((t[n-1][l-1] - delta[n-1] + delta[l-1]).real)              #diagonal elements of J2
                    J44 += Vm[l-1] * Ym[n-1][l-1] * math.sin((t[n-1][l-1] - delta[n-1] + delta[l-1]).real)              #diagonal elements of J4
                else:
                    pass
                if kb[n-1] != 3 and kb[l-1] != 3:
                    lk = int(nbus + l - ngs[l-1] - nss[l-1] - ns)
                    ll = int(l - nss[l-1])
                    # off diagonal elements of J1
                    A[nn-1][ll-1] = -Vm[n-1] * Vm[l-1] * Ym[n-1][l-1] * math.sin((t[n-1][l-1] - delta[n-1] + delta[l-1]).real)
                
                    if kb[l-1] == 1:
                        # off diagonal elements of J2
                         A[nn-1][lk-1] = Vm[n-1] * Ym[n-1][l-1] * math.cos((t[n-1][l-1] - delta[n-1] + delta[l-1]).real)
                
                    if kb[n-1] == 1:
                        # off diagonal elements of J3
                        A[lm-1][ll-1] = -Vm[n-1] * Vm[l-1] * Ym[n-1][l-1] * math.cos((t[n-1][l-1] - delta[n-1] + delta[l-1]).real)
                
                    if kb[n-1] == 1 and kb[l-1] == 1:
                        # off diagonal elements of J4
                        A[lm-1][lk-1] = -Vm[n-1] * Ym[n-1][l-1] * math.sin((t[n-1][l-1] - delta[n-1] + delta[l-1]).real)
                else:
                    pass
            else:
                pass
        Pk = Vm[n-1]**2 * Ym[n-1][n-1] * math.cos((t[n-1][n-1])) + J33
        Qk = -Vm[n-1]**2 * Ym[n-1][n-1] * math.sin((t[n-1][n-1])) - J11
        if kb[n-1] == 3:
            # swing bus 
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
        if kb[n-1] != 3:
            #diagonal elements of J1
            A[nn-1][nn-1] = J11        
            DC[nn-1] = P[n-1] - Pk
        if kb[n-1] == 1:
            #diagonal elements of J2
            A[nn-1][lm-1] = 2 * Vm[n-1] * Ym[n-1][n-1] * math.cos(t[n-1][n-1]) + J22    
            #diagonal elements of J3
            A[lm-1][nn-1] = J33
            #diagonal elements of J4         
            A[lm-1][lm-1] = -2 * Vm[n-1] * Ym[n-1][n-1] * math.sin(t[n-1][n-1]) - J44   
            DC[lm-1] = Q[n-1] - Qk
    #matrix A left division Matrix DC
    DX = np.linalg.solve(A, DC.T) 
    for n in range(1,nbus+1):
        nn = int(n - nss[n-1])
        lm = int(nbus + n - ngs[n-1] - nss[n-1] - ns)
        if kb[n-1] != 3:
            delta[n-1] +=  DX[nn-1]
        if kb[n-1] == 1:
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
    if kb[n] == 3:
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
SLT = 0
#####################################################################################
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
###################################################################################
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
