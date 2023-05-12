#power summation with first node is a gen in first node
import pandas as pd
import numpy as np
import copy
import cmath
#########################################################
def backwardsweep(dictionary):
#calculate from endnode Snode and sweep until dictionary empty
    for key in dictionary:
        leng=len(dictionary[key])
        for indice in range(leng-1,-1,-1):  
            last_item=dictionary[key][indice]-1
            Snode[last_item]=Sd[last_item]+sum(Sline[last_item][:])+np.conj(B[last_item])*(V[last_item]**2)
            print('S2phay(%i)'%last_item,"=Sd(%i)"%last_item,'Sline[%i][:]'%last_item,)
            if indice!=0:
                second_last=dictionary[key][indice-1]-1
                Sline[second_last][last_item]=Zbus[second_last][last_item]*(Snode[last_item]/V[last_item])**2 + Snode[last_item]
                print('Sphay(%i)='%second_last,'(%i)'%last_item,'Zbus[%i]'%second_last,'[%i]'%last_item,'*(Snode[%i]'%last_item,'/U[%i])**2'%last_item, "+ Snode[%i]"%last_item  )
    return
    

#################################################################
def forwardsweep(my_dict):
#calculate forward sweep
    for key in my_dict:
        for value in my_dict[key]:
            if value==my_dict[key][0]:
                prevalue=value
            else:
                V[value-1]=V[prevalue-1]-np.conj(Sline[prevalue-1][value-1]/V[prevalue-1])*Zbus[prevalue-1][value-1]
                #print("U(%i)" %value,"=U(%i)" %prevalue,"-Sline(%i)" %prevalue,"(%i)" %value)
                prevalue=value
    return
###################################################################################
def accur(U1,U2):
#calculate voltage mismatch and return biggest value
    epsilon=0
    for i in range(len(U1)):
        if i!=0:
            if epsilon < abs(U1[i].real-U2[i].real):
                epsilon=abs(U1[i].real-U2[i].real)
            if epsilon < abs(U1[i].imag-U2[i].imag):
                epsilon=abs(U1[i].imag-U2[i].imag)
    return epsilon
#####################################################################
#calculate power summation necessary parameter
j = 1j   
linedata=pd.read_excel("G:\do_an_tot_nghiep\input3bus.xlsx",sheet_name = 1,skiprows=1)
nl = linedata['FROMBUS'].values
nr = linedata['TOBUS'].values
R = linedata['R(Ohm)'].values
X = linedata['X(Ohm)'].values
Bc = j*linedata['B(microSiemens)'].values
basemva=100
basevolt=12
zbase=((basevolt**2))/basemva
nbr = len(linedata['FROMBUS'])
print(nbr)
nbus = max(max(nl),max(nr))
B=np.zeros(nbus,dtype=complex)
Zbus = np.zeros((int(nbus),int(nbus)), dtype=complex)
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
for i in range(0,nbr):
    #convert to per unit system
    Zbus[nl[i]-1][nr[i]-1]= (R[i] + j*X[i])/zbase
    Zbus[nr[i]-1][nl[i]-1]= (R[i] + j*X[i])/zbase
    # defining branch admittance as complex
    Bc[i]=Bc[i]/2*zbase*10**(-6)
#formation of the B array:
for n in range(nbus):
    for k in range(nbr):
        if nl[k] == n+1:
            B[n] +=Bc[k]
        elif nr[k] == n+1:
            B[n] += Bc[k]
        else:
            pass
dictionary={1:[1, 2, 3, 4, 5, 6],2:[1, 2, 3, 4, 9],3:[1, 2, 7, 8]}
Snode=np.zeros(nbus,dtype=complex)
Sline=np.zeros((nbus,nbus),dtype=complex)
busdata=pd.read_excel("G:\do_an_tot_nghiep\input3bus.xlsx",sheet_name = 0,skiprows=1)
kb=busdata['CODE'].values
Pd=np.zeros(nbus,dtype=float)
Qd=np.zeros(nbus,dtype=float)
Pg=np.zeros(nbus)
Qg=np.zeros(nbus)
Qsh=np.zeros(nbus,dtype=float)
V = np.zeros(nbus,dtype=complex)
P=np.zeros(nbus,dtype=complex)
Q=np.zeros(nbus,dtype=complex)
for k in range(nbus):
    n=int(busdata.iloc[k][0])
    n-=1
    V[n]=busdata.iloc[k][8]
    Pd[n]=busdata.iloc[k][4]/(10**3)    #P load
    Qsh[n]=busdata.iloc[k][6]/(10**3)
    Qd[n]=busdata.iloc[k][5]/(10**3)    #Q load
    Qsh[n]=busdata.iloc[k][6]/(10**3)  #injected Q from shunt capacitor
Sd=np.zeros(nbus,dtype=complex)
basemva=100
j=1j
for i in range(len(Pd)):
    Sd[i]=(Pd[i]+Qd[i]*j-Qsh[i]*j)/(basemva)

#in backward sweep, it will delete global dictionary so it's necessary to save and return dict
maxiter=100
accuracy=0.00001
converge=1
maxerror=10
iter=0
while maxerror >= accuracy and iter <= maxiter:
    iter+=1
    backwardsweep(dictionary)
    #copy value of U bc backward sweep will change value of U
    Vstor=copy.deepcopy(V)
    forwardsweep(dictionary)
    maxerror=accur(V,Vstor)
    Sline=np.zeros((nbus,nbus),dtype=complex)
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
    tech = '                   Power Flow Solution by Power Summation Method'
print(Snode)
Vm=abs(V)
deltad=np.zeros(nbus)
S=Snode
k = 0
for n in range(nbus):
    deltad[n] = cmath.phase(V[n]) * 180 / cmath.pi
    if kb[n] == 3:
        k += 1
        P[n]=S[n].real
        Q[n]=S[n].imag
        Pg[n] = P[n] * basemva + Pd[n]
        Qg[n] = Q[n] * basemva + Qd[n] - Qsh[n]
        
Pgt = sum(Pg)
Qgt = sum(Qg)
Pdt = sum(Pd)
Qdt = sum(Qd)
Qsht = sum(Qsh)
#S loss total
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