import sparse_matrix_lfybus as lfybus
import sparse_matrix_lfnewton as lfnewton
import numpy as np
import math
basemva=lfnewton.basemva
P=lfnewton.P
Q=lfnewton.Q
S=lfnewton.S
nbus=lfybus.nbus
nbr=lfybus.nbr
nl=lfybus.nl
nr=lfybus.nr
a=lfybus.a
V=lfnewton.V
y=lfybus.y
Bc=lfybus.Bc
""""""
for i in range(nbr):
    #return bus value of nr and nl when nr and nl previously were indices
    nr[i]+=1
    nl[i]+=1
SLT = 0
print('\n')
print('                           Line Flow and Losses \n\n')
print('     --Line--  Power at bus & line flow    --Line loss--  Transformer\n')
print('     from  to    MW      Mvar     MVA       MW      Mvar      tap\n')

for n in range(nbus):
    busprt = 0
    for L in range(nbr):
        if busprt == 0:
            print('   \n'), 
            print('%6g' % n, end='')
            print('      %9.3f' % (P[n]*basemva), end='')
            print('%9.3f' % (Q[n]*basemva), end='')
            print('%9.3f\n' % (abs(S[n]*basemva)), end='')

            busprt = 1
        else:
            pass

        if nl[L] == n+1:
            k = nr[L]
            In = (V[n] - a[L]*V[k])*y[L]/a[L]**2 + Bc[L]/a[L]**2*V[n]
            Ik = (V[k] - V[n]/a[L])*y[L] + Bc[L]*V[k]
            Snk = V[n]*np.conj(In)*basemva
            Skn = V[k]*np.conj(Ik)*basemva
            SL  = Snk + Skn
            SLT = SLT + SL
        elif nr[L] == n+1:
            k = nl[L]
            In = (V[n] - V[k]/a[L])*y[L] + Bc[L]*V[n]
            Ik = (V[k] - a[L]*V[n])*y[L]/a[L]**2 + Bc[L]/a[L]**2*V[k]
            Snk = V[n]*np.conj(In)*basemva
            Skn = V[k]*np.conj(Ik)*basemva
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

            if nl[L] == n+1 and a[L] != 1:
                print('%9.3f' % SL.imag, end='')
                print('%9.3f\n' % a[L], end='')
            else:
                print('%9.3f\n' % SL.imag, end='')
        else:
            pass

SLT = SLT/2
print('   \n'), 
print('    Total loss                         ', end='')
print('%9.3f' % SLT.real, end='')
print('%9.3f\n' % SLT.imag, end='')
