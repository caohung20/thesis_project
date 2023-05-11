#power summation with first node is a gen in first node
import pandas as pd
import numpy as np
import numpy as np
import copy
def firstEndnode(my_dict):
    #find first endnode
    firstE=set()
    for key in my_dict:
        leng=len(my_dict[key])
        firstE.add(my_dict[key][leng-1])
    return firstE
#########################################################
def backwardsweep(endnode,dictionary):
#calculate from endnode Snode and sweep until dictionary empty
    if  not endnode:
        return 
    else:
        for i in endnode: 
            for key in dictionary:
                leng=len(dictionary[key])
                if leng!=0:
                #check if leng !=0 continue process
                    if i==dictionary[key][leng-1]:
                        #check if dictionary[key] contains i
                        Snode[i-1]=Sd[i-1]+sum(Sline[i-1][:])+np.conj(B[i-1]*j)*(U[i-1]**2)
                        if leng>1:
                            Sline[dictionary[key][leng-2]-1][i-1]=Zbus[dictionary[key][leng-2]-1][i-1]*(Snode[i-1]/U[i-1])**2 + Snode[i-1]
                        if leng!=0:   
                            dictionary[key].remove(i)
        endnode.clear()
        for key in dictionary:
            leng=len(dictionary[key])
            if leng!=0:
                #add new endnode to set
                endnode.add(dictionary[key][leng-1])
        backwardsweep(endnode,dictionary)
#################################################################
def forwardsweep(my_dict):
#calculate forward sweep
    for key in my_dict:
        for value in my_dict[key]:
            if value==my_dict[key][0]:
                prevalue=value
            else:
                U[value-1]=U[prevalue-1]-np.conj(Sline[prevalue-1][value-1]/U[prevalue-1])*Zbus[prevalue-1][value-1]
                prevalue=value
    return
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
j = 1j   
linedata=pd.read_excel("G:\do_an_tot_nghiep\input_test.xlsx",sheet_name = 1,skiprows=1)
nl = linedata['FROMBUS'].values
nr = linedata['TOBUS'].values
R = linedata['R(Ohm)'].values
X = linedata['X(Ohm)'].values
Bc = j*linedata['B(microSiemens)'].values
basemva=100
basevolt=110
zbase=basevolt**2/basemva
nbr = len(linedata['FROMBUS'])
print(nbr)
nbus = max(max(nl), max(nr))
B=np.zeros(nbus,dtype=complex)
Zbus = np.zeros((int(nbus),int(nbus)), dtype=complex) 
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
endnode=firstEndnode(dictionary)
firstEndnode_store=copy.deepcopy(endnode)
Snode=np.zeros(nbus,dtype=complex)
Sline=np.zeros((nbus,nbus),dtype=complex)
busdata=pd.read_excel("G:\do_an_tot_nghiep\input_test.xlsx",sheet_name = 0,skiprows=1)
Pd=busdata['PLOAD[kw]'].values
Qd=busdata['QLOAD[kvar]'].values
U = np.zeros(nbus,dtype=complex)
for k in range(nbus):
    n=int(busdata.iloc[k][0])
    n-=1
    U[n]=busdata.iloc[k][8]

Sd=np.zeros(nbus,dtype=complex)
basemva=100
j=1j
for i in range(len(Pd)):
    Sd[i]=(Pd[i]+Qd[i]*j)/(basemva*10**(3))

#in backward sweep, it will delete global dictionary so it's necessary to save and return dict
stor_dict=copy.deepcopy(dictionary)
maxiter=100
accuracy=0.0001
converge=1
maxerror=10
iter=0
while maxerror >= accuracy and iter <= maxiter:
    iter+=1
    backwardsweep(endnode,dictionary)
    #copy value of U bc backward sweep will change value of U
    Ustor=copy.deepcopy(U)
    dictionary=copy.deepcopy(stor_dict)
    endnode=copy.deepcopy(firstEndnode_store)
    forwardsweep(dictionary)
    maxerror=accur(U,Ustor)
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
print(abs(U))
###############################################

def createLayer(my_dic):
    #create a dictionary of layer from the dfs search dictionary
    max_len = max(len(s) for s in my_dic.values())
    layer=dict()
    for i in range(max_len-1,-1,-1):
        layer[i]=list()
        for k in my_dic.keys():
            if len(my_dic[k])==(i+1):
                #add last item of dictionary[k] to layer[i]
                last_item = my_dic[k].pop()
                layer[i].append(last_item)
    return layer

    