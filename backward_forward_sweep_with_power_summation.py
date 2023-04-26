import pandas as pd
import numpy as np
import lfybus
import lfybus
import numpy as np
import pandas as pd
nbus=lfybus.nbus
Ybus=lfybus.Ybus
nbr=lfybus.nbr
nl=lfybus.nl
nr=lfybus.nr
dictionary={1:[1,2,7,8],2:[1,2,3,4,5,6],3:[1,2,3,4,5,11],4:[1,2,3,4,5,10,12]}
layer={6: [12], 5: [6, 11, 10], 4: [5, 5, 5], 3: [8, 4, 4, 4], 2: [7, 3, 3, 3], 1: [2, 2, 2, 2], 0: [1, 1, 1, 1]}
endnode={12,11,6,8}
Snode=np.zeros(nbus,dtype=complex)
Sline=np.zeros((nbus,nbus),dtype=complex)
busdata=pd.read_excel("G:\do_an_tot_nghiep\Inputs12_2.xlsx",sheet_name = 1,skiprows=1)
Pd=busdata['PLOAD[kw]'].values
Qd=busdata['QLOAD[kvar]'].values
U = busdata['Vscheduled[pu]'].values
for i in range(len(U)):
    U[i]=1
Sd=np.zeros(nbus,dtype=complex)
basemva=100
basevolt=11
for i in range(len(Pd)):
    Sd[i]=(Pd[i]+Qd[i])/basemva
def backwardsweep(endnode,dictionary):
#calculate from endnode Snode and sweep until dictionary empty
    if  not endnode:
        return 
    else:
        for i in endnode: 
            for key in dictionary:
                leng=len(dictionary[key])
                #check if dictionary[key] contains i
                if i==dictionary[key][leng-1]:
                    Snode[i-1]=Sd[i-1]+sum(Sline[i-1][:])
                    Sline[dictionary[key][leng-2]-1][i-1]=np.conj(Ybus[dictionary[key][leng-2]-1][i-1])*(Snode[i-1]/U[i-1])**2
                    if leng!=0:   
                        dictionary[key].remove(i)
        endnode.clear()
        for key in dictionary:
            leng=len(dictionary[key])
            if leng!=0:
                #add new endnode to set
                endnode.add(dictionary[key][leng-1])
        backwardsweep(endnode,dictionary)
backwardsweep(endnode,dictionary)
def forwardsweep(my_dict):
    for key in my_dict:
        for value in my_dict[key]:
            if value==my_dict[key][0]:
                prevalue=value
            else:
                U[value-1]=U[prevalue-1]-np.conj(Snode[prevalue-1]/U[prevalue-1])
                prevalue=value
while maxerror >= accuracy and iter <= maxiter:
    
print(Snode)             
        
###############################################

dictionary={1:[1,2,7,8],2:[1,2,3,4,5,6],3:[1,2,3,4,5,11],4:[1,2,3,4,5,10,12]}
layer={1:[1,2],2:[7,3],3:[8,4],4:[5],5:[6,11,10],6:[12]}
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

def backwardSweep(my_dic,my_layer):
    
#base power to convert to per unit
basemva=100 
accuracy=0.001
#maximum of iteration
maxiter=100 
busdata=pd.read_excel("G:\do an tot nghiep\myExample.xlsx",sheet_name = 0)
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