import time
time_start=time.time()
import pandas as pd
import numpy as np
#count branch loop(countBL)
countBL=0
def searchdel(arr1,arr2,k):
    return
"""
loop=[1,2,3,4,5,6,8]
def breakloop(loop,source):
    leng=len(loop)
    if source in loop:
        for i in range(1,leng+1):
            config1=[]
            config2=[]
            config1.append(source)
            config2.append(source)
            for k in range(1,i):
                config1.append(loop[k])
            for k in range (leng-1,i-1,-1):
                config2.append(loop[k])
            if len(config1)==1:
                config1=[]
            else:
                print('config1:',config1) 
            if len(config2)==1:
                config2=[]
            else:   
                print('config2:',config2)
                
    else:
        for i in range(0,leng):
            config1=[]
            config2=[]
            for k in range(0,i):
                config1.append(loop[k])
            for k in range (leng-1,i-1,-1):
                config2.append(loop[k])
            if len(config1)!=0:
                print('config1:',config1) 
            if len(config2)!=0:
                print('config2:',config2)    
    return
breakloop(loop,7)
"""
def searchList(arr1, arr2, root):
#find the next branch of the root, it will modify 
#the arr1 and arr2 to delete where it has searched
# if there is a branch loop, it will return error and 
# value of the loop
    k=0
    while k<len(arr1):
        if k==len(arr1):
            break
        if root== arr1[k] or root == arr2[k]:
            if root== arr1[k]:
                root = arr2[k]
                arr1 = np.delete(arr1, k)
                arr2 = np.delete(arr2, k)
                if root not in visited:
                    visited.append(root)
                else:
                    branchL=set()
                    iError=visited.index(root)
                    #while specify iError in visited, from visited[iError] to len(visited) is a loop
                    for i in range(iError,len(visited)):
                        branchL.add(visited[i])
                    branchLoop[len(branchLoop)]=branchL
                    blErr=True
                    break
            #check if k == leng arr1 because when it delete in the above statement, 
            # it indice may surpass len (arr)
            if k==len(arr1):
                break
            if root == arr2[k]:
                root = arr1[k]
                arr1 = np.delete(arr1, k)
                arr2 = np.delete(arr2, k)
                if root not in visited:
                    visited.append(root)
                else:
                    branchL=set()
                    iError=visited.index(root)
                    #while specify iError in visited, from visited[iError] to len(visited) is a loop
                    for i in range(iError,len(visited)):
                        branchL.add(visited[i])
                    branchLoop[len(branchLoop)]=branchL
                    blErr=True
                    break
        else:
            k+=1
        
    return arr1,arr2,root
def checkisland(island,my_dict):
    for key,value in enumerate(my_dict):
        if value in island:
            island.remove(value)
    if not island:
        return
    else:
        print(island,'is/are island(s)')



def dfs(arr1, arr2, root, check,count):
        if check == False:
            #variable count to eliminate cases when it retreats branch
            if count==0:
                print(visited)
                for i,value in enumerate(visited):
                    if i!=0:
                        for k in range(len(gen)):
                            if gen[k]==value:
                                genLoop.add(gen[k])
                                genLoop.add(visited[0])
            visited.remove(visited[len(visited)-1])
            count+=1
            if not visited:
                return
            root = visited[len(visited)-1]
            dfs(arr1, arr2, root, True,count)
        else:
            if root in arr1 or root in arr2:
                (arr1,arr2,root)=searchList(arr1, arr2, root)   
                count=0
                dfs(arr1, arr2, root, True,count)
            else:
                dfs(arr1, arr2, root, False,count)
busdata=pd.read_excel("G:\do_an_tot_nghiep\input_test.xlsx",sheet_name = 0,skiprows=1)
linedata=pd.read_excel("G:\do_an_tot_nghiep\input_test.xlsx",sheet_name = 1,skiprows=1)
nl = linedata['FROMBUS'].values
nr = linedata['TOBUS'].values
flag = linedata['FLAG'].values
Vsch = busdata['Vscheduled[pu]'].values
nbus = busdata['NO'].values
#Check active line
i=0
while i<len(flag):
    if flag[i]==0:
        nr = np.delete(nr, i)
        nl = np.delete(nl, i)
        flag = np.delete(flag, i)
    else:
        i+=1
    if i==len(flag):
        break
#specify generator
gen=[]
for i in range(len(nbus)):
    if Vsch[i]>0:
        gen.append(nbus[i])

island =  {n+1 for n in range(len(nbus))}

for k,value in enumerate(gen):
    visited=[]
    visited.append(value)
    genLoop=set()
    branchLoop=dict()
    dfs(nl,nr,value,True,0)
    if bool(genLoop) or bool(branchLoop):
        if bool(genLoop) and bool(branchLoop)!=0:
            print(' exist generation loop in bus: ',genLoop,'and branchLoop among bus',branchLoop)
        elif bool(genLoop):
            print(' exist generation loop in bus: ',genLoop)
        else:
            print('exist branch loop in bus',branchLoop)
        break
time_end=time.time()
if bool(island):
    print(island,'are/is island bus(es)')
print(time_end-time_start)
