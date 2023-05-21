import pandas as pd
import numpy as np
import math
import copy
#class line:

def returnactiveline(arr1,arr2,func_flag):
    #eliminate off line to 
    flagoff=0
    while flagoff in func_flag:
        for i in range(len(func_flag)):
            if i>=len(func_flag):
                break
            if func_flag[i]==0:
                arr1=np.delete(arr1,i)
                arr2=np.delete(arr2,i)
                func_flag=np.delete(func_flag,i)
    return arr1, arr2

def stateofbus(num):
    #convert num from decimal to binary and assign these value to flag in order to 
    # return all possible configuration of the radical distribution system.
    # This function also return count to assign 0 to the rest of changeableline 
    global count
    if num>=1:
        flag[changeableline[count]]=num%2
        stateofbus(num//2)
    count+=1
    return

def searchList(arr1, arr2, root):
#find the next branch of the root, it will modify 
#the arr1 and arr2 to delete where it has searched
# if there is a branch loop, it will return error and 
# value of the loop
    global blErr
    k=0
    while k<len(arr1) and not blErr:
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
                    blErr=True
                    break
        else:
            k+=1
        
    return arr1,arr2,root

def dfs(arr1, arr2, root, check, count):
    #depth first search
    if check == False:
        #variable count to eliminate cases when it retreats branch
        if count==0:
            #print(visited)
            iter=len(config)+1
            config[iter]=visited.copy()
        visited.remove(visited[len(visited)-1])
        count+=1
        if not visited:
            return
        root = visited[len(visited)-1]
        dfs(arr1, arr2, root, True,count)
    else:
        if root in arr1 or root in arr2:
            (arr1,arr2,root)=searchList(arr1, arr2, root)
            if blErr is True:
                return    
            count=0
            dfs(arr1, arr2, root, True,count)
        else:
            dfs(arr1, arr2, root, False,count)
        
def checkisland(island,my_dict):
    check=False
    for key in my_dict:
        for value in my_dict[key]:
            if value in island:
                island.discard(value)
    if not island:
        check=False
    else:
        check=True
    return check

def checkgenloop(func_gen,mydic):
    iEgen=False
    for key in mydic:
        nogen=0
        for i in range(len(func_gen)):
            if func_gen[i] in mydic[key]:
                nogen+=1
            if nogen>1:
                iEgen=True
                break

    return iEgen
        
if __name__ == "__main__":
    busdata=pd.read_excel("G:\do_an_tot_nghiep\Inputs12_2.xlsx",sheet_name = 1,skiprows=1)
    linedata=pd.read_excel("G:\do_an_tot_nghiep\Inputs12_2.xlsx",sheet_name = 2,skiprows=1)
    nl = linedata['FROMBUS'].values
    nbr=len(nl)
    nr = linedata['TOBUS'].values
    nbus=(max(nl))
    id=np.zeros((nbus,nbus),dtype=int)
    for i in range (nbr):
        id[nl[i]-1][nr[i]-1]=i+1
        id[nr[i]-1][nl[i]-1]=i+1
    flag = linedata['FLAG'].values
    flag3= linedata['FLAG3'].values
    Qmin= busdata['QgenMin[kvar]'].values
    #specify if a line can switch or not
    changeableline=[]
    for i in range(len(flag)):
        if flag3[i]!=0:
            changeableline.append(i)
        else:
            flag[i]=1
    #declare and specify the number of generators in bus
    gen=[]
    for i in range(nbus):
        if Qmin[i]!=0 and not math.isnan(Qmin[i]):
            gen.append(i+1)
    
    iter=0  
    for deci in range(2**(len(changeableline))):
        #initialized island
        island =  {n+1 for n in range(nbus)} 
        count=0
        stateofbus(deci)
        if count<=len(changeableline):
            for i in range(count,len(changeableline)):
                flag[changeableline[i]]=0
        (nl_new,nr_new)=returnactiveline(nl,nr,flag)
        config=dict()
        blErr=False
        for k,value in enumerate(gen):
            visited=[]
            visited.append(value)
            dfs(nl_new,nr_new,value,True,0)
            if blErr is True:
                break
        #print(config)
        if blErr is True:
            pass
        else:
            checkIsland=checkisland(island,config)
            if checkIsland is True:
                pass
            else:
                checkloop=checkgenloop(gen,config)
                if checkloop is True:
                    pass
                else: 
                    print('this config can calculate:',config)

