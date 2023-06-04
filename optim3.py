__author__    = "Dr. Pham Quang Phuong"
__copyright__ = "Copyright 2023"
__license__   = "All rights reserved"
__email__     = "phuong.phamquang@hust.edu.vn"
__status__    = "Released"
__version__   = "1.1.4"
"""
about: ....
"""
from KERNELPF3 import add2CSV,POWERFLOW
import os,sys,time
import argparse
import random
import math
from pyswarms.utils.plotters import plot_cost_history
import matplotlib.pyplot as plt

PARSER_INPUTS = argparse.ArgumentParser(epilog= "")
PARSER_INPUTS.usage = 'Distribution network analysis Tools'
PARSER_INPUTS.add_argument('-fi' , help = '*(str) Input file path .xlsx' , default = '',type=str,metavar='')
PARSER_INPUTS.add_argument('-fo' , help = ' (str) Output file path .csv' , default = '',type=str,metavar='')
ARGVS = PARSER_INPUTS.parse_known_args()[0]
#
def func1(x,pf):
    if type(x).__name__=='ndarray':
        x = x.tolist()[0]
    return pf.run1Config_WithObjective(varFlag=x)['Objective']
#
def monteCarlo(nIter,lineOff0=[],shuntOff0=[]):
    nIter = int(nIter)
    print('Running monteCarlo nIter=%i'%nIter)
    pf = POWERFLOW(ARGVS.fi)
    rs = [[],[time.ctime(),'MonteCarlo init_pos lineOff0=%s shuntOff0=%s'%(str(lineOff0),str(shuntOff0))],['iter','Objective','DeltaA','RateMax[%]','Umax[pu]','Umin[pu]','LineOff','ShuntOff'] ]
    add2CSV(ARGVS.fo,rs,',')
    r0 = {'Objective':math.inf,'LineOff':lineOff0,'ShuntOff':shuntOff0}
    if lineOff0 or shuntOff0:
        varFlag = pf.getVarFlag(lineOff0,shuntOff0)
        r1 = pf.run1Config_WithObjective(varFlag=varFlag)
        rs = ['init','%.5f'%r1['Objective'],'%.5f'%r1['DeltaA'].real,'%.3f'%r1['RateMax[%]'],'%.3f'%r1['Umax[pu]'],'%.3f'%r1['Umin[pu]'],str(r1['LineOff']),str(r1['ShuntOff'])]
        add2CSV(ARGVS.fo,[rs],',')
        r0 = r1
    #
    for i in range(nIter):
        x = [random.randint(0,1) for _ in range(pf.nVar)]
        r1 = pf.run1Config_WithObjective(varFlag=x)
        if r1['FLAG']=='CONVERGENCE':
            if r1['Objective']<r0['Objective']:
                rs = [str(i),'%.5f'%r1['Objective'],'%.5f'%r1['DeltaA'].real,'%.3f'%r1['RateMax[%]'],'%.3f'%r1['Umax[pu]'],'%.3f'%r1['Umin[pu]'],str(r1['LineOff']),str(r1['ShuntOff'])]
                add2CSV(ARGVS.fo,[rs],',')
                r0 = r1
    #
    s1 = ['MonteCarlo','nIter',nIter,' time[s]','%.2f'%(time.time()-pf.t0),str(r0['LineOff']),str(r0['ShuntOff'])]
    add2CSV(ARGVS.fo,[s1],',')
    print('\nOutFile: '+os.path.abspath(ARGVS.fo))
    print(s1)
    print('%.5f'%pf.tcheck)
##    print(pf.nn)
#
def pso(nIter,lineOff0=[],shuntOff0=[]):
    nIter = int(nIter)
    print('Running PSO nIter=%i'%nIter)
    import pyswarms as ps
    import numpy as np
    #
    pf = POWERFLOW(ARGVS.fi)
    #
    options = {'c1': 0.5, 'c2': 0.3, 'w': 0.9,'k':1,'p':1}
    #
    pos0 = None
    if lineOff0 or shuntOff0 :
        pos0 = np.array([ pf.getVarFlag(lineOff0,shuntOff0) ])
    #
    op = ps.discrete.binary.BinaryPSO(n_particles=1,dimensions=pf.nVar,options=options,init_pos=pos0)
    cost_history = op.cost_history
    #
    cost, pos = op.optimize(func1, iters=nIter,pf=pf)
    #
    plot_cost_history(cost_history)
    plt.show()
    #
    r1 = pf.run1Config_WithObjective(varFlag=pos)
    rs = [[],[time.ctime(),'PSO init_pos (lineOff)',str(lineOff0)],['Objective','DeltaA','RateMax[%]','Umax[pu]','Umin[pu]','LineOff','ShuntOff'] ]
    rs.append( ['%.5f'%r1['Objective'],'%.5f'%r1['DeltaA'].real,'%.3f'%r1['RateMax[%]'],'%.3f'%r1['Umax[pu]'],'%.3f'%r1['Umin[pu]'],str(r1['LineOff']),str(r1['ShuntOff'])])
    add2CSV(ARGVS.fo,rs,',')
    #
    s1 = ['PSO','nIter',nIter,' time[s]','%.2f'%(time.time()-pf.t0),str(r1['LineOff']),str(r1['ShuntOff'])]
    add2CSV(ARGVS.fo,[s1],',')
    print('\nOutFile: '+os.path.abspath(ARGVS.fo))
    print(s1)
#
# "D:\\00_BK\\optimpy\\Python310_64\\python.exe" optim2.py
if __name__ == '__main__':
    #ARGVS.fi = 'Inputs33bus.xlsx'
    ARGVS.fi = 'Inputs12.xlsx'
    ARGVS.fi = 'Inputs12_2.xlsx'
    #ARGVS.fi = 'Inputs190shunt.xlsx'
    #
    ARGVS.fo = 'res\\resOptim12.csv'
    #
##    lineOff0 = [7, 8, 9, 11, 15] # init, =[] if no init
    lineOff0 = [] # no init
    shuntOff0 = []
##    lineOff0 = [66, 103, 110, 169, 191]
##    shuntOff0 = [47, 66, 80, 130]
    #
    monteCarlo(2e4,lineOff0,shuntOff0)
    #
    #pso(1e4,lineOff0,shuntOff0)


