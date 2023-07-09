__author__    = "Dr Pham Quang Phuong","Cao Anh Quoc Hung"
__copyright__ = "Copyright 2023"
__license__   = "All rights reserved"
__email__     = "phuong.phamquang@hust.edu.vn"
__status__    = "Released"
__version__   = "1.1.5"
"""
about: ....
"""
from new_kernel import add2CSV,PowerFlow,Parameter,Configuration
import os,sys,time
import argparse
import random
import math
from pyswarms.utils.plotters import plot_cost_history
import pygad
import matplotlib.pyplot as plt

PARSER_INPUTS = argparse.ArgumentParser(epilog= "")
PARSER_INPUTS.usage = 'Distribution network analysis Tools'
PARSER_INPUTS.add_argument('-fi' , help = '*(str) Input file path .xlsx' , default = '',type=str,metavar='')
PARSER_INPUTS.add_argument('-fo' , help = ' (str) Output file path .csv' , default = '',type=str,metavar='')
ARGVS = PARSER_INPUTS.parse_known_args()[0]
#
def func1(x,param):
    if type(x).__name__=='ndarray':
        x = x.tolist()[0]
    config = Configuration(param=param,varFlag=x)
    pf = PowerFlow(config)
    return pf.run1Config_WithObjective()['Objective']
#
def monteCarlo(nIter,lineOff0=[],shuntOff0=[],dgOff0=[]):
    nIter = int(nIter)
    print('Running monteCarlo nIter=%i'%nIter)
    param = Parameter(ARGVS.fi)
    rs = [[],[time.ctime(),'MonteCarlo init_pos lineOff0=%s shuntOff0=%s dgOff0=%s'%(str(lineOff0),str(shuntOff0),str(dgOff0))],['iter','Objective','DeltaA','RateMax[%]','Umax[pu]','Umin[pu]','LineOff','ShuntOff','DgOff'] ]
    add2CSV(ARGVS.fo,rs,',')
    r0 = {'Objective':math.inf,'LineOff':lineOff0,'ShuntOff':shuntOff0,'DgOff':dgOff0}
    if lineOff0 or shuntOff0 or dgOff0:
        config = Configuration(param=param,lineOff=lineOff0,shuntOff=shuntOff0,dgOff=dgOff0)
        pf = PowerFlow(config)
        r1 = pf.run1Config_WithObjective()
        rs = ['init','%.5f'%r1['Objective'],'%.5f'%r1['DeltaA'].real,'%.3f'%r1['RateMax[%]'],'%.3f'%r1['Umax[pu]'],'%.3f'%r1['Umin[pu]'],str(r1['LineOff']),str(r1['ShuntOff']),str(r1['DgOff'])]
        add2CSV(ARGVS.fo,[rs],',')
        r0 = r1
    #
    for i in range(nIter):
        x = [random.randint(0,1) for _ in range(param.nVar)]
        config = Configuration(param=param,varFlag=x)
        pf = PowerFlow(config)
        r1 = pf.run1Config_WithObjective()
        if r1['FLAG']=='CONVERGENCE':
            if r1['Objective']<r0['Objective']:
                rs = [str(i),'%.5f'%r1['Objective'],'%.5f'%r1['DeltaA'].real,'%.3f'%r1['RateMax[%]'],'%.3f'%r1['Umax[pu]'],'%.3f'%r1['Umin[pu]'],str(r1['LineOff']),str(r1['ShuntOff'])]
                add2CSV(ARGVS.fo,[rs],',')
                r0 = r1
    #
    s1 = ['MonteCarlo','nIter',nIter,' time[s]',str(r0['LineOff']),str(r0['ShuntOff']),str(r0['DgOff'])]
    add2CSV(ARGVS.fo,[s1],',')
    print('\nOutFile: '+os.path.abspath(ARGVS.fo))
    print(s1)
##    print(pf.nn)
#
def pso(nIter,lineOff0=[],shuntOff0=[],dgOff0=[]):
    nIter = int(nIter)
    print('Running PSO nIter=%i'%nIter)
    import pyswarms as ps
    import numpy as np
    #
    param = Parameter(ARGVS.fi)
    #

    options = {'c1': 0.5, 'c2': 0.3, 'w': 0.9,'k':1,'p':1}
    #
    pos0 = None
    if lineOff0 or shuntOff0 or dgOff0 :
        pos0 = np.array([config.getVarFlag()])
    #
    op = ps.discrete.binary.BinaryPSO(n_particles=1,dimensions=param.nVar,options=options,init_pos=pos0)
    cost_history = op.cost_history
    #
    cost, pos = op.optimize(func1, iters=nIter,param=param)
    #
    plot_cost_history(cost_history)
    plt.show()
    #
    config = Configuration(param=param,varFlag=pos)
    pf = PowerFlow(config)
    r1 = pf.run1Config_WithObjective()
    rs = [[],[time.ctime(),'PSO init_pos (lineOff)',str(lineOff0)],['Objective','DeltaA','RateMax[%]','Umax[pu]','Umin[pu]','LineOff','ShuntOff'] ]
    rs.append( ['%.5f'%r1['Objective'],'%.5f'%r1['DeltaA'].real,'%.3f'%r1['RateMax[%]'],'%.3f'%r1['Umax[pu]'],'%.3f'%r1['Umin[pu]'],str(r1['LineOff']),str(r1['ShuntOff'])])
    add2CSV(ARGVS.fo,rs,',')
    #
    s1 = ['PSO','nIter',nIter,str(r1['LineOff']),str(r1['ShuntOff']),str(r1['DgOff'])]
    add2CSV(ARGVS.fo,[s1],',')
    print('\nOutFile: '+os.path.abspath(ARGVS.fo))
    print(s1)
#
# "D:\\00_BK\\optimpy\\Python310_64\\python.exe" optim2.py
"""
def fitness_func(reconfiguration_ga,solution,solution_idx):

    r1 = pf.run1Config_WithObjective(varFlag=solution)
    #
    fitness = -r1['Objective']
    return fitness
"""
#
class GA():
    def __init__(self):
        self.param = Parameter(ARGVS.fi)
    def fitness_func(self,reconfiguration,solution,solution_idx):
        config = Configuration(param=self.param,varFlag=solution)
        pf = PowerFlow(config=config)
        r1 = pf.run1Config_WithObjective()
        fitness = -r1['Objective']
        return fitness
    def reconfiguration_ga(self,nIter,lineOff0=[],shuntOff0=[],dgOff0=[]):
        nIter = int(nIter)
        print('Running Genetic Algorithm nIter=%i'%nIter)
        param = Parameter(ARGVS.fi)
        rs = [[],[time.ctime(),'GeneticAlgorithm init_pos lineOff0=%s shuntOff0=%s dgOff0=%s'%(str(lineOff0),str(shuntOff0),str(dgOff0))],['iter','Objective','DeltaA','RateMax[%]','Umax[pu]','Umin[pu]','LineOff','ShuntOff','DgOff'] ]
        add2CSV(ARGVS.fo,rs,',')
        r0 = {'Objective':math.inf,'LineOff':lineOff0,'ShuntOff':shuntOff0,'DgOff':dgOff0}
        if lineOff0 or shuntOff0 or dgOff0:
            config = Configuration(param=param,lineOff=lineOff0,shuntOff=shuntOff0,dgOff=dgOff0)
            pf = PowerFlow(config)
            r1 = pf.run1Config_WithObjective()
            rs = ['init','%.5f'%r1['Objective'],'%.5f'%r1['DeltaA'].real,'%.3f'%r1['RateMax[%]'],'%.3f'%r1['Umax[pu]'],'%.3f'%r1['Umin[pu]'],str(r1['LineOff']),str(r1['ShuntOff']),str(r1['DgOff'])]
            add2CSV(ARGVS.fo,[rs],',')
            r0 = r1
        gene_space = [[0,1] for _ in range(param.nVar)]
        reconfiguration = pygad.GA(num_generations=nIter,
                        num_parents_mating=10,
                        fitness_func=self.fitness_func,
                        sol_per_pop=100,mutation_num_genes=5,
                        num_genes=param.nVar,gene_space=gene_space,
                        gene_type=int,keep_elitism=1,
                        crossover_probability=0.6,
                        mutation_probability=0.008,)
        reconfiguration.run()

        # Retrieve the best solution and its fitness value
        best_solution = reconfiguration.best_solution()

        # Print the best solution and its fitness value
        print("Best Solution: ", best_solution)
        #
        config = Configuration(param=param,varFlag=best_solution[0])
        lineOff = config.lineOff
        shuntOff = config.shuntOff
        dgOff = config.dgOff
        s1 = ['GA','nIter',nIter,str(lineOff),str(shuntOff),str(dgOff)]

        print(s1)
        #
        reconfiguration.plot_fitness()
        #
        return 
if __name__ == '__main__':
    """
    #ARGVS.fi = 'Inputs33bc_shunt100.xlsx'
    #ARGVS.fi = 'Inputs12_2-low2.xlsx'
    
    ARGVS.fi = 'Inputs190month.xlsx'
    #ARGVS.fi = 'Inputs102.xlsx'

    #
    
    #
##    lineOff0 = [7, 8, 9, 11, 15] # init, =[] if no init
    lineOff0 = [] # no init
    shuntOff0 = []
##    lineOff0 = [66, 103, 110, 169, 191]
##    shuntOff0 = [47, 66, 80, 130]
    #
    #
    #
    pso(1e6,lineOff0,shuntOff0,dgOff0)
    """
    ARGVS.fi = 'Inputs12_2-low.xlsx'
    lineOff0 = [] # no init
    shuntOff0 = []
    dgOff0 = []
    ARGVS.fo = 'res\\resOptim12.csv'
    pso(1e4,lineOff0,shuntOff0,dgOff0)
    #GA().reconfiguration_ga(1000,lineOff0,shuntOff0,dgOff0)
    
    
    
