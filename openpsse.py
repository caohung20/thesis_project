import os
import sys

sys_paths = [
        r'C:\Program Files\PTI\PSSE35\35.3\PSSPY39'
        ]

env_paths = [
        r'C:\Program Files\PTI\PSSE35\35.3\PSSBIN', 
        r'C:\Program Files\PTI\PSSE35\35.3\PSSLIB', 
        ]
for path in sys_paths:
    sys.path.append(path)

for path in env_paths:
    os.environ['PATH'] = os.environ['PATH'] + ';' +  path

import psse35
import psspy

# Define the name and location of the new case file
case_name = 'new_case.sav'
case_path = os.path.join(os.getcwd(), case_name)
psspy.psseinit(80000)
# Create a new PSSE case
ierr = psspy.newcas(basemva=100)
#Integer IBUS bus number (input; no default allowed).
"""
Integer INTGAR (4) array of 4 elements specifying (input).
INTGAR(1) IDE, bus type code (1 by default)
INTGAR(2) AREA, area number (1 by default)
INTGAR(3) ZONE, zone number (1 by default)
INTGAR(4) OWNER, owner number (1 by default)
Real REALAR (3) array of 3 elements specifying (input).
REALAR(1) BASKV, bus base voltage in kV (0.0 by default)
REALAR(2) VM, bus voltage magnitude in pu (1.0 by default)
REALAR(3) VA, bus voltage phase angle (0.0 by default)"""
ierr = psspy.bus_data_2(1, intgar=[3, 1, 1, 1,], realar=[100.0, 0.0, 0.0],name='Bus 1')
ierr = psspy.bus_data_2(2, intgar=[2, 1, 1, 1,], realar=[100.0, 0.0, 0.0],name='Bus 2')
ierr = psspy.bus_data_2(3, intgar=[2, 1, 1, 1,], realar=[100.0, 0.0, 0.0],name='Bus 3')
"""
Integer IBUS bus number of from bus (input; no default allowed).
Integer JBUS bus number of to bus (input; no default allowed).
Character*2 CKT circuit identifier (input; '1' by default).
Integer INTGAR (6) array of 6 elements specifying (input).
INTGAR(1) ST, branch status (alias is ST) (1 by default)
INTGAR(2) METBUS, metered end bus number (IBUS or
JBUS) (alias is METBUS) (IBUS by default)
INTGAR(3) O1, first owner number (alias is O1) (owner
of bus IBUS by default)
INTGAR(4) O2, second owner number (alias is O2) (0 by
default)
INTGAR(5) O3, third owner number (alias is O3) (0 by
default)
INTGAR(6) O4, fourth owner number (alias is O4) (0 by
default)
Real REALAR (15) array of 15 elements specifying (input).
REALAR(1) R, nominal branch resistance (alias is R) (0.0
by default)
REALAR(2) X, nominal branch reactance (alias is X)
(THRSHZ by default; 0.0001 if THRSHZ = 0.0)
REALAR(3) B, total line charging (alias is B) (0.0 by default)
REALAR(4) GI, real line shunt at bus IBUS end (alias is GI)
(0.0 by default)
REALAR(5) BI, reactive line shunt at bus IBUS end (alias
is BI) (0.0 by default)
RATINGS(6) RATE6, rating set 6 line rating (alias is RATE6)
(0.0 by default)
RATINGS(7) RATE7, rating set 7 line rating (alias is RATE7)
(0.0 by default)
RATINGS(8) RATE8, rating set 8 line rating (alias is RATE8)
(0.0 by default)
RATINGS(9) RATE9, rating set 9 line rating (alias is RATE9)
(0.0 by default)
RATINGS(10) RATE10, rating set 10 line rating (alias is
RATE10) (0.0 by default)
RATINGS(11) RATE11, rating set 11 line rating (alias is
RATE11) (0.0 by default)
RATINGS(12) RATE12, rating set 12 line rating (alias is
RATE12) (0.0 by default)
"""
ierr = psspy.branch_data_3(ibus=1, jbus=2, ckt='1', intgar=[1,1,1,0,0,0], realar=[0.63,1.23,0,300,0,0,0,3,1,1,1,1])
ierr = psspy.branch_data_3(ibus=2, jbus=3, ckt='1', intgar=[1,2,1,0,0,0], realar=[0.735,1.435,0,200,0,0,0,3.5,1,1,1,1])
"""
Integer IBUS bus number (input; no default allowed).
Character*2 ID load identifier (input; '1' by default).
Integer INTGAR (4) array of 4 elements specifying (input).
INTGAR(1) STATUS, load status (1 by default)
INTGAR(2) AREA, area number (area of bus IBUS by default)
INTGAR(3) ZONE, zone number (zone of bus IBUS by default)
INTGAR(4) OWNER, owner number (owner of bus IBUS
by default)
Real REALAR (6) array of 6 elements specifying (input).
REALAR(1) PL, constant power active load (0.0 by default)
REALAR(2) QL, constant power reactive load (0.0 by default)
REALAR(3) IP, constant current active load (0.0 by default)
REALAR(4) IQ, constant current reactive load (0.0 by default)
REALAR(5) YP, constant admittance active load (0.0 by
default)
REALAR(6) YQ, constant admittance reactive load (0.0
by default)
"""
ierr = psspy.load_data(ibus=2, id='1', intgar=[1,1,1,1], realar=[320,160,0,0,0,0])
ierr = psspy.load_data(ibus=3, id='1', intgar=[1,1,1,1], realar=[320,160,0,0,0,0])
ierr = psspy.fnsl([0,2,1,0,0,1,0,1])

# Save the case to a file
ierr = psspy.save(case_path)

# Close the PSSE case
ierr = psspy.pssehalt()