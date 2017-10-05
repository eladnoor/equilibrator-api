import os
import inspect
from numpy import log

SCRIPT_DIR = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
DATA_DIR = os.path.join(os.path.split(SCRIPT_DIR)[0], 'data')
TEST_DIR = os.path.join(os.path.split(SCRIPT_DIR)[0], 'tests')

POSSIBLE_REACTION_ARROWS = ['<=>', '<=', '=>', '->', '<->', '=']

R = 8.31e-3   # kJ/(K*mol)
DEFAULT_TEMP = 298.15  # K
DEFAULT_IONIC_STRENGTH = 0.1  # mM
DEFAULT_PH = 7.0
DEFAULT_PMG = 14.0
DEFAULT_PHASE = 'aqueous'
RT = R * DEFAULT_TEMP
RTlog10 = RT * log(10)

DEBYE_HUECKLE_A = 2.91482
DEBYE_HUECKLE_B = 1.6
MG_FORMATION_ENERGY = -455.3  # kJ/mol, formation energy of Mg2+
FARADAY = 96.485  # kC/mol

DEFAULT_CONC_LB = 1e-6
DEFAULT_CONC_UB = 1e-2