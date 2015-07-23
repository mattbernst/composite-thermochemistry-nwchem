'''
G3(MP2) THERMOCHEMICAL METHOD FOR NWCHEM
Daniel R. Haney 2015
  use proper HLC coefficients for CCSDT and QCISDT 6/23/2015

g3mp2.py implements the G3(MP2) composite thermochemical method
in Python 2.7 for inclusion in an NWChem6.5 input file.

It requires:
    NWChem6.5 Zero Point Energy fix
    g3mp2large basis set
    revised 6-31gs basis set, if 3rd row element energies are desired

The original G3(MP2) method is described in:
  "Gaussian-3 theory using reduced Moller-Plesset order"
  Curtiss,Redfern,Raghavachari,Rassolov, and Pople
  J.Chem.Phys. 110 4703 (1999)

As in the GAMESS ab initio package, it has been modified to
use CCSD(T) instead of QCISD(T) in the base energy term.
This results in slightly higher total energies than reported by
Curtiss in the G2/97 test set, faster run times, and lower mean
average deviations in isodesmic reaction energy calculations.

In use, this python file must reside in the same directory
as the NWChem6.5+ input file.

Usage example:

### BEGIN EXAMPLE ###
start ammonia

print low

title "ammonia NH3, G3(MP2)"

geometry units angstroms print xyz
  symmetry C1
  N       -0.01598        0.04072        0.02285
  H        0.99885       -0.01406       -0.00792
  H       -0.33501       -0.01406       -0.94103
  H       -0.33493       -0.81030        0.47878
end

python noprint
import os
import sys
sys.path.append(os.getcwd())

import g3mp2

# invoke it
g3mp2.G3MP2()
end

task python
### END EXAMPLE ###


Invoked as:
    <nwchem executable> NH3.nw >NH3.nwo

Console output:

 ammonia -- NWChem G3(MP2,CCSD(T)) Composite Method
HF(zpe).MP2 optimization.MP2(fc).CCSD(T).GMP2large.HLC.done.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           SUMMARY OF G3(MP2,CCSD(T)) CALCULATIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MP2/6-31G(d)    =   -56.354212   CCSD(T)/6-31G(d) =   -56.372038
MP2/G3MP2large  =   -56.448179   delta(MP2)       =    -0.093967
ZPE(HF/6-31G(d))=     0.033024   ZPE Scale Factor =     0.892900
HLC             =    -0.036680   Free Energy      =     0.000000
Thermal Energy  =     0.039834   Thermal Enthalpy =     0.040778
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
E(G3(MP2)) @ 0K =   -56.469661   E(G3(MP2)) @298K =   -56.466812
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TO DO:
  Free Energy values

'''
import sys
import math
import time

# 'nwchem.py' is an intrinsic module
import nwchem


#______________________________________________________
#
#________________  PHYSICAL CONSTANTS _________________
#______________________________________________________
kCalPerHartree = 627.509451
Boltzmann = 1.3806488E-23
Avogadro = 6.02214129E+23
JoulePerKcal = 4.184E+03
T298 = 298.15

kT_298_perMol = (Boltzmann * T298 * Avogadro) / JoulePerKcal / kCalPerHartree


#______________________________________________________
#
#_________________  GLOBAL VARIABLES  _________________
#______________________________________________________

# Tracing = True          # console stderr output
Tracing = False        # console stderr output


# if False: use CCSD(T), if True: use QCISD(T)
QCISDT_flag = False


def use_qcisdt():
    global QCISDT_flag
    QCISDT_flag = True


def use_ccsdt():
    global QCISDT_flag
    QCISDT_flag = False


''''
    Zero Point Energy scale factor for NWChem HF/6-31G(d)
    Curtiss(original) = 0.8929
    GAMESS-US    = 0.9092
    NWChem6.5    = 0.9232
'''
ZPEScaleFactor = 0.8929


Ezpe = 0.0
Emp2full = 0.0
Emp2frozen = 0.0
Ecc = 0.0
Eg3mp2large = 0.0
Ehlc = 0.0
Ethermal = 0.0
Hthermal = 0.0
E0 = 0.0
E298 = 0.0
H298 = 0.0

dhf0 = 0.0
dhf298 = 0.0

# valence electron variables
nAlpha = 0
nBeta = 0
nFrozen = 0

Charge = 0

#______________________________________________________
# HFtype  normal values are 1 and "RHF", respectively.
# HFtype may be RHF or UHF

HFtype = 'R'    # or 'U'


def get_HFtype():
    '''return HFtype character = '[RU]' + \x00
    '''
    global HFtype
    return HFtype


def set_HFtype(HFstring='RHF'):
    '''set HF theory as RHF or UHF
    default to RHF

    param theory:   RHF,UHF
    theory type :   string
    return      :   nothing
    '''

    global HFtype
    hf = HFstring.upper().strip()[0]
    if hf == 'U':
        HFtype = hf
    else:
        HFtype = 'R'

#______________________________________________________
# spin multiplicity

Multiplicity = 0
multiplets = ["(null)", "singlet", "doublet",
              "triplet", "quartet", "quintet",
              "hextet", "septet", "octet"]


def get_multiplicity_str(mult=1):
    global Multiplicity

    if Multiplicity < 1 and Multiplicity >= range(multiplets):
        Multiplicity = 1

    return multiplets[Multiplicity]


def get_multiplicity():
    '''return multiplicity int
    '''
    global Multiplicity
    return Multiplicity


def set_multiplicity(mult):
    '''determine multiplicity and save globally

    arg mult :  singlet, doublet, ... octet
    mult type:  string
    return   :  nothing
    '''
    global Multiplicity

    mult = mult.lower().strip()
    if mult in multiplets:
        Multiplicity = multiplets.index(mult)
        if Multiplicity < 1:
            Multiplicity = 1
    else:
        Multiplicity = 1
        report("Defaulting to singlet spin multiplicity.")
        report("MULTIPLICITY argument must a string.")
        report("Valid args are singlet, doublet, triplet,")
        report("quartet, quintet, hextet, septet, octet,")
        report("and are cAsE INsenSITive.")

#    debug("Multiplicity is %d\n" % (Multiplicity))

#______________________________________________________
# list of atoms
#  e.g., CH3OH gives ['C','H','H','H','O''H']

AtomsList = []
NumAtoms = 0


def get_atoms_list():
    global AtomsList
    return AtomsList


def set_atoms_list():
    global AtomsList
    global NumAtoms

    # rtdb_get() returns a string if only one atom
    #       or a list of atoms (i.e., tags).
    # Absent type handling,
    # len('CU') is the same as len(['C','U'])
    #
    tags = nwchem.rtdb_get("geometry:geometry:tags")
    tag_type = type(tags).__name__

    if type(tags).__name__ == 'str':
        AtomsList.append(tags)
        debug('AtomsList: %s\n' % tags)

    elif tag_type == 'list':
        AtomsList.extend(tags)
        if debug_flag:
            atmstr = ''
            for atm in tags:
                atmstr += ' ' + atm
            debug('AtomsList: %s\n' % atmstr)

    else:
        AtomsList = []

    NumAtoms = len(AtomsList)
    debug('NumAtoms=%d\n' % NumAtoms)


def is_molecule():
    ''' True if number of atoms > 1
        Needs AtomsList populated from Ezpe
    '''
    return (NumAtoms > 1)


def is_atom():
    ''' True if number of atoms == 1
        Needs AtomsList populated from Ezpe
    '''
    return (NumAtoms == 1)

#______________________________________________________
# return True if molecule is H2
#


def is_H2():

    if (NumAtoms == 2) and \
            (atomic_number(AtomsList[0]) ==  1) and \
            (atomic_number(AtomsList[1]) == 1):
        return True
    else:
        return False

#______________________________________________________
# several console log utilities
##


def say(s):
    '''write to stderr console
       no implicit newline "\n"
    '''
    if (nwchem.ga_nodeid() == 0):
        sys.stderr.write(s)


def log(s):
    '''write to stdout console
    '''
    if (nwchem.ga_nodeid() == 0):
        print(s)


def report(s):
    '''write to stderr,stdout
       add newline to stderr
    '''
    say(s + '\n')
    log(s)


# debug_flag = True          # console stderr output
debug_flag = False        # console stderr output

def set_debug (onoff=0):
    global debug_flag
    debug_flag = (onoff != 0)

def debug(s):
    if debug_flag is True:
        say('DEBUG: %s\n' % s)

#______________________________________________________


def report_summary():
    '''Report results in GAMESS G3(MP2) output format
       for easy comparison.
       log() normally redirects to log file
       say() appears in terminal session
    '''
    global Ezpe
    global Emp2full
    global Emp2frozen
    global Ecc
    global Eg3mp2large
    global Ehlc
    global Ethermal
    global Hthermal
    global E0
    global E298
    global H298
    global QCISDT_flag

    Szpe = Ezpe * ZPEScaleFactor
    deltaMP2 = Eg3mp2large - Emp2frozen

    if QCISDT_flag:
        ccstring = 'QCISDT '
    else:
        ccstring = 'CCSD(T)'

    summary = [
        ("\n"),
        ("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~NWChem6.5"),
        ("               SUMMARY OF G3(MP2,%s) CALCULATIONS              " %
         ccstring),
        ("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"),
        ("    MP2/6-31G(d)    = % 12.6f   %s/6-31G(d) = % 12.6f" %
         (Emp2frozen, ccstring, Ecc)),
        ("    MP2/G3MP2large  = % 12.6f   delta(MP2)       = % 12.6f" %
         (Eg3mp2large, deltaMP2)),
        ("    ZPE(HF/6-31G(d))= % 12.6f   ZPE Scale Factor = % 12.6f" %
         (Szpe, ZPEScaleFactor)),
        ("    HLC             = % 12.6f   Free Energy      = % 12.6f" %
         (Ehlc, 0.0)),
        ("    Thermal Energy  = % 12.6f   Thermal Enthalpy = % 12.6f" %
         (Ethermal, Hthermal)),
        ("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"),
        ("    E(G3(MP2)) @ 0K = % 12.6f   E(G3(MP2)) @298K = % 12.6f" %
         (E0, E298)),
        ("    H(G3(MP2))      = % 12.6f   G(G3(MP2))       = % 12.6f" %
         (H298, 0.0)),
        ("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    ]

    for line in summary:
        report(line)

#______________________________________________________


def report_dHf():
    global dhf0
    global dhf298

    heatsOfFormation = [
        #("\n"),
        ("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"),
        ("          HEAT OF FORMATION   (0K): % 10.2f kCal/mol" % dhf0),
        ("          HEAT OF FORMATION (298K): % 10.2f kCal/mol" % dhf298),
        ("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    ]

    for line in heatsOfFormation:
        report(line)

#______________________________________________________


def reportAll():
    report_summary()
    report_dHf()

#______________________________________________________


def send_nwchem_cmd(s):
    nwchem.input_parse(s)
    debug("\ncmd: [%s]" % s)

#______________________________________________________


def set_charge(chg=0):
    global Charge
    Charge = chg
    send_nwchem_cmd("charge %s" % chg)

#______________________________________________________
'''Look up atomic number from element name or symbol
    since the atoms list may contain either.
'''

def atomic_number(element=''):

    element_dict= {
        'H':1,'HE':2,'LI':3,'BE':4,'B':5,'C':6,'N':7,'O':8,'F':9,'NE':10,
        'NA':11,'MG':12,'AL':13,'SI':14,'P':15,'S':16,'CL':17,'AR':18,
        'K':19,'CA':20,
        'SC':21,'TI':22,'V':23,'CR':24,'MN':25,
        'FE':26,'CO':27,'NI':28,'CU':29,'ZN':30,
        'GA':31,'GE':32,'AS':33,'SE':34,'BR':35,'KR':36,
        'HYDROGEN':1,'HELIUM':2,'LITHIUM':3,'BERYLLIUM':4,'BORON':5,
        'CARBON':6,'NITROGEN':7,'OXYGEN':8,'FLUORINE':9,'NEON':10,
        'SODIUM':11,'MAGNESIUM':12,'ALUMINIUM':13,'SILICON':14,
        'PHOSPHORUS':15,'SULFUR':16,'CHLORINE':17,'ARGON':18,
        'POTASSIUM':19,'CALCIUM':20,'SCANDIUM':21,'TITANIUM':22,
        'VANADIUM':23,'CHROMIUM':24,'MANGANESE':25,'IRON':26,
        'COBALT':27,'NICKEL':28,'COPPER':29,'ZINC':30,'GALLIUM':31,
        'GERMANIUM':32,'ARSENIC':33,'SELENIUM':34,'BROMINE':35,'KRYPTON':36
        }

    try:
        atmnum = element_dict[element.upper()]
    except:
        atmnum = 0

    return atmnum

#______________________________________________________


def atom_core_orbitals(atomicNumber):
    '''This replicates the core electron pair lookup table
        in ./src/geom/geom_core.F
    '''
    # NWChem version 6.5
    '''
    nCoreOrbitals = [0,        # zero index place holder
        0,                                                 0,
        1, 1,                               1, 1, 1, 1, 1, 1,
        5, 5,                               5, 5, 5, 5, 5, 5,
        9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
        18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,
        27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,
        27,27,27,27,27,27,27,27,27,27,27,27,27,27,43,43,43,43,
        43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,
        43,43,43,43,]
    '''

    # GAMESS version 12
    nCoreOrbitals = [
        0,        # zero index place holder
        0, 0,
        1, 1, 1, 1, 1, 1, 1, 1,
        5, 5, 5, 5, 5, 5, 5, 5,
        9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 14, 14, 14, 14, 14, 14,
        23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
        27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27,
        39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39,
        34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 39, 39, 39, 39, 39, 39,
        43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 50]

    if atomicNumber <= len(nCoreOrbitals):
        return nCoreOrbitals[atomicNumber]
    else:
        return 0


####
'''
 SUM the total number of frozen core orbitals in a molecule.

    nFrozen isn't consistently logged to RTDB by
    Tensor Contraction Engine methods.
'''


def sum_core_orbitals():
    global AtomsList

    ncore = 0
    for atom in AtomsList:
        ncore += atom_core_orbitals(atomic_number(atom))

    return ncore

#_______________________________________________________


def ESO(atomic_number=0, charge=0):
    ''' EspinOrbit tabulates the spin orbit energies
    of atomic species in the first three rows as listed in

    Gaussian-4 theory
    Larry A. Curtiss,Paul C. Redfern,Krishnan Raghavachari
    JOURNAL OF CHEMICAL PHYSICS 126, 084108 (2007)
    DOI: 10.1063/1.2436888

    This table contains lists of spin orbit energies for
    [neutral,positive,negative] species.
    Where Curtiss lists no values, 0.0 is returned.

    Values may not agree with current NIST listings.

    Although table values are in milli-Hartrees,
    function ESO returns values in Hartrees
    '''
    #   [neutral, Z+,     Z- ]
    ESpinOrbit = [
        [ 0.0,   0.0,    0.0 ],     # 00 zero index place holder
        [ 0.0,   0.0,    0.0 ],     # 01 H     Hydrogen
        [ 0.0,   0.0,    0.0 ],     # 02 He    Helium
        [ 0.0,   0.0,    0.0 ],     # 03 Li    Lithium
        [ 0.0,   0.0,    0.0 ],     # 04 Be    Beryllium
        [-0.05,  0.0,   -0.03],     # 05 B     Boron
        [-0.14, -0.2,    0.0 ],     # 06 C     Carbon
        [ 0.0,  -0.43,   0.0 ],     # 07 N     Nitrogen
        [-0.36,  0.0,   -0.26],     # 08 O     Oxygen
        [-0.61, -0.67,   0.0 ],     # 09 F     Fluorine
        [ 0.0,  -1.19,   0.0 ],     # 10 Ne    Neon
        [ 0.0,   0.0,    0.0 ],     # 11 Na    Sodium
        [ 0.0,   0.0,    0.0 ],     # 12 Mg    Magnesium
        [-0.34,  0.0,   -0.28],     # 13 Al    Aluminum
        [-0.68, -0.93,   0.0 ],     # 14 Si    Silicon
        [ 0.0,  -1.43,  -0.45],     # 15 P     Phosphorus
        [-0.89,  0.0,   -0.88],     # 16 S     Sulfur
        [-1.34, -1.68,   0.0 ],     # 17 Cl    Chlorine
        [ 0.0,  -2.18,   0.0 ],     # 18 Ar    Argon
        [ 0.0,   0.0,    0.0 ],     # 19 K     Potassium
        [ 0.0,   0.0,    0.0 ],     # 20 Ca    Calcium
        [ 0.0,   0.0,    0.0 ],     # 21 Sc    Scandium
        [ 0.0,   0.0,    0.0 ],     # 22 Ti    Titanium
        [ 0.0,   0.0,    0.0 ],     # 23 V     Vanadium
        [ 0.0,   0.0,    0.0 ],     # 24 Cr    Chromium
        [ 0.0,   0.0,    0.0 ],     # 25 Mn    Manganese
        [ 0.0,   0.0,    0.0 ],     # 26 Fe    Iron
        [ 0.0,   0.0,    0.0 ],     # 27 Co    Cobalt
        [ 0.0,   0.0,    0.0 ],     # 28 Ni    Nickel
        [ 0.0,   0.0,    0.0 ],     # 29 Cu    Copper
        [ 0.0,   0.0,    0.0 ],     # 30 Zn    Zinc
        [-2.51,  0.0,    0.0 ],     # 31 Ga    Gallium
        [-4.41, -5.37,   0.0 ],     # 32 Ge    Germanium
        [ 0.0,  -8.04,   0.0 ],     # 33 As    Arsenic
        [-4.3,   0.0,    0.0 ],     # 34 Se    Selenium
        [-5.6,  -6.71,   0.0 ],     # 35 Br    Bromine
        [ 0.0,  -8.16,   0.0 ]      # 36 Kr    Krypton
    ]

    if not atomic_number in range(len(ESpinOrbit)):
        return 0.0

    if charge > 0:
        ion = 1
    elif charge < 0:
        ion = 2
    else:
        ion = 0

    milliHa_to_Ha = 0.001
    return (ESpinOrbit[atomic_number][ion] * milliHa_to_Ha)

####


def spin_orbit_energy():
    global AtomsList
    global Charge

    if is_molecule():   # no spin orbit corrections for molecules
        return 0.0
    else:       # Its an atom, so AtomsList has only one member
        atom = AtomsList[0]
        return ESO(atomic_number(atom), Charge)

#______________________________________________________


def atomic_DHF(elementNum=0):
    '''return atomic heats of formation
       returns:      tuple
       return type:     2 doubles
    '''
#      atom [dHf(0), dHf(298)]  in kcal/mol

    atomDHF = [
        [0.0,0.0],            # 00 zero placeholder
        [51.63,   52.103  ],  # 01  Hydrogen
        [0.00,    0.00    ],  # 02  Helium
        [37.70,   38.07   ],  # 03  Lithium
        [76.40,   77.40   ],  # 04  Beryllium
        [135.10,  136.30  ],  # 05  Boron
        [169.98,  171.29  ],  # 06  Carbon
        [112.53,  112.97  ],  # 07  Nitrogen
        [58.99,   59.56   ],  # 08  Oxygen
        [18.47,   18.97   ],  # 09  Fluorine
        [0.00,    0.00    ],  # 10  Neon
        [25.76,   25.69   ],  # 11  Sodium
        [34.87,   35.16   ],  # 12  Magnesium
        [80.20,   80.80   ],  # 13  Aluminum
        [107.20,  108.20  ],  # 14  Silicon
        [75.45,   75.65   ],  # 15  Phosphorus
        [65.71,   66.25   ],  # 16  Sulfur
        [28.59,   28.99   ],  # 17  Chlorine
        [0.00,    0.00    ],  # 18  Argon
        [21.27,   21.49   ],  # 19  Potassium
        [42.50,   42.29   ],  # 20  Calcium
        [0.0,0.0],[0.0,0.0],  # transition elements 21-30 Sc-Zn
        [0.0,0.0],[0.0,0.0],  # transition elements
        [0.0,0.0],[0.0,0.0],  # transition elements
        [0.0,0.0],[0.0,0.0],  # transition elements
        [0.0,0.0],[0.0,0.0],  # transition elements
        [65.00,   65.00   ],  # 31  Gallium
        [88.91,   88.91   ],  # 32  Germanium
        [73.90,   72.42   ],  # 33  Arsenic
        [55.76,   54.27   ],  # 34  Selenium
        [26.74,   28.18   ],  # 35  Bromine
        [0.0,     0.0     ],  # 36  Krypton
    ]

    debug('atomic_DHF: elementNum=%d' % elementNum)
    if elementNum < len(atomDHF):
        debug('atomic_DHF: E,H=%.2f,%.2f' %
              (atomDHF[elementNum][0], atomDHF[elementNum][1]))
        return atomDHF[elementNum][0], atomDHF[elementNum][1]
    else:
        debug('atomic_DHF: error: element %d not in table?' % elementNum)
        return 0.0, 0.0


#_______________________________________________________
#
# return the total G3(MP2) energy of an atom
# for use in Delta Hf calculations
#
def E0_atom_qcisdt(elementNum=0):
    e0_qcisdt = [
         0.0,  # 00 zero index place holder
        -0.501839,     # 01 H     Hydrogen
        -2.902543,     # 02 He    Helium
        -7.434048,     # 03 Li    Lithium
        -14.629262,     # 04 Be    Beryllium
        -24.607093,     # 05 B     Boron
        -37.789349,     # 06 C     Carbon
        -54.525198,     # 07 N     Nitrogen
        -74.989850,     # 08 O     Oxygen
        -99.641120,     # 09 F     Fluorine
        -128.828970,     # 10 Ne    Neon
        -161.848004,     # 11 Na    Sodium
        -199.650845,     # 12 Mg    Magnesium
        -241.936973,     # 13 Al    Aluminum
        -288.939460,     # 14 Si    Silicon
        -340.826670,     # 15 P     Phosphorus
        -397.663794,     # 16 S     Sulfur
        -459.687272,     # 17 Cl    Chlorine
        -527.060963,     # 18 Ar    Argon
        -599.160512,     # 19 K     Potassium
        -676.789424,     # 20 Ca    Calcium
                 # elements 21-30 Sc-Zn
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        -1923.538354,     # 31 Ga    Gallium
        -2075.639002,     # 32 Ge    Germanium
        -2234.516874,     # 33 As    Arsenic
        -2400.180197,     # 34 Se    Selenium
        -2572.784022,     # 35 Br    Bromine
        -2752.417979       # 36 Kr    Krypton
        ]

    # Ideal gas kinetic energy contribution

    Ethermal = (5.0 / 2) * kT_298_perMol

    if elementNum < len(e0_qcisdt) and e0_qcisdt[elementNum] < 0.0:
        e0 = e0_qcisdt[elementNum]
        e298 = e0 + Ethermal
        return (e0, e298)
    else:
        return (0.0, 0.0)

#_______________________________________________________


def E0_atom_ccsdt(elementNum=0):

    e0_ccsdt = [
        0.0,  # 00 zero index place holder
        -0.501765,     # 01 H     Hydrogen
        -2.902353,     # 02 He    Helium
        -7.433974,     # 03 Li    Lithium
        -14.629072,     # 04 Be    Beryllium
        -24.606789,     # 05 B     Boron
        -37.788989,     # 06 C     Carbon
        -54.524779,     # 07 N     Nitrogen
        -74.989201,     # 08 O     Oxygen
        -99.640198,     # 09 F     Fluorine
        -128.827752,     # 10 Ne    Neon
        -161.847930,     # 11 Na    Sodium
        -199.650655,     # 12 Mg    Magnesium
        -241.936660,     # 13 Al    Aluminum
        -288.939067,     # 14 Si    Silicon
        -340.826225,     # 15 P     Phosphorus
        -397.663215,     # 16 S     Sulfur
        -459.686583,     # 17 Cl    Chlorine
        -527.060194,     # 18 Ar    Argon
        -599.160438,     # 19 K     Potassium
        -676.789234,     # 20 Ca    Calcium
        # elements 21-30 Sc-Zn
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        -1923.536619,     # 31 Ga
        -2075.637891,     # 32 Ge
        -2234.516031,     # 33 As
        -2400.179359,     # 34 Se
        -2572.783154,     # 35 Br
        -2752.417075       # 36 Kr
        ]

    # Ideal gas kinetic energy contribution
    eth = (5.0 / 2) * kT_298_perMol

    if elementNum < len(e0_ccsdt) and e0_ccsdt[elementNum] < 0.0:
        e0 = e0_ccsdt[elementNum]
        e298 = e0 + eth
        return (e0, e298)
    else:
        return (0.0, 0.0)

#_______________________________________________________


def E0_atom(elementNum=0):

    if QCISDT_flag:
        return E0_atom_qcisdt(elementNum)
    else:
        return E0_atom_ccsdt(elementNum)


#_______________________________________________________
# calculate heat of formation at 0K and 298K
#
# ALWAYS returns False, silently fails so that
# the summary is printed.
#
def calc_deltaHf():
    global E0
    global H298
    global dhf0
    global dhf298

    sum_atoms_E0 = 0.0
    sum_atoms_H298 = 0.0
    sum_atoms_dhf0 = 0.0
    sum_atoms_dhf298 = 0.0

    for atom in AtomsList:
        e0, h298 = E0_atom(atomic_number(atom))
        sum_atoms_E0 += e0
        sum_atoms_H298 += h298

        d0, d298 = atomic_DHF(atomic_number(atom))
        sum_atoms_dhf0 += d0
        sum_atoms_dhf298 += d298
        debug('sumDHF0,sumDHF298 = %.2f,%.2f' %
            (sum_atoms_dhf0, sum_atoms_dhf298))

    dhf0 = (E0 - sum_atoms_E0) * kCalPerHartree + sum_atoms_dhf0

    dhf298 = (H298 - sum_atoms_H298) * kCalPerHartree + sum_atoms_dhf298

    debug('dhf0,dhf298 = %.2f,%.2f' % (dhf0, dhf298))

    return False

#______________________________________________________


def init_g3mp2(charge=0, mult='singlet', use_qcisdt_f=False):
    '''say hello
    '''

    set_charge(charge)
    set_multiplicity(mult)

    if get_multiplicity() > 1:
        set_HFtype('uhf')
    else:
        set_HFtype('rhf')

    if use_qcisdt_f:
        use_qcisdt()
    else:
        use_ccsdt()

    if QCISDT_flag:
        ccstring = 'QCISD(T)'
    else:
        ccstring = 'CCSD(T)'

    title = nwchem.rtdb_get("title")
    say(" %s -- NWChem G3(MP2,%s) Composite Method\n" % (title, ccstring))

    return False


def build_SCF_cmd():

    multstr = get_multiplicity_str()
    hftype = get_HFtype()
    return ("scf ; direct ; %sHF ; %s ; end" % (hftype, multstr))


#______________________________________________________

def HF_optimize():
    '''
       HF/6-31G(d) optimization
    '''
    say('optimize.')
    send_nwchem_cmd("basis noprint ; * library 6-31G* ; end")
    scfcmd = build_SCF_cmd()
    send_nwchem_cmd(scfcmd)

    # fetch and copy atom names list (tags) which enumerates atoms.
    # only available _after_ SCF statement

    set_atoms_list()

    # optimize the geometry, ignore energy and gradient results
    if is_atom():
        en = nwchem.task_energy("scf")
    else:
        en, grad = nwchem.task_optimize("scf")

    debug("HF/6-31G(d) SCF:energy = %.7f Ha" % (en))
    return False


#______________________________________________________
def HF_zpe():
    ''' run hessian on equilibrium geometry
        get zero point energy.
                                                                          n
       note: linear tri-atomics and larger
       give high ZPE values in NWChem @ HF/6-31G*
    '''

    global Ezpe
    global Ethermal
    global Hthermal
    global AtomsList

    #AUKCAL = 627.5093314   # bogus
    AUKCAL = kCalPerHartree
    c = 2.99792458E+10
    h = 6.62606957E-27
    kgas = 1.3806488E-16     # cgs units
    Rgas = 1.9872041 / 1000.0 / AUKCAL     # atomic units

    temperature = T298

    say("zpe.")

    if is_atom():
        Ezpe = 0.0
        Ethermal = 1.5 * Rgas * temperature
        Hthermal = Ethermal + (Rgas * temperature)
        # Ethermal = 0.001416      # 3/2 * RT
        # Hthermal = 0.002360      # Ethermal + kT
        return False

    # run hessian on equilibrium geometry
    # ignore ZPE, calculate it from vibrations list
    zpe, vibs, intens = nwchem.task_freq("scf")

    # Handroll the ZPE because NWChem's zpe accumulates
    # truncation error from 3 sigfig physical constants.

    vibsum = 0.0
    for freq in vibs:
        if (freq > 0.1):
            vibsum += freq

    cm2Ha = 219474.6    # cm-1 to Hartree conversion
    Ezpe = vibsum / (2.0 * cm2Ha)

    # get total thermal energy, enthalpy
    eth = nwchem.rtdb_get("vib:Ethermal")
    hth = nwchem.rtdb_get("vib:Hthermal")
    debug("NWChem zpe,E,H thermal= %.6f, %.6f, %.6f\n" % (Ezpe,eth, hth))

    eth = 0.0
    hth = 0.0
    xdum = 0.0

    for freq in vibs:
        if (freq > 0.1):
            # freqency temperature in Kelvin from cm-1
            thetav = freq * (h * c / kgas)
            if (temperature > 0.0):
                xdum = math.exp(-thetav / temperature)
            else:
                xdum = 0.0

            xdum = xdum / (1.0 - xdum)
            eth = eth + thetav * (0.5 + xdum)

    # linear boolean is available only after task_freq('scf') runs
    # NWChem only writes the flag if molecule is linear

    eth = eth * Rgas

    try:
        is_linear = nwchem.rtdb_get("vib:linear")
    except:
        is_linear = False

    if (is_linear and QCISDT_flag):
        # translational(3/2RT) and rotation(2/2RT) thermal corrections
        eth = eth + 2.5 * Rgas * temperature
    else:
        # translational(3/2RT) and rotation(3/2RT) thermal corrections
        eth = eth + 3.0 * Rgas * temperature

    # Hthermal = eth+pV=eth+RT, since pV=RT
    hth = eth + Rgas * temperature
    debug("RgasT= %.9f, kT_298_perMol= %.9f\n" %
        (Rgas * temperature, kT_298_perMol))
    debug("Handrolled E,H thermal= %.6f, %.6f\n" % (eth, hth))

    Ethermal = eth
    Hthermal = hth

    return False
#


#______________________________________________________
def MP2_optimize():
    '''Optimize geometry at MP2(full)/6-31G(d)
    '''

    global Emp2full
    global Emp2frozen

    say('MP2 optimize.')

    if Multiplicity > 1:
        send_nwchem_cmd("tce ; scf ; mp2 ; end")
        if is_atom():
            en = nwchem.task_energy("tce")
        else:
            en, grad = nwchem.task_optimize("tce")
    else:
        if is_atom():
            en = nwchem.task_energy("mp2")
        else:
            en, grad = nwchem.task_optimize("mp2")

    debug('optimize: MP(2,full)/6-31G*= %.6f\n' % en)

    Emp2full = en

    return False

#______________________________________________________


def MP2_frozen():
    # MP2(fc)/6-31G* single point energy

    global Emp2frozen

    say('MP2(frozen).')

    scfcmd = build_SCF_cmd()
    send_nwchem_cmd(scfcmd)

    send_nwchem_cmd("unset mp2:*")
    send_nwchem_cmd("mp2 ; freeze atomic ; end")

    if Multiplicity > 1:
        send_nwchem_cmd("unset tce:*")
        send_nwchem_cmd("tce ; scf ; mp2 ; freeze atomic ; end")

        en = nwchem.task_energy("tce")
    else:
        en = nwchem.task_energy("mp2")

    debug('MP2 frozen: en=%.6f\n' % en)
    Emp2frozen = en

    return False


def CCSDT_frozen():
    '''get CCSDT(fc)/6-31G(d) energy
    '''
    global Ecc

    say("CCSD(T).")

    # NWChem6.5 direct CCSD(T) at 6-31G* fails on H2
    if Multiplicity > 1 or is_atom() or is_H2():
        send_nwchem_cmd("unset tce:*")
        send_nwchem_cmd("tce ; ccsd(t) ; freeze atomic ; end")
        en = nwchem.task_energy("tce")
    else:
        send_nwchem_cmd("ccsd ; freeze atomic ; end")
        en = nwchem.task_energy("ccsd(t)")

    debug(' CCSD(T): en=%.6f\n' % en)

    Ecc = en

    return False

#______________________________________________________

def QCISDT_frozen():

# QCISDT(frozen core)/6-31G(d) single point energy

    global Ecc

    say("QCISD(T).")

    nwchem.input_parse("tce ; qcisd(t) ; freeze atomic ; end")

    en = nwchem.task_energy("tce")
    Ecc = en

    return False

#______________________________________________________

def ccsdt_qcisdt_frozen(flag=False):

    if QCISDT_flag:
        QCISDT_frozen()
    else:
        CCSDT_frozen()

#______________________________________________________

def MP2_g3mp2large():
    '''get MP2(fc)/G3MP2large single point energy
    '''
    global Eg3mp2large

    say('GMP2large.')

    send_nwchem_cmd('''
        basis spherical
          * library g3mp2large
        end''')

    send_nwchem_cmd("unset mp2:*")
    send_nwchem_cmd("mp2 ; freeze atomic ; end")

    if Multiplicity > 1 or is_atom():
        send_nwchem_cmd("unset tce:*")
        send_nwchem_cmd("tce ; mp2 ; freeze atomic ; end")
        en = nwchem.task_energy("tce")
    else:
        en = nwchem.task_energy("mp2")

    debug(' g3mp2large: en=%.6f\n' % en)

    Eg3mp2large = en

    return False


#______________________________________________________

def HLC_generic(A=0.009170, B=0.004455, C=0.009155, D=0.001947):
    '''calculate High Level Correction term
        from alpha and beta VALENCE electron count.
    '''

    say('HLC.')

    nClosed = nwchem.rtdb_get("scf:nclosed")
    nOpen = nwchem.rtdb_get("scf:nopen")
    nElec = nwchem.rtdb_get("scf:nelec")
    nFrozen = sum_core_orbitals()

    # According to Curtiss,
    # nBeta = num valence pairs
    # nAlpha = num unpaired or remaining valence electrons
    # subject to the constraint that nAlpha >= nBeta

    nBeta = nClosed - nFrozen

    nAlpha = nElec - (nFrozen * 2) - nBeta

    debug('\nclosed=%d open=%d frozen=%d nAlpha=%d nBeta=%d\n' %
        (nClosed,nOpen,nFrozen,nAlpha,nBeta))

    if nAlpha < nBeta:
        nAlpha, nBeta = nBeta, nAlpha

    if is_molecule():
        Ehlc = -(A * nBeta) - B * (nAlpha - nBeta)
    else:  # it's an atom, calc is different
        Ehlc = -(C * nBeta) - D * (nAlpha - nBeta)

    return False

#______________________________________________________


def HLC_qcisdt():
    ''' empirical correction coefficients
    for Curtiss original G3(MP2) method
    '''
    A = 0.009279
    B = 0.004471
    C = 0.009345
    D = 0.002021
    return HLC_generic(A, B, C, D)


#______________________________________________________
def HLC_ccsdt():
    ''' empirical correction coefficients
    for later G3(MP2,CCSDT) method
    '''
    A = 0.009170
    B = 0.004455
    C = 0.009155
    D = 0.001947
    return HLC_generic(A, B, C, D)


def HLC():
    global QCISDT_flag

    if QCISDT_flag:
        return HLC_qcisdt()
    else:
        return HLC_ccsdt()

#______________________________________________________


def calc_total_energies():
    '''G3MP2 step 6
        calculate the E(G3MP2)@0K, @298K energies
    '''

    global Ezpe
    global Ethermal
    global Hthermal
    global Emp2frozen
    global Ecc
    global Eg3mp2large
    global Ehlc
    global E0
    global E298
    global H298
    global ESO
    

    E0 = Ecc + \
        (Eg3mp2large - Emp2frozen) + \
        (Ezpe * ZPEScaleFactor) + \
        Ehlc

    if is_atom():
        E0 += spin_orbit_energy()

    E298 = E0 + (Ethermal - Ezpe)
    H298 = E0 + (Hthermal - Ezpe)  # kT_298_perMol

    return False


#______________________________________________________
#______________________________________________________

def G3MP2(charge=0, mult='singlet', use_qcisdt_f=False):

    g3mp2_function = [
        HF_optimize,
        HF_zpe,
        MP2_optimize,
        MP2_frozen,
        ccsdt_qcisdt_frozen,
        MP2_g3mp2large,
        HLC,
        calc_total_energies,
        calc_deltaHf
    ]

    if init_g3mp2(charge, mult, use_qcisdt_f) == True:
        sys.exit(0)

    abnormal_end = False

    t0=time.time()

    for i in range(len(g3mp2_function)):
        if g3mp2_function[i]():
            abnormal_end = True
            break

    et=time.time()-t0
    report("\nWall: %.2f seconds" % et)

    if abnormal_end is False:
        reportAll()

######################### END G3MP2 MODULE #########################
