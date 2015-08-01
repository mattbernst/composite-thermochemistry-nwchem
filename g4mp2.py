'''
 G4(MP2) composite method for Python under NWChem 6.5
 Implementation by Matt B. Ernst and Daniel R. Haney
 7/18/2015

 Gaussian-4 theory using reduced order perturbation theory
 Larry A. Curtiss,Paul C. Redfern,Krishnan Raghavachari
 THE JOURNAL OF CHEMICAL PHYSICS 127, 124105 2007


 1 optimize     @ B3LYP/6-31G(2df,p)
 2 Ezpe         = zpe at B3LYP/6-31G(2df,p)
 3 E(MP2)       = MP2(fc)/6-31G(d)

 4 E(ccsd(t))   = CCSD(fc,T)/6-31G(d)
 5 E(HF/G3LXP)  = HF/G3LargeXP
 6 E(G3LargeXP) = MP2(fc)/G3LargeXP

 7 E(HF1)       = HF/g4mp2-aug-cc-pVTZ
 8 E(HF2)       = HF/g4mp2-aug-cc-pVQZ
   E(HFlimit)   = extrapolated HF limit, =CBS

   delta(HF)    = E(HFlimit) - E(HF/G3LargeXP)
   E(SO)        = spin orbit energy
   Ehlc         = High Level Correction

   E(G4(MP2)) = E(CCSD(T)) +
                E(G3LargeXP) - E(MP2) +
                Delta(HFlimit) +
                E(SO) +
                E(HLC) +
                Ezpe * scale_factor
'''
#______________________________________________________

import sys
import math
import time

# 'nwchem.py' is an intrinsic module
import nwchem

#______________________________________________________
#
#________________  PHYSICAL CONSTANTS _________________
#______________________________________________________

kCalPerHartree  = 627.509451
Boltzmann       = 1.3806488E-23
Avogadro        = 6.02214129E+23
JoulePerKcal    = 4.184E+03
T298            = 298.15

kT_298_perMol   = (Boltzmann * T298 * Avogadro) / JoulePerKcal / kCalPerHartree


#______________________________________________________
#
#_________________  GLOBAL VARIABLES  _________________
#______________________________________________________

#debug_flag = True          # console stderr output
debug_flag = False        # console stderr output

def set_debug(onoff=0):
    global debug_flag
    debug_flag = (onoff != 0)

''''
Zero Point Energy scale factor for B3LYP/6-31G(2df,p)
'''
ZPEScaleFactor = 0.9854    # Curtiss scale factor for Gaussian 09
#ZPEScaleFactor = 0.9798     # Truhlar scale Factor for NWChem 6.5

Ezpe        = 0.0
Emp2        = 0.0
Eccsdt      = 0.0
Ehfg3lxp    = 0.0
Emp2g3lxp   = 0.0
Ehf1        = 0.0
Ehf2        = 0.0
Ecbs        = 0.0
Ehlc        = 0.0
Ethermal    = 0.0
Hthermal    = 0.0
ESO         = 0.0
E0          = 0.0
E298        = 0.0
H298        = 0.0

dhf0        = 0.0
dhf298      = 0.0

Charge      = 0

#______________________________________________________
# HFtype  normal values are 1 and "RHF", respectively.
# HFtype may be RHF or UHF

HFtype = 'R'    # or 'U'

def get_HFtype():
    '''return HFtype character = '[RU]' + \x00
    '''
    global HFtype
    return HFtype


def set_HFtype (HFstring='RHF'):
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
multiplets  = ["(null)","singlet","doublet",
                "triplet","quartet","quintet",
                "hextet","septet", "octet"]


def get_multiplicity_str (mult=1):
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



    debug("Multiplicity is %d\n" % (Multiplicity))

#______________________________________________________
# list of atoms
#  e.g., CH3OH gives ['C','H','H','H','O''H']

AtomsList   = []
NumAtoms    = 0

def get_atoms_list():
    global AtomsList
    return AtomsList

def set_atoms_list ():
    global AtomsList
    global NumAtoms

    # rtdb_get() returns EITHER
    #       a string if only one atom, or
    #       a list of atoms (i.e., tags).
    #
    # Without type handling,
    # len('CU') is the same as len(['C','U'])
    #
    tags = nwchem.rtdb_get("geometry:geometry:tags")
    tag_type = type(tags).__name__

    if tag_type == 'str':
        AtomsList.append(tags)
        debug('AtomsList: %s\n' % tags)

    elif tag_type == 'list':
        AtomsList.extend(tags)
        if debug_flag:
            atmstr=''
            for atm in tags:
                atmstr += ' '+atm
            debug('AtomsList: %s\n' % atmstr)
    else:
        AtomsList = []

    NumAtoms = len(AtomsList)
    debug('NumAtoms=%d\n' % NumAtoms)

def is_molecule ():
    ''' True if number of atoms > 1
    '''
    return (NumAtoms > 1)

def is_atom ():
    ''' True if number of atoms == 1
    '''
    return (NumAtoms == 1)

# detect special case for when CCSD(T)/6-31G* fails
def is_H2 ():
    if NumAtoms == 2 and \
            atomic_number(AtomsList[0]) == 1 and \
            atomic_number(AtomsList[1]) == 1 :
        return True
    else:
        return False
#______________________________________________________
## several console log utilities
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
    #say(s+'\n')
    log(s)


def debug (s):
    if debug_flag:
        say('DEBUG: %s\n' % s)


def report_summary ():
    '''Report results in GAMESS G3(MP2) output format
       for easy comparison.
       log() normally redirects to log file
       say() appears in terminal session
    '''
    global Eb3lyp
    global Ezpe
    global Emp2
    global Eccsdt
    global Ehfg3lxp
    global Emp2g3lxp
    global Ehf1
    global Ehf2
    global Ecbs
    global Ehlc
    global Ethermal
    global Hthermal
    global E0
    global E298
    global H298

    Szpe = Ezpe * ZPEScaleFactor
    dMP2 = Emp2g3lxp - Emp2
    dHF  = Ecbs - Ehfg3lxp

    summary = [
        ("\n    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~NWChem6.5"),
        ("                  SUMMARY OF G4(MP2) CALCULATIONS"),
        ("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"),
        ("    B3LYP/6-31G(2df,p)= % 12.6f   HF/maug-cc-p(T+d)Z= % 12.6f" % (Eb3lyp,Ehf1)),
        ("    HF/CBS            = % 12.6f   HF/maug-cc-p(Q+d)Z= % 12.6f" % (Ecbs,Ehf2)),
        ("    MP2/6-31G(d)      = % 12.6f   CCSD(T)/6-31G(d)  = % 12.6f" % (Emp2,Eccsdt)),
        ("    HF/G3MP2LARGEXP   = % 12.6f   MP2/G3MP2LARGEXP  = % 12.6f" % (Ehfg3lxp,Emp2g3lxp)),
        ("    DE(MP2)           = % 12.6f   DE(HF)            = % 12.6f" % (dMP2,dHF)),
        ("    ZPE(B3LYP)        = % 12.6f   ZPE SCALE FACTOR  = % 12.6f" % (Szpe,ZPEScaleFactor)),
        ("    HLC               = % 12.6f   FREE ENERGY       = % 12.6f" % (Ehlc,0.0)),
        ("    THERMAL ENERGY    = % 12.6f   THERMAL ENTHALPY  = % 12.6f" % (Ethermal,Hthermal)),
        ("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"),
        ("    E(G4(MP2)) @ 0K   = % 12.6f   E(G4(MP2)) @298K  = % 12.6f" % (E0,E298)),
        ("    H(G4(MP2))        = % 12.6f   G(G4(MP2))        = % 12.6f" % (H298,0.0)),
        ("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"),
    ]

    for line in summary:
        report(line)

#______________________________________________________
def report_dHf ():
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

def report_all ():
    report_summary ()
    report_dHf ()

#______________________________________________________

def send_nwchem_cmd (s):
    nwchem.input_parse(s)
    debug("cmd: [%s]" % s)

#______________________________________________________

def set_charge (chg=0):
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

def atom_core_orbitals (atomicNumber):
    '''This replicates the core electron pair lookup table
        in ./src/geom/geom_core.F
    '''
    #NWChem version 6.5
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
    nCoreOrbitals = [0,        # zero index place holder
        0,                                                  0,
        1,  1,                               1, 1, 1, 1, 1, 1,
        5,  5,                               5, 5, 5, 5, 5, 5,
        9,  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,14,14,14,14,14,14,
        23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,
        27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,
        39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,
        34,34,34,34,34,34,34,34,34,34,39,39,39,39,39,39,
        43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,50]

    if atomicNumber <= len(nCoreOrbitals):
        return nCoreOrbitals[atomicNumber]
    else:
        return 0


#______________________________________________________
'''
 SUM the total number of frozen core orbitals in a molecule.

    nFrozen isn't consistently logged to RTDB by
    Tensor Contraction Engine methods.
'''
def sum_core_orbitals ():
    global AtomsList

    ncore = 0
    for atom in AtomsList:
        ncore += atom_core_orbitals(atomic_number(atom))

    return ncore

#_______________________________________________________

def E_spin_orbit (atomic_number=0, charge=0):
    ''' EspinOrbit tabulates the spin orbit energies
    of atomic species in the first three rows as listed in

    Gaussian-4 theory
    Larry A. Curtiss,Paul C. Redfern,Krishnan Raghavachari
    JOURNAL OF CHEMICAL PHYSICS 126, 084108 (2007)
    DOI: 10.1063/1.2436888

    This table contains lists of spin orbit energies for
    [neutral,positive,negative] species.
    When Curtiss lists no values, 0.0 is returned.

    Values may not agree with current NIST listings.

    Although table values are in milli-Hartrees,
    function E_spin_orbit returns values in Hartrees
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

    if charge >0:
        ion = 1
    elif charge <0:
        ion = 2
    else:
        ion =0

    milliHa_to_Ha = 0.001
    return (ESpinOrbit[atomic_number][ion] * milliHa_to_Ha)

####

def spin_orbit_energy ():
    global AtomsList
    global Charge

    if is_molecule():   # no spin orbit corrections for molecules
        eso = 0.0
    else:       # Its an atom, so AtomsList has only one member
        atom=AtomsList[0]
        eso = E_spin_orbit(atomic_number(atom),Charge)
        
    ESO = eso
    return False

#______________________________________________________

def atomic_DHF (elementNum=0):
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

    if elementNum <len(atomDHF):
        debug('atomic_DHF: E0,E298=%.2f,%.2f' %
            (atomDHF[elementNum][0], atomDHF[elementNum][1]))
        return atomDHF[elementNum][0], atomDHF[elementNum][1]
    else:
        report('atomic_DHF: error: element %d not in table?' % elementNum)
        return 0.0,0.0


#_______________________________________________________
def E0_atom (elementNum=0):
    '''List of precalculated atomic G4(MP2) energies at 0K
        returns E(0K), E(298.15K) tuple
    '''

    e0_g4mp2 = [ 0.0 ,  # 00 zero index place holder
           -0.502094 ,  # 01    H   Hydrogen
           -2.892437 ,  # 02    He  Helium
           -7.434837 ,  # 03    Li  Lithium
          -14.618701 ,  # 04    Be  Beryllium
          -24.610037 ,  # 05    B   Boron
          -37.794204 ,  # 06    C   Carbon
          -54.532825 ,  # 07    N   Nitrogen
          -75.002483 ,  # 08    O   Oxygen
          -99.659686 ,  # 09    F   Fluorine
         -128.854769 ,  # 10    Ne  Neon
         -161.860999 ,  # 11    Na  Sodium
         -199.646948 ,  # 12    Mg  Magnesium
         -241.944728 ,  # 13    Al  Aluminum
         -288.947800 ,  # 14    Si  Silicon
         -340.837016 ,  # 15    P   Phosphorus
         -397.676523 ,  # 16    S   Sulfur
         -459.703691 ,  # 17    Cl  Chlorine
         -527.083295 ,  # 18    Ar  Argon
         -599.166975 ,  # 19    K   Potassium
         -676.784184 ,  # 20    Ca  Calcium
         0.0,0.0,0.0,   # 21-23 transition metals
         0.0,0.0,0.0,   # 24-26 transition metals
         0.0,0.0,0.0,   # 27-29 transition metals
         0.0,           # 30    transition metals
        -1923.601298 ,  # 31    Ga  Gallium
        -2075.700329 ,  # 32    Ge  Germanium
        -2234.578295 ,  # 33    As  Arsenic
        -2400.243694 ,  # 34    Se  Selenium
        -2572.850476 ,  # 35    Br  Bromine
        -2752.487773 ,  # 36    Kr  Krypton
    ]


    # Ideal gas kinetic energy contribution

    eThermal = (5.0/2) * kT_298_perMol

    if elementNum <len(e0_g4mp2) and e0_g4mp2[elementNum]<0.0:
        e0 = e0_g4mp2[elementNum]
        e298 = e0 + eThermal
        return (e0,e298)

    else:
        return (0.0,0.0)

#_______________________________________________________
# calculate heat of formation at 0K and 298K
#
# ALWAYS returns False, silently fails so that
# the summary is printed.
#
def calc_deltaHf ():
    global E0
    global H298
    global dhf0
    global dhf298

    sum_atoms_E0     = 0.0
    sum_atoms_H298   = 0.0
    sum_atoms_dhf0   = 0.0
    sum_atoms_dhf298 = 0.0

    for atom in AtomsList:
        e0,h298 = E0_atom(atomic_number(atom))
        sum_atoms_E0 += e0
        sum_atoms_H298 += h298

        d0,d298 = atomic_DHF(atomic_number(atom))
        sum_atoms_dhf0       += d0
        sum_atoms_dhf298     += d298

        debug('sumDHF0,sumDHF298 = %.2f,%.2f' %
            (sum_atoms_dhf0,sum_atoms_dhf298))

    dhf0 = (E0 - sum_atoms_E0) * kCalPerHartree + sum_atoms_dhf0

    dhf298 = (H298 - sum_atoms_H298) * kCalPerHartree + sum_atoms_dhf298

    debug('dhf0,dhf298 = %.2f,%.2f' % (dhf0,dhf298))

    return False

#_______________________________________________________

def init_g4mp2(charge=0,mult='singlet'):
    '''say hello
    '''

    set_charge(charge)
    set_multiplicity(mult)

    if get_multiplicity() > 1:
        set_HFtype('uhf')
    else:
        set_HFtype('rhf')

    title = nwchem.rtdb_get("title")
    say(" %s -- NWChem G4(MP2) Composite Method\n" % (title))

    return False

#_______________________________________________________

def build_SCF_cmd ():
    global Multiplicity

    if Multiplicity != 1:
        multstr = get_multiplicity_str()
        hftype = get_HFtype()
        return ("scf ; %sHF ; %s ; end" % (hftype,multstr))
    else:
        return "scf ; direct ; singlet ; end"

#_______________________________________________________

def limits_high():
    """Increase iterations allowed for geometry optimization and electronic
    convergence.
    """

    #send_nwchem_cmd("scf; maxiter 999; end")
    send_nwchem_cmd("driver; maxiter 999; end")

def optimize():
# 1 optimize  B3LYP/6-31G(2df,p)

    global Eb3lyp
    global Multiplicity

    say('optimize.')
    limits_high()

    send_nwchem_cmd("basis noprint ; * library 6-31G(2df,p) ; end")

    # canonical B3LYP spec
    # GAMESS-US b3lyp uses VWN_5
    # Gaussian uses VWN_3
    # NWChem uses something entirely different

    b3lyp_GAMESS = 'xc HFexch 0.2 slater 0.8 becke88 nonlocal 0.72 vwn_5 0.19 lyp 0.81'
    b3lyp_Gaussian = 'xc HFexch 0.2 slater 0.8 becke88 nonlocal 0.72 vwn_3 0.19 lyp 0.81'
    b3lyp_NWChem = 'xc b3lyp'

    #blips = b3lyp_GAMESS
    #blips = b3lyp_NWChem
    blips = b3lyp_Gaussian

    if Multiplicity>1:
        send_nwchem_cmd('dft ; odft ; mult %d ; %s ; end' % (Multiplicity,blips))
    else:
        send_nwchem_cmd('dft ; %s ; end' % blips)

    # fetch and copy atom names list (tags) which enumerates atoms.
    # only available _after_ SCF statement

    set_atoms_list()

    if is_atom():
        en = nwchem.task_energy("dft")
    else:
        en,grad = nwchem.task_optimize("dft")

    Eb3lyp = en

    debug('B3LYP/6-31G(2df,p) energy = %f Ha' % (en))
    return False

#_______________________________________________________

def E_zpe ():
    say('ZPE.')
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

    # these identifiers replicate those used in NWChem freq calculations
    #AUKCAL  = 627.5093314  # WRONG!
    AUKCAL  = kCalPerHartree
    c       = 2.99792458E+10
    h       = 6.62606957E-27
    kgas    = 1.3806488E-16     # cgs units
    Rgas    = 1.9872041/1000.0/AUKCAL     # atomic units

    temperature = 298.15

    if is_atom():
        Ezpe = 0.0
        Ethermal = 1.5 * Rgas * temperature         # 3/2 * RT
        Hthermal = Ethermal + (Rgas * temperature)
        return False

    # run hessian on equilibrium geometry
    # ignore ZPE, calculate it from vibrations list
    try:
        zpe,vibs,intens = nwchem.task_freq("dft")
    except NWChemError,message:
        report("NWChem error: %s\n" % message)
        report("FAILED: Zero Point Energy calculation")
        return True

    # Handroll the ZPE because NWChem's zpe accumulates
    # truncation error from 3 sigfig physical constants.

    vibsum =0.0
    for freq in vibs:
        if (freq > 0.1):
            vibsum += freq

    cm2Ha = 219474.6    # cm-1 to Hartree conversion
    Ezpe    = vibsum / (2.0 * cm2Ha)

    # shamelessly swipe code from NWCHEM/src/vib_wrtFreq.F
    eth = 0.0
    hth = 0.0
    xdum = 0.0

    for freq in vibs:
        if (freq > 0.1):
            thetav = freq * (h * c / kgas)    #freqency temperature in Kelvin from cm-1
            if (temperature > 0.0):
                xdum   = math.exp(-thetav/temperature)
            else:
                xdum = 0.0

            xdum = xdum / (1.0 - xdum)
            eth = eth + thetav * (0.5 + xdum)

    eth = eth * Rgas

    # linear boolean is available only after task_freq('???') runs
    # NWChem only writes the flag if molecule is linear

    try:
        is_linear = nwchem.rtdb_get("vib:linear")
    except:
        is_linear = False

    if (is_linear):
        # translational(3/2RT) and rotation(2/2RT) thermal corrections
        eth = eth + 2.5 * Rgas * temperature
    else:
        # translational(3/2RT) and rotation(3/2RT) thermal corrections
        eth = eth + 3.0 * Rgas * temperature

    # Hthermal = eth+pV=eth+RT, since pV=RT
    hth = eth + Rgas * temperature

    debug("ZPE,E,H thermal= %.6f, %.6f, %.6f\n" % (Ezpe,eth,hth))

    Ethermal= eth
    Hthermal= hth

    return False
#

#_______________________________________________________
def E_mp2 ():
# 3 E_(MP2) =  MP2(fc)/6-31G(d)//B3LYP/6-31G(2df,p)

    global Multiplicity
    global Emp2

    say('MP2(fc).')

    send_nwchem_cmd("unset basis:*")
    send_nwchem_cmd("basis noprint ; * library 6-31G* ; end")

    scfcmd = build_SCF_cmd()
    send_nwchem_cmd(scfcmd)

    send_nwchem_cmd("unset mp2:*")
    send_nwchem_cmd("mp2 ; freeze atomic ; end")

    if Multiplicity>1:
        send_nwchem_cmd("unset tce:*")
        send_nwchem_cmd("tce ; scf ; mp2 ; freeze atomic ; end")
        en=nwchem.task_energy("tce")

    else:
        en=nwchem.task_energy("mp2")

    debug('MP2 frozen: en=%.6f\n' % en)

    Emp2 = en

    return False


#_______________________________________________________
def E_ccsdt ():
    '''get CCSDT(fc)/6-31G(d) energy
    '''
    global Multiplicity
    global Eccsdt

    say("CCSD(T).")

    if Multiplicity>1 or is_atom() or is_H2():
        send_nwchem_cmd("unset tce:*")
        send_nwchem_cmd("tce ; ccsd(t) ; freeze atomic ; end")
        en=nwchem.task_energy("tce")
    else:
        send_nwchem_cmd("ccsd ; freeze atomic ; end")
        en=nwchem.task_energy("ccsd(t)")

    debug('CCSD(T): en=%.6f\n' % en)

    Eccsdt = en

    return False

#_______________________________________________________

def E_hf_g3lxp ():
    global Ehfg3lxp

# energy @ HF/G3LargeXP

    send_nwchem_cmd("unset basis:*")
    send_nwchem_cmd('''
        basis spherical
          * library g3mp2largexp
        end''')

    en = nwchem.task_energy("scf")

    Ehfg3lxp = en

    debug("HF/G3LargeXP energy= %f Ha" % (en))

    return False

#_______________________________________________________

def E_mp2_g3lxp ():

    global Emp2g3lxp
# 5 E_(HF/G3LXP) = MP2fc/G3LargeXP

    send_nwchem_cmd("unset basis:*")
    send_nwchem_cmd('''
        basis spherical
          * library g3mp2largexp
        end''')

    send_nwchem_cmd("unset mp2:*")
    send_nwchem_cmd("mp2 ; freeze atomic ; end")

    if Multiplicity>1 or is_atom():
        send_nwchem_cmd("unset tce:*")
        send_nwchem_cmd("tce ; mp2 ; freeze atomic ; end")
        en=nwchem.task_energy("tce")
    else:
        en=nwchem.task_energy("mp2")

    debug('MP(2,fc)/g3mp2large: en=%.6f\n' % en)

    Emp2g3lxp = en

    return False

#_______________________________________________________
def E_hf_augccpvnz (basisset=''):

    basis_cmd = ('basis spherical ; * library %s ; end' % basisset)
    send_nwchem_cmd("unset basis:*")
    send_nwchem_cmd(basis_cmd)

    en = nwchem.task_energy("scf")

    return en

#_______________________________________________________
def E_hf1 ():

    global Ehf1
    say("HF1.")

    basisset = 'g4mp2-aug-cc-pvtz'

    en = E_hf_augccpvnz(basisset)
    Ehf1 = en

    debug("HF/%s energy= %f Ha" % (basisset,en))

    return False

#_______________________________________________________
def E_hf2 ():

    global Ehf2
    say("HF2.")

    basisset = 'g4mp2-aug-cc-pvqz'

    en = E_hf_augccpvnz(basisset)
    Ehf2 = en

    debug("HF/%s energy= %f Ha" % (basisset,en))

    return False

#_______________________________________________________
def E_cbs ():
#   E_(HFlimit)=   extrapolated HF limit
    global Ehf1
    global Ehf2
    global Ecbs

    say('CBS.')

    if abs(Ehf1) > abs(Ehf2):
        Ehf1,Ehf2 = Ehf2,Ehf1

    # Petersen CBS extrapolation
    a = -1.63
    cbs = (Ehf2 - Ehf1 * math.exp(a)) / (1 - math.exp(a))

    Ecbs = cbs

    debug('CBS energy = %.6f' % cbs)

    return False

#______________________________________________________

#_______________________________________________________
def E_hlc ():
#   E_hlc =    High Level Correction

    global Ehlc

    # Correction coefficients in milli Hartrees
    # closed shell
    A   = 0.009472

    # open shell
    Ap  = 0.009769
    B   = 0.003179

    # atoms and atom ions
    C   = 0.009741
    D   = 0.002115

    # single electron pair species
    # e.g.: Li2, Na2, LiNa, BeH2, BeH+
    E   = 0.002379

    say('HLC.')

    nClosed = nwchem.rtdb_get("scf:nclosed")
    nOpen   = nwchem.rtdb_get("scf:nopen")
    nElec   = nwchem.rtdb_get("scf:nelec")
    nFrozen = sum_core_orbitals()

    # According to Curtiss,
    # nBeta = num valence pairs
    # nAlpha = num unpaired or remaining valence electrons
    # subject to the constraint that nAlpha >= nBeta

    nBeta   = nClosed - nFrozen
    nAlpha  = nElec - nFrozen - nClosed

    debug ('\nclosed=%d open=%d frozen=%d nAlpha=%d nBeta=%d\n' % \
            (nClosed,nOpen,nFrozen,nAlpha,nBeta))

    if nAlpha < nBeta:
        nAlpha,nBeta = nBeta,nAlpha

    # test for single (valence) electron pair species
    if (nOpen == 0) and (nClosed - nFrozen) == 1 and \
            (nAlpha == 1) and (nBeta == 1):
        hlc = E

    elif is_atom():
        hlc = -C * (nBeta) - D * (nAlpha-nBeta)

    elif Multiplicity > 1:
        hlc = -Ap * (nBeta) - B * (nAlpha-nBeta)

    else:                   # USUAL CASE: singlet, closed shell
        hlc = -A * (nBeta)

    Ehlc = hlc

    return False

#_______________________________________________________
def E_g4mp2 ():
#   E_(G3(MP2)) =  E_(CCSD(T)) +
#         E_(G3LargeXP) - E_(MP2) +
#         Delta(HFlimit) +
#         E_(SO) +
#         E_(HLC) +
#         E_zpe * scale_factor
    global Ezpe
    global Emp2
    global Eccsdt
    global Ehfg3lxp
    global Emp2g3lxp
    global Ehf1
    global Ehf2
    global Ecbs
    global Ehlc
    global Ethermal
    global Hthermal
    global ESO
    global E0
    global E298
    global H298

    scaled_ZPE  = Ezpe * ZPEScaleFactor

    E0 =    Eccsdt + \
            (Emp2g3lxp - Emp2) + \
            (Ecbs - Ehfg3lxp) + \
            ESO + \
            Ehlc + \
            scaled_ZPE

    E298 =  E0 + (Ethermal - Ezpe)

    H298 =  E0 + (Hthermal - Ezpe)

    return False

#_______________________________________________________
def G4MP2 (charge=0, mult='singlet'):

    Edone = False

    g4mp2_function = [
        optimize,
        E_zpe,
        E_mp2,
        E_ccsdt,
        E_hf_g3lxp,
        E_mp2_g3lxp,
        E_hf1,
        E_hf2,
        E_cbs,
        E_hlc,
        spin_orbit_energy,
        E_g4mp2,
        calc_deltaHf
    ]

    if init_g4mp2(charge,mult) == True:
        sys.exit(0)

    abnormal_end = False

    t0=time.time()

    for i in range(len(g4mp2_function)):
        if g4mp2_function[i]():
            abnormal_end = True
            break

    et=time.time()-t0
    report("\nWall: %.2f seconds" % et)

    if abnormal_end is False:
        report_all()

#_______________________________________________________
#__________________  END G4MP2  ________________________
#_______________________________________________________
