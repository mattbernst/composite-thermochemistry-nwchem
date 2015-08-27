import glob
import hashlib
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
AUKCAL = kCalPerHartree
Rgas = 1.9872041 / 1000.0 / AUKCAL     # atomic units

kT_298_perMol   = (Boltzmann * T298 * Avogadro) / JoulePerKcal / kCalPerHartree

class Gn_common(object):
    def __init__(self, charge=0, multiplicity="singlet", tracing=False,
                 debug=False, integral_memory_cache=3500000000,
                 integral_disk_cache=0, force_c1_symmetry=False):

        self.dhf298      = 0.0
        self.dhf0        = 0.0
        self.force_c1_symmetry = force_c1_symmetry
        multiplets  = ["(null)", "singlet", "doublet", "triplet", "quartet",
                       "quintet", "hextet","septet", "octet"]
        self.integral_memory_cache = integral_memory_cache
        self.integral_disk_cache = integral_disk_cache
        self.nOpen = None
        self.multiplicity_numeric = multiplets.index(multiplicity.lower())
        self.charge = charge
        self.multiplicity = multiplicity
        if multiplicity != "singlet":
            self.hftype = "uhf"
        else:
            self.hftype = "rhf"

        self.tracing = tracing
        self.debug_flag = debug

        self.geohash = self.geometry_hash()
        self.atoms = []

    def say(self, s):
        """Write to stderr console. No implicit newline "\n".

        :param s: message to write
        :type s : str
        """

        if nwchem.ga_nodeid() == 0:
            sys.stderr.write(s)

    def log(self, s):
        """Write to stdout console.

        :param s: message to write
        :type s : str
        """

        if nwchem.ga_nodeid() == 0:
            sys.stdout.write(s + "\n")


    def report(self, s):
        """Write to stderr, stdout.
           Add newline to stderr.

           :param s: message to write
           :type s : str
        """

        #self.say(s + '\n')
        self.log(s)

    def debug(self, s):
        """Write message to stderr if debug_flag is on.

        :param s: message to write
        :type s : str
        """

        if self.debug_flag:
            self.say("DEBUG: {0}\n".format(s))

    def vib_thermo(self, vibs):
        """Handroll the ZPE because NWChem's zpe accumulates
        truncation error from 3 sigfig physical constants.
        """

        AUKCAL  = 627.5093314
        c       = 2.99792458E+10
        h       = 6.62606957E-27
        kgas    = 1.3806488E-16     # cgs units
        Rgas    = 1.9872041/1000.0/AUKCAL     # atomic units

        temperature = 298.15

        vibsum = 0.0
        for freq in vibs:
            if (freq > 0.1):
                vibsum += freq

        cm2Ha = 219474.6    # cm-1 to Hartree conversion
        self.Ezpe = vibsum / (2.0 * cm2Ha)

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

        # linear boolean is available only after task_freq('scf') runs
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
        self.debug("Handrolled E,H thermal= %.6f, %.6f\n" % (eth,hth))

        self.Ethermal = eth
        self.Hthermal = hth


    def geometry_hash(self):
        """Produce a hashed geometry identifier from the geometry in the
        RTDB. This is useful to generate file names for writing and reading
        geometry to/from disk.

        :return: sha1 hex digest of geometry
        :rtype : str
        """

        keys = [nwchem.rtdb_first()]
        while True:
            try:
                keys.append(nwchem.rtdb_next())
            except nwchem.NWChemError:
                break

        ckey = [k for k in keys if "coords" in k and "geometry" in k][0]
        tkey = [k for k in keys if "tags" in k and "geometry" in k][0]
        coords = nwchem.rtdb_get(ckey)
        tags = nwchem.rtdb_get(tkey)
        fused = " ".join([str(c) for c in coords]) + " ".join([str(t) for t in tags])
        result = hashlib.sha1(fused).hexdigest()

        return result

    def is_molecule(self):
        """Determine if this is a molecular system (more than 1 atom)

        :return: True if more than 1 atom, else False
        :rtype : bool
        """

        return len(self.atoms) > 1

    def is_atom(self):
        """Determine if this is an atomic system (just 1 atom)

        :return: True if exactly 1 atom, else False
        :rtype : bool
        """

        return len(self.atoms) == 1

    def send_nwchem_cmd(self, s):
        """Send a command to be parsed as NWChem job input language.

        :param s: command to sent
        :type s : str
        """

        nwchem.input_parse(s)
        self.debug("cmd: [%s]" % s)

    def set_charge(self, charge=0):
        """Set NWChem system charge

        :param charge: total system charge
        :type charge : int
        """

        self.charge = charge
        self.send_nwchem_cmd("charge %s" % charge)

    def element_number(self, element):
        """Get the atomic number associated with a full element name,
        like "lithium". Return 0 if lookup fails.

        :param element: element name
        :type element : str
        :return: atomic number
        :rtype : int
        """

        elementNames = [ 'zero',                        # placeholder
            'HYDROGEN','HELIUM','LITHIUM','BERYLLIUM','BORON','CARBON',
            'NITROGEN','OXYGEN','FLUORINE','NEON','SODIUM','MAGNESIUM',
            'ALUMINIUM','SILICON','PHOSPHORUS','SULFUR','CHLORINE','ARGON',
            'POTASSIUM','CALCIUM','SCANDIUM','TITANIUM','VANADIUM','CHROMIUM',
            'MANGANESE','IRON','COBALT','NICKEL','COPPER','ZINC','GALLIUM',
            'GERMANIUM','ARSENIC','SELENIUM','BROMINE','KRYPTON'
            ]

        number = elementNames.index(element.upper())
        if number == -1:
            number = 0
        return number

    def symbol_number(self, symbol):
        """Get the atomic number associated with an element symbol, like "Li".
        Return 0 if lookup fails.

        :param symbol: element symbol
        :type symbol : str
        :return: atomic number
        :rtype : int
        """

        atomicSymbols = [ 'zero',                       # placeholder
            'H',                               'HE',
            'LI','BE','B' ,'C' ,'N' ,'O' ,'F' ,'NE',
            'NA','MG','AL','SI','P' ,'S' ,'CL','AR',
            'K' ,'CA',
                 'SC','TI','V' ,'CR','MN','FE','CO','NI','CU','ZN',
                      'GA','GE','AS','SE','BR','KR'
            ]

        number = atomicSymbols.index(symbol.upper())
        if number == -1:
            number = 0
        return number

    def atomic_number(self, s):
        """Get the atomic number of an element symbol or name. Try to treat
        the input as a symbol first, then as an element if that fails.

        :param s: element symbol or name
        :type s : str
        :return: atomic number
        :rtype : int
        """

        return self.symbol_number(s) or self.element_number(s)

    def basis_prepare(self, basis, input="", output="",
                      coordinates="spherical", context="scf"):
        """Set up commands to store vectors to a file and/or project or
        load stored vectors for use as initial guess. Also handles
        switching between basis sets.

        :param basis: name of current basis set
        :type basis :str
        :param input: optional name of basis set for vector input
        :type input : str
        :param output: optional name of basis set for vector output
        :type output : str
        :param coordinates: "cartesian" or "spherical", for current basis
        :type coordinates : str
        :param context: "scf" or "dft" vectors setup context
        :type context : str
        """

        def simplename(basis_name):
            name = basis_name[:]
            t = {"-" : "_", "*" : "star", "+" : "plus",
                 "(" : "", ")" : "", "," : "_"}
            for key, value in t.items():
                name = name.replace(key, value)

            return name

        sn = simplename(basis)
        sn_input = simplename(input)
        sn_output = simplename(output)
        basis_cmd = "basis {0} {1} ; * library {2} ; end".format(sn, coordinates, basis)
        self.send_nwchem_cmd(basis_cmd)
        self.send_nwchem_cmd('set "ao basis" {0}'.format(sn))

        if input and output:
            t = "{context}; vectors input project {small} {small}.movecs output {large}.movecs; end"
            vectors = t.format(context=context, small=sn_input, large=sn_output)
        elif input:
            t = "{context}; vectors input {small}.movecs; end"
            vectors = t.format(context=context, small=sn_input)
        elif output:
            t = "{context}; vectors input atomic output {large}.movecs; end"
            vectors = t.format(context=context, large=sn_output)
        #no vector projection or storage; just wanted to set a different basis
        else:
            vectors = ""

        if vectors:
            self.send_nwchem_cmd(vectors)

    def initialize_atoms_list(self):
        """Use NWChem's RTDB geometry to initialize self.atoms,
        e.g. CH3OH gives ['C','H','H','H','O''H']
        """

        # rtdb_get() returns a string if only one atom
        #       or a list of atoms (i.e., tags).
        # Absent type handling,
        # len('CU') is the same as len(['C','U'])
        tags = nwchem.rtdb_get("geometry:geometry:tags")
        if type(tags) == str:
            self.atoms.append(tags)
            self.debug('AtomsList: %s\n' % tags)
        else:
            self.atoms.extend(tags)
            if self.debug_flag:
                atmstr=''
                for atm in tags:
                    atmstr += ' '+atm
                self.debug('AtomsList: %s\n' % atmstr)

        self.debug('NumAtoms={0}\n'.format(len(self.atoms)))

    def build_SCF_cmd(self):
        """Prepare SCF block with multiplicity-appropriate choice of
        HF form.

        :return: SCF control block
        :rtype : str
        """

        memory_cache_words = self.integral_memory_cache / 8
        disk_cache_words = self.integral_disk_cache / 8

        tpl = "scf ; semidirect memsize {0} filesize {1}; {2} ; {3} ; end"
        block = tpl.format(memory_cache_words, disk_cache_words,
                           self.multiplicity, self.hftype)
        return block

    def reset_symmetry(self):
        """Reload geometry and force symmetry down if TCE must be used.
        """

        if self.force_c1_symmetry:
            return

        xyzs = glob.glob(self.geohash + "*.xyz")
        xyzs.sort()
        geofile = xyzs[-1]

        #TODO: use largest Abelian subgroup symmetries instead of c1
        if self.multiplicity != "singlet":
            symmetry_block = "symmetry c1;"
            geoblock = "geometry units angstroms print xyz; {0} load {1}; end"
            self.send_nwchem_cmd(geoblock.format(symmetry_block, geofile))

    def report_dHf(self):
        """Report change in heat of formation going from 0 K to 298 K.
        """

        heatsOfFormation = [
        #("\n"),
        ("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"),
        ("          HEAT OF FORMATION   (0K): % 10.2f kCal/mol" % self.dhf0),
        ("          HEAT OF FORMATION (298K): % 10.2f kCal/mol" % self.dhf298),
        ("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        ]


        for line in heatsOfFormation:
            self.report(line)

    def spin_orbit_energy(self):
        """Get spin orbit energy correction according to nature of system and charge.

        :return: spin orbit energy correction
        :rtype : float
        """

        if self.is_molecule():   # no spin orbit corrections for molecules
            correction = 0.0
        else:       # It's an atom
            atom = self.atoms[0]
            correction = self.E_spin_orbit(self.atomic_number(atom),
                                           self.charge)

        return correction

    def E_spin_orbit(self, atomic_number, charge):
        """EspinOrbit tabulates the spin orbit energies
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
        function E_spin_orbit returns values in Hartrees.

        :param atomic_number: atomic number of atomic species
        :type atomic_number : int
        :param charge: charge on atomic species
        :type charge : int
        :return: energy in Hartrees
        :rtype : float
        """

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
        espin = ESpinOrbit[atomic_number][ion] * milliHa_to_Ha
        return espin

    def atomic_DHF (self, elementNum):
        """Get atomic heats of formation at 0 K and 298 K, in
        kcal/mol.

        :param elementNum: atomic number
        :type elementNum : int
        :return: atomic heats of formation
        :rtype : tuple
        """

        #atom [dHf(0), dHf(298)]  in kcal/mol
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
            [0.0,0.0],[0.0,0.0],
            [0.0,0.0],[0.0,0.0],
            [0.0,0.0],[0.0,0.0],
            [0.0,0.0],[0.0,0.0],
            [65.00,   65.00   ],  # 31  Gallium
            [88.91,   88.91   ],  # 32  Germanium
            [73.90,   72.42   ],  # 33  Arsenic
            [55.76,   54.27   ],  # 34  Selenium
            [26.74,   28.18   ],  # 35  Bromine
            [0.0,     0.0     ],  # 36  Krypton
        ]

        self.debug('atomic_DHF: elementNum=%d' % elementNum)
        try:
            result = atomDHF[elementNum]
            self.debug('atomic_DHF: E,H=%.2f,%.2f' % (result[0], result[1]))
        except IndexError:
            self.debug('atomic_DHF: error: element %d not in table?' % elementNum)
            result = (0.0, 0.0)

        return result

    def atom_core_orbitals(self, atomicNumber, convention="gamess"):
        """This replicates the core electron pair lookup table in
        src/geom/geom_core.F

        :param atomicNumber: atomic number of an atom
        :type atomicNumber : int
        :param convention: "gamess" or "nwchem" table
        :type convention : str
        :return: number of core orbitals for atom
        :rtype : int
        """

        #NWChem version 6.5
        nwchemCoreOrbitals = [0,        # zero index place holder
            0,                                                 0,
            1, 1,                               1, 1, 1, 1, 1, 1,
            5, 5,                               5, 5, 5, 5, 5, 5,
            9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
            18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,
            27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,
            27,27,27,27,27,27,27,27,27,27,27,27,27,27,43,43,43,43,
            43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,
            43,43,43,43,]

        # GAMESS version 12
        gamessCoreOrbitals = [0,        # zero index place holder
            0,                                                  0,
            1,  1,                               1, 1, 1, 1, 1, 1,
            5,  5,                               5, 5, 5, 5, 5, 5,
            9,  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,14,14,14,14,14,14,
            23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,
            27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,
            39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,
            34,34,34,34,34,34,34,34,34,34,39,39,39,39,39,39,
            43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,50]

        cmap = {"gamess" : gamessCoreOrbitals, "nwchem" : nwchemCoreOrbitals}
        try:
            nCoreOrbitals = cmap[convention.lower()]
        except KeyError:
            raise ValueError("Uknown core orbital convention " + repr(convention))

        if atomicNumber <= len(nCoreOrbitals):
            n = nCoreOrbitals[atomicNumber]
        else:
            n = 0

        return n

    def sum_core_orbitals(self, convention="gamess"):
        """Sum the total number of frozen core orbitals in a system.
        nFrozen isn't consistently logged to RTDB by
        Tensor Contraction Engine methods, so do the work ourselves.

        :param convention: orbital freeze convention, "gamess" or "nwchem"
        :type convention : str
        :return: sum of frozen core orbitals
        :rtype : int
        """

        total = sum([self.atom_core_orbitals(self.atomic_number(a),
                                             convention=convention)
                     for a in self.atoms])
        return total


class G4_mp2(Gn_common):
    """
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
    """

    def __init__(self, *args, **kw):
        super(G4_mp2, self).__init__(*args, **kw)

        self.correlated_basis = [("6-31G*", "cartesian"),
                                 ("g3mp2largexp", "spherical")]
        self.cbs_basis = [("g4mp2-aug-cc-pvtz", "spherical"),
                          ("g4mp2-aug-cc-pvqz", "spherical")]

        #Zero Point Energy scale factor for B3LYP/6-31G(2df,p)
        self.ZPEScaleFactor = 0.9854    # Curtiss scale factor for Gaussian 09
        #self.ZPEScaleFactor = 0.9798     # Truhlar scale Factor for NWChem 6.5

        self.Ezpe        = 0.0
        self.Emp2        = 0.0
        self.Eccsdt      = 0.0
        self.Ehfg3lxp    = 0.0
        self.Emp2g3lxp   = 0.0
        self.Ehf1        = 0.0
        self.Ehf2        = 0.0
        self.Ecbs        = 0.0
        self.Ehlc        = 0.0
        self.Ethermal    = 0.0
        self.Hthermal    = 0.0
        self.E0          = 0.0
        self.E298        = 0.0
        self.H298        = 0.0

        # valence electron variables
        self.nAlpha      = 0
        self.nBeta       = 0
        self.nFrozen     = 0

    def report_summary(self):
        """Report results in GAMESS G3(MP2) output format for easy comparison.
        """

        Szpe = self.Ezpe * self.ZPEScaleFactor
        dMP2 = self.Emp2g3lxp - self.Emp2
        dHF  = self.Ecbs - self.Ehfg3lxp

        summary = [
            ("\n    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~NWChem6.5"),
            ("                  SUMMARY OF G4(MP2) CALCULATIONS"),
            ("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"),
            ("    B3LYP/6-31G(2df,p)= % 12.6f   HF/maug-cc-p(T+d)Z= % 12.6f" % (self.Eb3lyp, self.Ehf1)),
            ("    HF/CBS            = % 12.6f   HF/maug-cc-p(Q+d)Z= % 12.6f" % (self.Ecbs, self.Ehf2)),
            ("    MP2/6-31G(d)      = % 12.6f   CCSD(T)/6-31G(d)  = % 12.6f" % (self.Emp2, self.Eccsdt)),
            ("    HF/G3MP2LARGEXP   = % 12.6f   MP2/G3MP2LARGEXP  = % 12.6f" % (self.Ehfg3lxp, self.Emp2g3lxp)),
            ("    DE(MP2)           = % 12.6f   DE(HF)            = % 12.6f" % (dMP2, dHF)),
            ("    ZPE(B3LYP)        = % 12.6f   ZPE SCALE FACTOR  = % 12.6f" % (Szpe, self.ZPEScaleFactor)),
            ("    HLC               = % 12.6f   FREE ENERGY       = % 12.6f" % (self.Ehlc, 0.0)),
            ("    THERMAL ENERGY    = % 12.6f   THERMAL ENTHALPY  = % 12.6f" % (self.Ethermal, self.Hthermal)),
            ("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"),
            ("    E(G4(MP2)) @ 0K   = % 12.6f   E(G4(MP2)) @298K  = % 12.6f" % (self.E0, self.E298)),
            ("    H(G4(MP2))        = % 12.6f   G(G4(MP2))        = % 12.6f" % (self.H298, 0.0)),
            ("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"),
        ]

        for line in summary:
            self.report(line)

    def report_all(self):
        self.report_summary ()
        self.report_dHf ()

    def E0_atom(self, elementNum):
        """List of precalculated atomic G4(MP2) energies at 0K
        Returns E(0K), E(298.15K) tuple.

        :param elementNum: atomic number
        :type elementNum : int
        :return: energy at 0 K and 298 K
        :rtype : tuple
        """

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
        eth = (5.0/2) * kT_298_perMol

        try:
            e0 = e0_g4mp2[elementNum]
            e298 = e0 + eth
            result = (e0, e298)
        except IndexError:
            result = (0.0, 0.0)

        return result

    def calc_deltaHf(self):
        """Calculate heat of formation at 0K and 298K.
        """

        sum_atoms_E0     = 0.0
        sum_atoms_E298   = 0.0
        sum_atoms_dhf0   = 0.0
        sum_atoms_dhf298 = 0.0

        if self.is_molecule():
            for atom in self.atoms:
                e0, e298 = self.E0_atom(self.atomic_number(atom))
                sum_atoms_E0 += e0
                sum_atoms_E298 += e298

                d0, d298 = self.atomic_DHF(self.atomic_number(atom))
                sum_atoms_dhf0       += d0
                sum_atoms_dhf298     += d298
                self.debug('sumDHF0,sumDHF298 = %.2f,%.2f' % (sum_atoms_dhf0, sum_atoms_dhf298))
        else:
            return False

        self.dhf0 = (self.E0 - sum_atoms_E0) * kCalPerHartree + sum_atoms_dhf0
        self.dhf298 = (self.H298 - sum_atoms_E298) * kCalPerHartree + sum_atoms_dhf298
        self.debug('dhf0,dhf298 = %.2f,%.2f' % (self.dhf0, self.dhf298))

    def init_g4mp2(self):
        """Say hello.
        """

        title = nwchem.rtdb_get("title")
        self.say(" %s -- NWChem G4(MP2) Composite Method\n" % (title))

    def prepare_scf_vectors(self):
        """Set up converged scf vectors before first correlated steps so they
        can be used to initialize later calculations.
        """

        self.basis_prepare(self.correlated_basis[0][0],
                           output=self.correlated_basis[0][0],
                           coordinates=self.correlated_basis[0][1])

        self.send_nwchem_cmd(self.build_SCF_cmd())
        nwchem.task_energy("scf")

    def optimize(self):
        """# 1 optimize  B3LYP/6-31G(2df,p)

        """

        self.say('optimize.')

        self.send_nwchem_cmd("basis noprint ; * library 6-31G(2df,p) ; end")
        scfcmd = self.build_SCF_cmd()
        self.send_nwchem_cmd(scfcmd)

        # canonical B3LYP spec
        # GAMESS-US b3lyp uses VWN_5
        # Gaussian uses VWN_3
        # NWChem uses something entirely different

        b3lyp_GAMESS = 'xc HFexch 0.2 slater 0.8 becke88 nonlocal 0.72 vwn_5 0.19 lyp 0.81'
        b3lyp_Gaussian = 'xc HFexch 0.2 slater 0.8 becke88 nonlocal 0.72 vwn_3 0.19 lyp 0.81'
        b3lyp_NWChem = 'xc b3lyp'

        blips = b3lyp_Gaussian
        memory_cache_words = self.integral_memory_cache / 8
        disk_cache_words = self.integral_disk_cache / 8
        mem = "semidirect memsize {0} filesize {1}".format(memory_cache_words,
                                                           disk_cache_words)

        if self.multiplicity != "singlet":
            self.send_nwchem_cmd('dft ; odft ; mult %d ; %s ; %s ; end' % (self.multiplicity_numeric, blips, mem))
        else:
            self.send_nwchem_cmd('dft ; %s ; %s ; end' % (blips, mem))

        # fetch and copy atom names list (tags) which enumerates atoms.
        # only available _after_ SCF statement

        self.initialize_atoms_list()
        self.send_nwchem_cmd("driver; maxiter 99; xyz {0}; end".format(self.geohash))
        self.send_nwchem_cmd("scf; maxiter 99; end")

        # optimize the geometry, ignore energy and gradient results
        if self.is_atom():
            en = nwchem.task_energy("dft")

        else:
            self.debug("task_optimize(dft)")
            en, grad = nwchem.task_optimize("dft")

        self.Eb3lyp = en

        #self.report("debug: HF/6-31G(2df,p) SCF:energy = %f Ha" % (en))

    def E_zpe(self):
        """Run hessian on equilibrium geometry, get zero point energy.
                                                                              n
        Note: linear tri-atomics and larger give high ZPE values in
        NWChem @ HF/6-31G*
        """

        self.say('ZPE.')

        temperature = T298

        if self.is_atom():
            self.Ezpe = 0.0
            self.Ethermal = 1.5 * Rgas * temperature # 3/2 * RT
            self.Hthermal = self.Ethermal + (Rgas * temperature)
            return False

        # run hessian on equilibrium geometry
        # ignore ZPE, calculate it from vibrations list
        zpe, vibs, intens = nwchem.task_freq("dft")
        self.vib_thermo(vibs)

    def E_mp2(self):
        """Calculate the MP2 energy at the B3LYP-optimized geometry.

        # 3 E_(MP2) =  MP2(fc)/6-31G(d)//B3LYP/6-31G(2df,p)

        :return: failure code (True for failure, False for success)
        :rtype : bool
        """

        self.say('MP2(fc).')

        self.basis_prepare(self.correlated_basis[0][0],
                           input=self.correlated_basis[0][0],
                           coordinates=self.correlated_basis[0][1])

        scfcmd = self.build_SCF_cmd()
        self.send_nwchem_cmd(scfcmd)

        self.send_nwchem_cmd("unset mp2:*")
        self.send_nwchem_cmd("mp2 ; freeze atomic ; end")

        try:
            if self.multiplicity != "singlet":
                self.send_nwchem_cmd("unset tce:*")
                self.send_nwchem_cmd("tce ; scf ; mp2 ; freeze atomic ; end")

                en = nwchem.task_energy("tce")
            else:
                en = nwchem.task_energy("mp2")

            self.debug('MP2 frozen: en=%.6f\n' % en)
        except:
            self.report("FAILED: MP2(fc)/6-31G(2df,p) energy")
            return True
        else:
            self.Emp2 = en

        return False

    def E_ccsdt(self):
        """# 4 E_(ccsd(t)) =  CCSD(fc,T)/6-31G(d)

        Get CCSDT(fc)/6-31G(d) energy
        """

        self.say("CCSD(T).")

        if self.multiplicity != "singlet" or self.is_atom():
            self.send_nwchem_cmd("unset tce:*")
            self.send_nwchem_cmd("tce ; ccsd(t) ; freeze atomic ; end")
            en = nwchem.task_energy("tce")
        else:
            self.send_nwchem_cmd("ccsd ; freeze atomic ; end")
            en = nwchem.task_energy("ccsd(t)")

            self.debug('CCSD(T): en=%.6f\n' % en)

        self.Eccsdt = en

    def E_hf_g3lxp(self):
        """# 5 E_(HF/G3LXP) = HF/G3LargeXP
        """

        self.basis_prepare(self.correlated_basis[1][0],
                           input=self.correlated_basis[0][0],
                           output=self.correlated_basis[1][0],
                           coordinates=self.correlated_basis[1][1])

        en = nwchem.task_energy("scf")
        self.Ehfg3lxp = en
        self.debug("HF/G3LargeXP SCF:energy = %f Ha" % (en))

    def E_mp2_g3lxp(self):
        """# 5 E_(HF/G3LXP) = MP2fc/G3LargeXP
        """

        self.basis_prepare(self.correlated_basis[1][0],
                           coordinates=self.correlated_basis[1][1])

        self.send_nwchem_cmd("unset mp2:*")
        self.send_nwchem_cmd("mp2 ; freeze atomic ; end")

        if self.multiplicity != "singlet" or self.is_atom():
            self.send_nwchem_cmd("unset tce:*")
            self.send_nwchem_cmd("tce ; mp2 ; freeze atomic ; end")
            en = nwchem.task_energy("tce")
        else:
            en = nwchem.task_energy("mp2")

        self.debug('MP(2,fc)/g3mp2large: en=%.6f\n' % en)

        self.Emp2g3lxp = en

    def E_hf1(self):
        """Use g4mp2-aug-cc-pvtz basis set to get first HF energy.
        """

        self.say("HF1.")

        self.basis_prepare(self.cbs_basis[0][0],
                           input=self.correlated_basis[0][0],
                           output=self.cbs_basis[0][0],
                           coordinates=self.cbs_basis[0][1])

        en = nwchem.task_energy("scf")

        self.Ehf1 = en
        self.debug("HF/%s: energy = %f Ha" % (self.cbs_basis[0][0], en))

    def E_hf2(self):
        """Use g4mp2-aug-cc-pvqz to get second
        HF energy.

        :return: failure code (True for failure, False for success)
        :rtype : bool
        """

        self.say("HF2.")

        self.basis_prepare(self.cbs_basis[1][0],
                           input=self.cbs_basis[0][0],
                           output=self.cbs_basis[1][0],
                           coordinates=self.cbs_basis[1][1])

        en = nwchem.task_energy("scf")

        self.Ehf2 = en
        self.debug("HF/%s: energy = %f Ha" % (self.cbs_basis[0][0], en))

    def E_cbs(self):
        """E_(HFlimit) = extrapolated HF limit

        :return: failure code (True for failure, False for success)
        :rtype : bool
        """

        self.say('CBS.')

        use_Petersen = True    # use Petersen CBS extrapolation

        #TODO: why abs() here?
        if abs(self.Ehf1) > abs(self.Ehf2):
            self.Ehf1, self.Ehf2 = self.Ehf2, self.Ehf1

        if use_Petersen:    # petersen CBS extrapolation
            a = -1.63
            cbs = (self.Ehf2 - self.Ehf1 * math.exp(a)) / (1 - math.exp(a))
        else:               # Truhlar CBS extrapolation
            a = 3.4
            k1 = math.pow(3,a) / (math.pow(3,a) - math.pow(2,a))
            k2 = math.pow(2,a) / (math.pow(3,a) - math.pow(2,a))
            cbs = k1 * self.Ehf2 - k2 * self.Ehf1

        self.Ecbs = cbs

        self.debug('CBS energy = %.6f' % cbs)

        return False

    def E_hlc(self):
        """E_hlc =    High Level Correction
        """

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

        self.say('HLC.')

        nClosed = nwchem.rtdb_get("scf:nclosed")
        nOpen   = nwchem.rtdb_get("scf:nopen")
        nElec   = nwchem.rtdb_get("scf:nelec")
        self.nFrozen = self.sum_core_orbitals()

        # According to Curtiss,
        # nBeta = num valence pairs
        # nAlpha = num unpaired or remaining valence electrons
        # subject to the constraint that nAlpha >= nBeta

        self.nBeta   = nClosed - self.nFrozen
        self.nAlpha  = nElec - self.nFrozen - nClosed

        self.debug('\nclosed=%d open=%d frozen=%d nAlpha=%d nBeta=%d\n' % \
                   (nClosed, nOpen, self.nFrozen, self.nAlpha, self.nBeta))

        if self.nAlpha < self.nBeta:
            self.nAlpha, self.nBeta = self.nBeta, self.nAlpha
            #self.say('** a<->b swap: nAlpha=%d nBeta=%d\n' % (self.nAlpha,self.nBeta))

        # test for single (valence) electron pair species
        if (nOpen == 0) and (nClosed - self.nFrozen) == 1 and \
                (self.nAlpha == 1) and (self.nBeta == 1):
            hlc = E

        elif self.is_atom():
            hlc = -C * (self.nBeta) - D * (self.nAlpha - self.nBeta)

        elif self.multiplicity != "singlet":
            hlc = -Ap * (self.nBeta) - B * (self.nAlpha - self.nBeta)

        else:                   # USUAL CASE: singlet, closed shell
            hlc = -A * (self.nBeta)

        self.Ehlc = hlc

    def E_g4mp2(self):
        """
        E_(G3(MP2)) =  E_(CCSD(T)) +
        E_(G3LargeXP) - E_(MP2) +
        Delta(HFlimit) +
        E_(SO) +
        E_(HLC) +
        E_zpe * scale_factor

        :return: failure code (True for failure, False for success)
        :rtype : bool
        """

        scaled_ZPE  = self.Ezpe * self.ZPEScaleFactor

        self.E0 = self.Eccsdt + \
                  (self.Emp2g3lxp - self.Emp2) + \
                  (self.Ecbs - self.Ehfg3lxp) + \
                  self.spin_orbit_energy() + \
                  self.Ehlc + \
                  scaled_ZPE

        self.E298 =  self.E0 + (self.Ethermal - self.Ezpe)

        self.H298 =  self.E0 + (self.Hthermal - self.Ezpe)

        return False

    def run(self):
        """Calculate G4MP2 energy for a system that has already been prepared.

        """

        self.Edone = False

        g4mp2_function = [
            self.optimize,
            self.E_zpe,
            self.prepare_scf_vectors,
            self.E_hf_g3lxp,
            self.E_hf1,
            self.E_hf2,
            self.reset_symmetry,
            self.E_mp2,
            self.E_ccsdt,
            self.E_mp2_g3lxp,
            self.E_cbs,
            self.E_hlc,
            self.spin_orbit_energy,
            self.E_g4mp2,
            self.calc_deltaHf
        ]

        self.init_g4mp2()
        t0=time.time()

        for i in range(len(g4mp2_function)):
            g4mp2_function[i]()

        et=time.time()-t0
        self.report("\nWall: %.2f seconds" % et)

        self.report_all()

class G3_mp2(Gn_common):
    """
     G3(MP2) THERMOCHEMICAL METHOD FOR NWCHEM
    Daniel R. Haney 2015

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
    permit CCSD(T) instead of QCISD(T) in the base energy term.
    This results in slightly higher total energies than reported by
    Curtiss in the G2/97 test set, faster run times, and lower mean
    average deviations in isodesmic reaction energy calculations.
    """

    def __init__(self, *args, **kw):
        if kw.pop("use_qcisdt", False):
            self.highest_correlated = "qcisd(t)"
        else:
            self.highest_correlated = "ccsd(t)"

        super(G3_mp2, self).__init__(*args, **kw)

        self.ZPEScaleFactor = 0.8929


        self.Ezpe = 0.0
        self.Emp2full = 0.0
        self.Emp2frozen = 0.0
        self.Ecc = 0.0
        self.Eg3mp2large = 0.0
        self.Ehlc = 0.0
        self.Ethermal = 0.0
        self.Hthermal = 0.0
        self.E0 = 0.0
        self.E298 = 0.0
        self.H298 = 0.0

        # valence electron variables
        self.nAlpha = 0
        self.nBeta = 0
        self.nFrozen = 0

    def HF_optimize(self):
        '''
           HF/6-31G(d) optimization
        '''

        self.say('optimize.')
        self.send_nwchem_cmd("scf; maxiter 99; end")
        self.send_nwchem_cmd("driver; maxiter 99; xyz {0}; end".format(self.geohash))

        #self.send_nwchem_cmd("basis noprint ; * library 6-31G* ; end")
        self.basis_prepare("6-31G*", output="6-31G*", coordinates="cartesian")
        scfcmd = self.build_SCF_cmd()
        self.send_nwchem_cmd(scfcmd)

        # fetch and copy atom names list (tags) which enumerates atoms.
        # only available _after_ SCF statement

        self.initialize_atoms_list()

        # optimize the geometry, ignore energy and gradient results
        if self.is_atom():
            en = nwchem.task_energy("scf")
        else:
            en, grad = nwchem.task_optimize("scf")

    def HF_zpe(self):
        ''' run hessian on equilibrium geometry
            get zero point energy.
                                                                              n
           note: linear tri-atomics and larger
           give high ZPE values in NWChem @ HF/6-31G*
        '''

        temperature = T298

        self.say("zpe.")

        if self.is_atom():
            self.Ezpe = 0.0
            self.Ethermal = 1.5 * Rgas * temperature
            self.Hthermal = self.Ethermal + (Rgas * temperature)
            return

        # run hessian on equilibrium geometry
        # ignore ZPE, calculate it from vibrations list
        zpe, vibs, intens = nwchem.task_freq("scf")
        self.vib_thermo(vibs)

    def MP2_optimize(self):
        '''Optimize geometry at MP2(full)/6-31G(d)
        '''

        self.say('MP2 optimize.')

        if self.multiplicity != "singlet":
            self.send_nwchem_cmd("tce ; scf ; mp2 ; end")
            if self.is_atom():
                en = nwchem.task_energy("tce")
            else:
                en, grad = nwchem.task_optimize("tce")
        else:
            if self.is_atom():
                en = nwchem.task_energy("mp2")
            else:
                en, grad = nwchem.task_optimize("mp2")

        self.debug('optimize: MP(2,full)/6-31G*= %.6f\n' % en)

        self.Emp2full = en

    def MP2_frozen(self):
        # MP2(fc)/6-31G* single point energy

        self.say('MP2(frozen).')

        scfcmd = self.build_SCF_cmd()
        self.send_nwchem_cmd(scfcmd)

        self.send_nwchem_cmd("unset mp2:*")
        self.send_nwchem_cmd("mp2 ; freeze atomic ; end")

        if self.multiplicity != "singlet":
            self.send_nwchem_cmd("unset tce:*")
            self.send_nwchem_cmd("tce ; scf ; mp2 ; freeze atomic ; end")
            en = nwchem.task_energy("tce")

        else:
            en = nwchem.task_energy("mp2")

        self.debug('MP2 frozen: en=%.6f\n' % en)
        self.Emp2frozen = en

    def ccsdt_qcisdt_frozen(self):
        #highest level correlated calculation
        if self.highest_correlated == "qcisd(t)":
            self.say("QCISD(T).")
            tce = "tce ; qcisd(t) ; freeze atomic ; end"
            nwchem.input_parse(tce)
            en = nwchem.task_energy("tce")
            self.Ecc = en

        elif self.highest_correlated == "ccsd(t)":
            self.say("CCSD(T).")
            if self.multiplicity != "singlet" or self.is_atom():
                self.send_nwchem_cmd("unset tce:*")
                self.send_nwchem_cmd("tce ; ccsd(t) ; freeze atomic ; end")
                en = nwchem.task_energy("tce")
            else:
                self.send_nwchem_cmd("ccsd ; freeze atomic ; end")
                en = nwchem.task_energy("ccsd(t)")

            self.debug(' CCSD(T): en=%.6f\n' % en)
            self.Ecc = en

        else:
            raise ValueError("Unknown correlation treatment {}".format(repr(self.highest_correlated)))

    def MP2_g3mp2large(self):
        '''get MP2(fc)/G3MP2large single point energy
        '''

        self.say('GMP2large.')

        self.basis_prepare("g3mp2large", input="6-31G*", coordinates="spherical")

        self.send_nwchem_cmd("unset mp2:*")

        if self.multiplicity != "singlet" or self.is_atom():
            self.send_nwchem_cmd("unset tce:*")
            self.send_nwchem_cmd("tce ; mp2 ; end")
            en = nwchem.task_energy("tce")
        else:
            en = nwchem.task_energy("mp2")

        self.debug(' g3mp2large: en=%.6f\n' % en)

        self.Eg3mp2large = en

    def HLC_generic(self, A, B, C, D):
        '''calculate High Level Correction term
            from alpha and beta VALENCE electron count.
        '''

        self.say('HLC.')

        nClosed = nwchem.rtdb_get("scf:nclosed")
        nOpen = nwchem.rtdb_get("scf:nopen")
        nElec = nwchem.rtdb_get("scf:nelec")
        nFrozen = self.sum_core_orbitals()

        # According to Curtiss,
        # nBeta = num valence pairs
        # nAlpha = num unpaired or remaining valence electrons
        # subject to the constraint that nAlpha >= nBeta

        nBeta = nClosed - nFrozen

        nAlpha = nElec - (nFrozen * 2) - nBeta

        self.debug('\nclosed=%d open=%d frozen=%d nAlpha=%d nBeta=%d\n' %
            (nClosed,nOpen,nFrozen,nAlpha,nBeta))

        if nAlpha < nBeta:
            nAlpha, nBeta = nBeta, nAlpha

        if self.is_molecule():
            self.Ehlc = -(A * nBeta) - B * (nAlpha - nBeta)
        else:  # it's an atom, calc is different
            self.Ehlc = -(C * nBeta) - D * (nAlpha - nBeta)

    def HLC_qcisdt(self):
        ''' empirical correction coefficients
        for Curtiss original G3(MP2) method
        '''

        A = 0.009279
        B = 0.004471
        C = 0.009345
        D = 0.002021

        return self.HLC_generic(A, B, C, D)

    def HLC_ccsdt(self):
        ''' empirical correction coefficients
        for later G3(MP2,CCSDT) method
        '''

        A = 0.009170
        B = 0.004455
        C = 0.009155
        D = 0.001947

        return self.HLC_generic(A, B, C, D)

    def choose_correlated(self, key, map):
        #return setting based on correlation scheme, alert on mismatch
        try:
            return map[key]
        except KeyError:
            raise ValueError("Unknown correlation treatment {}".format(repr(key)))

    def HLC(self):
        hlc = self.choose_correlated(self.highest_correlated,
                                     {"qcisd(t)" : self.HLC_qcisdt,
                                      "ccsd(t)" : self.HLC_ccsdt})

        v = hlc()
        return v

    def init_g3mp2(self):
        '''say hello
        '''

        ccstring = self.choose_correlated(self.highest_correlated,
                                          {"qcisd(t)" : "QCISD(T)",
                                           "ccsd(t)" : "CCSD(T)"})

        title = nwchem.rtdb_get("title")
        self.say(" %s -- NWChem G3(MP2,%s) Composite Method\n" % (title, ccstring))

    def calc_total_energies(self):
        '''G3MP2 step 6
            calculate the E(G3MP2)@0K, @298K energies
        '''


        self.E0 = self.Ecc + \
            (self.Eg3mp2large - self.Emp2frozen) + \
            (self.Ezpe * self.ZPEScaleFactor) + \
            self.Ehlc

        if self.is_atom():
            self.E0 += self.spin_orbit_energy()

        self.E298 = self.E0 + (self.Ethermal - self.Ezpe)
        self.H298 = self.E0 + (self.Hthermal - self.Ezpe)  # kT_298_perMol

    def E0_atom_qcisdt(self, elementNum):
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


    def E0_atom_ccsdt(self, elementNum):

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


    def E0_atom(self, elementNum):
        fn = self.choose_correlated(self.highest_correlated,
                                    {"qcisd(t)" : self.E0_atom_qcisdt,
                                     "ccsd(t)" : self.E0_atom_ccsdt})

        return fn(elementNum)



    def calc_deltaHf(self):
        """calculate heat of formation at 0K and 298K
        """

        sum_atoms_E0 = 0.0
        sum_atoms_H298 = 0.0
        sum_atoms_dhf0 = 0.0
        sum_atoms_dhf298 = 0.0

        for atom in self.atoms:
            e0, h298 = self.E0_atom(self.atomic_number(atom))
            sum_atoms_E0 += e0
            sum_atoms_H298 += h298

            d0, d298 = self.atomic_DHF(self.atomic_number(atom))
            sum_atoms_dhf0 += d0
            sum_atoms_dhf298 += d298
            self.debug('sumDHF0,sumDHF298 = %.2f,%.2f' %
                (sum_atoms_dhf0, sum_atoms_dhf298))

        self.dhf0 = (self.E0 - sum_atoms_E0) * kCalPerHartree + sum_atoms_dhf0

        self.dhf298 = (self.H298 - sum_atoms_H298) * kCalPerHartree + sum_atoms_dhf298

        self.debug('dhf0,dhf298 = %.2f,%.2f' % (self.dhf0, self.dhf298))

    def report_summary(self):
        '''Report results in GAMESS G3(MP2) output format
           for easy comparison.
           log() normally redirects to log file
           say() appears in terminal session
        '''

        Szpe = self.Ezpe * self.ZPEScaleFactor
        deltaMP2 = self.Eg3mp2large - self.Emp2frozen

        ccstring = self.choose_correlated(self.highest_correlated,
                                          {"qcisd(t)" : "QCISDT ",
                                           "ccsd(t)" : "CCSD(T)"})


        summary = [
            ("\n"),
            ("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~NWChem6.5"),
            ("               SUMMARY OF G3(MP2,%s) CALCULATIONS              " %
             ccstring),
            ("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"),
            ("    MP2/6-31G(d)    = % 12.6f   %s/6-31G(d) = % 12.6f" %
             (self.Emp2frozen, ccstring, self.Ecc)),
            ("    MP2/G3MP2large  = % 12.6f   delta(MP2)       = % 12.6f" %
             (self.Eg3mp2large, deltaMP2)),
            ("    ZPE(HF/6-31G(d))= % 12.6f   ZPE Scale Factor = % 12.6f" %
             (Szpe, self.ZPEScaleFactor)),
            ("    HLC             = % 12.6f   Free Energy      = % 12.6f" %
             (self.Ehlc, 0.0)),
            ("    Thermal Energy  = % 12.6f   Thermal Enthalpy = % 12.6f" %
             (self.Ethermal, self.Hthermal)),
            ("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"),
            ("    E(G3(MP2)) @ 0K = % 12.6f   E(G3(MP2)) @298K = % 12.6f" %
             (self.E0, self.E298)),
            ("    H(G3(MP2))      = % 12.6f   G(G3(MP2))       = % 12.6f" %
             (self.H298, 0.0)),
            ("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        ]

        for line in summary:
            self.report(line)


    def reportAll(self):
        self.report_summary()
        self.report_dHf()


    def run(self):
        g3mp2_function = [
            self.HF_optimize,
            self.HF_zpe,
            self.reset_symmetry,
            self.MP2_optimize,
            self.MP2_frozen,
            self.ccsdt_qcisdt_frozen,
            self.MP2_g3mp2large,
            self.HLC,
            self.calc_total_energies,
            self.calc_deltaHf
        ]

        self.init_g3mp2()

        t0=time.time()

        for i in range(len(g3mp2_function)):
            g3mp2_function[i]()

        et=time.time()-t0
        self.report("\nWall: %.2f seconds" % et)

        self.reportAll()
