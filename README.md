# composite-thermochemistry-nwchem
Composite thermochemical models implemented with NWChem's Python functionality

##Installation

###Requirements

* A supported NWChem platform: this code has been tested under Linux and OS X. In theory, NWChem runs under Windows but the authors have no experience with it. The Windows port appears to be minimally supported. Windows users without a lot of cross-platform compilation experience should consider a Linux virtual machine instead.
* A _recent_, _Python enabled_ build of NWChem: ideally after svn revision 27160 but at least following 27110. Use any release after 6.6 from the [NWChem download page](http://www.nwchem-sw.org/index.php/Download).
* The basis set data for the models (included), which must be copied into NWChem's library directory.

###Installing

After you have compiled NWChem with embedded Python support and installed it, copy everything from `basis_data` to NWChem's basis library directory. The library directory is `$NWCHEM_TOP/src/basis/libraries/`.

If you find your jobs crashing with ga errors before physical RAM is exhausted, make sure you have set `memory_total` to a sufficiently high value in your `~/.nwchemrc`.
 
Here is an example `.nwchemrc`. Replace `/opt/science/nwchem/current/` with the path to your actual NWChem build.

```
memory_total 1000000000
nwchem_basis_library /opt/science/nwchem/current/src/basis/libraries/
nwchem_nwpw_library /opt/science/nwchem/current/src/nwpw/libraryps/
ffield amber
amber_1 /opt/science/nwchem/current/src/data/amber_s/
amber_2 /opt/science/nwchem/current/src/data/amber_q/
amber_3 /opt/science/nwchem/current/src/data/amber_x/
amber_4 /opt/science/nwchem/current/src/data/amber_u/
spce    /opt/science/nwchem/current/src/data/solvents/spce.rst
charmm_s /opt/science/nwchem/current/src/data/charmm_s/
charmm_x /opt/science/nwchem/current/src/data/charmm_x/
```

The memory total is specified as a number of 64 bit doubles. These settings allow NWChem to use a maximum of 1 billion doubles, 8 billion bytes, of memory. The memory total is specified per-core. In a compute cluster this means that NWChem expects to find 8 billion bytes of RAM for every processor on every compute node. If you are running on a single machine, in some corner cases you may need to adjust number of processors used or memory limits to avoid having the N processors * M bytes-per-processor exceed available physical memory.

##Running

###Basic

The driver program `ctc.py` is the easiest way to run one or many composite calculations with a local installation of NWChem. It can run calculations on individual systems from the command line or run many calculations in succession from a .csv spreadsheet such as those included under `test_data`. It will automatically detect the number of CPU cores and available memory and run calculations with available resources. It can also automatically detect and recover from many common error conditions.

Basic example:

    ./ctc.py -g test_data/G2-97/ammonia.xyz
    Running: g3mp2-ccsdt test_data/G2-97/ammonia.xyz 0 singlet
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~NWChem6.5
                   SUMMARY OF G3(MP2,CCSD(T)) CALCULATIONS              
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MP2/6-31G(d)    =   -56.354212   CCSD(T)/6-31G(d) =   -56.372038
        MP2/G3MP2large  =   -56.448179   delta(MP2)       =    -0.093967
        ZPE(HF/6-31G(d))=     0.033039   ZPE Scale Factor =     0.892900
        HLC             =    -0.036680   Free Energy      =     0.000000
        Thermal Energy  =     0.039853   Thermal Enthalpy =     0.040798
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        E(G3(MP2)) @ 0K =   -56.469645   E(G3(MP2)) @298K =   -56.466794
        H(G3(MP2))      =   -56.465850   G(G3(MP2))       =     0.000000
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              HEAT OF FORMATION   (0K):      -8.42 kCal/mol
              HEAT OF FORMATION (298K):     -10.10 kCal/mol
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
     Task  times  cpu:        3.6s     wall:        3.7s

To run many calculations from a CSV:

    ./ctc.py --csv test_data/G2_97-inputs.csv
    
    ...
    
    (many minutes later)
    
    113 of 113 jobs successfully ran

This will generate an output CSV file collecting the results of every calculation, ending in `$MODELNAME-out.csv`. The preceding example generates `test_data/G2_97-inputs-g3mp2-ccsdt-out.csv`, for instance.

Whether run individually or from a CSV file, each completed calculation will also produce a corresponding .js file, e.g. `ammonia_g3mp2-ccsdt_singlet_0.js`. These files contain a parsed JSON representation of the calculation results.  

Run `./ctc.py` without arguments to see command line options and available thermochemical models. 

###Advanced

To run without using `ctc.py`, e.g. on a remote computing facility, ensure that the `.py` file implementing your model is copied to the calculation's working directory. You will also need to add basis set libraries to the global installation or copy them along with your input deck. Inside your `.nw` input file, invoke the model as shown in `examples/NH3.nw`. If you have trouble, make sure this ammonia example works correctly.  

##Background

The _Gaussian-N_ methods developed by Larry A. Curtiss and coworkers combine multiple quantum chemical calculations with empirical fitting and correction factors to closely reproduce experimental enthalpies of formation for a variety of organic molecules. The combination of several modest quantum calculations yields better accuracy per CPU-hour expended than e.g. one heroic calculation with high level correlation treatment and an enormous basis set. The third major iteration of the theory, Gaussian-3, was introduced in Curtiss, L. A., Raghavachari, K., Redfern, P. C., Rassolov, V., & Pople, J. A. (1998). [Gaussian-3 (G3) theory for molecules containing first and second-row atoms.](http://scitation.aip.org/content/aip/journal/jcp/109/18/10.1063/1.477422) _The Journal of chemical physics, 109_(18), 7764-7776. The particularly cost effective G3(MP2) variation was published later that year as Curtiss, L. A., Redfern, P. C., Raghavachari, K., Rassolov, V., & Pople, J. A. (1999). [Gaussian-3 theory using reduced Møller-Plesset order.](http://scitation.aip.org/content/aip/journal/jcp/110/10/10.1063/1.478385) _The Journal of chemical physics, 110_(10), 4703-4709.

This software currently implements G3(MP2) in two "flavors": one that optionally uses a QCISD(T) single point calculation as its single most expensive step, as in the original Gaussian-3 publication, and a default using a CCSD(T) calculation instead. The CCSD(T) variant was first published in Curtiss, L. A., Raghavachari, K., Redfern, P. C., Baboul, A. G., & Pople, J. A. (1999). [Gaussian-3 theory using coupled cluster energies.](http://www.sciencedirect.com/science/article/pii/S0009261499011264) _Chemical physics letters, 314_(1), 101-107. The CCSD(T) default is usually expected to produce slightly better accuracy than the QCISD(T) variant, and to run noticeably faster on closed shell systems. GAMESS-US implements _only_ the CCSD(T) version. The QCISD(T) version is included for compatability and novelty value.

It also implements G4(MP2).

The "gn-" variants are somewhat refactored versions of Daniel Haney's model implementations. They use more tricks to run quickly but should produce the exact same results as his originals.

Available models:

* g3mp2-ccsdt: Like the G3MP2 implementation in GAMESS-US
* g3mp2-qcisdt: Uses the older QCISD(T) model instead of CCSD(T) as highest correlated calculation, as in the original Curtiss G3MP2 model
* g4mp2: Implements the G4MP2 model
* gn-g3mp2-ccsdt: A refactored, speed-optimized implementation of g3mp2-ccsdt. This is the default model.
* gn-g3mp2-qcisdt: A refactored, speed-optimized implementation of g3mp2-qcisdt.
* gn-g4mp2: A refactored, speed-optimized implementation of g4mp2.

##Limitations

* The models are not yet fully calibrated. Agreement with experimental data and prior published results of Curtiss are reasonably good, but some deviations remain.
* This implementation currently handles only first and second row elements. Handling third row elements requires yet-to-be-written additional control over orbital freezing schemes during the correlated wavefunction calculations. Here's the paper with a roadmap for third row calculations: Curtiss, L. A., Redfern, P. C., Rassolov, V., Kedziora, G., & Pople, J. A. (2001). [Extension of Gaussian-3 theory to molecules containing third-row atoms K, Ca, Ga–Kr.](http://scitation.aip.org/content/aip/journal/jcp/114/21/10.1063/1.1366337) _The Journal of Chemical Physics, 114_(21), 9287-9295.
* Calculations that force use of QCISD(T) instead of CCSD(T), and all calculations on open shell systems, run significantly slower than calculations on closed shell systems with CCSD(T). NWChem includes hand written routines for coupled cluster calculations on closed shell systems. QCISD(T) calculations or open shell calculations force correlated calculations via NWChem's Tensor Contraction Engine, which implements autogenerated, very flexible, but slower (at least for small machines) correlated wavefunction calculations.
* Calculations that demand use of the TCE (see previous point) typically run into memory limitations on small machines before they become computationally unbearable. [Contrary to the NWChem manual](http://www.nwchem-sw.org/index.php/Release65:TCE#IO_--_parallel_I.2FO_scheme), TCE calculations _only_ work with the all-in-RAM `ga` I/O scheme. Nonetheless, all the systems in the G2 or G3 test sets should be tractable with any combination of settings on a single modern laptop or desktop computer. Closed-shell-only calculations that use CCSD(T) can use disk for temporary storage, as the hand written routines permit it.
* The reported enthalpy of formation at 298 K is corrected for the constant volume heat capacity, but it _is not_ the Gibbs energy calculated at 298 K. Calculating the Gibbs energy requires incorporation of entropic terms that is not yet automated. The rotational entropy is dependent on correct symmetry determination.
 
