# composite-thermochemistry-nwchem
Composite thermochemical models implemented with NWChem's Python functionality

##Installation

###Requirements

* A supported NWChem platform: this code has been tested under Linux and OS X. In theory, NWChem runs under Windows but the authors have no experience with it. The Windows port appears to be minimally supported. Windows users without a lot of cross-platform compilation experience should consider a Linux virtual machine instead.
* A _recent_, _Python enabled_ build of NWChem: ideally after svn revision 27160 but at least following 27110. Use the latest source code development snapshot from the [NWChem download page](http://www.nwchem-sw.org/index.php/Download).
* The g3mp2large basis set (included), which must be copied into NWChem's library directory.

###Installing

After you have compiled NWChem with embedded Python support and installed it, copy the `g3mp2large` file from `basis_data` to NWChem's basis library directory. The library directory is `$NWCHEM_TOP/src/basis/libraries/`. Place the `g3mp2.py` file in a location where NWChem can pick it up from the `PYTHONPATH`. If you are unsure, copy `g3mp2.py` into the working directory where you are setting up a calculation. Try running the included ammonia examples. You should get results very similar to the example ammonia output.

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

##Background

The _Gaussian-N_ methods developed by Larry A. Curtiss and coworkers combine multiple quantum chemical calculations with semi-empirical fitting factors to closely reproduce experimental enthalpies of formation for a variety of organic molecules. The combination of several modest quantum calculations yields better accuracy per CPU-hour expended than e.g. one heroic calculation with high level correlation treatment and an enormous basis set. The third major iteration of the theory, Gaussian-3, was introduced in Curtiss, L. A., Raghavachari, K., Redfern, P. C., Rassolov, V., & Pople, J. A. (1998). [Gaussian-3 (G3) theory for molecules containing first and second-row atoms.](http://scitation.aip.org/content/aip/journal/jcp/109/18/10.1063/1.477422) _The Journal of chemical physics, 109_(18), 7764-7776. The particularly cost effective G3(MP2) variation was published later that year as Curtiss, L. A., Redfern, P. C., Raghavachari, K., Rassolov, V., & Pople, J. A. (1999). [Gaussian-3 theory using reduced Møller-Plesset order.](http://scitation.aip.org/content/aip/journal/jcp/110/10/10.1063/1.478385) _The Journal of chemical physics, 110_(10), 4703-4709.

This software currently implements G3(MP2) in two "flavors": one that optionally uses a QCISD(T) single point calculation as its single most expensive step, as in the original Gaussian-3 publication, and a default using a CCSD(T) calculation instead. The CCSD(T) variant was first published in Curtiss, L. A., Raghavachari, K., Redfern, P. C., Baboul, A. G., & Pople, J. A. (1999). [Gaussian-3 theory using coupled cluster energies.](http://www.sciencedirect.com/science/article/pii/S0009261499011264) _Chemical physics letters, 314_(1), 101-107. The CCSD(T) default is usually expected to produce slightly better accuracy than the QCISD(T) variant, and to run noticeably faster on closed shell systems. GAMESS-US implements _only_ the CCSD(T) version. The QCISD(T) version is included for compatability and novelty value.

##Limitations

* This implementation currently handles only first and second row elements. Handling third row elements requires yet-to-be-written additional control over orbital freezing schemes during the correlated wavefunction calculations. Here's the paper with a roadmap for third row calculations: Curtiss, L. A., Redfern, P. C., Rassolov, V., Kedziora, G., & Pople, J. A. (2001). [Extension of Gaussian-3 theory to molecules containing third-row atoms K, Ca, Ga–Kr.](http://scitation.aip.org/content/aip/journal/jcp/114/21/10.1063/1.1366337) _The Journal of Chemical Physics, 114_(21), 9287-9295.
* Calculations that force use of QCISD(T) instead of CCSD(T), and all calculations on open shell systems, run significantly slower than calculations on closed shell systems with CCSD(T). NWChem includes hand written routines for MP2 and coupled cluster calculations on closed shell systems. QCISD(T) calculations or open shell calculations force correlated calculations via NWChem's Tensor Contraction Engine, which implements autogenerated, very flexible, but slower (at least for small machines) correlated wavefunction calculations.
* Calculations that demand use of the TCE (see previous point) typically run into memory limitations on small machines before they become computationally unbearable. Contrary to the NWChem manual, TCE calculations _only_ work with the all-in-RAM `ga` I/O scheme. Nonetheless, all the systems in the G2 or G3 test sets should be tractable with any combination of settings on a single modern laptop or desktop computer. Closed-shell-only calculations that use CCSD(T) can use disk for temporary storage, as the hand written MP2 and CCSD(T) routines both permit it.
* Example calculations currently force C1 symmetry, e.g. do not exploit symmetry groups. Calculations on symmetric molecules thus run slower than they could. NWChem's automatic symmetry detection is fairly reliable, but it can find non-Abelian symmetries that the TCE chokes on. The TCE must be used for open shell systems and for QCISD(T)-forced runs. In the interest of unexpected-failure-free operation, our examples thus show consistent use of C1 symmetry. If you know your molecule's Abelian symmetry group you can input it by hand instead of "C1". If you are doing only closed-shell calculations with CCSD(T) you can remove symmetry directives and let NWChem find your symmetry automatically.
* The reported enthalpy of formation at 298 K is corrected for the constant volume heat capacity, but it _is not_ the Gibbs energy calculated at 298 K. Calculating the Gibbs energy requires incorporation of entropic terms. The rotational entropy is dependent on correct symmetry determination. We cannot leave NWChem automatic symmetry determination on in the general case (see previous point). Therefore, the 298 K results reported here are more akin to those reported by GAMESS-US than those from Gaussian 09.
 
