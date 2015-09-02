#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
import sys
import argparse
import os
import shutil
import subprocess
import shlex
import json
import time
import csv
import pprint
import fcntl

def _make_fd_blocking(file_obj):
    """
    Updates the flags of the file descriptor to make it blocking.
    file_obj is a `file` object, which has a `.fileno()` method.
    """
    
    fd = file_obj.fileno()
    flags = fcntl.fcntl(fd, fcntl.F_GETFL)
    if flags & os.O_NONBLOCK:
        fcntl.fcntl(fd, fcntl.F_SETFL, flags & ~os.O_NONBLOCK)

def make_stdio_blocking():
    """
    Makes stdout and stderr blocking.
    This prevents resource contention issues with subprocesses.
    See https://github.com/cobrateam/splinter/issues/257
    (This is needed to prevent IOErrors under OS X)
    """

    _make_fd_blocking(sys.__stderr__)
    _make_fd_blocking(sys.__stdout__)

tpl = """start {startname}

memory {memory} mb
title {jobname}

geometry units angstroms print xyz
 {symmetry}
 load {structure}
end

python noprint
import os
import sys
sys.path.append(os.getcwd())
{composite}
end

task python
"""

class Runner(object):
    models = ["g3mp2-ccsdt", "g3mp2-qcisdt", "g4mp2", "gn-g3mp2-ccsdt",
              "gn-g3mp2-qcisdt", "gn-g4mp2"]
    def __init__(self, model, geofile, charge, multiplicity, nproc, memory,
                 tmpdir, verbose, noclean, force):
        self.model = model
        self.geofile = geofile
        self.charge = charge
        self.multiplicity = self.get_multiplicity(multiplicity)
        self.nproc = nproc or self.get_nproc()
        self.memory = memory or self.get_memory()
        self.verbose = verbose
        self.tmpdir = tmpdir
        self.noclean = noclean
        self.force = force
        self.start_time = time.time()

    def get_multiplicity(self, mult):
        """Validate or translate multiplicity.
        """

        m = mult.lower()
        multiplets = ["(null)", "singlet", "doublet",
                      "triplet", "quartet", "quintet",
                      "hextet", "septet", "octet"]
        if m in multiplets:
            result = m
        else:
            try:
                result = multiplets[int(m)]
            except:
                raise ValueError("Invalid multiplicity {0}".format(repr(m)))

        return result

    def get_deck(self, force_c1_symmetry=False):
        """Create a complete job deck for execution. Also return data needed
        to set up job execution.
        """

        memory_per_core = self.memory / self.nproc
        startname = os.path.basename(self.geofile).split(".xyz")[0]
        jobname = "{}_{}_{}_{}".format(startname, self.model,
                                       self.multiplicity, self.charge)

        #symmetry auto-detection will activate if no explicit symmetry set
        symmetry = {True : "symmetry c1", False : ""}[force_c1_symmetry]

        #G3 (MP2, CCSDT)
        if self.model == "g3mp2-ccsdt":
            pymodel = "g3mp2.py"
            m = """import g3mp2
g3mp2.G3MP2(charge={charge}, mult={mult})""".format(charge=self.charge, mult=repr(self.multiplicity))
            if self.multiplicity != "singlet":
                symmetry = "symmetry c1"

        #G3 (MP2, QCISDT)
        elif self.model == "g3mp2-qcisdt":
            pymodel = "g3mp2.py"
            m = """import g3mp2
g3mp2.G3MP2(charge={charge}, mult={mult}, use_qcisdt_f=True)""".format(charge=self.charge, mult=repr(self.multiplicity))
            symmetry = "symmetry c1"

        #G4 (MP2)
        elif self.model == "g4mp2":
            pymodel = "g4mp2.py"
            m = """import g4mp2
g4mp2.G4MP2(charge={charge}, mult={mult})""".format(charge=self.charge, mult=repr(self.multiplicity))
            if self.multiplicity != "singlet":
                symmetry = "symmetry c1"

        #G4 (MP2), alternative implementation
        elif self.model == "gn-g4mp2":
            #allow up to 40% of memory to be used for SCF integral caching
            #value is in bytes rather than megabytes
            #N.B.: default memory partitioning allocates 50% of total
            #memory to global data, and integral cache must fit within
            #global memory section
            integral_cache = int(memory_per_core * 2 ** 20 * 0.4)
            pymodel = "Gn.py"
            m = """import Gn
model=Gn.G4_mp2(charge={charge}, multiplicity={mult}, integral_memory_cache={cache}, force_c1_symmetry={force})
model.run()""".format(charge=self.charge, mult=repr(self.multiplicity), cache=integral_cache, force=force_c1_symmetry)

        #G3 (MP2) CCSD(T), alternative implementation
        elif self.model == "gn-g3mp2-ccsdt":
            integral_cache = int(memory_per_core * 2 ** 20 * 0.4)
            pymodel = "Gn.py"
            m = """import Gn
model=Gn.G3_mp2(charge={charge}, multiplicity={mult}, integral_memory_cache={cache}, force_c1_symmetry={force})
model.run()""".format(charge=self.charge, mult=repr(self.multiplicity), cache=integral_cache, force=force_c1_symmetry)

        #G3 (MP2) QCISD(T), alternative implementation
        elif self.model == "gn-g3mp2-qcisdt":
            integral_cache = int(memory_per_core * 2 ** 20 * 0.4)
            pymodel = "Gn.py"
            m = """import Gn
model=Gn.G3_mp2(charge={charge}, multiplicity={mult}, integral_memory_cache={cache}, use_qcisdt=True, force_c1_symmetry={force})
model.run()""".format(charge=self.charge, mult=repr(self.multiplicity), cache=integral_cache, force=force_c1_symmetry)

        deck = tpl.format(startname=startname, memory=memory_per_core,
                          jobname=jobname,
                          structure=os.path.basename(self.geofile),
                          composite=m, symmetry=symmetry)
        
        tmpdir = self.tmpdir + jobname
        deckfile = jobname + ".nw"
        jsfile = deckfile[:-3] + ".js"
        logfile = deckfile[:-3] + ".log"
        log_location = tmpdir + "/" + logfile

        return {"deck" : deck, "pymodel" : pymodel, "geometry" : self.geofile,
                "jobname" : jobname, "jsfile" : jsfile, "tmpdir" : tmpdir,
                "deckfile" : deckfile, "logfile" : logfile,
                "log_location" : log_location,
                "force_c1_symmetry" : force_c1_symmetry}

    def get_banner(self, prefix):
        """Get a small banner to describe the job currently being run.
        """
        
        banner = " ".join([prefix, self.model, self.geofile,
                           str(self.charge), self.multiplicity])
        return banner

    def run(self, jobdata):
        """Run NWChem for a given deck.
        """

        tmpdir = jobdata["tmpdir"]
        deckfile = jobdata["deckfile"]
        logfile = jobdata["logfile"]
        
        t = jobdata["jobname"]
        os.system("rm -rf {}".format(tmpdir))
        os.makedirs(tmpdir)
        
        with open(tmpdir + "/" + deckfile, "w") as outfile:
            outfile.write(jobdata["deck"])

        shutil.copy(jobdata["pymodel"], tmpdir)
        shutil.copy(jobdata["geometry"], tmpdir)

        if self.verbose:
            redirector = "| tee "
        else:
            redirector = "&> "

        if self.nproc == 1:
            runner = "cd {0} && nwchem {1} {2} {3}".format(tmpdir, deckfile, redirector, logfile)
        else:
            runner = "cd {0} && mpirun -np {1} nwchem {2} {3} {4}".format(tmpdir, self.nproc, deckfile, redirector, logfile)

        banner = self.get_banner("Running:")
        print(banner)

        if not self.verbose:
            cmd = shlex.split(runner)
            command = ["/bin/bash", "-i", "-c"] + [" ".join(cmd)]
            p = subprocess.Popen(command, stdout=subprocess.PIPE,
                                 stdin=subprocess.PIPE)
            output = p.communicate()[0]

        else:
            os.system(runner)

        make_stdio_blocking()

    def parse_summary(self, summary):
        """Parse the summary banner from NWChem log file to get individual
        energy components from the job.

        :param summary: block of summary text to parse
        :type summary : str
        :return: parsed components
        :rytpe : dict
        """

        D = {}
        for line in summary.split("\n"):
            #ZPE(HF/6-31G(d))=     0.033039   ZPE Scale Factor =     0.892900
            if "=" in line:
                pieces = []
                for chunk in line.split():
                    if chunk.endswith("="):
                        pieces.append(chunk[:-1])
                        pieces.append("=")
                    else:
                        pieces.append(chunk)

                buffer = []
                for piece in pieces:
                    try:
                        fv = float(piece)
                    except:
                        fv = None
                    if fv is not None:
                        key = " ".join(buffer).strip()
                        D[key] = fv
                        buffer = []
                    elif piece != "=":
                        buffer.append(piece)

            #HEAT OF FORMATION   (0K):     -21.70 kCal/mol
            elif ":" in line:
                try:
                    a, b = line.split(":")
                    a = a.strip()
                    f = float(b.split()[0])
                    D[a] = f
                except ValueError:
                    pass

        return D

    def run_and_extract(self, jobdata):
        """Run calculation job and handle logged output, including helpful
        error messages.

        Testing error cases:
         -For unparameterized elements, try hydrogen selenide
         -For bad multiplicity, try neutral ammonia, triplet
         -For geometry optimization failure, try +4 ammonia, triplet
          (it is surprisingly hard to make geometry optimization fail)
         -For symmetry failure, try g3mp2-ccsdt h2 (with NWChem < r27285)

        :param jobdata: data to control current job
        :type jobdata : dict
        :return: extracted data or error
        :rtype : dict
        """

        log_location = jobdata["log_location"]

        if not self.force and os.path.exists(jobdata["jsfile"]):
            rerun = True
            banner = self.get_banner("Skipping, already ran:")
            print(banner)
            with open(jobdata["jsfile"]) as jsin:
                serialized = jsin.read()
            deserialized = json.loads(serialized)
            print(deserialized["summary"])
            return deserialized

        else:
            rerun = False
            self.run(jobdata)
            elapsed = time.time() - self.start_time
            
        with open(log_location) as lf:
            log = lf.readlines()

        extracting = False
        extracted = []
        for line in log:
            if "~~~" in line:
                extracting = True

            if extracting:
                extracted.append(line)
                
            if line.strip().split()[:2] == ["Task", "times"]:
                extracting = False

        summary = "".join(extracted)
        print(summary)

        #Job ran to expected completion
        if summary:
            parsed = self.parse_summary(summary)
            records = {"summary" : summary, "multiplicity" : self.multiplicity,
                       "nproc" : self.nproc, "memory" : self.memory,
                       "geofile" : self.geofile, "model" : self.model,
                       "charge" : self.charge, "elapsed" : elapsed,
                       "components" : parsed}
        
            with open(jobdata["jsfile"], "w") as jshandle:
                json.dump(records, jshandle, sort_keys=True, indent=2)

            if not self.noclean:
                os.system("rm -rf {0}".format(jobdata["tmpdir"]))

            return records

        #Job failed somehow
        else:
            logdata = "".join(log)
            errors = {"no. of electrons and multiplicity not compatible" :
                      "The multiplicity appears to be incorrect for the given system and charge.",
                      "no. of closed-shell electrons is not even" :
                      "The multiplicity appears to be incorrect for the given system and charge.",
                      "bas_tag_lib: no such basis available" :
                      "Input contains unparameterized elements. These methods are tested only for main group elements through the second row.",
                      "driver_energy_step: energy failed" :
                      "Geometry optimization failure. You may need to provide an input geometry that is closer to equilibrium. See details in " + log_location,
                      "driver: task_gradient failed" :
                      "Geometry optimization failure. You may need to provide an input geometry that is closer to equilibrium. See details in " + log_location,
                      "Failed to converge in maximum number of steps" :
                      "Geometry optimization failure. You may need to provide an input geometry that is closer to equilibrium. See details in " + log_location,
                      "sym_center_map is inconsistent" :
                      "Symmetry problems with geometry. Forcing C1.",
                      "non-Abelian symmetry not permitted" :
                      "Symmetry problems with geometry. Forcing C1.",
                      "AUTOZ failed" :
                      "Symmetry problems with geometry. Forcing C1.",
                      "ran out of MA memory" :
                      "Memory shortage",
                      "MA error" :
                      "Memory shortage",
                      "dgesv failed" :
                      "Numerical failure",
                      "moints_epair_eval: zero denominator" :
                      "Numerical failure",
                      "maximum iterations exceeded" :
                      "Convergence failure",
                      "calculations not reaching convergence" :
                      "Convergence failure"}

            cause = ""
            for k, v in sorted(errors.items()):
                if k in logdata:
                    sys.stderr.write(v + "\n")
                    cause = v

            #Sometimes autosym fails, but forcing C1 allows job to run
            #Also sometimes random numerical problems are fixed this way
            symmetry_problems = ["Symmetry problems", "Numerical failure",
                                 "Geometry optimization failure", "Convergence failure"]
            symfail = [cause.startswith(p) for p in symmetry_problems]
            if True in symfail and not jobdata["force_c1_symmetry"]:
                sys.stderr.write("Trying c1 symmetry to fix previous failure\n")
                deck = self.get_deck(force_c1_symmetry=True)
                return self.run_and_extract(deck)

            #If memory runs short, try reducing nproc to give more memory per core
            elif cause.startswith("Memory shortage") and self.nproc > 1:
                self.nproc /= 2
                deck = self.get_deck()
                return self.run_and_extract(deck)

            if not cause:
                sys.stderr.write("Unknown error. See details in {}\n".format(log_location))
                cause = "unknown - {}".format(log_location)

            return {"error" : cause}

    def get_memory(self):
        """Automatically get available memory (Linux only)
        """

        megabytes = 0
        #default memory allocation scheme gives too little
        #to global memory, so pad total memory value to get
        #closer to real practical RAM limits
        padding = 1.5
        
        try:
            with open("/proc/meminfo") as infile:
                data = infile.read()
            for line in data.split("\n"):
                if "MemTotal" in line:
                    kilobytes = int(line.strip().split()[1])
                    megabytes = int((kilobytes * padding) / 1024)
        except IOError:
            megabytes = 1000

        return max(megabytes, 1000)

    def get_nproc(self):
        """Automatically get number of processors (Linux only)
        """

        nproc = 0
        amd = True

        try:
            with open("/proc/cpuinfo") as infile:
                data = infile.read()
            for line in data.split("\n"):
                if "GenuineIntel" in line:
                    amd = False
                elif "processor" in line:
                    try:
                        nproc = int(line.strip().split()[-1]) + 1
                    except ValueError:
                        pass
            #assume hyperthreading if intel processor, use only real cores
            if not amd:
                nproc /= 2
        #can't read cpuinfo, so default to 1
        except IOError:
            nproc = 1

        return max(1, nproc)

def main(args):
    try:
        m = Runner(args.model, args.xyz, args.charge, args.multiplicity,
                   args.nproc, args.memory, args.tmpdir, args.verbose,
                   args.noclean, args.force)
        deck = m.get_deck()
    except:
        return True

    m.run_and_extract(deck)

def csvwrite(outname, rows, header):
    with open(outname, "w") as outfile:
        d = csv.DictWriter(outfile, fieldnames=header)
        d.writeheader()
        for row in rows:
            d.writerow(row)
    
def csvmain(args):
    failures = []
    successes = []
    header = []

    columns = ["elapsed", "memory", "nproc"]
    
    with open(args.csv) as csvfile:
        reader = csv.DictReader(csvfile)
        rows = [r for r in reader]

    for j, row in enumerate(rows):
        msg = "Running {} of {}".format(j + 1, len(rows))
        print(msg)
        m = Runner(args.model, row["System"], row["Charge"],
                   row["Multiplicity"], args.nproc, args.memory, args.tmpdir,
                   args.verbose, args.noclean, args.force)
        deck = m.get_deck()
        result = m.run_and_extract(deck)
        if result.get("error"):
            failures.append((row["System"], row["Charge"], row["Multiplicity"],
                             result.get("error")))

        else:
            if not header:
                extras = result["components"].keys() + columns
                header = reader.fieldnames + sorted(extras)
                model = result["model"]
                outname = args.csv.replace(".csv", "-{model}-out.csv".format(model=model))
            s = row
            
            for col in columns:
                s[col] = result[col]
                
            s.update(result["components"])
            successes.append(s)
            csvwrite(outname, successes, header)

    print "{} of {} jobs successfully ran".format(len(rows) - len(failures),
                                                  len(rows))
    if failures:
        print("Failures:")
        pprint.pprint(failures)

if __name__ == "__main__":
    available_models = ", ".join(Runner.models)
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="Treat a chemical system with one of the following composite thermochemical models: " + available_models + ". An .xyz file or appropriate .csv file is required as input.")
    parser.add_argument("-n", "--nproc", help="Number of processor cores to use (auto-assigned if not chosen)", type=int,default=0)
    parser.add_argument("--memory", help="Maximum memory to use, in megabytes (auto-assigned if not chosen)", type=int, default=0)
    parser.add_argument("--multiplicity", help="System spin multiplicity", default="singlet")
    parser.add_argument("-m", "--model", help="Thermochemical model to use", default="g3mp2-ccsdt")
    parser.add_argument("-c", "--charge", help="System charge", type=int, default=0)
    parser.add_argument("-g", "--xyz", help="XYZ geometry file", default="")
    parser.add_argument("--csv", help="CSV file describing a batch of systems to run", default="")
    parser.add_argument("-v", "--verbose", help="If activated, show job output as it executes", action="store_true", default=False)
    parser.add_argument("--noclean", help="If active, don't clean up temporary files after calculation", action="store_true", default=False)
    parser.add_argument("--force", help="If active, re-run a calculation even when output file already exists", action="store_true", default=False)
    parser.add_argument("--tmpdir", help="Temporary directory", default="/tmp/")
    args = parser.parse_args()

    if args.model not in Runner.models:
        sys.stderr.write("Unknown model. Available models are: {}\n".format(available_models))
        sys.exit(1)
    
    if not (args.xyz or args.csv):
        parser.print_help()
        sys.stderr.write("\nYou must supply an .xyz system geometry file or a .csv file describing multiple systems.\n")
        
    elif args.xyz:
        error = main(args)
        if error:
            parser.print_help()

    elif args.csv:
        csvmain(args)
