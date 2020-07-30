import subprocess 
import sys
        
#======================================= FILE PREPARATION ==============================================

class ionsMDP:
    def __call__(self):
        self.paramList = []
        self.paramAllName = ['integrator','emtol','emsteps','nsteps','nstlist','cutoffscheme','ns_type','coulombtype','rcoulomb','rvdw','pbc']
        self.integrator = 'integrator = steep'
        self.emtol = 'emtol = 1000.0'
        self.emsteps = 'emstep = 0.01'
        self.nsteps = 'nsteps = 50000'
        self.nstlist = 'nstlist = 1'
        self.cutoffscheme = 'cutoff-scheme = Verlet'
        self.ns_type = 'ns_type = grid'
        self.coulombtype = 'coulombtype = cutoff'
        self.rcoulomb = 'rcoulomb = 1.0'
        self.rvdw = 'rvdw = 1.0'
        self.pbc = 'pbc = xyz'

    def param(self):
        print('Here are parameters of adding ions step:')
        for name in self.paramAllName:
            print(eval('self.'+str(name)))
        print('\n')
        print('Please specify the parameter that you would like to change or type n if you do not want to: ')
        paramName = input()
        if str(paramName) == str('n'):
            print('\n')
        elif str(paramName) in self.paramAllName:
            print('Specify the new value of this parameter: ')
            new = input()
            setattr(self,str(paramName),str(paramName) + ' = ' + str(new))
            print('\n')
            self.param()
        elif str(paramName) == str('cutoff-scheme'):
            print('Specify the new value of this parameter: ')
            new = input()
            self.cutoffscheme = 'cutoff-scheme = %s' %(new)
            print('\n')
            self.param()
        else:
            print('This parameter does not exist!')
            print('\n')
            self.param()
    
    def paramAppend(self):    
        for name in self.paramAllName:
          self.paramList.append(eval('self.'+str(name)))
          
        return
        
    def toFile(self):
        sys.stdout = open('ions.mdp', 'w')
        for i in range(len(self.paramList)):
            print(self.paramList[i])            
        sys.stdout.flush()
        sys.stdout.close()
        sys.stdout = sys.__stdout__

            
def IONSMDP():
    mdp = ionsMDP()
    mdp()
    mdp.param()
    mdp.paramAppend()
    mdp.toFile()    

class minimMDP:
    def __call__(self):
        self.paramList = []
        self.paramAllName = ['integrator','emtol','emsteps','nsteps','nstlist','cutoffscheme','ns_type','coulombtype','rcoulomb','rvdw','pbc']
        self.integrator = 'integrator = steep'
        self.emtol = 'emtol = 1000.0'
        self.emsteps = 'emstep = 0.01'
        self.nsteps = 'nsteps = 50000'
        self.nstlist = 'nstlist = 1'
        self.cutoffscheme = 'cutoff-scheme = Verlet'
        self.ns_type = 'ns_type = grid'
        self.coulombtype = 'coulombtype = PME'
        self.rcoulomb = 'rcoulomb = 1.0'
        self.rvdw = 'rvdw = 1.0'
        self.pbc = 'pbc = xyz'

    def param(self):
        print('Here are parameters of energy minimization step:')
        for name in self.paramAllName:
            print(eval('self.'+str(name)))
        print('\n')
        print('Please specify the parameter that you would like to change or type n if you do not want to: ')
        paramName = input()
        if str(paramName) == str('n'):
            print('\n')
        elif str(paramName) in self.paramAllName:
            print('Specify the new value of this parameter: ')
            new = input()
            setattr(self,str(paramName),str(paramName) + ' = ' + str(new))
            print('\n')
            self.param()
        elif str(paramName) == str('cutoff-scheme'):
            print('Specify the new value of this parameter: ')
            new = input()
            self.cutoffscheme = 'cutoff-scheme = %s' %(new)
            print('\n')
            self.param()
        else:
            print('This parameter does not exist!')
            print('\n')
            self.param()
    
    def paramAppend(self):                
        for name in self.paramAllName:
          self.paramList.append(eval('self.'+str(name)))
          
        return
        
    def toFile(self):
        sys.stdout = open('minim.mdp', 'w')
        for i in range(len(self.paramList)):
            print(self.paramList[i])            
        sys.stdout.flush()
        sys.stdout.close()
        sys.stdout = sys.__stdout__

            
def MINIMMDP():
    mdp = minimMDP()
    mdp()
    mdp.param()
    mdp.paramAppend()
    mdp.toFile()

class nvtMDP:
    def __call__(self):
        self.paramList = []
        self.paramAllName = ['title','define','integrator','nsteps','dt','nstxout','nstvout','nstenergy','nstlog', \
                         'continuation','constraint_algorithm','constraints','lincs_iter','lincs_order','cutoffscheme','ns_type','nstlist','rcoulomb', \
                         'rvdw','DispCorr','coulombtype','pme_order','fourierspacing','tcoupl','tcgrps','tau_t','ref_t','pcoupl', \
                         'pbc','gen_vel','gen_temp','gen_seed']
        self.title = 'title = NVT equilibration'
        self.define = 'define = -DPOSRES'
        self.integrator = 'integrator = md'
        self.nsteps = 'nsteps = 50000'
        self.dt = 'dt = 0.002'
        self.nstxout = 'nstxout = 500' 
        self.nstvout = 'nstvout = 500'
        self.nstenergy = 'nstenergy = 500'
        self.nstlog = 'nstlog = 500'
        self.continuation = 'continuation = no' 
        self.constraint_algorithm = 'constraint_algorithm = lincs' 
        self.constraints = 'constraints = h-bonds'
        self.lincs_iter = 'lincs_iter = 1'
        self.lincs_order = 'lincs_order = 4'
        self.cutoffscheme = 'cutoff-scheme = Verlet'
        self.ns_type = 'ns_type = grid'
        self.nstlist = 'nstlist = 10'
        self.rcoulomb = 'rcoulomb = 1.0'
        self.rvdw = 'rvdw = 1.0'
        self.DispCorr = 'DispCorr = EnerPres'
        self.coulombtype = 'coulombtype = PME'
        self.pme_order = 'pme_order = 4'
        self.fourierspacing = 'fourierspacing = 0.16'
        self.tcoupl = 'tcoupl = V-rescale'
        self.tcgrps = 'tc-grps = Protein Non-Protein'
        self.tau_t = 'tau_t = 0.1     0.1'
        self.ref_t = 'ref_t = 300     300'
        self.pcoupl = 'pcoupl = no'
        self.pbc = 'pbc = xyz'
        self.gen_vel = 'gen_vel = yes'
        self.gen_temp = 'gen_temp = 300'
        self.gen_seed = 'gen_seed = -1'

    def param(self):
        print('Here are parameters of NVT equilibration step:')
        for name in self.paramAllName:
            print(eval('self.'+str(name)))
        print('\n')
        print('Please specify the parameter that you would like to change or type n if you do not want to: ')
        paramName = input()
        if str(paramName) == str('n'):
            print('\n')
        elif str(paramName) in self.paramAllName:
            print('Specify the new value of this parameter: ')
            new = input()
            setattr(self,str(paramName),str(paramName) + ' = ' + str(new))
            print('\n')
            self.param()
        elif str(paramName) == str('cutoff-scheme'):
            print('Specify the new value of this parameter: ')
            new = input()
            self.cutoffscheme = 'cutoff-scheme = %s' %(new)
            print('\n')
            self.param()
        elif str(paramName) == str('tc-grps'):
            print('Specify the new value of this parameter: ')
            new = input()
            self.tcgrps = 'tc-grps = %s' %(new)
            print('\n')
            self.param()
        else:
            print('This parameter does not exist!')
            print('\n')
            self.param()
    
    def paramAppend(self):        
        for name in self.paramAllName:
          self.paramList.append(eval('self.'+str(name)))
          
        return
        
    def toFile(self):
        sys.stdout = open('nvt.mdp', 'w')
        for i in range(len(self.paramList)):
            print(self.paramList[i])            
        sys.stdout.flush()
        sys.stdout.close()
        sys.stdout = sys.__stdout__

            
def NVTMDP():
    mdp = nvtMDP()
    mdp()
    mdp.param()
    mdp.paramAppend()
    mdp.toFile()

class nptMDP:
    def __call__(self):
        self.paramList = []
        self.paramAllName = ['title','define','integrator','nsteps','dt','nstxout','nstvout','nstenergy','nstlog', \
                         'continuation','constraint_algorithm','constraints','lincs_iter','lincs_order','cutoffscheme','ns_type','nstlist','rcoulomb', \
                         'rvdw','DispCorr','coulombtype','pme_order','fourierspacing','tcoupl','tcgrps','tau_t','ref_t','pcoupl','pcoupltype', \
                         'tau_p','ref_p','compressibility','refcoord_scaling','pbc','gen_vel']
        self.title = 'title = NPT equilibration'
        self.define = 'define = -DPOSRES'
        self.integrator = 'integrator = md'
        self.nsteps = 'nsteps = 50000'
        self.dt = 'dt = 0.002'
        self.nstxout = 'nstxout = 500' 
        self.nstvout = 'nstvout = 500'
        self.nstenergy = 'nstenergy = 500'
        self.nstlog = 'nstlog = 500'
        self.continuation = 'continuation = yes' 
        self.constraint_algorithm = 'constraint_algorithm = lincs' 
        self.constraints = 'constraints = h-bonds'
        self.lincs_iter = 'lincs_iter = 1'
        self.lincs_order = 'lincs_order = 4'
        self.cutoffscheme = 'cutoff-scheme = Verlet'
        self.ns_type = 'ns_type = grid'
        self.nstlist = 'nstlist = 10'
        self.rcoulomb = 'rcoulomb = 1.0'
        self.rvdw = 'rvdw = 1.0'
        self.DispCorr = 'DispCorr = EnerPres'
        self.coulombtype = 'coulombtype = PME'
        self.pme_order = 'pme_order = 4'
        self.fourierspacing = 'fourierspacing = 0.16'
        self.tcoupl = 'tcoupl = V-rescale'
        self.tcgrps = 'tc-grps = Protein Non-Protein'
        self.tau_t = 'tau_t = 0.1     0.1'
        self.ref_t = 'ref_t = 300     300'
        self.pcoupl = 'pcoupl = Parrinello-Rahman'
        self.pcoupltype = 'pcoupltype = isotropic'
        self.tau_p = 'tau_p = 2.0'
        self.ref_p = 'ref_p = 1.0'
        self.compressibility = 'compressibility = 4.5e-5'
        self.refcoord_scaling = 'refcoord_scaling = com'
        self.pbc = 'pbc = xyz'
        self.gen_vel = 'gen_vel = no'

    def param(self):
        print('Here are parameters of NPT equilibration step:')
        for name in self.paramAllName:
            print(eval('self.'+str(name)))
        print('\n')
        print('Please specify the parameter that you would like to change or type n if you do not want to: ')
        paramName = input()
        if str(paramName) == str('n'):
            print('\n')
        elif str(paramName) in self.paramAllName:
            print('Specify the new value of this parameter: ')
            new = input()
            setattr(self,str(paramName),str(paramName) + ' = ' + str(new))
            print('\n')
            self.param()
        elif str(paramName) == str('cutoff-scheme'):
            print('Specify the new value of this parameter: ')
            new = input()
            self.cutoffscheme = 'cutoff-scheme = %s' %(new)
            print('\n')
            self.param()
        elif str(paramName) == str('tc-grps'):
            print('Specify the new value of this parameter: ')
            new = input()
            self.tcgrps = 'tc-grps = %s' %(new)
            print('\n')
            self.param()
        else:
            print('This parameter does not exist!')
            print('\n')
            self.param()
    
    def paramAppend(self):        
        for name in self.paramAllName:
          self.paramList.append(eval('self.'+str(name)))
          
        return
        
    def toFile(self):
        sys.stdout = open('npt.mdp', 'w')
        for i in range(len(self.paramList)):
            print(self.paramList[i])            
        sys.stdout.flush()
        sys.stdout.close()
        sys.stdout = sys.__stdout__

            
def NPTMDP():
    mdp = nptMDP()
    mdp()
    mdp.param()
    mdp.paramAppend()
    mdp.toFile()

class mdMDP:
    def __call__(self):
        self.paramList = []
        self.paramAllName = ['title','integrator','nsteps','dt','nstxout','nstvout','nstfout','nstenergy','nstlog','nstxoutcompressed', \
                         'compressedxgrps','continuation','constraint_algorithm','constraints','lincs_iter','lincs_order','cutoffscheme','ns_type','nstlist','rcoulomb', \
                         'rvdw','coulombtype','pme_order','fourierspacing','tcoupl','tcgrps','tau_t','ref_t','pcoupl','pcoupltype', \
                         'tau_p','ref_p','compressibility','pbc','DispCorr','gen_vel']
        self.title = 'title = Molecular dynamics (MD)'
        self.integrator = 'integrator = md'
        self.nsteps = 'nsteps = 500000'
        self.dt = 'dt = 0.002'
        self.nstxout = 'nstxout = 0' 
        self.nstvout = 'nstvout = 0'
        self.nstfout = 'nstfout = 0'
        self.nstenergy = 'nstenergy = 5000'
        self.nstlog = 'nstlog = 5000'
        self.nstxoutcompressed = 'nstxout-compressed = 5000'
        self.compressedxgrps = 'compressed-x-grps = System'
        self.continuation = 'continuation = yes' 
        self.constraint_algorithm = 'constraint_algorithm = lincs' 
        self.constraints = 'constraints = h-bonds'
        self.lincs_iter = 'lincs_iter = 1'
        self.lincs_order = 'lincs_order = 4'
        self.cutoffscheme = 'cutoff-scheme = Verlet'
        self.ns_type = 'ns_type = grid'
        self.nstlist = 'nstlist = 10'
        self.rcoulomb = 'rcoulomb = 1.0'
        self.rvdw = 'rvdw = 1.0'
        self.coulombtype = 'coulombtype = PME'
        self.pme_order = 'pme_order = 4'
        self.fourierspacing = 'fourierspacing = 0.16'
        self.tcoupl = 'tcoupl = V-rescale'
        self.tcgrps = 'tc-grps = Protein Non-Protein'
        self.tau_t = 'tau_t = 0.1     0.1'
        self.ref_t = 'ref_t = 300     300'
        self.pcoupl = 'pcoupl = Parrinello-Rahman'
        self.pcoupltype = 'pcoupltype = isotropic'
        self.tau_p = 'tau_p = 2.0'
        self.ref_p = 'ref_p = 1.0'
        self.compressibility = 'compressibility = 4.5e-5'
        self.pbc = 'pbc = xyz'
        self.DispCorr = 'DispCorr = EnerPres'
        self.gen_vel = 'gen_vel = no'

    def param(self):
        print('Here are parameters of molecular dynamics (MD) step:')
        for name in self.paramAllName:
            print(eval('self.'+str(name)))
        print('\n')
        print('Please specify the parameter that you would like to change or type n if you do not want to: ')
        paramName = input()
        if str(paramName) == str('n'):
            print('\n')
        elif str(paramName) in self.paramAllName:
            print('Specify the new value of this parameter: ')
            new = input()
            setattr(self,str(paramName),str(paramName) + ' = ' + str(new))
            print('\n')
            self.param()
        elif str(paramName) == str('nstxout-compressed'):
            print('Specify the new value of this parameter: ')
            new = input()
            self.nstxoutcompressed = 'nstxout-compressed = %s' %(new)
            print('\n')
            self.param()
        elif str(paramName) == str('compressed-x-grps'):
            print('Specify the new value of this parameter: ')
            new = input()
            self.compressedxgrps = 'compressed-x-grps = %s' %(new)
            print('\n')
            self.param()
        elif str(paramName) == str('cutoff-scheme'):
            print('Specify the new value of this parameter: ')
            new = input()
            self.cutoffscheme = 'cutoff-scheme = %s' %(new)
            print('\n')
            self.param()
        elif str(paramName) == str('tc-grps'):
            print('Specify the new value of this parameter: ')
            new = input()
            self.tcgrps = 'tc-grps = %s' %(new)
            print('\n')
            self.param()
        else:
            print('This parameter does not exist!')
            print('\n')
            self.param()
    
    def paramAppend(self):        
        for name in self.paramAllName:
          self.paramList.append(eval('self.'+str(name)))
          
        return
        
    def toFile(self):
        sys.stdout = open('md.mdp', 'w')
        for i in range(len(self.paramList)):
            print(self.paramList[i])            
        sys.stdout.flush()
        sys.stdout.close()
        sys.stdout = sys.__stdout__

            
def MDMDP():
    mdp = mdMDP()
    mdp()
    mdp.param()
    mdp.paramAppend()
    mdp.toFile()

def removeNonProtein():
    structure = sys.argv[1]

    print("Welcome to AutoSolPro")
    print("Automated molecular dynamics protocols for soluble proteins")
    print("By Khattiya Pongsirijinda")
    print("Skolkovo Institute of Science and Technology")
    print("Khattiya.Pongsirijinda@skoltech.ru")
    print("============================")
    
    nonProtein = ['ANISOU','HETATM']
    with open(structure) as oldfile, open('clean.pdb', 'w') as newfile:
        for line in oldfile:
            if not any(s in line for s in nonProtein):
                newfile.write(line)
    
def forceField():
    pdbToGmx = subprocess.Popen(['gmx','pdb2gmx','-f','clean.pdb','-o','struct.gro','-water','spce'], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    stdin_data = pdbToGmx.communicate(input=b'15 \n')[0]

def solv():
    subprocess.run(['gmx','editconf','-f','struct.gro','-o','box.gro','-c','-d','1.0','-bt','cubic'])
    subprocess.run(['gmx','solvate','-cp','box.gro','-cs','spc216.gro','-o','solv.gro','-p','topol.top'])
    
class ions:
    def saltConc(self):
        print('Default concentration during adding ions step is 0 mol/liter. Would you like to change? y/n: ')
        choice = input()
        if str(choice) == str('y'):
            print('Specify the desired concentration (only the value): ')
            self.conc = input()
            return
        if str(choice) == str('n'):
            self.conc = '0'
            return
        else:   
            print('Please choose y or n.')
            self.saltConc()    

    def addIons(self):
        subprocess.run(['gmx','grompp','-f','ions.mdp','-c','solv.gro','-p','topol.top','-o','ions.tpr'])
        AddIons = subprocess.Popen(['gmx','genion','-s','ions.tpr','-o','solv_ions.gro','-p','topol.top','-pname','NA','-nname','CL','-conc',self.conc,'-neutral'], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        stdin_data = AddIons.communicate(input=b'13 \n')[0]


#====================================== EM PART =================================================================


def gromppEm():
    subprocess.run(['gmx','grompp','-f','minim.mdp','-c','solv_ions.gro','-p','topol.top','-o','em.tpr'])

def mdpEmLauncher():
    subprocess.run(['gmx','mdrun','-v','-deffnm','em'])
    

#====================================== EQ PART =================================================================


def gromppEqNVT():
    subprocess.run(['gmx','grompp','-f','nvt.mdp','-c','em.gro','-r','em.gro','-p','topol.top','-o','nvt.tpr'])

def mdpEqNVTLauncher():
    subprocess.run(['gmx','mdrun','-deffnm','nvt'])
    
def gromppEqNPT():
    subprocess.run(['gmx','grompp','-f','npt.mdp','-c','nvt.gro','-r','nvt.gro','-t','nvt.cpt','-p','topol.top','-o','npt.tpr'])

def mdpEqNPTLauncher():
    subprocess.run(['gmx','mdrun','-deffnm','npt'])
    

#====================================== MD PART =================================================================


def gromppMD():
    subprocess.run(['gmx','grompp','-f','md.mdp','-c','npt.gro','-t','npt.cpt','-p','topol.top','-o','md_0_1.tpr'])

def mdpMDLauncher():
    subprocess.run(['gmx','mdrun','-deffnm','md_0_1'])
    

#===================================== RUN PART ================================================================


def runner():
    removeNonProtein()
    i = ions()
    i.saltConc()
    IONSMDP()
    MINIMMDP()
    NVTMDP()
    NPTMDP()
    MDMDP()
    forceField()
    solv()
    i.addIons()
    print('Files are ready.')
    gromppEm()
    print('EM is ready.')
    mdpEmLauncher()
    print('EM is done.')
    gromppEqNVT()
    print('EQ-NVT is ready.')
    mdpEqNVTLauncher()
    print('EQ-NVT is done.')
    gromppEqNPT()
    print('EQ-NPT is ready.')
    mdpEqNPTLauncher()
    print('EQ-NPT is done.')
    gromppMD()
    print('MD is ready.')
    mdpMDLauncher()
    print('MD is done.')
    

#===================================== EoF PART ================================================================


if __name__ == '__main__': 
    runner()
    