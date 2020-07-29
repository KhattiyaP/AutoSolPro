import subprocess 
import sys
        
#======================================= FILE PREPARATION ==============================================


class mdpFile(object):
    
    def __init__(self):
        self.paramList = []
        
    def __call__(self):
        self.paramList.append('title                   = Molecular dynamics (MD)')
        self.paramList.append('; Run parameters')
        self.paramList.append('integrator              = md        ; leap-frog integrator')
        
    def numberOfSteps(self): 
        print('Default number of steps is 500000 (1 ns). Would you like to change? y/n: ')
        default = input()
        if str(default) == str('y'):
            print('Specify the desired number of steps: ')
            new = input()
            new = int(new)
            nsteps = 'nsteps                  = %d' %(new)
            self.paramList.append(nsteps)   
            return
        if str(default) == str('n'):
            self.paramList.append('nsteps                  = 500000')
            return
        else:   
            print('Please choose y or n.')
            self.numberOfSteps()
            
    def param(self):
        self.paramList.append('dt                      = 0.002     ; 2 fs')
        self.paramList.append('; Output control')
        self.paramList.append('nstxout                 = 0         ; suppress bulky .trr file by specifying') 
        self.paramList.append('nstvout                 = 0         ; 0 for output frequency of nstxout,')
        self.paramList.append('nstfout                 = 0         ; nstvout, and nstfout')
        self.paramList.append('nstenergy               = 5000      ; save energies every 10.0 ps')
        self.paramList.append('nstlog                  = 5000      ; update log file every 10.0 ps')
        self.paramList.append('nstxout-compressed      = 5000      ; save compressed coordinates every 10.0 ps')
        self.paramList.append('compressed-x-grps       = System    ; save the whole system')
        self.paramList.append('; Bond parameters')
        self.paramList.append('continuation            = yes       ; Restarting after NPT') 
        self.paramList.append('constraint_algorithm    = lincs     ; holonomic constraints') 
        self.paramList.append('constraints             = h-bonds   ; bonds involving H are constrained')
        self.paramList.append('lincs_iter              = 1         ; accuracy of LINCS')
        self.paramList.append('lincs_order             = 4         ; also related to accuracy')
        self.paramList.append('; Neighborsearching')
        self.paramList.append('cutoff-scheme           = Verlet    ; Buffered neighbor searching')
        self.paramList.append('ns_type                 = grid      ; search neighboring grid cells')
        self.paramList.append('nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet scheme')
        self.paramList.append('rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)')
        self.paramList.append('rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)')
        self.paramList.append('; Electrostatics')
        self.paramList.append('coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics')
        self.paramList.append('pme_order               = 4         ; cubic interpolation')
        self.paramList.append('fourierspacing          = 0.16      ; grid spacing for FFT')
        self.paramList.append('; Temperature coupling is on')
        self.paramList.append('tcoupl                  = V-rescale             ; modified Berendsen thermostat')
        self.paramList.append('tc-grps                 = Protein Non-Protein   ; two coupling groups - more accurate')
        self.paramList.append('tau_t                   = 0.1     0.1           ; time constant, in ps')
        self.paramList.append('ref_t                   = 300     300           ; reference temperature, one for each group, in K')
        self.paramList.append('; Pressure coupling is on')
        self.paramList.append('pcoupl                  = Parrinello-Rahman     ; Pressure coupling on in NPT')
        self.paramList.append('pcoupltype              = isotropic             ; uniform scaling of box vectors')
        self.paramList.append('tau_p                   = 2.0                   ; time constant, in ps')
        self.paramList.append('ref_p                   = 1.0                   ; reference pressure, in bar')
        self.paramList.append('compressibility         = 4.5e-5                ; isothermal compressibility of water, bar^-1')
        self.paramList.append('; Periodic boundary conditions')
        self.paramList.append('pbc                     = xyz       ; 3-D PBC')
        self.paramList.append('; Dispersion correction')
        self.paramList.append('DispCorr                = EnerPres  ; account for cut-off vdW scheme')
        self.paramList.append('; Velocity generation')
        self.paramList.append('gen_vel                 = no        ; Velocity generation is off')
        
        return
        
    def toFile(self):
        sys.stdout = open('md.mdp', 'w')
        for i in range(len(self.paramList)):
            print(self.paramList[i])            
        sys.stdout.flush()
        sys.stdout.close()
        sys.stdout = sys.__stdout__

            
def MDmdp():
    mdp = mdpFile()
    mdp()
    mdp.numberOfSteps()
    mdp.param()
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

def addIons():
    subprocess.run(['gmx','grompp','-f','ions.mdp','-c','solv.gro','-p','topol.top','-o','ions.tpr'])
    AddIons = subprocess.Popen(['gmx','genion','-s','ions.tpr','-o','solv_ions.gro','-p','topol.top','-pname','NA','-nname','CL','-neutral'], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
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


def mdRun():
    removeNonProtein()
    MDmdp()
    forceField()
    solv()
    addIons()
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
    mdRun()
    