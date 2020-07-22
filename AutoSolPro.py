import subprocess 
import sys
        
#======================================= FILE PREPARATION ==============================================


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
    