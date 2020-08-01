import subprocess 
import sys
import argparse
        
#======================================= PARAM INPUT ==============================================


class paramInput:
    def __call__(self):
        print("Welcome to AutoSolPro")
        print("Automated molecular dynamics protocols for soluble proteins")
        print("By Khattiya Pongsirijinda")
        print("Skolkovo Institute of Science and Technology")
        print("Khattiya.Pongsirijinda@skoltech.ru")
        print("============================")
        
    def param(self):
        param_ions = ['integrator_ions','emtol_ions','emstep_ions','nsteps_ions','nstlist_ions','cutoffscheme_ions','ns_type_ions','coulombtype_ions','rcoulomb_ions','rvdw_ions','pbc_ions']
        param_minim = ['integrator_minim','emtol_minim','emstep_minim','nsteps_minim','nstlist_minim','cutoffscheme_minim','ns_type_minim','coulombtype_minim','rcoulomb_minim','rvdw_minim','pbc_minim']
        param_nvt = ['title_nvt', 'define_nvt', 'integrator_nvt', 'nsteps_nvt', 'dt_nvt', 'nstxout_nvt', 'nstvout_nvt', 'nstenergy_nvt', 'nstlog_nvt', 'continuation_nvt', 'constraint_algorithm_nvt', 'constraints_nvt', 'lincs_iter_nvt', 'lincs_order_nvt', 'cutoffscheme_nvt', 'ns_type_nvt', 'nstlist_nvt', 'rcoulomb_nvt', 'rvdw_nvt', 'DispCorr_nvt', 'coulombtype_nvt', 'pme_order_nvt', 'fourierspacing_nvt', 'tcoupl_nvt', 'tcgrps_nvt', 'tau_t_nvt', 'ref_t_nvt', 'pcoupl_nvt', 'pbc_nvt', 'gen_vel_nvt', 'gen_temp_nvt', 'gen_seed_nvt']
        param_npt = ['title_npt', 'define_npt', 'integrator_npt', 'nsteps_npt', 'dt_npt', 'nstxout_npt', 'nstvout_npt', 'nstenergy_npt', 'nstlog_npt', 'continuation_npt', 'constraint_algorithm_npt', 'constraints_npt', 'lincs_iter_npt', 'lincs_order_npt', 'cutoffscheme_npt', 'ns_type_npt', 'nstlist_npt', 'rcoulomb_npt', 'rvdw_npt', 'DispCorr_npt', 'coulombtype_npt', 'pme_order_npt', 'fourierspacing_npt', 'tcoupl_npt', 'tcgrps_npt', 'tau_t_npt', 'ref_t_npt', 'pcoupl_npt', 'pcoupltype_npt', 'tau_p_npt', 'ref_p_npt', 'compressibility_npt', 'refcoord_scaling_npt', 'pbc_npt', 'gen_vel_npt']
        param_md = ['title_md', 'integrator_md', 'nsteps_md', 'dt_md', 'nstxout_md', 'nstvout_md', 'nstfout_md', 'nstenergy_md', 'nstlog_md', 'nstxoutcompressed_md', 'compressedxgrps_md', 'continuation_md', 'constraint_algorithm_md', 'constraints_md', 'lincs_iter_md', 'lincs_order_md', 'cutoffscheme_md', 'ns_type_md', 'nstlist_md', 'rcoulomb_md', 'rvdw_md', 'coulombtype_md', 'pme_order_md', 'fourierspacing_md', 'tcoupl_md', 'tcgrps_md', 'tau_t_md', 'ref_t_md', 'pcoupl_md', 'pcoupltype_md', 'tau_p_md', 'ref_p_md', 'compressibility_md', 'pbc_md', 'DispCorr_md', 'gen_vel_md']
        paramMDP = param_ions + param_minim + param_nvt + param_npt + param_md
        
        param_file = ['file']
        param_ff = ['water_ff']
        param_ed = ['bt_ed','d_ed']
        param_gen = ['pname_gen','nname_gen','conc_gen']
        paramRest = param_file + param_ff + param_ed + param_gen
        
        parser = argparse.ArgumentParser()
        #pdb file
        parser.add_argument('--file',help='Input pdb file name',required=True)
        
        #ions.mdp
        parser.add_argument('--integrator_ions',help='Algorithm (steep = steepest descent minimization)',default='steep')
        parser.add_argument('--emtol_ions',help='Stop minimization when the maximum force < 1000.0 kJ/mol/nm',default='1000.0')
        parser.add_argument('--emstep_ions',help='Minimization step size',default='0.01')
        parser.add_argument('--nsteps_ions',help='Maximum number of (minimization) steps to perform',default='50000')
        parser.add_argument('--nstlist_ions',help='Frequency to update the neighbor list and long range forces',default='1')
        parser.add_argument('--cutoffscheme_ions',help='Buffered neighbor searching',default='Verlet')
        parser.add_argument('--ns_type_ions',help='Method to determine neighbor list (simple, grid)',default='grid')
        parser.add_argument('--coulombtype_ions',help='Treatment of long range electrostatic interactions',default='cutoff')
        parser.add_argument('--rcoulomb_ions',help='Short-range electrostatic cut-off',default='1.0')
        parser.add_argument('--rvdw_ions',help='Short-range Van der Waals cut-off',default='1.0')
        parser.add_argument('--pbc_ions',help='Periodic Boundary Conditions in all 3 dimensions',default='xyz')
        
        #minim.mdp
        parser.add_argument('--integrator_minim',help='Algorithm (steep = steepest descent minimization)',default='steep')
        parser.add_argument('--emtol_minim',help='Stop minimization when the maximum force < 1000.0 kJ/mol/nm',default='1000.0')
        parser.add_argument('--emstep_minim',help='Minimization step size',default='0.01')
        parser.add_argument('--nsteps_minim',help='Maximum number of (minimization) steps to perform',default='50000')
        parser.add_argument('--nstlist_minim',help='Frequency to update the neighbor list and long range forces',default='1')
        parser.add_argument('--cutoffscheme_minim',help='Buffered neighbor searching',default='Verlet')
        parser.add_argument('--ns_type_minim',help='Method to determine neighbor list (simple, grid)',default='grid')
        parser.add_argument('--coulombtype_minim',help='Treatment of long range electrostatic interactions',default='PME')
        parser.add_argument('--rcoulomb_minim',help='Short-range electrostatic cut-off',default='1.0')
        parser.add_argument('--rvdw_minim',help='Short-range Van der Waals cut-off',default='1.0')
        parser.add_argument('--pbc_minim',help='Periodic Boundary Conditions in all 3 dimensions',default='xyz')
        
        #nvt.mdp
        parser.add_argument('--title_nvt',help='Title',default='NVT equilibration')
        parser.add_argument('--define_nvt',help='position restrain the protein',default='-DPOSRES')
        parser.add_argument('--integrator_nvt',help='leap-frog integrator',default='md')
        parser.add_argument('--nsteps_nvt',help='2 * 50000 = 100 ps',default='50000')
        parser.add_argument('--dt_nvt',help='2 fs',default='0.002')
        parser.add_argument('--nstxout_nvt',help='save coordinates every 1.0 ps',default='500')
        parser.add_argument('--nstvout_nvt',help='save velocities every 1.0 ps',default='500')
        parser.add_argument('--nstenergy_nvt',help='save energies every 1.0 ps',default='500')
        parser.add_argument('--nstlog_nvt',help='update log file every 1.0 ps',default='500')
        parser.add_argument('--continuation_nvt',help='first dynamics run',default='no')
        parser.add_argument('--constraint_algorithm_nvt',help='holonomic constraints',default='lincs')
        parser.add_argument('--constraints_nvt',help='bonds involving H are constrained',default='h-bonds')
        parser.add_argument('--lincs_iter_nvt',help='accuracy of LINCS',default='1')
        parser.add_argument('--lincs_order_nvt',help='also related to accuracy',default='4')
        parser.add_argument('--cutoffscheme_nvt',help='Buffered neighbor searching',default='Verlet')
        parser.add_argument('--ns_type_nvt',help='search neighboring grid cells',default='grid')
        parser.add_argument('--nstlist_nvt',help='20 fs, largely irrelevant with Verlet',default='10')
        parser.add_argument('--rcoulomb_nvt',help='short-range electrostatic cutoff (in nm)',default='1.0')
        parser.add_argument('--rvdw_nvt',help='short-range van der Waals cutoff (in nm)',default='1.0')
        parser.add_argument('--DispCorr_nvt',help='account for cut-off vdW scheme',default='EnerPres')
        parser.add_argument('--coulombtype_nvt',help='Particle Mesh Ewald for long-range electrostatics',default='PME')
        parser.add_argument('--pme_order_nvt',help='cubic interpolation',default='4')
        parser.add_argument('--fourierspacing_nvt',help='grid spacing for FFT',default='0.16')
        parser.add_argument('--tcoupl_nvt',help='modified Berendsen thermostat',default='V-rescale')
        parser.add_argument('--tcgrps_nvt',help='two coupling groups - more accurate',default='Protein Non-Protein')
        parser.add_argument('--tau_t_nvt',help='time constant, in ps',default='0.1     0.1')
        parser.add_argument('--ref_t_nvt',help='reference temperature, one for each group, in K',default='300     300')
        parser.add_argument('--pcoupl_nvt',help='no pressure coupling in NVT',default='no')
        parser.add_argument('--pbc_nvt',help='3-D PBC',default='xyz')
        parser.add_argument('--gen_vel_nvt',help='assign velocities from Maxwell distribution',default='yes')
        parser.add_argument('--gen_temp_nvt',help='temperature for Maxwell distribution',default='300')
        parser.add_argument('--gen_seed_nvt',help='generate a random seed',default='-1')
        
        #npt.mdp
        parser.add_argument('--title_npt',help='Title',default='NPT equilibration')
        parser.add_argument('--define_npt',help='position restrain the protein',default='-DPOSRES')
        parser.add_argument('--integrator_npt',help='leap-frog integrator',default='md')
        parser.add_argument('--nsteps_npt',help='2 * 50000 = 100 ps',default='50000')
        parser.add_argument('--dt_npt',help='2 fs',default='0.002')
        parser.add_argument('--nstxout_npt',help='save coordinates every 1.0 ps',default='500')
        parser.add_argument('--nstvout_npt',help='save velocities every 1.0 ps',default='500')
        parser.add_argument('--nstenergy_npt',help='save energies every 1.0 ps',default='500')
        parser.add_argument('--nstlog_npt',help='update log file every 1.0 ps',default='500')
        parser.add_argument('--continuation_npt',help='Restarting after NVT',default='yes')
        parser.add_argument('--constraint_algorithm_npt',help='holonomic constraints',default='lincs')
        parser.add_argument('--constraints_npt',help='bonds involving H are constrained',default='h-bonds')
        parser.add_argument('--lincs_iter_npt',help='accuracy of LINCS',default='1')
        parser.add_argument('--lincs_order_npt',help='also related to accuracy',default='4')
        parser.add_argument('--cutoffscheme_npt',help='Buffered neighbor searching',default='Verlet')
        parser.add_argument('--ns_type_npt',help='search neighboring grid cells',default='grid')
        parser.add_argument('--nstlist_npt',help='20 fs, largely irrelevant with Verlet',default='10')
        parser.add_argument('--rcoulomb_npt',help='short-range electrostatic cutoff (in nm)',default='1.0')
        parser.add_argument('--rvdw_npt',help='short-range van der Waals cutoff (in nm)',default='1.0')
        parser.add_argument('--DispCorr_npt',help='account for cut-off vdW scheme',default='EnerPres')
        parser.add_argument('--coulombtype_npt',help='Particle Mesh Ewald for long-range electrostatics',default='PME')
        parser.add_argument('--pme_order_npt',help='cubic interpolation',default='4')
        parser.add_argument('--fourierspacing_npt',help='grid spacing for FFT',default='0.16')
        parser.add_argument('--tcoupl_npt',help='modified Berendsen thermostat',default='V-rescale')
        parser.add_argument('--tcgrps_npt',help='two coupling groups - more accurate',default='Protein Non-Protein')
        parser.add_argument('--tau_t_npt',help='time constant, in ps',default='0.1     0.1')
        parser.add_argument('--ref_t_npt',help='reference temperature, one for each group, in K',default='300     300')
        parser.add_argument('--pcoupl_npt',help='Pressure coupling on in NPT',default='Parrinello-Rahman')
        parser.add_argument('--pcoupltype_npt',help='uniform scaling of box vectors',default='isotropic')
        parser.add_argument('--tau_p_npt',help='time constant, in ps',default='2.0')
        parser.add_argument('--ref_p_npt',help='reference pressure, in bar',default='1.0')
        parser.add_argument('--compressibility_npt',help='isothermal compressibility of water, bar^-1',default='4.5e-5')
        parser.add_argument('--refcoord_scaling_npt',help='refcoord scaling',default='com')
        parser.add_argument('--pbc_npt',help='3-D PBC',default='xyz')
        parser.add_argument('--gen_vel_npt',help='Velocity generation is off',default='no')
        
        #md.mdp
        parser.add_argument('--title_md',help='Title',default='Molecular dynamics (MD)')
        parser.add_argument('--integrator_md',help='leap-frog integrator',default='md')
        parser.add_argument('--nsteps_md',help='2 * 500000 = 1000 ps (1 ns)',default='500000')
        parser.add_argument('--dt_md',help='2 fs',default='0.002')
        parser.add_argument('--nstxout_md',help='suppress bulky .trr file by specifying 0',default='0')
        parser.add_argument('--nstvout_md',help='suppress bulky .trr file by specifying 0',default='0')
        parser.add_argument('--nstfout_md',help='suppress bulky .trr file by specifying 0',default='0')
        parser.add_argument('--nstenergy_md',help='save energies every 10.0 ps',default='5000')
        parser.add_argument('--nstlog_md',help='update log file every 10.0 ps',default='5000')
        parser.add_argument('--nstxoutcompressed_md',help='save compressed coordinates every 10.0 ps',default='5000')
        parser.add_argument('--compressedxgrps_md',help='save the whole system',default='System')
        parser.add_argument('--continuation_md',help='Restarting after NVT',default='yes')
        parser.add_argument('--constraint_algorithm_md',help='holonomic constraints',default='lincs')
        parser.add_argument('--constraints_md',help='bonds involving H are constrained',default='h-bonds')
        parser.add_argument('--lincs_iter_md',help='accuracy of LINCS',default='1')
        parser.add_argument('--lincs_order_md',help='also related to accuracy',default='4')
        parser.add_argument('--cutoffscheme_md',help='Buffered neighbor searching',default='Verlet')
        parser.add_argument('--ns_type_md',help='search neighboring grid cells',default='grid')
        parser.add_argument('--nstlist_md',help='20 fs, largely irrelevant with Verlet',default='10')
        parser.add_argument('--rcoulomb_md',help='short-range electrostatic cutoff (in nm)',default='1.0')
        parser.add_argument('--rvdw_md',help='short-range van der Waals cutoff (in nm)',default='1.0')
        parser.add_argument('--coulombtype_md',help='Particle Mesh Ewald for long-range electrostatics',default='PME')
        parser.add_argument('--pme_order_md',help='cubic interpolation',default='4')
        parser.add_argument('--fourierspacing_md',help='grid spacing for FFT',default='0.16')
        parser.add_argument('--tcoupl_md',help='modified Berendsen thermostat',default='V-rescale')
        parser.add_argument('--tcgrps_md',help='two coupling groups - more accurate',default='Protein Non-Protein')
        parser.add_argument('--tau_t_md',help='time constant, in ps',default='0.1     0.1')
        parser.add_argument('--ref_t_md',help='reference temperature, one for each group, in K',default='300     300')
        parser.add_argument('--pcoupl_md',help='Pressure coupling on in NPT',default='Parrinello-Rahman')
        parser.add_argument('--pcoupltype_md',help='uniform scaling of box vectors',default='isotropic')
        parser.add_argument('--tau_p_md',help='time constant, in ps',default='2.0')
        parser.add_argument('--ref_p_md',help='reference pressure, in bar',default='1.0')
        parser.add_argument('--compressibility_md',help='isothermal compressibility of water, bar^-1',default='4.5e-5')
        parser.add_argument('--pbc_md',help='3-D PBC',default='xyz')
        parser.add_argument('--DispCorr_md',help='account for cut-off vdW scheme',default='EnerPres')
        parser.add_argument('--gen_vel_md',help='Velocity generation is off',default='no')
        
        #forceField
        parser.add_argument('--water_ff',help='Water model to use: select, none, spc, spce, tip3p, tip4p, tip5p, tips3p',default='spce')
        
        #editconf
        parser.add_argument('--bt_ed',help='Box type for -box and -d: triclinic, cubic, dodecahedron, octahedron',default='cubic')
        parser.add_argument('--d_ed',help='Distance between the solute and the box',default='1.0')
        
        #genion
        parser.add_argument('--pname_gen',help='Name of the positive ion',default='NA')
        parser.add_argument('--nname_gen',help='Name of the negative ion',default='CL')
        parser.add_argument('--conc_gen',help='Specify salt concentration (mol/liter).',default='0')
        
        args = parser.parse_args()
        
        for name in paramMDP:
            vars()[name] = getattr(args, name)
            
        for name in paramRest:
            setattr(paramInput,name,getattr(args, name))
        
        f_ions = open('ions.mdp', 'w')
        for name in param_ions:
            if name == 'cutoffscheme_ions':
                f_ions.write('cutoff-scheme = {}\n'.format(eval(name)))
            else:
                f_ions.write(name[:-5] + ' = {}\n'.format(eval(name)))
        f_ions.close()
        
        f_minim = open('minim.mdp', 'w')    
        for name in param_minim:
            if name == 'cutoffscheme_minim':
                f_minim.write('cutoff-scheme = {}\n'.format(eval(name)))
            else:
                f_minim.write(name[:-6] + ' = {}\n'.format(eval(name)))
        f_minim.close()
            
        f_nvt = open('nvt.mdp', 'w')    
        for name in param_nvt:
            if name == 'cutoffscheme_nvt':
                f_nvt.write('cutoff-scheme = {}\n'.format(eval(name)))
            elif name == 'tcgrps_nvt':
                f_nvt.write('tc-grps = {}\n'.format(eval(name)))
            else:
                f_nvt.write(name[:-4] + ' = {}\n'.format(eval(name)))
        f_nvt.close()
                
        f_npt = open('npt.mdp', 'w')    
        for name in param_npt:
            if name == 'cutoffscheme_npt':
                f_npt.write('cutoff-scheme = {}\n'.format(eval(name)))
            elif name == 'tcgrps_npt':
                f_npt.write('tc-grps = {}\n'.format(eval(name)))
            else:
                f_npt.write(name[:-4] + ' = {}\n'.format(eval(name)))
        f_npt.close()
    
        f_md = open('md.mdp', 'w')    
        for name in param_md:
            if name == 'nstxoutcompressed_md':
                f_md.write('nstxout-compressed = {}\n'.format(eval(name)))
            elif name == 'compressedxgrps_md':
                f_md.write('compressed-x-grps = {}\n'.format(eval(name)))
            elif name == 'cutoffscheme_md':
                f_md.write('cutoff-scheme = {}\n'.format(eval(name)))
            elif name == 'tcgrps_md':
                f_md.write('tc-grps = {}\n'.format(eval(name)))
            else:
                f_md.write(name[:-3] + ' = {}\n'.format(eval(name)))
        f_md.close()
        
def Param():
    p = paramInput()
    p()
    p.param()         

#======================================= FILE PREPARATION ==============================================

def removeNonProtein():
    nonProtein = ['ANISOU','HETATM']
    with open(paramInput.file) as oldfile, open('clean.pdb', 'w') as newfile:
        for line in oldfile:
            if not any(s in line for s in nonProtein):
                newfile.write(line)
    
def forceField():
    pdbToGmx = subprocess.Popen(['gmx','pdb2gmx','-f','clean.pdb','-o','struct.gro','-water',paramInput.water_ff], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    stdin_data = pdbToGmx.communicate(input=b'15 \n')[0]

def solv():
    subprocess.run(['gmx','editconf','-f','struct.gro','-o','box.gro','-c','-d',paramInput.d_ed,'-bt',paramInput.bt_ed])
    subprocess.run(['gmx','solvate','-cp','box.gro','-cs','spc216.gro','-o','solv.gro','-p','topol.top'])
     
def addIons():
    subprocess.run(['gmx','grompp','-f','ions.mdp','-c','solv.gro','-p','topol.top','-o','ions.tpr'])
    AddIons = subprocess.Popen(['gmx','genion','-s','ions.tpr','-o','solv_ions.gro','-p','topol.top','-pname',paramInput.pname_gen,'-nname',paramInput.nname_gen,'-conc',paramInput.conc_gen,'-neutral'], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
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
    Param()
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
    runner()
    