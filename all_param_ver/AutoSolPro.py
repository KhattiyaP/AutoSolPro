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
        param_ff = ['chainsep_ff','merge_ff','water_ff','angle_ff','dist_ff','posrefc_ff','vsite_ff']
        param_ed = ['bt_ed','d_ed','density_ed','resnr_ed','rvdw_ed','label_ed']
        param_sol = ['radius_sol','scale_sol','shell_sol','maxsol_sol']
        param_gai = ['time_gai','maxwarn_gai']
        param_gem = ['time_gem','maxwarn_gem']
        param_gnvt = ['time_gnvt','maxwarn_gnvt']
        param_gnpt = ['time_gnpt','maxwarn_gnpt']
        param_gmd = ['time_gmd','maxwarn_gmd']
        param_gen = ['np_gen','pname_gen','pq_gen','nn_gen','nname_gen','nq_gen','rmin_gen','seed_gen','conc_gen']
        param_mem = ['npme_mem','nt_mem','ntmpi_mem','ntomp_mem','ntomp_pme_mem','pinoffset_mem','pinstride_mem','rdd_mem','rcon_mem','dds_mem','nstlist_mem','pforce_mem','cpt_mem','nsteps_mem','maxh_mem','replex_mem','nex_mem','reseed_mem']
        param_mnvt = ['npme_mnvt','nt_mnvt','ntmpi_mnvt','ntomp_mnvt','ntomp_pme_mnvt','pinoffset_mnvt','pinstride_mnvt','rdd_mnvt','rcon_mnvt','dds_mnvt','nstlist_mnvt','pforce_mnvt','cpt_mnvt','nsteps_mnvt','maxh_mnvt','replex_mnvt','nex_mnvt','reseed_mnvt']
        param_mnpt = ['npme_mnpt','nt_mnpt','ntmpi_mnpt','ntomp_mnpt','ntomp_pme_mnpt','pinoffset_mnpt','pinstride_mnpt','rdd_mnpt','rcon_mnpt','dds_mnpt','nstlist_mnpt','pforce_mnpt','cpt_mnpt','nsteps_mnpt','maxh_mnpt','replex_mnpt','nex_mnpt','reseed_mnpt']
        param_mmd = ['npme_mmd','nt_mmd','ntmpi_mmd','ntomp_mmd','ntomp_pme_mmd','pinoffset_mmd','pinstride_mmd','rdd_mmd','rcon_mmd','dds_mmd','nstlist_mmd','pforce_mmd','cpt_mmd','nsteps_mmd','maxh_mmd','replex_mmd','nex_mmd','reseed_mmd']
        paramRest = param_file + param_ff + param_ed + param_sol + param_gai + param_gem + param_gnvt + param_gnpt + param_gmd + param_gen + param_mem + param_mnvt + param_mnpt + param_mmd
        
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
        parser.add_argument('--chainsep_ff',help='Condition in PDB files when a new chain should be started (adding termini): id_or_ter, id_and_ter, ter, id, interactive',default='id_or_ter')
        parser.add_argument('--merge_ff',help='Merge multiple chains into a single [moleculetype]: no, all, interactive',default='no')
        parser.add_argument('--water_ff',help='Water model to use: select, none, spc, spce, tip3p, tip4p, tip5p, tips3p',default='spce')
        parser.add_argument('--angle_ff',help='Minimum hydrogen-donor-acceptor angle for a H-bond (degrees)',default='135')
        parser.add_argument('--dist_ff',help='Maximum donor-acceptor distance for a H-bond (nm)',default='0.3')
        parser.add_argument('--posrefc_ff',help='Force constant for position restraints',default='1000')
        parser.add_argument('--vsite_ff',help='Convert atoms to virtual sites: none, hydrogens, aromatics',default='none')
        
        #editconf
        parser.add_argument('--bt_ed',help='Box type for -box and -d: triclinic, cubic, dodecahedron, octahedron',default='cubic')
        parser.add_argument('--d_ed',help='Distance between the solute and the box',default='1.0')
        parser.add_argument('--density_ed',help='Renumber residues starting from resnr',default='1000')
        parser.add_argument('--resnr_ed',help='',default='-1')
        parser.add_argument('--rvdw_ed',help='Default Van der Waals radius (in nm) if one can not be found in the database or if no parameters are present in the topology file',default='0.12')
        parser.add_argument('--label_ed',help='Add chain label for all residues',default='A')
        
        #solvate
        parser.add_argument('--radius_sol',help='Default van der Waals distance',default='0.105')
        parser.add_argument('--scale_sol',help='Scale factor to multiply Van der Waals radii',default='0.57')
        parser.add_argument('--shell_sol',help='Thickness of optional water layer around solute',default='0')
        parser.add_argument('--maxsol_sol',help='Maximum number of solvent molecules to add if they fit in the box',default='0')
        
        #grompp addIons
        parser.add_argument('--time_gai',help='Take frame at or first after this time.',default='-1')
        parser.add_argument('--maxwarn_gai',help='Number of allowed warnings during input processing.',default='0')
        
        #grompp EM
        parser.add_argument('--time_gem',help='Take frame at or first after this time.',default='-1')
        parser.add_argument('--maxwarn_gem',help='Number of allowed warnings during input processing.',default='0')
        
        #grompp NVT
        parser.add_argument('--time_gnvt',help='Take frame at or first after this time.',default='-1')
        parser.add_argument('--maxwarn_gnvt',help='Number of allowed warnings during input processing.',default='0')
        
        #grompp NPT
        parser.add_argument('--time_gnpt',help='Take frame at or first after this time.',default='-1')
        parser.add_argument('--maxwarn_gnpt',help='Number of allowed warnings during input processing.',default='0')
        
        #grompp MD
        parser.add_argument('--time_gmd',help='Take frame at or first after this time.',default='-1')
        parser.add_argument('--maxwarn_gmd',help='Number of allowed warnings during input processing.',default='0')
        
        #genion
        parser.add_argument('--np_gen',help='Number of positive ions',default='0')
        parser.add_argument('--pname_gen',help='Name of the positive ion',default='NA')
        parser.add_argument('--pq_gen',help='Charge of the positive ion',default='1')
        parser.add_argument('--nn_gen',help='Number of negative ions',default='0')
        parser.add_argument('--nname_gen',help='Name of the negative ion',default='CL')
        parser.add_argument('--nq_gen',help='Charge of the negative ion',default='-1')
        parser.add_argument('--rmin_gen',help='Minimum distance between ions',default='0.6')
        parser.add_argument('--seed_gen',help='Seed for random number generator',default='1993')
        parser.add_argument('--conc_gen',help='Specify salt concentration (mol/liter).',default='0')
        
        #MD EM
        parser.add_argument('--npme_mem',help='Number of separate ranks to be used for PME, -1 is guess',default='-1')
        parser.add_argument('--nt_mem',help='Total number of threads to start (0 is guess)',default='0')
        parser.add_argument('--ntmpi_mem',help='Number of thread-MPI threads to start (0 is guess)',default='0')
        parser.add_argument('--ntomp_mem',help='Number of OpenMP threads per MPI rank to start (0 is guess)',default='0')
        parser.add_argument('--ntomp_pme_mem',help='Number of OpenMP threads per MPI rank to start (0 is -ntomp)',default='0')
        parser.add_argument('--pinoffset_mem',help='The lowest logical core number to which mdrun should pin the first thread',default='0')
        parser.add_argument('--pinstride_mem',help='Pinning distance in logical cores for threads, use 0 to minimize the number of threads per physical core',default='0')
        parser.add_argument('--rdd_mem',help='The maximum distance for bonded interactions with DD (nm), 0 is determine from initial coordinates',default='0')
        parser.add_argument('--rcon_mem',help='Maximum distance for P-LINCS (nm), 0 is estimate',default='0')
        parser.add_argument('--dds_mem',help='Fraction in (0,1) by whose reciprocal the initial DD cell size will be increased',default='0.8')
        parser.add_argument('--nstlist_mem',help='Set nstlist when using a Verlet buffer tolerance (0 is guess)',default='0')
        parser.add_argument('--pforce_mem',help='Print all forces larger than this (kJ/mol nm)',default='-1')
        parser.add_argument('--cpt_mem',help='Checkpoint interval (minutes)',default='15')
        parser.add_argument('--nsteps_mem',help='Run this number of steps, overrides .mdp file option (-1 means infinite, -2 means use mdp option, smaller is invalid)',default='-2')
        parser.add_argument('--maxh_mem',help='Terminate after 0.99 times this time (hours)',default='-1')
        parser.add_argument('--replex_mem',help='Attempt replica exchange periodically with this period (steps)',default='0')
        parser.add_argument('--nex_mem',help='Number of random exchanges to carry out each exchange interval',default='0')
        parser.add_argument('--reseed_mem',help='Seed for replica exchange, -1 is generate a seed',default='-1')
        
        #MD NVT
        parser.add_argument('--npme_mnvt',help='Number of separate ranks to be used for PME, -1 is guess',default='-1')
        parser.add_argument('--nt_mnvt',help='Total number of threads to start (0 is guess)',default='0')
        parser.add_argument('--ntmpi_mnvt',help='Number of thread-MPI threads to start (0 is guess)',default='0')
        parser.add_argument('--ntomp_mnvt',help='Number of OpenMP threads per MPI rank to start (0 is guess)',default='0')
        parser.add_argument('--ntomp_pme_mnvt',help='Number of OpenMP threads per MPI rank to start (0 is -ntomp)',default='0')
        parser.add_argument('--pinoffset_mnvt',help='The lowest logical core number to which mdrun should pin the first thread',default='0')
        parser.add_argument('--pinstride_mnvt',help='Pinning distance in logical cores for threads, use 0 to minimize the number of threads per physical core',default='0')
        parser.add_argument('--rdd_mnvt',help='The maximum distance for bonded interactions with DD (nm), 0 is determine from initial coordinates',default='0')
        parser.add_argument('--rcon_mnvt',help='Maximum distance for P-LINCS (nm), 0 is estimate',default='0')
        parser.add_argument('--dds_mnvt',help='Fraction in (0,1) by whose reciprocal the initial DD cell size will be increased',default='0.8')
        parser.add_argument('--nstlist_mnvt',help='Set nstlist when using a Verlet buffer tolerance (0 is guess)',default='0')
        parser.add_argument('--pforce_mnvt',help='Print all forces larger than this (kJ/mol nm)',default='-1')
        parser.add_argument('--cpt_mnvt',help='Checkpoint interval (minutes)',default='15')
        parser.add_argument('--nsteps_mnvt',help='Run this number of steps, overrides .mdp file option (-1 means infinite, -2 means use mdp option, smaller is invalid)',default='-2')
        parser.add_argument('--maxh_mnvt',help='Terminate after 0.99 times this time (hours)',default='-1')
        parser.add_argument('--replex_mnvt',help='Attempt replica exchange periodically with this period (steps)',default='0')
        parser.add_argument('--nex_mnvt',help='Number of random exchanges to carry out each exchange interval',default='0')
        parser.add_argument('--reseed_mnvt',help='Seed for replica exchange, -1 is generate a seed',default='-1')
        
        #MD NPT
        parser.add_argument('--npme_mnpt',help='Number of separate ranks to be used for PME, -1 is guess',default='-1')
        parser.add_argument('--nt_mnpt',help='Total number of threads to start (0 is guess)',default='0')
        parser.add_argument('--ntmpi_mnpt',help='Number of thread-MPI threads to start (0 is guess)',default='0')
        parser.add_argument('--ntomp_mnpt',help='Number of OpenMP threads per MPI rank to start (0 is guess)',default='0')
        parser.add_argument('--ntomp_pme_mnpt',help='Number of OpenMP threads per MPI rank to start (0 is -ntomp)',default='0')
        parser.add_argument('--pinoffset_mnpt',help='The lowest logical core number to which mdrun should pin the first thread',default='0')
        parser.add_argument('--pinstride_mnpt',help='Pinning distance in logical cores for threads, use 0 to minimize the number of threads per physical core',default='0')
        parser.add_argument('--rdd_mnpt',help='The maximum distance for bonded interactions with DD (nm), 0 is determine from initial coordinates',default='0')
        parser.add_argument('--rcon_mnpt',help='Maximum distance for P-LINCS (nm), 0 is estimate',default='0')
        parser.add_argument('--dds_mnpt',help='Fraction in (0,1) by whose reciprocal the initial DD cell size will be increased',default='0.8')
        parser.add_argument('--nstlist_mnpt',help='Set nstlist when using a Verlet buffer tolerance (0 is guess)',default='0')
        parser.add_argument('--pforce_mnpt',help='Print all forces larger than this (kJ/mol nm)',default='-1')
        parser.add_argument('--cpt_mnpt',help='Checkpoint interval (minutes)',default='15')
        parser.add_argument('--nsteps_mnpt',help='Run this number of steps, overrides .mdp file option (-1 means infinite, -2 means use mdp option, smaller is invalid)',default='-2')
        parser.add_argument('--maxh_mnpt',help='Terminate after 0.99 times this time (hours)',default='-1')
        parser.add_argument('--replex_mnpt',help='Attempt replica exchange periodically with this period (steps)',default='0')
        parser.add_argument('--nex_mnpt',help='Number of random exchanges to carry out each exchange interval',default='0')
        parser.add_argument('--reseed_mnpt',help='Seed for replica exchange, -1 is generate a seed',default='-1')
        
        #MD MD
        parser.add_argument('--npme_mmd',help='Number of separate ranks to be used for PME, -1 is guess',default='-1')
        parser.add_argument('--nt_mmd',help='Total number of threads to start (0 is guess)',default='0')
        parser.add_argument('--ntmpi_mmd',help='Number of thread-MPI threads to start (0 is guess)',default='0')
        parser.add_argument('--ntomp_mmd',help='Number of OpenMP threads per MPI rank to start (0 is guess)',default='0')
        parser.add_argument('--ntomp_pme_mmd',help='Number of OpenMP threads per MPI rank to start (0 is -ntomp)',default='0')
        parser.add_argument('--pinoffset_mmd',help='The lowest logical core number to which mdrun should pin the first thread',default='0')
        parser.add_argument('--pinstride_mmd',help='Pinning distance in logical cores for threads, use 0 to minimize the number of threads per physical core',default='0')
        parser.add_argument('--rdd_mmd',help='The maximum distance for bonded interactions with DD (nm), 0 is determine from initial coordinates',default='0')
        parser.add_argument('--rcon_mmd',help='Maximum distance for P-LINCS (nm), 0 is estimate',default='0')
        parser.add_argument('--dds_mmd',help='Fraction in (0,1) by whose reciprocal the initial DD cell size will be increased',default='0.8')
        parser.add_argument('--nstlist_mmd',help='Set nstlist when using a Verlet buffer tolerance (0 is guess)',default='0')
        parser.add_argument('--pforce_mmd',help='Print all forces larger than this (kJ/mol nm)',default='-1')
        parser.add_argument('--cpt_mmd',help='Checkpoint interval (minutes)',default='15')
        parser.add_argument('--nsteps_mmd',help='Run this number of steps, overrides .mdp file option (-1 means infinite, -2 means use mdp option, smaller is invalid)',default='-2')
        parser.add_argument('--maxh_mmd',help='Terminate after 0.99 times this time (hours)',default='-1')
        parser.add_argument('--replex_mmd',help='Attempt replica exchange periodically with this period (steps)',default='0')
        parser.add_argument('--nex_mmd',help='Number of random exchanges to carry out each exchange interval',default='0')
        parser.add_argument('--reseed_mmd',help='Seed for replica exchange, -1 is generate a seed',default='-1')
        
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
    pdbToGmx = subprocess.Popen(['gmx','pdb2gmx','-f','clean.pdb','-o','struct.gro','-chainsep',paramInput.chainsep_ff,'-merge',paramInput.merge_ff,'-water',paramInput.water_ff,'-angle',paramInput.angle_ff,'-dist',paramInput.dist_ff,'-posrefc',paramInput.posrefc_ff,'-vsite',paramInput.vsite_ff], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    stdin_data = pdbToGmx.communicate(input=b'15 \n')[0]

def solv():
    subprocess.run(['gmx','editconf','-f','struct.gro','-o','box.gro','-c','-bt',paramInput.bt_ed,'-d',paramInput.d_ed,'-density',paramInput.density_ed,'-resnr',paramInput.resnr_ed,'-rvdw',paramInput.rvdw_ed,'-label',paramInput.label_ed])
    subprocess.run(['gmx','solvate','-cp','box.gro','-cs','spc216.gro','-o','solv.gro','-p','topol.top','-radius',paramInput.radius_sol,'-scale',paramInput.scale_sol,'-shell',paramInput.shell_sol,'-maxsol',paramInput.maxsol_sol])
     
def addIons():
    subprocess.run(['gmx','grompp','-f','ions.mdp','-c','solv.gro','-p','topol.top','-o','ions.tpr','-time',paramInput.time_gai,'-maxwarn',paramInput.maxwarn_gai])
    AddIons = subprocess.Popen(['gmx','genion','-s','ions.tpr','-o','solv_ions.gro','-p','topol.top','-np',paramInput.np_gen,'-pname',paramInput.pname_gen,'-pq',paramInput.pq_gen,'-nn',paramInput.nn_gen,'-nname',paramInput.nname_gen,'-nq',paramInput.nq_gen,'-rmin',paramInput.rmin_gen,'-seed',paramInput.seed_gen,'-conc',paramInput.conc_gen,'-neutral'], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    stdin_data = AddIons.communicate(input=b'13 \n')[0]


#====================================== EM PART =================================================================


def gromppEm():
    subprocess.run(['gmx','grompp','-f','minim.mdp','-c','solv_ions.gro','-p','topol.top','-o','em.tpr','-time',paramInput.time_gem,'-maxwarn',paramInput.maxwarn_gem])

def mdpEmLauncher():
    subprocess.run(['gmx','mdrun','-v','-deffnm','em','-npme',paramInput.npme_mem,'-nt',paramInput.nt_mem,'-ntmpi',paramInput.ntmpi_mem,'-ntomp',paramInput.ntomp_mem,'-ntomp_pme',paramInput.ntomp_pme_mem,'-pinoffset',paramInput.pinoffset_mem,'-pinstride',paramInput.pinstride_mem,'-rdd',paramInput.rdd_mem,'-rcon',paramInput.rcon_mem,'-dds',paramInput.dds_mem,'-nstlist',paramInput.nstlist_mem,'-pforce',paramInput.pforce_mem,'-cpt',paramInput.cpt_mem,'-nsteps',paramInput.nsteps_mem,'-maxh',paramInput.maxh_mem,'-replex',paramInput.replex_mem,'-nex',paramInput.nex_mem,'-reseed',paramInput.reseed_mem])
    

#====================================== EQ PART =================================================================


def gromppEqNVT():
    subprocess.run(['gmx','grompp','-f','nvt.mdp','-c','em.gro','-r','em.gro','-p','topol.top','-o','nvt.tpr','-time',paramInput.time_gnvt,'-maxwarn',paramInput.maxwarn_gnvt])

def mdpEqNVTLauncher():
    subprocess.run(['gmx','mdrun','-deffnm','nvt','-npme',paramInput.npme_mnvt,'-nt',paramInput.nt_mnvt,'-ntmpi',paramInput.ntmpi_mnvt,'-ntomp',paramInput.ntomp_mnvt,'-ntomp_pme',paramInput.ntomp_pme_mnvt,'-pinoffset',paramInput.pinoffset_mnvt,'-pinstride',paramInput.pinstride_mnvt,'-rdd',paramInput.rdd_mnvt,'-rcon',paramInput.rcon_mnvt,'-dds',paramInput.dds_mnvt,'-nstlist',paramInput.nstlist_mnvt,'-pforce',paramInput.pforce_mnvt,'-cpt',paramInput.cpt_mnvt,'-nsteps',paramInput.nsteps_mnvt,'-maxh',paramInput.maxh_mnvt,'-replex',paramInput.replex_mnvt,'-nex',paramInput.nex_mnvt,'-reseed',paramInput.reseed_mnvt])
    
def gromppEqNPT():
    subprocess.run(['gmx','grompp','-f','npt.mdp','-c','nvt.gro','-r','nvt.gro','-t','nvt.cpt','-p','topol.top','-o','npt.tpr','-time',paramInput.time_gnpt,'-maxwarn',paramInput.maxwarn_gnpt])

def mdpEqNPTLauncher():
    subprocess.run(['gmx','mdrun','-deffnm','npt','-npme',paramInput.npme_mnpt,'-nt',paramInput.nt_mnpt,'-ntmpi',paramInput.ntmpi_mnpt,'-ntomp',paramInput.ntomp_mnpt,'-ntomp_pme',paramInput.ntomp_pme_mnpt,'-pinoffset',paramInput.pinoffset_mnpt,'-pinstride',paramInput.pinstride_mnpt,'-rdd',paramInput.rdd_mnpt,'-rcon',paramInput.rcon_mnpt,'-dds',paramInput.dds_mnpt,'-nstlist',paramInput.nstlist_mnpt,'-pforce',paramInput.pforce_mnpt,'-cpt',paramInput.cpt_mnpt,'-nsteps',paramInput.nsteps_mnpt,'-maxh',paramInput.maxh_mnpt,'-replex',paramInput.replex_mnpt,'-nex',paramInput.nex_mnpt,'-reseed',paramInput.reseed_mnpt])
    

#====================================== MD PART =================================================================


def gromppMD():
    subprocess.run(['gmx','grompp','-f','md.mdp','-c','npt.gro','-t','npt.cpt','-p','topol.top','-o','md_0_1.tpr','-time',paramInput.time_gmd,'-maxwarn',paramInput.maxwarn_gmd])

def mdpMDLauncher():
    subprocess.run(['gmx','mdrun','-deffnm','md_0_1','-npme',paramInput.npme_mmd,'-nt',paramInput.nt_mmd,'-ntmpi',paramInput.ntmpi_mmd,'-ntomp',paramInput.ntomp_mmd,'-ntomp_pme',paramInput.ntomp_pme_mmd,'-pinoffset',paramInput.pinoffset_mmd,'-pinstride',paramInput.pinstride_mmd,'-rdd',paramInput.rdd_mmd,'-rcon',paramInput.rcon_mmd,'-dds',paramInput.dds_mmd,'-nstlist',paramInput.nstlist_mmd,'-pforce',paramInput.pforce_mmd,'-cpt',paramInput.cpt_mmd,'-nsteps',paramInput.nsteps_mmd,'-maxh',paramInput.maxh_mmd,'-replex',paramInput.replex_mmd,'-nex',paramInput.nex_mmd,'-reseed',paramInput.reseed_mmd])
    

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
    
