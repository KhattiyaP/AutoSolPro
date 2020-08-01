# AutoSolPro
### Automated molecular dynamics protocols for soluble proteins
## Getting Started
### Prerequisites
* [GROMACS](http://www.gromacs.org/) needs to be installed.
* If you use Windows, you may run by Ubuntu in Windows Subsystem for Linux.
### How to Use
1. Create a directory (no space) and put there AutoSolPro.py and yourfile.pdb (your protein file).
    
2. Change the directory in Command Prompt to your folder and type
```
python3 AutoSolPro.py --file 'yourfile.pdb'
```
and if you would like to change the parameters from the default, they can use the command, for example,
```
python3 AutoSolPro.py --file 'yourfile.pdb' --nsteps_md '5000'
```
You can find out more about the parameters by the command,
```
python3 AutoSolPro.py --help'
```

## Author
**[Khattiya Pongsirijinda](mailto:Khattiya.Pongsirijinda@skoltech.ru)** - Skolkovo Institute of Science and Technology

## References
* [GROMACS Tutorials](http://www.mdtutorials.com/gmx/)
* [MYGENerator (Polish)](https://mygen.wbbib.uj.edu.pl/mygenerator)
