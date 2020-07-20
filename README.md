# AutoSolPro
### Automated molecular dynamics protocols for soluble proteins
## Getting Started
### Prerequisites
* [GROMACS](http://www.gromacs.org/) needs be installed.
* If you use Windows, you may run by Ubuntu in Windows Subsystem for Linux.
### How to Use
1. Create a directory (no space) and put there all the files:
    * AutoSolPro.py
    * ions.mdp
    * md.mdp
    * minim.mdp
    * npt.mdp
    * nvt.mdp
    * yourfile.pdb (your .pdb file)
    
2. Change the directory in Command Prompt to your folder and type
```
python3 AutoSolPro.py yourfile.pdb
```

## Authors
**[Khattiya Pongsirijinda](mailto:Khattiya.Pongsirijinda@skoltech.ru)** - Skolkovo Institute of Science and Technology

## References
* [GROMACS Tutorials](http://www.mdtutorials.com/gmx/)
* [MYGENerator (Polish)](https://mygen.wbbib.uj.edu.pl/mygenerator)
