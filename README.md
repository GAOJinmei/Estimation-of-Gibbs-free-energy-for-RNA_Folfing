# Creation of an objective function for the RNA folding problem
## Authors - Jinmei GAO


For a given ribonucleotide chain, the RNA folding problem consists in finding the native fold among the astronomically large number of possible conformations. The native fold being the one with the lowest Gibbs free energy, the objective function should be an estimator of this energy.

<sub>School project M2 GENIOMHE paris saclay - Structural Bioinformatics - speaker Guillaume Postic</sub>


## Goals 

 3 Python scripts will be implemented :
 
* train the objective function, using interatomic distance distributions that are computed from a dataset of known (10 RNA structures) (i.e., experimentally determined) 3D structures;
the 10 RNA structures are 10 pdb files which you can find in the pdb folder

* plot the scoring profiles, i.e. the score (or estimated Gibbs free energy) as a function of the interatomic distance;

* use the objective function to evaluate predicted structures from the RNA-Puzzles dataset.
 

## Overall Constraint 

* Only C3' atoms are taken into account

* Only consider intrachain distances are considered

* Only consider residues separed by 3 position on the sequence

* Compute observed frequencies ( range of 0 to 20 Ångström of observing two residues i and j separared by a distance bin r 


 ## Training script : Compute the interactomic distance from a given set of PDB files : 
                     
* Compute observed frequencies, range of 0 to 20 Ångström of observing two residues i and j separared by a distance bin r is calculated as follows :

$$ f_{i,j} ^{OBS}(r) = { N_{i,j}(r) \over N_{i,j} } $$
                      
 Where $N_{i,j}(r)$ is the count of i and j within the distance bin r and $N_{i,j}$ is the count of i and j for all distance bins 
 
                      
* Compute reference frequency , different residue types are  indistinct X                    
 
 $$ f_{X,X} ^{REF}(r) = { N_{X,X}(r) \over N_{X,X}} $$
 

* Compute log-ratio of the two frequency (score of pseudo energy)

$$ u_{i,j}(r) = { -log \left( f _{i,j} ^{OBS}(r) \over f_{i,j} ^{REF}(r) \right) } $$


Training script should therefore generate 10 files of 20 lines. Maximum scoring value will be arbitrarily set to 10 


## Testing script : Evaluate predicted structures : 

It will compute all the distances for a given structure using the same treshold as the training script. For each distance, ascoring value will be computed using a linear interpolation as follow :

$$ Linear Interpolation(Y) = { y1 + (x - x1) \left (y2 - y1 \over x2 - x1 \right) } $$


## Package required 
These  packages are required to run the script 

```javascript
conda install -c anaconda numpy 
```
or

```bash
pip install numpy 
pip install pandas

```
and 
```javascript
conda install -c anaconda numpy 
```

## Clone the repository :

Clone the repository :
```bash
$ git clone https://github.com/hamidsta/RNA-folding-problem.gitjinmei
```


## Running training script1.py
check how to use this script

```bash
  python script1.py --help
```
using the script by using the directory where keeps all the pdb files
```bash
  python script1.py -d 'pdb/'
```

## Running ploting script2.py

To run script2, run the following command

```bash
  python script2.py
```

## Running testing script3.py
check how to use this script

```bash
  python script3.py --help
```
using the script by using a RNA structure as input (a pdb file)
```bash
  python script2.py -pdb '2lx1.pdb'
```

## Data availability 

Data are present inside the !!!!! put folder name !!!!! 
And can be download through the following link :  https://www.rcsb.org/ 





