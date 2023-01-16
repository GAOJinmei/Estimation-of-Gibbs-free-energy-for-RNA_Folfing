
# Creation of an objective function for the RNA folding problem


For a given ribonucleotide chain, the RNA folding problem consists in finding the native fold among the astronomically large number of possible conformations. The native fold being the one with the lowest Gibbs free energy, the objective function should be an estimator of this energy.




## Authors

- Jinmei GAO



# Goals
I will implement 3 Python scripts to:
- train the objective function, using interatomic distance distributions that are computed from a dataset of known (10 RNA structures) (i.e., experimentally determined) 3D structures;
the 10 RNA structures are 10 pdb files which you can find in the pdb folder

- plot the scoring profiles, i.e. the score (or estimated Gibbs free energy) as a function of the interatomic distance;

- use the objective function to evaluate predicted structures from the RNA-Puzzles dataset.





## Usage/Examples
First these packages are required


```javascript
conda install -c anaconda numpy 
or
pip install numpy


pip install pandas

conda install -c conda-forge matplotlib
```


## script1.py
check how to use this script

```bash
  python script1.py --help
```
using the script by using the directory where keeps all the pdb files
```bash
  python script1.py -d 'pdb/'
```
```


## script2.py

To run script2, run the following command

```bash
  python script2.py
```


## script3.py
check how to use this script

```bash
  python script3.py --help
```
using the script by using a RNA structure as input (a pdb file)
```bash
  python script2.py -pdb '2lx1.pdb'
```
```

