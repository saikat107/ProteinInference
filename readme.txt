ProteinInference
-------------------------------------------------------

Open the project in nettbeans.
sigmaTestData Contails several files for running the test
Test file name can be changed from the Configuration.java file.

Detailed Description:
1. protein_peptide_map.csv and protein_peptide_map_new.csv these 3 files represent the mapping between protein and peptides which I produced from parsing Sigma49 dataset. You will see the file in the following format
a,b ; where both a and b are integers. a represents a protein and b represents peptide. and each (a,b) pair is a relation between protein and peptide.
2. The above mentioned two files acts as a database in my code.
3. In SigmaTestData folder there are 18 test files. Test file format is.
pr1,pr2,pr3................... : pep1, pep2, pep3,..................
where the code works pep1,pep2,pep3,....... is the peptide set and pr1,pr2,pr3,....... is the expected protein set. In the code, we tried to infer proteins from the peptide set and evaluate the result with the expected protein set.
4. In the code, I used memtic algorithm. which first builds a candidate protein set from the above mentioned database. 
5. Then I generated initial random population population.
6. I started memetic algorithm then. Parent selection procedures are a) RouletteWheel Selection, b) Fitness uniform selction(FUSS). Here is a short introduction of fuss (http://www.hutter1.de/ai/sfuss.pdf)
7. Mutation was random.
8. I implemented both Simulated Annealing and HillClimbing as local improvement procedure. Simulated annealing seems to be working slower. 
9. For evaluation function, please notice the MAgPI.pdf in this folder.