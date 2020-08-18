## Authors

**Laura González Antiga** - [Linkedin](https://www.linkedin.com/in/laura-g-antiga/)
**Olfat Khannous Lleiffe** - [Linkedin](https://www.linkedin.com/in/olfat-khannous-lleiffe-155a22122/)
**Marina Murillo Recio** - [Linkedin](https://www.linkedin.com/in/marina-murillo-recio-986168b2/)

# PYSBI

Given a set of interacting pairs (prot-prot) reconstruct the complete macrocomplex.

The following file contains a description of the standalone application (PYSBI) that solve the specific problem proposed in the Python and Structural Bioinformatics’s subjects of the Master in Bioinformatics for Health Sciences at Pompeu Fabra University.

This information will help you to install (installation) our program and to understand how to use it properly (tutorial).

## INSTALLATION

1-Download the macrocomplex_rec.tar.gz folder

2-Unzip the folder:  tar -xvzf  macrocomplex_rec.tar.gz

```
tar -xvzf  macrocomplex_rec.tar.gz
```

## FOLDER CONTENT

•    modules_arg.py

•    maincompl.py

•    functions.py

Several folders in order to try the program:
•    His3

•    1FN3

Brief explanation of each of the files.py:

-modules_arg.py: contains the modules that have to be imported and the argparser that makes the program easy to write user-friendly command-line interfaces.

-maincompl.py: contains the main program.

-functions.py: all the functions that are implemented to solve the problem.


## COMMAND-LINE ARGUMENTS

--input (-i): Name of the directory that contains the files.

--output (-o): Directory that contains the set of outputs that are obtained (pdb files with the resulting macrocomplexes).

--verbose (-v): Using this option the user can follow step by step what the program is doing (It’s printed in the Standard Error).

--distance (-d): Interaction distance of the chains of the input pdb files. Used to filter the input. The default value is 6.


## COMMAND LINE APPLICATION

The following command line is an example of the simplest way to run our program correctly:

```
python3 maincompl.py –i example1 –o OUTPUT
```

If the user wants to see brief explanations of the steps that the program follows to obtain the output can add the verbose argument (-v) and these descriptions will appear in the standard error:
```
python3 maincompl.py –i example1 –o OUTPUT –v
```
Another optional argument is the distance (-d). Our program filters and uses the input pdbs that have 2 interacting chains using a default distance of 6. If the user wants to use another one can use this argument to change it:
```
python3 maincompl.py –i example1 –o OUTPUT –v –d 8
```







