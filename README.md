# Master-thesis
All work related to my master thesis will be uploaded here

## Structure
Everything is found in the documantation subfolder *doc*, which again contains the folders *data*, *plots*, *scripts*, *src* and *text*. 

- *data* will typically contain text files with values, written by the program. An example will be energy values from VMC runnings.
- *plots* is the folder where one can find plots produced by the Python scripts. 
- *scripts* is where all the additional scripts are placed. Those are typically Python scripts which are not communicating with the main code.
- *src* is the place where the main source code can be found. This will most likely be written in C++.
- *text* is the final subfolder, where the Master Thesis PDF and TEX-file gonna be. 

## Requirements
The source code was created considering legibility and to be easy to run. The only package that was used which is not a standard C++ package is the Eigen package, but everything you need can be found in the source folder (everything is already setup). 

When it comes to the scripts, you gonna need to install
- Numpy
- Matplotlib.pyplot

## How to run?
All parameters that should be changed are found in *main.cpp*. Basically run *main.py* and the magic will happen.
