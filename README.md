# Repository Aim and Description
This repository provides downloadable R code to train a Bayes Linear Emulator of a complex computer model.

## Main Content
Only the first file is needed for you to train your own emulators using the provided code. The other two files are however key to give context to Bayes Linear Emulation (if you are new to that) and explain how to use the provided code.
* **`Emulation.R`**: Defines the function `BL.Emul`, which builds a Bayes Linear Emulator of an unknown function f. It will have to be sourced in your R session.
* **`Documentation.pdf`**: Explains the setting of Bayes Linear Emulation and provides details about the arguments of `BL.Emul`. Read this!
* **`Emulation_Example.R`**: This script steps you through the process of building an emulator of a function f over a 2D domain, using `BL.Emul`. Two emulators of the same function are built: one leaves all optional arguments of `BL.Emul` to their defaults, the other makes more tailored choices. The script is densely commented.

## What is an Emulator?
Suppose to have a function f expensive to evaluate: for example, f(x) is the output of 
a simulator, modelling the dynamics of a complex physical system when the system parameters take the values specified by x (a vector). 

An **emulator** is a fast statistical surrogate of f: once the value of f is observed at a small number of inputs
<img src="http://latex.codecogs.com/svg.latex?x_1,&space;\dots,&space;x_n" title="http://latex.codecogs.com/svg.latex?x_1, \dots, x_n" />,
the emulator provides a prediction of f(x) at any other input x, associating with the prediction a measure of its uncertainty. 
The emulator implemented in this repository relies on Bayes Linear updates to accomplish the task.

## How to Use the Provided Code to Train Your Own Emulator
Download the scripts `Emulation.R` and `Corr_fun.R` and `source` them in `R`. The function you need to use is `BL.Emul()`. 
You will need to provide the set of design points
<img src="http://latex.codecogs.com/svg.latex?x_1,&space;\dots,&space;x_n" title="http://latex.codecogs.com/svg.latex?x_1, \dots, x_n" />,
the associated observed values 
<img src="http://latex.codecogs.com/svg.latex?f(x_1),&space;\dots,&space;f(x_n)" title="http://latex.codecogs.com/svg.latex?x_1, \dots, x_n" />,
and the set of test inputs at which the emulator prediction is requested.

The function has several optional arguments, which allow to specify *e.g.* regression terms, correlation lengths, prior variance, nugget, etc. If any is not provided, a default choice is made. The script `Emulation_Example.R` illustrates use of the function on a toy problem, with both default and tailored choices: Plots of the emulator performance are below.

# More Details on the Emulator and Meaning of all Optional Arguments

An overview of the role of each optional argument is provided inside 



## Computational Note
but outperforms a "more natural" nested-loop structure by orders of magnitude. The  

# A bit more detail
