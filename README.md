# Short Description
This repository provides downloadable R code to train a Bayes Linear emulator of a complex computer model.

## What is an Emulator?
Suppose to have a function f expensive to evaluate: for example, f(x) is the output of 
a simulator, modelling the dynamics of a complex physical system when the system parameters take the values specified by x (a vector). 

An **emulator** is a fast statistical surrogate of f: once the value of f is observed at a small number of inputs
<img src="http://latex.codecogs.com/svg.latex?x_1,&space;\dots,&space;x_n" title="http://latex.codecogs.com/svg.latex?x_1, \dots, x_n" />,
the emulator provides a prediction of f(x) at any other input x, associating with the prediction a measure of its uncertainty. 
The emulator implemented in this repository relies on Bayes Linear updates to accomplish the task.

## How to Train Your Own Using the Provided Code
Download the scripts `emulation.R` and `corr_function.R` and `source` them in `R`. The function you need to use is `BL.Emul()`. 
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
