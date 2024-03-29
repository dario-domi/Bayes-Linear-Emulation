# Repository Description
This repository provides downloadable R code to train a Bayes Linear Emulator of a complex computer model.

## What's an Emulator? Brief version
Suppose to have a function *f* expensive to evaluate: for example, *f(x)* is the output of 
a simulator, modelling the dynamics of a complex physical system when the system parameters take the values specified by *x* (a vector). 

An **emulator** is a fast statistical surrogate of *f*: once the value of *f* is observed at a small number of inputs
<img src="http://latex.codecogs.com/svg.latex?x_1,&space;\dots,&space;x_n" title="http://latex.codecogs.com/svg.latex?x_1, \dots, x_n" />,
the emulator provides a prediction of *f(x)* at any other input *x*, quantifying the uncertainty associated with the prediction. 
The emulator implemented in this repository relies on Bayes Linear updates to accomplish the task. More details at the end of this README or in `Documentation.pdf` :)

## Main Files
You only need the first of the following three files to train your own emulators using the provided code. The other two files are key to give context to Bayes Linear Emulation (if you are new to that) and for illustration of how to use the provided code.
* **`Emulation.R`**: Defines the R function `BL.Emul`, which builds a Bayes Linear Emulator of an unknown function *f*. It will have to be sourced in your R session.
* **`Documentation.pdf`**: Explains the setting of Bayes Linear Emulation and provides details about the arguments of `BL.Emul`. Read this!
* **`Emulation_Example.R`**: This R script steps you through the process of building emulators using `BL.Emul`. Two emulators of the same function _f_ over a 2D domain are built and compared: one leaves all optional arguments of `BL.Emul` to their defaults, the other makes more tailored choices. The script is densely commented.

## How to Use the Provided Code to Train Your Own Emulator
1. Download the repo: Code -> Download ZIP *(or type "git clone https://github.com/dario-domi/Bayes-Linear-Emulation.git" from terminal)*
2. Copy the file Emulation.R into the folder you launch R from.
3. Type `source("Emulation.R")` in your R session.

The function you need to use is `BL.Emul()`. 
You will need to provide the set of design points
<img src="http://latex.codecogs.com/svg.latex?x_1,&space;\dots,&space;x_n" title="http://latex.codecogs.com/svg.latex?x_1, \dots, x_n" />,
the associated observed values 
<img src="http://latex.codecogs.com/svg.latex?f(x_1),&space;\dots,&space;f(x_n)" title="http://latex.codecogs.com/svg.latex?x_1, \dots, x_n" />,
and the set of test inputs at which the emulator prediction is requested.

The function has several optional arguments, which allow to specify *e.g.* regression terms, correlation lengths, prior variance, nugget, etc. If any is not provided, a default choice is made. The script `Emulation_Example.R` illustrates use of the function on a toy problem, with both default and tailored choices.

# More Details on the Emulator and Meaning of all Optional Arguments
<!--
An overview of the role of each optional argument is provided inside 
-->
This part of the README is being written, please be patient :)

<img width="500" src="https://github.com/dario-domi/Bayes-Linear-Emulation/blob/main/Pictures/Emul_Var.png">
<img width="500" src="https://github.com/dario-domi/Bayes-Linear-Emulation/blob/main/Pictures/Function_f.png">


<!--
## Computational Note
but outperforms a "more natural" nested-loop structure by orders of magnitude. The  
# A bit more detail
-->
