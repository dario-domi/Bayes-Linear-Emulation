# Repository Description
This repository provides downloadable R code to train a Bayes Linear emulator of a complex computer model.
## What is an emulator? - In short
Suppose to have a function f expensive to evaluate: for example, f(x) is the output of 
a simulator, modelling the dynamics of a complex physical system when the system parameters take the values specified by x (a vector). 
An emulator is a fast statistical surrogate of f: once the value of f is observed at a small number of inputs
<img src="http://latex.codecogs.com/svg.latex?x_1,&space;\dots,&space;x_n" title="http://latex.codecogs.com/svg.latex?x_1, \dots, x_n" />,
the emulator provides a prediction of f(x) at any other input x, together with a measure of uncertainty of the prediction.
