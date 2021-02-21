# Short Overview
This repository provides downloadable R code to train a Bayes Linear emulator of a complex computer model.
## What is an emulator?
Suppose to have a function f(x) which is expensive to evaluate: for example, f is represented by computer code
modelling the dynamics of complex physical systems according to the values of system parameters x (a vector). 
An emulator is a fast statistical surrogate of f: once the value of f is observed at a small number of inputs
<img src="http://latex.codecogs.com/svg.latex?x_1,&space;\dots,&space;x_n" title="http://latex.codecogs.com/svg.latex?x_1, \dots, x_n" />,
the emulator provides a prediction of f(x) at any other input x, associating with theprediction a measure of its uncertainty.
