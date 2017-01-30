# morphogen-curve-fitting
Matlab functions for curve-fitting to morphogen gradient data.

Description
===========

This repository contains the Matlab code that was used to perform curve fitting on morphogen gradients. Specifically, three types of functions were fitted to the data:
- decaying exponential
- two-domain model
- two-domain-gradual-sink model

Descriptions of these functions, and the rationale behind fitting them, are described in the document 'supplemental_modeling_notes.pdf. The models were developed in the Mathematica notebook 'dpp_gradient_modeling.nb' (developed in Mathematica 10.0.0.0). 

The main function is 'morphogenGradientCurveFitting.m'. Documentation on function usage is provided in the comments at the top of the function source code.

'curveFittingExamples.m' provides examples demonstrating function usage.

The functions were tested in Matlab R2012b.
