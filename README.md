# 3D_Synthesis_based_on_Fine-grained_Parts
Code for synthesizing 3D Man-made Shapes using fine-grained segments. This is the code for the paper: 
[
Gonzalez D. and van Kaick Oliver. "3D Synthesis of Man-made Objects based on Fine-grained Parts". Computers & Graphics, (Shape Modeling International - SMI). 2018. In Press.
]

*********
GENERAL INFORMATION:
Our code requires pre-processed data to work. Just like our paper, the code in this repository is focused on using a set of fine-grained segments for allowing a user to synthesize different 3D shapes of 3 datasets: Chairs, Tables, and a Hybrid set.

Therefore, we provide the code for synthesizing 3D models, employing pre-processed data. That is, our code uses segments previously computed to generate new 3D shapes.

1) We provide the pre-processed data in the form of several matlab files (.mat files). 
NOTE: some of these .mat files are large (more than 50 Mb).

2) We provide also the Matlab and C++ code used by our synthesis application.

For questions about the code, pleas contact Diego Gonzalez: diego.gonzalez@carleton.ca
*********

*********
//---
USAGE INSTRUCTIONS FOR THE CODE IN THIS REPOSITORY:
1) Requirements:
(*) Platform: Windows. We have tested the code on Windows 10.

Required Software:
1.1) MATLAB: we have tested the code using Matlab R2016 (R).
1.2) C++: a part of the code uses C++.
Our code to compute and Object Oriented Bounding Box (OBB) for each one of our pre-computed fine-grained segments is based on the method by Fish et al.: 
" Fish, N, Averkiou, M, Van Kaick, O, Sorkine-Hornung, O, Cohen-Or, D, Mitra, NJ. Meta-representation of shape families. ACM Trans on
Graphics 2014;33(4):34:1–34:11. "

WE PROVIDE PRECOMPILED CODE
//---

//---
USE CASES FOR OUR CODE:
1) Replicate the results
//---
