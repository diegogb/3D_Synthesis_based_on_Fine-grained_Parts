# 3D_Synthesis_based_on_Fine-grained_Parts
Code for synthesizing 3D Man-made Shapes using fine-grained segments. This is the code for the paper: 
[
Gonzalez D. and van Kaick Oliver. "3D Synthesis of Man-made Objects based on Fine-grained Parts". Computers & Graphics, Special Issue on Shape Modeling International SMI. 2018. In Press.
]

Please remember to cite our paper if you find useful our code and/or data for your own project.

*********
GENERAL INFORMATION:
The code in this repository is focused on using a set of fine-grained segments for allowing a user to synthesize different 3D shapes of 3 datasets: Chairs, Tables, and a Hybrid set.

Therefore, we provide the code for synthesizing 3D models using our approach, and employing pre-processed data. That is, our code uses fine-grained segments previously computed to generate new 3D shapes.

1) We provide the pre-processed data in the form of several matlab files (.mat files). 
NOTE: some of these .mat files are large (more than 50 Mb).

2) We provide also the Matlab and C++ code used by our synthesis application.

3) To synthesize 3D shapes using our method you require: all the files in the root of this repository AND the files contained in the folder called "OBB"

4) We also include in this repository the shapes of our experimental input sets, as well as the synthesized shapes that we present in our paper. All the data is inside the folder "data".
NOTE: some of the files contained in the data folder are large (more than 50 Mb)
Please check the file "readme.md" inside the folder "data" for details about the
NOTE: the files contained in the "data" folder are NOT required to execute our code.


For questions about the code, you may contact Diego Gonzalez: diego.gonzalez@carleton.ca (or diegogb2007@gmail.com)
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

The C++ code for computing an OBB is within the folder OBB
//---

//---
EXECUTING THE CODE:

USE CASES FOR OUR CODE:
1) Replicate the results
//---
