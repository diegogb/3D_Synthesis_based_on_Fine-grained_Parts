# 3D_Synthesis_based_on_Fine-grained_Parts
Code for synthesizing 3D Man-made Shapes using fine-grained segments. This is the code for the work: 

**Gonzalez D. and van Kaick Oliver. "3D Synthesis of Man-made Objects based on Fine-grained Parts". Computers & Graphics, Special Issue on Shape Modeling International SMI. 2018. In Press.**

Please cite our publication if you use our code and/or for your own project.

*********
**GENERAL INFORMATION:**
We provide the code for synthesizing 3D models using our approach and employing pre-processed data. That is, our code uses fine-grained segments previously computed to generate new 3D shapes.

1) We provide the pre-processed data in the form of several matlab files (.mat files). 
NOTE: some of these .mat files are large (more than 50 Mb in some cases).

2) We provide also the Matlab and C++ code used by our synthesis application.

3) To synthesize 3D shapes using our method you require: all the files in the root of this repository AND the files contained in the folder called "OBB"

4) This repository also includes a folder called "data", which contains the shapes of our experimental input sets, as well as the synthesized shapes presented in our paper. Please check the file "readme.md" inside the folder "data/input" for details about our input sets.

NOTES: some of the files contained in the "data" folder are large ZIP files (more than 20Mb in some cases).

The files contained in the "data" folder are NOT required to execute our code.


- For questions about the code, you may contact Diego Gonzalez: diego.gonzalez@carleton.ca
*********

*********
//---
** USAGE INSTRUCTIONS FOR THE CODE IN THIS REPOSITORY: **
1) Requirements:

- Platform: Windows. We have tested the code on Windows 10.

Required Software:

1.1) MATLAB: we have tested our code using Matlab R2016a (R).

1.2) Our code to compute an Object Oriented Bounding Box (OBB) is writen in C++ and is based on the method by Fish et al.: 

Fish, N, Averkiou, M, Van Kaick, O, Sorkine-Hornung, O, Cohen-Or, D, Mitra, NJ. Meta-representation of shape families. ACM Trans on
Graphics 2014;33(4):34:1–34:11.

The C++ code for computing an OBB is inside the folder "OBB".
//---

//---
** EXECUTING THE CODE:**
1) Download all the code and pre-processed data (.mat files) from this repository into any folder in your machine.

2) Open Matlab in your computer and set the Matlab's path to the folder containing our code.

3) Compiling the C++ code:

We provide pre-compiled MEX files for our C++ code, which we compiled on Matlab R2016a. (We employed "Microsoft Visual C++ 2015 comunity edition" as C++ compiler).

You may need to recompile the code and generate new .MEX files for your machine.

To generate new .MEX files in your machine:

3.1) You must have installed a C++ compiler compatible with Matlab. 

There are several free C++ compiler that you can use. For a complete updated list of C++ compilers compatible with Matlab, you can check the following website:
https://www.mathworks.com/support/compilers.html

3.2) Open and run the script make.m that is inside the folder "OBB".

---
**
4) RUNNING THE SCRIPTS:

4.1) To replicate our paper's figures: on matlab, run the script files called "FigN.m" to generate the 3D shapes (point cloud and mesh) for the corresponding figure.

4.2) To generate timing statistics: on matlab, run the script files called: 

4.2.1) SynthesizeChairCollection_timings.m: to generate a .txt file with execution times for our input dataset of chairs. (This script will also synthesize a new collection of 3D chair models).

4.2.2) SynthesizeTableCollection_timings.m: to generate a .txt file with execution times for our input dataset of tables. (This script will also synthesize a new collection of 3D table models).

4.2.3) SynthesizeHybridCollection_timings.m: to generate a .txt file with execution times for our hybrid input dataset. (This script will also synthesize a new collection of 3D shapes using our hybrid set).


4.3) To generate new collections of 3D shapes (without recording times):

4.3.1) SynthesizeChairCollection.m: synthesize a new collection of 3D chair models.

4.2.2) SynthesizeTableCollection.m: synthesize a new collection of 3D table models.

4.2.3) SynthesizeHybridCollection.m: synthesize a new collection of 3D shapes using our hybrid set.
**

Each script contains detailed comments and usage instructions.
//---
*********
