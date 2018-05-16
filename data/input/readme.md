This folder contains all the input shapes employed for the experiments described in our paper.

% ------- %

IMPORTANT NOTE:
These 3D shapes are included here for completeness and in case you want to use the same input datasets that we have employed. 

**For synthesizing 3D shapes using our method, you can use the pre-processed data and the scripts provided in this repository: see the Readme file on the root of the repository for details.**

% ------- %

This folder is divided into two subfolders that contain the 3D models of the datasets described in "Section 6 - Experimental Results" of our paper:

1) Chairs subfolder: a family of 30 chair 3D models.

2) Tables subfolder: a family of 15 table 3D models.

In each of these two subfolders you can find two types of files:

(a) .obj files: containing each input 3D mesh.

(b) .ply files: containing each fine-grained segmentation that we obtained from the corresponding 3D mesh

NOTE: The hybrid set mentioned in our paper combines all the table models with a subset of 15 chair models, which correspond to the following .obj files: ChairArms1.obj, ChairArms10.obj, ChairArms2.obj, ChairArms4.obj, ChairArms9.obj, ChairNoArms1.obj, ChairNoArms12.obj, ChairNoArms17.obj, ChairNoArms18.obj, ChairNoArms19.obj, ChairNoArms2.obj, ChairNoArms4.obj, ChairNoArms5.obj, ChairNoArms8.obj, ChairNoArms9.obj

% ------- %
NOTES ABOUT THE FIGURES IN OUR PAPER:

(1) The input 3D shapes shown in the figures of our paper in gray color were generated from to the following .obj files:

  - Figure 2 (inset), and Figure 9 (Template-Bottom): Table11.obj.
  
  - Figure 6(a): Table7.obj; Figure 6(b): Table4.obj; Figure 6(c): Table15.obj; Figure 6(d): Table3.obj; Figure 6(e): Table9.obj; Figure 6(f): Table9.obj
  
  - Figure 7(a): ChairArms1.obj; Figure 7(b): ChairNoArms12.obj; Figure 7(c): ChairNoArms6.obj; Figure 7(d): ChairNoArms13.obj; Figure 7(e): ChairArms10.obj; Figure 7(f): ChairNoArms18.obj
  
  - Figure 8 (Template-Bottom): ChairNoArms19.obj
  
  - Figure 10(a): ChairArms1.obj; Figure 10(b): Table6.obj; Figure 10(c): Table15.obj; Figure 10(d): ChairArms4.obj; Figure 10(e): ChairArms10.obj; Figure 10(f): Table10.obj
  
  - Figure 11(a): Table8.obj; Figure 11(b): Table6.obj; Figure 11(c): ChairArms3.obj; Figure 11(d): ChairArms9.obj; Figure 11(e): Table10.obj; Figure 11(f): ChairArms2.obj
  
(2) The fine-grained segmentations shown in the figures of our paper were generated from to the .ply files corresponding to each 3D mesh mentioned above. For example:

  - Figure 2 (right), and Figure 9 (Template-Top): seg-Table11.ply
  - Figure 6 (a): seg-Table7.ply; Figure 6(b): seg-Table4.ply; and so on.
  - The fine-grained segmentation shown in Figure 4 corresponds to the file seg-ChairArms8.ply
% ------- %

