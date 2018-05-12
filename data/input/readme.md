This folder contains all the input shapes employed in the experiments of our paper. 
% -------
IMPORTANT NOTE:
This input 3D shapes are included for convenience only.
FOR SYNTHESIZING 3D SHAPES USING OUR APPROACH AND OUR INPUT 3D SHAPES, YOU CAN LOAD THE .mat FILES INCLUDED IN THE ROOT OF THE data FOLDER:
SEE THE Readme.md FILE ON THE ROOT OF THIS REPOSITORY FOR DETAILS
% -------

The folder is divided into two subfolders for two of the three datasets described in "Section 6 - Experimental Results" of our paper:
1) Chairs subfolder: a family of 30 chair 3D models.
2) Tables subfolder: a family of 15 table 3D models

Note that the hybrid set mentioned in our paper combines all the table models with a subset of 15 chair models.

In each of these two subfolder you can find two types of files:
(*) .obj files: containing each input 3D mesh
(*) .ply files: containing each fine-grained segmentation that we obtained from the corresponding 3D mesh


% -------
NOTES:
(1) The input 3D shapes shown in the figures of our paper in gray color were generated from to the following .obj files:
  (*) Figure 2 (inset), and Figure 9 (Template-Bottom): Table11.obj 
  (*) Figure 6(a): Table7.obj; Figure 6(b): Table4.obj; Figure 6(c): Table15.obj; Figure 6(d): Table3.obj; Figure 6(e): Table9.obj; Figure 6(f): Table9.obj;
  (*) Figure 7(a): ChairArms1.obj; Figure 7(b): ChairNoArms12.obj; Figure 7(c): ChairNoArms6.obj; Figure 7(d): ChairNoArms13.obj; Figure 7(e): ChairArms10.obj; Figure 7(f): ChairNoArms18.obj; 
  (*) Figure 8 (Template-Bottom): ChairNoArms19.obj
  (*) Figure 10(a): ChairArms1.obj; Figure 10(b): Table6.obj; Figure 10(c): Table15.obj; Figure 10(d): ChairArms4.obj; Figure 10(e): ChairArms10.obj; Figure 10(f): Table10.obj;
  (*) Figure 11(a): Table8.obj; Figure 11(b): Table6.obj; Figure 11(c): ChairArms3.obj; Figure 11(d): ChairArms9.obj; Figure 11(e): Table10.obj; Figure 11(f): ChairArms2.obj; 
  
(2) The fine-grained segmentations shown in the figures of our paper were generated from to the .ply files corresponding to each 3D mesh mentioned above. For example:
  (*) Figure 2 (right), and Figure 9 (Template-Top): seg-Table11.ply  
  (*) Figure 6 (a): seg-Table7.ply; Figure 6(b): seg-Table4.ply;
  etc.
% -------
