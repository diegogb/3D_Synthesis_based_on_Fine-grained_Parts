%addpath('obb');
try
    currDir= pwd();
    newDir= fullfile(currDir, 'obb');
    if (exist(newDir, 'dir'))
        cd(newDir)
        mex -O -I. symOBB_PointSet.cpp Wm5ContBox3.cpp Wm5DistPoint3Triangle3.cpp Wm5MathematicsPCH.cpp Wm5Distance.cpp Wm5Vector3.cpp Wm5Math.cpp Wm5Assert.cpp Wm5ApprGaussPointsFit3.cpp Wm5Quaternion.cpp Wm5EigenDecomposition.cpp Wm5Matrix2.cpp Wm5Matrix3.cpp Wm5SingularValueDecomposition.cpp Wm5Memory.cpp Wm5CorePCH.cpp Wm5Mutex.cpp Wm5Matrix4.cpp Wm5Vector4.cpp Wm5Vector2.cpp Wm5Query.cpp Wm5ContMinBox2.cpp Wm5ContSymBox3.cpp Wm5ConvexHull2.cpp Wm5ConvexHull1.cpp Wm5ConvexHull.cpp Wm5ConvexHull3.cpp Wm5FileIO.cpp Wm5Endian.cpp;
        mex -O -I. sym_plane.cpp Wm5DistPoint3Triangle3.cpp Wm5MathematicsPCH.cpp Wm5Distance.cpp Wm5Vector3.cpp Wm5Math.cpp Wm5Assert.cpp Wm5FileIO.cpp Wm5Endian.cpp;
		cd(currDir);
    end
catch ME
    errordlg(['Error compiling the MEX function "symOBB_PointSet" ' ME.message]);
end