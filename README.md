This code is written in MATLAB. 

This model evaluates the bulk rheology of a heterogenous shear zone. The frictional viscous model (MT_frictional_viscous_StressBC.m) allows local frictional failure when the local stress in a clast and/or the matrix attains the yield strength. This package also includes the pure viscous model without allowing frictional failure (MT_viscous.m and MT_viscous_StressBC.m).

The executable file 'test. mlx' shows an example of the frictional-viscous model to calculate the bulk rheology of a shear zone containing 200 ellipsoidal amphibolite clasts within an interconnected phyllosilicate-dominated matrix.
Please include the 'Routine' folder and the clasts' shape and orientation installed in the 'ellipsoidals_a_q.mat' file before running the 'test. mlx' file.
