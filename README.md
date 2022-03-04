# GroundwaterRechargeParFlow
These routines were published by Bastian Waldowski on the 02.03.2022. 
They allow for postprocessing ParFlow binary files to postprocessed binary files. 
Use 'compileCode.bsh' to compile the Fortran scripts. Start the routines like this: './recharge Input_Script.txt'. 
The name of the input .txt file can be changed. 
The input .txt file must contain information about the model (paths to input files, grid information etc.). Please look at 'Input_Script.txt' as an example input .txt file.
In the input .txt file:
	- The grid is defined top to bottom
	- Write 'TIME' to refer to timestep in ParFlow index format (e.g. '00001')
	- Write 'ENSE' for ensemble member in ParFlow index format (e.g. '00001')
	- Write 'NUM' for straight number of ensemble member (e.g. '1')
'reshape_3dbin.m' is a function to read and reshape the postprocessed binary files with MATLAB. 
