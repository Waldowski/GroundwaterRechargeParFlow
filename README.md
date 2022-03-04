# GroundwaterRechargeParFlow
These routines were published by Bastian Waldowski on the 02.03.2022. 
They allow for postprocessing ParFlow binary files to postprocessed binary files. 
Use 'compileCode.bsh' to compile the Fortran scripts. Start the routines like this: './recharge Input_Script.txt'. 
The name of the input .txt file can be changed. 
The input .txt file must contain information about the model (paths to input files, grid information etc.). Please look at 'Input_Script.txt' as an example input .txt file.
'reshape_3dbin.m' is a function to read and reshape the postprocessed binary files with MATLAB. 
