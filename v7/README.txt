Hello and welcome to the DB GCS RS theta model

This should be version 7. You can tell because the folder (directory) is 'ThetaModel/v7' . Is this correct?

In this directory (...v7/) there should be:

6 .png files
3 .m files
1 .mat file


The .png files:
are exemplary data sets to try. If you don't know where to begin, try one of those. 
Note about datasetC... the conditions that made that simulation allowed n to go above 1 (y to go below 0). So it's useless, but also useful! You decide.


The .mat file:
Stores the last values of the previous simulation. Just ignore it. To use those values as the initial conditions of the next simulation, make sure you specifiy 'istate = 4' in your call to synctheta_v7


The .m files:
The main.m file is the file to run. It kicks off the simulation.
The synctheta_v7.m file is the simulation itself. Go here to change most parameters.
The DBPlot_v2.m file is a function to plot the results. If you are interested in making new plots, you can edit this script or write your own.


To run:
Open main.m
hit Run

To change params:
Open synctheta_v7
Scroll to User Params section (currently line 58)
play with params
Remember, run main to run simulation

The script will automatically save the simulation figure as a .png