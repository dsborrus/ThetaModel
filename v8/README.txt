Hello and welcome to the DB GCS RS theta model

This should be version 8. You can tell because the folder (directory) is 'ThetaModel/v8' . Is this correct?

In this directory (...v8/) there should be:

4 .m files
several folders

Each folder represents a different "study"

For new users, try opening JustForKicks. This is the sandbox folder.

Inside is a main.m file and a lastconditions.mat file.


The main.m file is the jumping off point for running the simulation. So open that up. 
Here we can play with params, and then run the simulation. 
All experiment folders have a main.m file, and they all call the same simulation function, that is v8/simulate_v*.m.


The .m files:
The main.m file is the file to run. It kicks off the simulation.
The simulate_v1.m file is the simulation itself.
The plot.m file is a function to plot the results. If you are interested in making new plots, you can edit this script or write your own.
updategif.m is another figure making file. It was for making a gif!


To run:
Open JustforKicks/main.m
hit Run

To change params:
Open JustforKicks/main
Scroll to User Params section
play with params
Remember, run main to run simulation

The script will automatically save the simulation figure as a .png


UPDATE HISTORY ( Not for the lay-user)
%
% Version 8 - Project CleanUp
%
% Version 7 - delaying synaptic depression
% 
% Version 6 - removing continous synaptic current, replacing it with
% discontinous jumps, super similar to synpatic depression
%
%
% Version 5 - updating pulses, no longe rpulse coupled, but upstream
% neurons create synaptic current in downstream neuron
% Notes on v5: it is trash. Could never get it to oscillate correctly :(
% v6 will drop the continous synaptic current and use discontinous
% like synaptic depression is
% I think the problem was in the histerises loop? See FindingBistability
% folder (that's actually the only folder that got changed along with main
% directoty.
%
% Version 4 - allowing for multiple spikes in one time "frame"
%
% Miniupdate C - Increase speed of simulation by removing the large matrix
% multiplication for connectivity matrix. Now picks and chooses only
% neurons it needs!
%
% Version 3-  Update name: Restore.
% Now, changing the theta pulse. Before it was D/n.
% But we must also accoount for the connectivity of the network
% D/((N-1)P).
% Also graphics are getting an update. Need raster plot (DB plot)
%
%
% Version 2 - Tried to speed up code, it failed. Restored back to version 1
%
%
% Version 1 - including synaptic depression