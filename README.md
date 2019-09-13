# iceberg_Fe_model
Iceberg Fe model

Readme for modified Moon et al. 2017 iceberg melt model
Original code by David Sutherland, modified by Dustin Carroll

README, iceberg melt code, Dave Sutherland, September 2018

These files were created by Dave Sutherland and used in Moon et al. 2017 paper on iceberg in Sermilik Fjord. The main wrapper code here is run_meltcases.m, which uses constant forcing (wind, sun, T/S profiles, velocity, etc.) and replicates the summer and winter cases described in the paper. To run the along-fjord variation or the monthly variation, you would need a different wrapper and some extra files that are available upon request. 

To get the summer time results we used in the paper, execute run_meltcases as written in the dropbox folder. It will spit out text going through 20 size classes of icebergs, from length L = 50 m to L = 1000 m in 50 m increments. Please read the Moon et al. methods section to get a feel for the flow of what happens when (i.e. initializing icebergs, calculating melt, adjusting size of iceberg, etc.).

The main inputs you can change to run_meltcases are:
- T/S data: currently this points to Aug 2010 CTD data and is in the structure serm loaded. Note it needs to be able to find temp and sal variables in whatever structure you give it. 
- casename: name your run here
- Velocity: this refers to both horizontal velocity as a function of depth and a constant vertical velocity (wvel). We initialize the horizontal velocity using the file indicated. 
- Tair: air temperature
- SWflux: shortwave radiation flux (W/m2)
- Winds: m/s wind speed
- IceC: sea ice concentration. This affects both the wave driven melt and if it equals 1, then icebergs are considered static and don’t move, which affects the relative velocity (forced convection melt)

You can also change the flags in do_melt to turn on/off certain melt terms. Additionally, you can change the size classes you do (though note that changing this means you’ll need to change the distribution vectors as well!), timespan run for, and dz (note that only dz = 10 m and dz = 5 m have been tested). 

And then it runs! It uses these inputs in the main function called iceberg_melt.m, which runs for each size class of iceberg. Then, it scales up the output in mberg by using distributions of icebergs from Sulak et al. 2017 (per and Nbergs). A lot of why all these files are here are so we could run the code over multiple size classes and calculate melt fluxes, not just get melt rates...

Output
The main output from iceberg_melt is a structure called mberg, which is 1 x # of size classes (in the original case, 20). 

Inside mberg, for a particular size class (i.e., try mberg(1)), you’ll find lots of variables. They are:
ice_init: initial geometry of this size class iceberg
M[****]: depending on what melt terms you called, you’ll get these variables that start with capital M for each melt process and total. These are time series (and depth if applicable) of the total melt flux (i.e., melt rate times a relevant area) for this size class in m3/sec. 
i_m[****]: In contrast to the M terms described above, these are the actual melt rates (in m/day)  in time for each process and total, again depending on what melt terms you called. 

The other output variables are how the geometry changed in time. 

Note the depth dependent terms are # vertical layers x 1 x # timesteps (in this case, 50 x 1 x 30), while the depth independent terms are just 1 x timesteps (1 x 30). 

Dependencies
The iceberg_melt code calls on a host of mfiles that all should be in the Icebergs_code folder. We also have a run_meltcases_seasonal.m that uses mooring data to get monthly timescales and run_meltcases_alongfjord.m that uses spatially varying conditions and ice distributions. You’d need extra files to run these. I have now put those files on dropbox. Basically, the along-fjord case uses 5 zones (~20km long each for Sermilik) to bin icebergs by, as well as T/S profiles and ice concentration. We have to use the same velocity info since we don’t have any idea how it varies along-fjord.

For the seasonal case, it utilizes Fiamma’s mooring info for T/S by month, then builds up a seasonal ADCP record based on Becca’s work. It also varies the atmospheric and sea ice concentrations by month. However, the iceberg distributions are assumed constant throughout the year. 


Added notes:
- note the rolling and slab breakoff processes do not function at the moment and are hard-coded to not be implemented. 
- we run the code for a month because of the time dependent horizontal velocity (i.e. in Sermilik Fjord there is a strong intermediary exchange flow). 
- the melt rates are most sensitive to the velocities. I was able to justify the horizontal velocities based on actual data, and the vertical velocities based on means from running MITgcm simulations of similar sized icebergs in summer and winter T/S conditions 
