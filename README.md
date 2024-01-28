# V2X-System-Level-Simulator
This repository provides the system-level implementation of the V2X communication.
# Runnig the Simulation

The main file of the LTE System Level Simulator (LTE SL Simulator) is LTE_sim_main.m, though you may normally
run the simulation through a batch file such as LTE_sim_launcher.m, which performs the following tasks:
- Loading a configuration file of choice. See Section III for a list of configurable parameters.
- Executing the LTE_sim_main.m main simulation file.
The simulation parameters are loaded from the LTE_load_params_*.m script or alternatively from similarly-named
files. Do note that the LTE_load_params_dependant.m script is used for automatically generating additional simulation
parameters from the base parameters specified. A basic example of a configuration file used with the LTE SL Simulator can
be found in LTE_load_params_hex_grid_tilted.m.

#  SIMULATION PARAMETERS
Below you can find a list of the parameters that can be configured in the configuration parameters files. You can find all of
them in the +simulation_config simulator subfolder. The files are loaded via the LTE_load_params function. Check
LTE_sim_main_launcher_examples and LTE_load_params function for a list of possible preconfigured simulation
files.
When called with no arguments, LTE_load_params loads the +simulation_config/hex_grid_tilted file. The
optional string argument allows you (among others) to load the following other preconfigured setups (see
LTE_sim_main_launcher_examples):
•  tri_sector
•  tri_sector_tilted, tri_sector_tilted_4x2 tri_sector_tilted_4x4.
•  tri_sector_plus_femtocells
•  six_sector_tilted
•  capesso_pathlossmaps
•  omnidirectional_eNodeBs

## A. Debug options
• LTE_config.debug_level: configures how much debug text output is shown. Options are:
– 0: no output.
– 1: basic output.
– 2: extended output.
## B. Plotting options
• LTE_config.show_network.: configures how much plots are shown. Options are:
– 0: no plots shown.
– 1: show some plots.
– 2: show all plots, which includes one showing the moving User Equipments (UEs), which may slow down simulations
significantly.
– 3: show also the plots of the generated microscale fading traces.
## C. General parameters
• LTE_config.frequency.: frequency in which the system is operating [Hz].
• LTE_config.bandwidth.: system bandwidth. Allowed values are 1.4 MHz, 3 MHz, 5 MHz, 10 MHz, 15 MHz, and
20 MHz. This bandwidths are equivalent to 6, 15, 25, 50, 75, and 100 Resource Blocks (RBs) respectively.
• LTE_config.nTX: number of transmit antennas. Used to generate the channel trace.
• LTE_config.nRX: number of receive antennas. Used to generate the channel trace.
• LTE_config.tx_mode: the transmission modes are defined in TS 36.213-820 Section 7.1, page 12 [3].
– 1: single antenna.
– 2: Transmission Diversity (TxD).
– 3: Open Loop Spatial Multiplexing (OLSM). Spatial multiplexing with Large Cyclic Delay Diversity (CDD).
– 4: Closed Loop Spatial Multiplexing (CLSM).
– 5: Multiuser MIMO (not yet implemented).

## D. Random number generation options

These options allow you to actually reproduce the same exact simulation by means of resetting the random number generator
to a known seed.
• LTE_config.seedRandStream: in order to allow repeatability, it is possible to seed MATLAB’s default random
number generator. Set it to either true or false.
• LTE_config.RandStreamSeed: if the above is set to true, it specifies the seed. Seeds must be an integer between
0 and 2^32.

## E. Simulation time
 LTE_config.simulation_time_tti: length of the simulation in Transmission Time Intervals (TTIs)

## F. Cache options
• LTE_config.cache_network: whether you want to save the generated eNodeBs, Pathloss map and Shadow fading
map to a .mat file. Either true or false. All cache options work in the following way:
– cache=true and file exists: read cache file.
– cache=true and file does not exist: create and then store data in cache file.
– cache=false: do not use cache at all.
• LTE_config.network_cache: the name of the cache file. set it to auto if you want the simulator to assign a name
automatically.
• LTE_config.delete_ff_trace_at_end: since the microscala fading trace takes up large amounts of space, when
doing the final save command, it is preferable to delete it, so as not to have too large result files.
• LTE_config.delete_pathloss_at_end: Further reduces the amount of space needed to store the traces by
deleting the pathloss maps from the results file.
• LTE_config.UE_cache: whether to save the user position to a file. Either true or false.
• LTE_config.UE_cache_file: the name of the cache file. set it to auto if you want the simulator to assign a name
automatically.
## G. Network layout and macroscopic pathloss parameters
These parameters specify how the network layout is created. However, if the map is loaded, these parameters will be
overridden by the loaded map.
• LTE_config.network_source: Available options
– generated: A hexagonal grid of equidistantly-spaced eNodeB sites with three sectors each will be created.
– capesso: eNodeB position, configuration, and pathloss data are read from data exported from and written from
the CapessoTMplanning tool (see Section X). When using this source, shadow fading data is not generated, as the
imported pathloss maps should already have it incorporated.
### 1) Generated network parameters:
• LTE_config.inter_eNodeB_distance: in meters. When the network is generated, this determines the distance
between the eNodeBs.
• LTE_config.map_resolution: in meters/pixel. Also the resolution used for initial user creation.
• LTE_config.nr_eNodeB_rings: number of eNodeB rings. 0 rings specifies that just a single eNodeB will be created.
• LTE_config.minimum_coupling_loss (optional): describes the minimum loss in signal [dB] between Base Station
(BS) and UE or UE and UE in the worst case and is defined as the minimum distance loss including antenna gains measured
between antenna connectors. Recommended values [5] are 70 dB for urban areas, 80 dB for rural.
• LTE_config.macroscopic_pathloss_model: sets what macroscopic pathloss model is to be used. Depending
on the choice, different choices are available for
LTE_config.macroscopic_pathloss_model_settings.environment. The available macroscopic pathloss
models are:
– free space: free space pathloss. More for testing purposes than for actual use with simulations. 

 
