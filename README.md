
# Dynamic Counter Propagating All Normald Dispersion Laser - Simulation
Emergence of counter-propagating single-cavity dual-comb lasers (SCDCLs) have revolutionized the field of dual-comb metrology due to its simpler architecture and free-running configuration. Conventional time-domain dual-comb technique suffers from the trade-off problem between temporal resolution and measurement frame rate requiring large acquisition time for high resolution measurements with significant dead time. To overcome this trade-off problem in SCDCL, here we introduce, for the first time, an all optical dynamic modulation technique for SCDCLs. We demonstrate the novel idea by modulating the repetition rate difference (Δfr) about zero using the newly developed CANDi fiber laser. Owing to the non-reciprocal longitudinal gain distribution in the gain fiber, we achieve dynamic tuning of Δfr by tuning the pump laser. Using dispersive Fourier transform and spectral interferometry, we further verify that the relative timing jitter is unaffected due to this dynamic tuning technique. We apply this technique to characterize the relaxation time constant of an off-shelf SESAM and achieve a framerate enhancement factor of 143, which can be further improved by modulating at a faster rate. This dynamic SCDCL opens many new exciting possibilities in nonlinear dual-comb metrology. In this repository supporting code is provided to better understand the physical process of the NPR mode-locking and dynamic modulation induced by the spatially varying YDF gain.


For further information you can reach Jonathan Musgrave @jomu3154@colorado.edu!
## Authors

- Jonathan Musgrave - Developed simulation
- Neeraj Prakash - Experimental and lead author
- Bowen Li - Experimental and corresponding author
- Shu-Wei Huang - Corresponding author

## Documentation

The simulation works by solving the IVP problem of the CNLSE coupled to the Rate_equations. This is completed by defining the parameters of the fiber in the fiber struct as well as pump and input signal. The sys struct defines the system variables important to the numerical sampling such as window size, number of temporal steps, and the frep. Note that the frep defines the pulse energy through the defined signal average power for both the CW and CCW and is a reference repetition rate at which the simulation time window moves at. Thus this is not the same as the variable reptition rate of the CW and CCW pulse directions. 

Each segement of fiber is treated as independent struct variables. For instance in CANDI_Model_Symmetric_Simple.m the system is seperated into the passive1, YDF, and passive2, with corresponding fields fiber.type defined passive, active, and passive respectively signifying if the fiber is doped or not. The simulation consists of solving the CNLSE and RE in a SSFM via the following shooting method:

pseudo code:
1. Define fibers, system, signals, and pump structs and relevant variables.
2. Pass variables to Fiber_Prop_Lz_v1.m and check for defaults
3. If fiber is passive propagate CNLSE without solving RE. Else if fiber is active solve fiber in forward direction initializing backward pump power as no absortion. Solve the IVP problem via the shooting method.
4. Propagate signal to the next segment of fiber via Fiber_Prop_Lz_v1.m 
5. Once all fibers have been traversed apply PBS, spectral filter,  and waveplates and save variables.
6. If modulate pump on this round trip than do. 
7. repeat 1-6 until completion.


