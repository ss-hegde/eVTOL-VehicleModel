import FLOWUnsteady as uns
import FLOWVLM as vlm

function run_single_simulation(run_name, airfoil_name, magVinf, AOA, b, ar, tr, lambda, gamma)

  run_name        = run_name                  # Name of this simulation

  save_path       = run_name                  # Where to save this simulation
  paraview        = false                      # Whether to visualize with Paraview

  # ----------------- SIMULATION PARAMETERS --------------------------------------
  AOA             = AOA                       # (deg) angle of attack
  magVinf         = magVinf                   # (m/s) freestream velocity
  rho             = 1.225                     # (kg/m^3) air density
  qinf            = 0.5*rho*magVinf^2         # (Pa) static pressure
  orientation(t)  = zeros(3)                  # Incidence angles

  Vinf(X, t)      = magVinf*[cosd(AOA), 0.0, sind(AOA)]  # Freestream function


  # ----------------- GEOMETRY PARAMETERS ----------------------------------------
  b               = b
  ar              = ar                       # Aspect ratio b/c_tip
  tr              = tr                       # Taper ratio c_tip/c_root
  twist_root      = 0.0                       # (deg) twist at root
  twist_tip       = 0.0                       # (deg) twist at tip
  lambda          = lambda                      # (deg) sweep
  gamma           = gamma                       # (deg) dihedral

  # Discretization
  n               = 50                        # Number of spanwise elements per side
  r               = 10.0                      # Geometric expansion of elements
  central         = false                     # Whether expansion is central

  # NOTE: A geometric expansion of 10 that is not central means that the last
  #       element is 10 times wider than the first element. If central, the
  #       middle element is 10 times wider than the peripheral elements.

  # ----------------- SOLVER PARAMETERS ------------------------------------------
  # Time parameters
  wakelength      = 2.75*b                    # (m) length of wake to be resolved
  ttot            = wakelength/magVinf        # (s) total simulation time
  nsteps          = 200                       # Number of time steps

  # VLM and VPM parameters
  p_per_step      = 1                         # Number of particle sheds per time step
  lambda_vpm      = 2.0                       # VPM core overlap
  sigma_vpm_overwrite = lambda_vpm * magVinf * (ttot/nsteps)/p_per_step # Smoothing core size
  sigma_vlm_solver= -1                        # VLM-on-VLM smoothing radius (deactivated with <0)
  sigma_vlm_surf  = 0.05*b                    # VLM-on-VPM smoothing radius

  shed_starting   = true                      # Whether to shed starting vortex
  vlm_rlx         = 0.7                       # VLM relaxation




  # ----------------- 1) VEHICLE DEFINITION --------------------------------------
  println("Generating geometry...")

  # Generate wing
  wing = vlm.simpleWing(b, ar, tr, twist_root, lambda, gamma;
                          twist_tip=twist_tip, n=n, r=r, central=central);


  println("Generating vehicle...")

  # Generate vehicle
  system = vlm.WingSystem()                   # System of all FLOWVLM objects
  vlm.addwing(system, "Wing", wing)

  vlm_system = system;                        # System solved through VLM solver
  wake_system = system;                       # System that will shed a VPM wake

  vehicle = uns.VLMVehicle(   system;
                              vlm_system=vlm_system,
                              wake_system=wake_system
                          );


  # ------------- 2) MANEUVER DEFINITION -----------------------------------------

  Vvehicle(t) = zeros(3)                      # Translational velocity of vehicle over time
  anglevehicle(t) = zeros(3)                  # Angle of the vehicle over time

  # Varying angle of the vehicle over time
  # frequency = 5
  # amplitude = 1
  # function anglevehicle(t)

  #   # NOTE: This function receives the non-dimensional time `t` and returns the
  #   #       attitude of the vehicle (a vector with inclination angles with
  #   #       respect to each axis of the global coordinate system)

  #   return [0, amplitude*sin(2*pi*frequency * t*ttot), 0]
  # end


  angle = ()                                  # Angle of each tilting system (none)
  RPM = ()                                    # RPM of each rotor system (none)

  maneuver = uns.KinematicManeuver(angle, RPM, Vvehicle, anglevehicle)
  uns.plot_maneuver(maneuver)

  # ------------- 3) SIMULATION DEFINITION ---------------------------------------

  Vref = 0.0                                  # Reference velocity to scale maneuver by
  RPMref = 0.0                                # Reference RPM to scale maneuver by
  Vinit = Vref*Vvehicle(0)                    # Initial vehicle velocity
  Winit = pi/180*(anglevehicle(1e-6) - anglevehicle(0))/(1e-6*ttot)  # Initial angular velocity

                                              # Maximum number of particles
  max_particles = (nsteps+1)*(vlm.get_m(vehicle.vlm_system)*(p_per_step+1) + p_per_step)

  simulation = uns.Simulation(vehicle, maneuver, Vref, RPMref, ttot;
                                                      Vinit=Vinit, Winit=Winit);




  # ------------- 4) MONITORS DEFINITIONS ----------------------------------------

  calc_aerodynamicforce_fun = uns.generate_calc_aerodynamicforce(;
                                      add_parasiticdrag=true,
                                      add_skinfriction=false,
                                      airfoilpolar="xf-naca6412-il-1000000.csv"
                                      )

  D_dir = [cosd(AOA), 0.0, sind(AOA)]         # Direction of drag
  L_dir = uns.cross(D_dir, [0,1,0])           # Direction of lift

  figs, figaxs = [], []                       # Figures generated by monitor

  # Generate wing monitor
  monitor_wing = uns.generate_monitor_wing(wing, Vinf, b, ar,
                                              rho, qinf, AOA, orientation, anglevehicle, Vvehicle, nsteps;
                                              calc_aerodynamicforce_fun=calc_aerodynamicforce_fun,
                                              L_dir=L_dir,
                                              D_dir=D_dir,
                                              out_figs=figs,
                                              out_figaxs=figaxs,
                                              save_path=save_path,
                                              run_name=run_name,
                                              figname="wing monitor",
                                              );




  # ------------- 5) RUN SIMULATION ----------------------------------------------
  println("Running simulation...")

  uns.run_simulation(simulation, nsteps;
                      # ----- SIMULATION OPTIONS -------------
                      Vinf=Vinf,
                      rho=rho,
                      # ----- SOLVERS OPTIONS ----------------
                      p_per_step=p_per_step,
                      max_particles=max_particles,
                      sigma_vlm_solver=sigma_vlm_solver,
                      sigma_vlm_surf=sigma_vlm_surf,
                      sigma_rotor_surf=sigma_vlm_surf,
                      sigma_vpm_overwrite=sigma_vpm_overwrite,
                      shed_starting=shed_starting,
                      vlm_rlx=vlm_rlx,
                      extra_runtime_function=monitor_wing,
                      # ----- OUTPUT OPTIONS ------------------
                      save_path=save_path,
                      run_name=run_name
                      );


  # ------------------ SAVE GEOMETRY DATA-----------------------------------------
  airfoil_name = airfoil_name
  geom_fname = "wing_geometry_data.csv"
  if save_path!=nothing
    fname = joinpath(save_path, geom_fname)
    f = open(fname, "w")
    print(f, "b_wing,ar_wing,tr,lambda,gamma,airfoil \n")
    print(f, b, ",", ar, ",",tr, ",", lambda, ",", gamma, ",", airfoil_name)
    close(f)
  end

  # ----------------- 6) VISUALIZATION -------------------------------------------
  if paraview
    println("Calling Paraview...")

    # Files to open in Paraview
    files = joinpath(save_path, run_name*"_Wing_vlm...vtk")
    files *= ";"*run_name*"_pfield...xmf;"

    # Call Paraview
    run(`paraview --data=$(files)`)

  end

end

function schedule_simulations()
  # List of simulation parameters for each job
  simulation_params = [
      
      # Training Cases:
      # ("wing_dataset_mw_aoa_-3.8_vinf_3.6",   "xf-naca4412-il-1000000", 3.6, -3.8, 1.077, 4.12, 0.88, -4.5, -3.6),
      # ("wing_dataset_mw_aoa_6.4_vinf_16.5",   "xf-naca4412-il-1000000", 16.5, 6.4, 0.933, 4.488, 0.84, 4.8, -4.9),
      # ("wing_dataset_mw_aoa_13.9_vinf_19.7",  "xf-naca4412-il-1000000", 19.7, 13.9, 1.136, 5.688, 0.55, -3.9, -1.6),
      # ("wing_dataset_mw_aoa_11.3_vinf_5.7",   "xf-naca4412-il-1000000", 5.7, 11.3, 1.045, 4.516, 0.63, 3.4, -3.9),
      # ("wing_dataset_mw_aoa_12.2_vinf_3.7",   "xf-naca4412-il-1000000", 3.7, 12.2, 1.13, 4.538, 0.87, -6.4, -3.8),
      # ("wing_dataset_mw_aoa_12.2_vinf_15.7",  "xf-naca4412-il-1000000", 15.7, 12.2, 1.187, 4.015, 0.86, 8.6, 0.3),
      # ("wing_dataset_mw_aoa_12.1_vinf_15.0",  "xf-naca4412-il-1000000", 15.0, 12.1, 1.11, 5.341, 0.93, 4.3, 0.7),
      # ("wing_dataset_mw_aoa_-0.8_vinf_8.9",   "xf-naca4412-il-1000000", 8.9, -0.8, 1.013, 5.798, 0.99, 9.7, 2.6),
      # ("wing_dataset_mw_aoa_13.6_vinf_8.3",   "xf-naca4412-il-1000000", 8.3, 13.6, 0.815, 4.632, 0.65, -6.7, 4.4),
      # ("wing_dataset_mw_aoa_6.0_vinf_10.1",   "xf-naca4412-il-1000000", 10.1, 6.0, 0.883, 4.123, 0.61, -10.0, 3.9),
      
      # ("wing_dataset_mw_aoa_14.4_vinf_17.8",  "xf-naca4412-il-1000000", 17.8, 14.4, 1.164, 5.295, 0.88, -0.2, -4.8) ,
      # ("wing_dataset_mw_aoa_13.7_vinf_12.0",  "xf-naca4412-il-1000000", 12.0, 13.7, 1.03, 4.057, 0.97, -4.2, -2.4) ,
      # ("wing_dataset_mw_aoa_13.6_vinf_5.3",   "xf-naca4412-il-1000000", 5.3, 13.6, 1.076, 4.284, 0.53, 7.7, -0.5) ,
      # ("wing_dataset_mw_aoa_8.8_vinf_17.4",   "xf-naca4412-il-1000000", 17.4, 8.8, 1.07, 5.721, 0.76, -1.6, 1.3) ,
      # ("wing_dataset_mw_aoa_1.9_vinf_5.7",    "xf-naca4412-il-1000000", 5.7, 1.9, 0.905, 5.625, 0.77, -7.7, 1.7) ,
      # ("wing_dataset_mw_aoa_10.4_vinf_24.9",  "xf-naca4412-il-1000000", 24.9, 10.4, 1.017, 4.325, 0.72, 8.4, 1.3) ,
      # ("wing_dataset_mw_aoa_8.3_vinf_12.7",   "xf-naca4412-il-1000000", 12.7, 8.3, 1.098, 5.423, 0.55, -4.6, 3.9) ,
      # ("wing_dataset_mw_aoa_8.7_vinf_14.0",   "xf-naca4412-il-1000000", 14.0, 8.7, 0.898, 5.962, 0.54, -2.4, -2.8) ,
      # ("wing_dataset_mw_aoa_-3.3_vinf_18.2",  "xf-naca4412-il-1000000", 18.2, -3.3, 1.165, 5.189, 0.67, 6.5, 2.7) ,
      # ("wing_dataset_mw_aoa_-4.5_vinf_20.5",  "xf-naca4412-il-1000000", 20.5, -4.5, 0.917, 5.476, 0.82, 0.2, -3.1) ,
      # ("wing_dataset_mw_aoa_4.1_vinf_24.6",   "xf-naca4412-il-1000000", 24.6, 4.1, 1.022, 4.892, 0.82, 3.0, -0.7) ,
      
      # ("wing_dataset_mw_aoa_0.5_vinf_23.2",  "xf-naca4412-il-1000000", 23.2, 0.5, 0.919, 5.924, 0.61, -6.8, -0.2) ,
      # ("wing_dataset_mw_aoa_0.1_vinf_9.7",   "xf-naca4412-il-1000000", 9.7, 0.1, 0.971, 4.292, 0.85, 8.9, 0.0) ,
      # ("wing_dataset_mw_aoa_1.9_vinf_5.5",   "xf-naca4412-il-1000000", 5.5, 1.9, 0.901, 4.917, 0.76, 6.5, -4.3) ,
      # ("wing_dataset_mw_aoa_5.7_vinf_12.4",  "xf-naca4412-il-1000000", 12.4, 5.7, 1.197, 4.09, 0.87, 9.6, -3.8) ,
      # ("wing_dataset_mw_aoa_10.6_vinf_5.2",  "xf-naca4412-il-1000000", 5.2, 10.6, 1.07, 5.541, 0.79, 9.2, -1.4) ,
      # ("wing_dataset_mw_aoa_3.2_vinf_6.2",   "xf-naca4412-il-1000000", 6.2, 3.2, 0.999, 5.505, 0.8, 5.1, -3.6) ,
      # ("wing_dataset_mw_aoa_8.1_vinf_22.6",  "xf-naca4412-il-1000000", 22.6, 8.1, 1.183, 5.654, 0.72, 3.6, 0.7) ,
      # ("wing_dataset_mw_aoa_14.9_vinf_6.1",  "xf-naca4412-il-1000000", 6.1, 14.9, 0.985, 5.435, 0.61, -6.8, 3.0) ,
      # ("wing_dataset_mw_aoa_13.0_vinf_6.4",  "xf-naca4412-il-1000000", 6.4, 13.0, 1.052, 5.038, 0.8, 3.5, -2.7) ,
      # ("wing_dataset_mw_aoa_2.1_vinf_9.0",   "xf-naca4412-il-1000000", 9.0, 2.1, 1.002, 5.195, 0.6, -9.5, 4.6) ,
      # ("wing_dataset_mw_aoa_14.0_vinf_11.6", "xf-naca4412-il-1000000", 11.6, 14.0, 1.046, 4.291, 0.66, -1.3, 0.6) ,
      
      # ("wing_dataset_mw_aoa_-0.7_vinf_7.5",  "xf-naca4412-il-1000000", 7.5, -0.7, 0.86, 4.745, 0.74, -1.5, -3.3) ,
      # ("wing_dataset_mw_aoa_14.7_vinf_4.0",  "xf-naca4412-il-1000000", 4.0, 14.7, 1.084, 4.526, 0.95, 8.5, -5.0) ,
      # ("wing_dataset_mw_aoa_14.2_vinf_8.5",  "xf-naca4412-il-1000000", 8.5, 14.2, 0.834, 4.277, 0.84, 5.2, 4.1) ,
      # ("wing_dataset_mw_aoa_8.4_vinf_3.3",   "xf-naca4412-il-1000000", 3.3, 8.4, 0.95, 5.979, 0.82, -9.1, -2.0) ,
      # ("wing_dataset_mw_aoa_14.8_vinf_4.6",  "xf-naca4412-il-1000000", 4.6, 14.8, 0.937, 5.579, 0.79, 8.8, 0.5) ,
      # ("wing_dataset_mw_aoa_2.5_vinf_5.4",   "xf-naca4412-il-1000000", 5.4, 2.5, 0.965, 4.333, 0.56, -6.2, -2.6) ,
      # ("wing_dataset_mw_aoa_3.7_vinf_4.4",   "xf-naca4412-il-1000000", 4.4, 3.7, 0.946, 5.624, 0.52, 6.3, 3.6) ,
      # ("wing_dataset_mw_aoa_11.1_vinf_6.7",  "xf-naca4412-il-1000000", 6.7, 11.1, 0.994, 4.566, 0.73, -4.4, -4.0) ,

      # # Testing Cases:
      # ("wing_dataset_mw_aoa_9.7_vinf_20.5",  "xf-naca4412-il-1000000", 20.5, 9.7, 0.875, 4.688, 0.79, -4.3, 0.3) ,
      # ("wing_dataset_mw_aoa_6.5_vinf_9.2",   "xf-naca4412-il-1000000", 9.2, 6.5, 1.043, 4.108, 0.88, -2.7, 1.2) ,
      # ("wing_dataset_mw_aoa_1.2_vinf_17.0",  "xf-naca4412-il-1000000", 17.0, 1.2, 1.063, 5.739, 0.88, -1.8, 1.6) ,
      ("wing_dataset_mw_aoa_14.9_vinf_17.1", "xf-naca4412-il-1000000", 17.1, 14.9, 1.164, 5.673, 0.59, 2.1, -3.5) ,
      ("wing_dataset_mw_aoa_0.4_vinf_16.4",  "xf-naca4412-il-1000000", 16.4, 0.4, 1.164, 4.89, 0.88, 6.8, 3.7) ,
      ("wing_dataset_mw_aoa_0.1_vinf_23.3",  "xf-naca4412-il-1000000", 23.3, 0.1, 0.988, 4.836, 0.9, -2.8, 2.6) ,
      ("wing_dataset_mw_aoa_-4.5_vinf_11.8", "xf-naca4412-il-1000000", 11.8, -4.5, 1.082, 4.596, 0.51, -1.6, 2.7) ,
      ("wing_dataset_mw_aoa_-0.9_vinf_11.4", "xf-naca4412-il-1000000", 11.4, -0.9, 1.073, 5.341, 0.92, -5.8, -3.1) ,
      ("wing_dataset_mw_aoa_4.6_vinf_23.1",  "xf-naca4412-il-1000000", 23.1, 4.6, 1.055, 4.005, 0.91, 9.8, -1.4) ,
      ("wing_dataset_mw_aoa_1.7_vinf_22.9",  "xf-naca4412-il-1000000", 22.9, 1.7, 1.056, 4.811, 0.74, -8.1, 1.6) ,

      

      
  ]


  # Run each simulation sequentially
  for (run_name, airfoil_name, magVinf, AOA, b, ar, tr, lambda, gamma) in simulation_params
    run_single_simulation(run_name, airfoil_name, magVinf, AOA, b, ar, tr, lambda, gamma)
  end

end

schedule_simulations()