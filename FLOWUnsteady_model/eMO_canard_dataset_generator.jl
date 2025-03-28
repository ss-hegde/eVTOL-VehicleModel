import FLOWUnsteady as uns
import FLOWVLM as vlm

function run_single_simulation(run_name, airfoil_name, magVinf, AOA, b, ar, tr, lambda, gamma, frequency)

  run_name        = run_name                  # Name of this simulation

  save_path       = run_name                  # Where to save this simulation
  paraview        = false                      # Whether to visualize with Paraview

  # ----------------- SIMULATION PARAMETERS --------------------------------------
  AOA             = AOA                       # (deg) angle of attack
  magVinf         = magVinf                   # (m/s) freestream velocity
  rho             = 1.225                     # (kg/m^3) air density
  qinf            = 0.5*rho*magVinf^2         # (Pa) static pressure
  orientation(t)  = [0.0, 1.0, 0.0]           # Incidence angles [Phi, theta, psi]

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

  wing_pos = [0.0, 0.0, 0.0]
  wing_axis = uns.gt.rotation_matrix(orientation(0)[1], -orientation(0)[2], orientation(0)[3])
  vlm.setcoordsystem(wing, wing_pos, wing_axis)


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

  # Vvehicle(t) = zeros(3)                      # Translational velocity of vehicle over time
  
  frequency = frequency
  amplitude = 0.02
  function Vvehicle(t)

    # NOTE: This function receives the non-dimensional time `t` and returns the
    #       attitude of the vehicle (a vector with inclination angles with
    #       respect to each axis of the global coordinate system)

    return [0, 0, amplitude*sin(2*pi*frequency * t)]
  end
  
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

  Vref = 1.0                                  # Reference velocity to scale maneuver by
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
      # ("wing_dataset_ca_aoa_1.0_vinf_5", "xf-naca6412-il-1000000", 5.0, 1.0, 7.684, 15.368, 1.0, 0.0, 0.0),
      # Training Cases:
      # ("wing_dataset_ca_aoa_13.7_vinf_4.3",  "xf-naca6412-il-1000000", 4.3,  13.7, 0.896, 13.561, 0.91, 7.2, 3.8),
      # ("wing_dataset_ca_aoa_9.4_vinf_5.4",   "xf-naca6412-il-1000000", 5.4,  9.4,  0.814, 12.782, 0.71, 5.6, 3.5),
      # ("wing_dataset_ca_aoa_14.5_vinf_10.6", "xf-naca6412-il-1000000", 10.6, 14.5, 0.814, 16.82, 0.78, 1.2, 1.0),
      # ("wing_dataset_ca_aoa_-4.4_vinf_6.7",  "xf-naca6412-il-1000000", 6.7,  -4.4, 1.162, 14.851, 0.9, 6.0, 0.9),
      # ("wing_dataset_ca_aoa_0.9_vinf_16.3",  "xf-naca6412-il-1000000", 16.3, 0.9,  1.088, 12.641, 0.52, 1.1, 3.9),
      # ("wing_dataset_ca_aoa_0.7_vinf_6.4",   "xf-naca6412-il-1000000", 6.4,  0.7,  1.041, 17.683, 0.83, 6.7, -2.8),
      # ("wing_dataset_ca_aoa_-1.1_vinf_16.4", "xf-naca6412-il-1000000", 16.4, -1.1, 1.135, 17.983, 0.69, 9.9, -2.3),
      # ("wing_dataset_ca_aoa_5.9_vinf_11.3",  "xf-naca6412-il-1000000", 11.3, 5.9,  1.057, 15.061, 0.56, 10.0, 2.1),
      # ("wing_dataset_ca_aoa_-1.4_vinf_19.2", "xf-naca6412-il-1000000", 19.2, -1.4, 1.071, 13.587, 0.81, 4.5, 2.3),
      # ("wing_dataset_ca_aoa_3.8_vinf_15.1",  "xf-naca6412-il-1000000", 15.1, 3.8,  1.08, 14.964, 0.55, 4.6, 3.9),

      # ("wing_dataset_ca_aoa_1.0_vinf_3.7",   "xf-naca6412-il-1000000", 3.7,  1.0,  0.919, 16.66, 0.86, 7.0, 3.9),
      # ("wing_dataset_ca_aoa_9.0_vinf_3.5",   "xf-naca6412-il-1000000", 3.5,  9.0,  0.806, 14.328, 0.67, 3.2, 1.2),
      # ("wing_dataset_ca_aoa_4.4_vinf_3.3",   "xf-naca6412-il-1000000", 3.3,  4.4,  1.119, 12.49, 0.77, 3.9, -0.7),
      # ("wing_dataset_ca_aoa_9.4_vinf_1.2",   "xf-naca6412-il-1000000", 1.2,  9.4,  0.922, 12.582, 0.74, 4.0, -3.8),
      # ("wing_dataset_ca_aoa_2.3_vinf_2.4",   "xf-naca6412-il-1000000", 2.4,  2.3,  1.056, 16.631, 0.51, 2.8, -1.7),
      # ("wing_dataset_ca_aoa_5.9_vinf_7.4",   "xf-naca6412-il-1000000", 7.4,  5.9,  1.004, 16.813, 0.59, 6.9, -4.1),
      # ("wing_dataset_ca_aoa_9.5_vinf_15.0",  "xf-naca6412-il-1000000", 15.0, 9.5,  0.89,  12.413, 0.89, 6.5, -1.7),
      # ("wing_dataset_ca_aoa_-4.6_vinf_14.4", "xf-naca6412-il-1000000", 14.4, -4.6, 0.937, 14.782, 0.78, 5.0, 4.4),
      # ("wing_dataset_ca_aoa_14.5_vinf_9.0",  "xf-naca6412-il-1000000", 9.0,  14.5, 0.913, 16.737, 0.68, 1.1, -3.2),
      # ("wing_dataset_ca_aoa_12.0_vinf_7.2",  "xf-naca6412-il-1000000", 7.2,  12.0, 1.159, 14.659, 0.7, 0.3, 2.4),

      # ("wing_dataset_ca_aoa_-1.8_vinf_11.1", "xf-naca6412-il-1000000", 11.1, -1.8, 1.117, 13.431, 0.65, 3.7, 3.3),
      # ("wing_dataset_ca_aoa_5.2_vinf_13.1",  "xf-naca6412-il-1000000", 13.1, 5.2,  0.937, 15.009, 0.89, 1.7, 0.3),
      # ("wing_dataset_ca_aoa_3.7_vinf_12.5",  "xf-naca6412-il-1000000", 12.5, 3.7,  1.03,  17.503, 0.72, 7.4, -4.2),
      # ("wing_dataset_ca_aoa_7.3_vinf_16.9",  "xf-naca6412-il-1000000", 16.9, 7.3,  0.989, 12.238, 0.87, 6.1, 3.1),
      # ("wing_dataset_ca_aoa_7.0_vinf_8.0",   "xf-naca6412-il-1000000", 8.0,  7.0,  0.977, 15.825, 0.81, 0.3, -0.3),
      # ("wing_dataset_ca_aoa_8.1_vinf_18.9",  "xf-naca6412-il-1000000", 18.9, 8.1,  1.101, 14.774, 0.59, 9.7, 4.4),
      # ("wing_dataset_ca_aoa_5.9_vinf_11.7",  "xf-naca6412-il-1000000", 11.7, 5.9,  0.897, 16.984, 0.55, 7.6, 1.3),
      # ("wing_dataset_ca_aoa_10.8_vinf_15.9", "xf-naca6412-il-1000000", 15.9, 10.8, 1.011, 12.058, 0.77, 2.3, -0.5),
      # ("wing_dataset_ca_aoa_10.8_vinf_7.5",  "xf-naca6412-il-1000000", 7.5,  10.8, 0.971, 12.593, 1.0, 5.4, 4.1),
      # ("wing_dataset_ca_aoa_-3.4_vinf_6.4",  "xf-naca6412-il-1000000", 6.4, -3.4,  0.839, 13.309, 0.57, 2.7, 1.0),

      # ("wing_dataset_ca_aoa_2.8_vinf_1.7",   "xf-naca6412-il-1000000", 1.7,  2.8,  0.991, 16.382, 0.97, 0.0, 0.1),
      # ("wing_dataset_ca_aoa_-2.0_vinf_8.5",  "xf-naca6412-il-1000000", 8.5,  -2.0, 1.084, 15.091, 0.81, 5.4, -2.5),
      # ("wing_dataset_ca_aoa_1.2_vinf_9.9",   "xf-naca6412-il-1000000", 9.9,  1.2,  0.926, 16.464, 0.84, 1.5, -3.7),
      # ("wing_dataset_ca_aoa_14.0_vinf_14.6", "xf-naca6412-il-1000000", 14.6, 14.0, 0.82,  13.322, 0.66, 6.0, 1.1),
      # ("wing_dataset_ca_aoa_-1.9_vinf_15.7", "xf-naca6412-il-1000000", 15.7, -1.9, 1.081, 13.901, 0.69, 7.5, 4.2),
      # ("wing_dataset_ca_aoa_6.1_vinf_1.2",   "xf-naca6412-il-1000000", 1.2,  6.1,  1.139, 17.306, 0.88, 8.1, -2.5),
      # ("wing_dataset_ca_aoa_1.5_vinf_3.9",   "xf-naca6412-il-1000000", 3.9,  1.5,  0.948, 16.463, 0.84, 9.3, -4.7),
      # ("wing_dataset_ca_aoa_4.4_vinf_8.6",   "xf-naca6412-il-1000000", 8.6,  4.4,  1.028, 12.779, 0.69, 7.0, 1.3),
      # ("wing_dataset_ca_aoa_13.9_vinf_13.8", "xf-naca6412-il-1000000", 13.8, 13.9, 0.984, 17.88, 0.68, 0.2, 1.8),
      # ("wing_dataset_ca_aoa_9.6_vinf_9.7",   "xf-naca6412-il-1000000", 9.7,  9.6,  0.882, 15.17, 0.69, 9.4, -4.8),
      
      # Unsteady dataset
      # ("wing_dataset_ca_aoa_4.9_vinf_3.5_freq_66.1",  "xf-naca6412-il-1000000", 3.5, 4.9, 1.131, 15.934, 0.93, 6.4, 1.0, 66.1) ,
      # ("wing_dataset_ca_aoa_0.9_vinf_1.3_freq_57.0",  "xf-naca6412-il-1000000", 1.3, 0.9, 0.904, 14.884, 0.68, 4.9, 1.0, 57.0) ,
      # ("wing_dataset_ca_aoa_0.7_vinf_6.9_freq_55.2",  "xf-naca6412-il-1000000", 6.9, 0.7, 0.858, 16.939, 0.54, 1.3, 2.5, 55.2) ,
      # ("wing_dataset_ca_aoa_3.7_vinf_8.0_freq_68.4",  "xf-naca6412-il-1000000", 8.0, 3.7, 0.931, 12.84, 0.53, 4.5, 1.2, 68.4) ,
      # ("wing_dataset_ca_aoa_4.4_vinf_7.1_freq_66.2",  "xf-naca6412-il-1000000", 7.1, 4.4, 1.184, 14.236, 0.59, 2.1, 4.3, 66.2) ,
      # ("wing_dataset_ca_aoa_2.4_vinf_11.4_freq_69.3", "xf-naca6412-il-1000000", 11.4, 2.4, 0.862, 15.846, 0.95, 1.9, 4.1, 69.3) ,
      # ("wing_dataset_ca_aoa_5.7_vinf_2.1_freq_63.7",  "xf-naca6412-il-1000000", 2.1, 5.7, 1.036, 15.909, 0.67, 5.7, 1.4, 63.7) ,
      # ("wing_dataset_ca_aoa_4.6_vinf_10.9_freq_68.6", "xf-naca6412-il-1000000", 10.9, 4.6, 1.061, 14.308, 0.83, 2.4, 0.4, 68.6) ,
      # ("wing_dataset_ca_aoa_4.0_vinf_11.3_freq_69.7", "xf-naca6412-il-1000000", 11.3, 4.0, 0.856, 16.585, 0.9, 8.5, 4.4, 69.7) ,
      # ("wing_dataset_ca_aoa_2.8_vinf_9.9_freq_51.3",  "xf-naca6412-il-1000000", 9.9, 2.8, 0.935, 16.847, 0.58, 3.6, 4.6, 51.3) ,
      
      # ("wing_dataset_ca_aoa_2.3_vinf_6.7_freq_68.3",  "xf-naca6412-il-1000000", 6.7, 2.3, 0.805, 13.743, 0.63, 2.0, 4.2, 68.3) ,
      # ("wing_dataset_ca_aoa_1.6_vinf_10.2_freq_58.1", "xf-naca6412-il-1000000", 10.2, 1.6, 1.1, 13.816, 0.5, 4.2, 3.5, 58.1) ,
      # ("wing_dataset_ca_aoa_3.6_vinf_3.3_freq_59.9",  "xf-naca6412-il-1000000", 3.3, 3.6, 1.018, 16.542, 0.55, 7.1, 3.8, 59.9) ,
      # ("wing_dataset_ca_aoa_0.3_vinf_2.2_freq_58.8",  "xf-naca6412-il-1000000", 2.2, 0.3, 1.169, 14.404, 0.51, 8.1, 4.3, 58.8) ,
      # ("wing_dataset_ca_aoa_2.4_vinf_6.9_freq_55.5",  "xf-naca6412-il-1000000", 6.9, 2.4, 0.957, 12.699, 0.64, 4.0, 1.4, 55.5) ,
      # ("wing_dataset_ca_aoa_3.8_vinf_7.1_freq_65.5",  "xf-naca6412-il-1000000", 7.1, 3.8, 1.009, 13.128, 0.66, 4.0, 1.7, 65.5) ,
      # ("wing_dataset_ca_aoa_5.2_vinf_10.9_freq_60.2", "xf-naca6412-il-1000000", 10.9, 5.2, 1.15, 16.907, 0.94, 9.3, 1.4, 60.2) ,
      # ("wing_dataset_ca_aoa_4.3_vinf_9.2_freq_57.0",  "xf-naca6412-il-1000000", 9.2, 4.3, 1.128, 14.515, 0.63, 4.9, 2.4, 57.0) ,
      # ("wing_dataset_ca_aoa_3.5_vinf_1.5_freq_55.2",  "xf-naca6412-il-1000000", 1.5, 3.5, 1.012, 12.693, 0.89, 8.1, 2.4, 55.2) ,
      
      # ("wing_dataset_ca_aoa_1.4_vinf_6.7_freq_68.2", "xf-naca6412-il-1000000", 6.7, 1.4, 0.824, 13.992, 0.69, 7.3, 4.2, 68.2) ,
      # ("wing_dataset_ca_aoa_1.9_vinf_9.5_freq_66.2", "xf-naca6412-il-1000000", 9.5, 1.9, 1.198, 14.083, 0.71, 8.4, 3.1, 66.2) ,
      # ("wing_dataset_ca_aoa_2.9_vinf_1.6_freq_65.6", "xf-naca6412-il-1000000", 1.6, 2.9, 1.13, 16.07, 0.57, 6.7, 0.1, 65.6) ,
      # ("wing_dataset_ca_aoa_4.3_vinf_5.9_freq_63.4", "xf-naca6412-il-1000000", 5.9, 4.3, 1.084, 17.894, 0.9, 4.7, 2.9, 63.4) ,
      # ("wing_dataset_ca_aoa_2.9_vinf_5.6_freq_69.2", "xf-naca6412-il-1000000", 5.6, 2.9, 1.143, 12.441, 0.64, 0.4, 4.8, 69.2) ,
      # ("wing_dataset_ca_aoa_5.2_vinf_4.7_freq_60.3", "xf-naca6412-il-1000000", 4.7, 5.2, 1.046, 14.595, 0.69, 0.3, 3.4, 60.3) ,
      # ("wing_dataset_ca_aoa_1.7_vinf_5.7_freq_66.0", "xf-naca6412-il-1000000", 5.7, 1.7, 0.837, 17.762, 0.81, 4.6, 2.2, 66.0) ,
      # ("wing_dataset_ca_aoa_3.6_vinf_9.1_freq_67.7", "xf-naca6412-il-1000000", 9.1, 3.6, 1.01, 12.43, 0.61, 7.4, 2.8, 67.7) ,
      # ("wing_dataset_ca_aoa_1.8_vinf_5.2_freq_70.0", "xf-naca6412-il-1000000", 5.2, 1.8, 0.83, 17.615, 0.7, 2.2, 2.3, 70.0) ,
      # ("wing_dataset_ca_aoa_5.8_vinf_4.4_freq_57.1", "xf-naca6412-il-1000000", 4.4, 5.8, 1.148, 17.299, 0.68, 0.5, 4.3, 57.1) ,
      
      # ("wing_dataset_ca_aoa_4.3_vinf_4.5_freq_55.4",  "xf-naca6412-il-1000000", 4.5, 4.3, 1.046, 12.798, 0.73, 7.5, 1.4, 55.4) ,
      # ("wing_dataset_ca_aoa_5.5_vinf_6.4_freq_58.9",  "xf-naca6412-il-1000000", 6.4, 5.5, 0.918, 12.519, 0.51, 5.8, 2.3, 58.9) ,
      # ("wing_dataset_ca_aoa_2.4_vinf_8.7_freq_60.4",  "xf-naca6412-il-1000000", 8.7, 2.4, 1.077, 15.615, 0.81, 8.5, 4.8, 60.4) ,
      # ("wing_dataset_ca_aoa_0.8_vinf_10.8_freq_59.4", "xf-naca6412-il-1000000", 10.8, 0.8, 0.862, 16.2, 0.58, 2.4, 0.5, 59.4) ,
      # ("wing_dataset_ca_aoa_2.4_vinf_5.0_freq_66.2",  "xf-naca6412-il-1000000", 5.0, 2.4, 0.817, 16.437, 0.57, 1.2, 0.9, 66.2) ,
      # ("wing_dataset_ca_aoa_5.7_vinf_8.7_freq_57.5",  "xf-naca6412-il-1000000", 8.7, 5.7, 0.92, 15.918, 0.89, 0.5, 0.0, 57.5) ,
      # ("wing_dataset_ca_aoa_2.7_vinf_10.8_freq_60.7", "xf-naca6412-il-1000000", 10.8, 2.7, 1.126, 16.72, 0.88, 7.6, 4.4, 60.7) ,
      # ("wing_dataset_ca_aoa_2.7_vinf_11.0_freq_66.6", "xf-naca6412-il-1000000", 11.0, 2.7, 1.191, 12.926, 0.78, 8.0, 4.3, 66.6) ,
      # ("wing_dataset_ca_aoa_3.3_vinf_3.3_freq_64.9",  "xf-naca6412-il-1000000", 3.3, 3.3, 1.01, 12.787, 0.56, 6.1, 2.1, 64.9) ,
      # ("wing_dataset_ca_aoa_4.0_vinf_3.4_freq_68.0",  "xf-naca6412-il-1000000", 3.4, 4.0, 1.093, 14.767, 0.8, 4.8, 2.7, 68.0) ,
      # ("wing_dataset_ca_aoa_4.4_vinf_6.9_freq_68.6",  "xf-naca6412-il-1000000", 6.9, 4.4, 1.095, 12.804, 0.58, 6.0, 4.3, 68.6) ,

      # # # Testing Cases
      # ("wing_dataset_ca_aoa_1.9_vinf_2.3",   "xf-naca6412-il-1000000", 2.3,  1.9,  0.904, 14.059, 0.68, 1.2, 0.4),
      # ("wing_dataset_ca_aoa_8.4_vinf_5.6",   "xf-naca6412-il-1000000", 5.6,  8.4,  1.174, 12.063, 0.85, 8.7, -4.0),
      # ("wing_dataset_ca_aoa_2.0_vinf_14.3",  "xf-naca6412-il-1000000", 14.3, 2.0,  1.05,  13.981, 0.56, 8.1, 3.3),
      # ("wing_dataset_ca_aoa_-2.7_vinf_1.3",  "xf-naca6412-il-1000000", 1.3,  -2.7, 1.107, 13.129, 0.92, 3.3, 3.8),
      # ("wing_dataset_ca_aoa_2.7_vinf_6.5",   "xf-naca6412-il-1000000", 6.5,  2.7,  1.082, 15.554, 0.99, 4.9, 4.7),
      # ("wing_dataset_ca_aoa_6.9_vinf_15.0",  "xf-naca6412-il-1000000", 15.0, 6.9,  1.141, 17.56, 0.65, 3.2, 3.1),
      # ("wing_dataset_ca_aoa_5.3_vinf_1.4",   "xf-naca6412-il-1000000", 1.4,  5.3,  1.046, 17.169, 0.52, 9.6, -2.7),
      # ("wing_dataset_ca_aoa_7.3_vinf_15.9",  "xf-naca6412-il-1000000", 15.9, 7.3,  0.989, 14.504, 0.9, 0.7, 4.9),
      # ("wing_dataset_ca_aoa_-1.2_vinf_13.9", "xf-naca6412-il-1000000", 13.9, -1.2, 0.835, 16.734, 0.78, 2.5, -4.4),
      # ("wing_dataset_ca_aoa_10.7_vinf_18.8", "xf-naca6412-il-1000000", 18.8, 10.7, 1.092, 12.81, 0.64, 7.1, 0.9),

      # Unsteady dataset
      # ("wing_dataset_ca_aoa_4.6_vinf_9.2_freq_52.4",  "xf-naca6412-il-1000000", 9.2, 4.6, 1.002, 13.74, 0.71, 5.1, 3.7, 52.4) ,
      # ("wing_dataset_ca_aoa_3.6_vinf_10.9_freq_58.4", "xf-naca6412-il-1000000", 10.9, 3.6, 1.038, 16.055, 0.74, 7.1, 4.3, 58.4) ,
      # ("wing_dataset_ca_aoa_4.5_vinf_2.0_freq_54.6",  "xf-naca6412-il-1000000", 2.0, 4.5, 0.808, 16.906, 0.63, 1.8, 3.0, 54.6) ,
      # ("wing_dataset_ca_aoa_1.0_vinf_4.2_freq_51.9",  "xf-naca6412-il-1000000", 4.2, 1.0, 1.198, 17.705, 0.7, 10.0, 2.1, 51.9) ,
      ("wing_dataset_ca_aoa_5.8_vinf_5.7_freq_52.8",  "xf-naca6412-il-1000000", 5.7, 5.8, 1.137, 16.56, 0.69, 6.6, 2.8, 52.8) ,
      ("wing_dataset_ca_aoa_0.5_vinf_7.7_freq_63.5",  "xf-naca6412-il-1000000", 7.7, 0.5, 1.004, 17.032, 0.57, 1.9, 0.7, 63.5) ,
      ("wing_dataset_ca_aoa_3.8_vinf_7.5_freq_64.5",  "xf-naca6412-il-1000000", 7.5, 3.8, 1.027, 14.431, 0.59, 8.7, 2.6, 64.5) ,
      ("wing_dataset_ca_aoa_3.5_vinf_9.8_freq_55.6",  "xf-naca6412-il-1000000", 9.8, 3.5, 1.115, 17.897, 0.78, 1.2, 0.8, 55.6) ,
      ("wing_dataset_ca_aoa_2.9_vinf_10.8_freq_68.8", "xf-naca6412-il-1000000", 10.8, 2.9, 1.102, 13.777, 0.73, 8.3, 0.7, 68.8) ,
      ("wing_dataset_ca_aoa_5.5_vinf_11.1_freq_65.4", "xf-naca6412-il-1000000", 11.1, 5.5, 1.171, 16.935, 0.72, 1.7, 2.8, 65.4) ,



  ]


  # Run each simulation sequentially
  for (run_name, airfoil_name, magVinf, AOA, b, ar, tr, lambda, gamma, frequency) in simulation_params
    run_single_simulation(run_name, airfoil_name, magVinf, AOA, b, ar, tr, lambda, gamma, frequency)
  end

end

schedule_simulations()