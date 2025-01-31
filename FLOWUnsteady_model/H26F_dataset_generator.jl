# Import necessary modules
import FLOWUnsteady as uns
import FLOWVLM as vlm
import FLOWVPM as vpm

# Function to run a single simulation
function run_single_simulation(run_name, rotor_file, RPM, magVinf, J, AOA, pitch, tilt, CW)
    
    save_path       = run_name           # Where to save this simulation
    paraview        = false               # Whether to visualize with Paraview
    
    # Rotor geometry
    rotor_data_path = joinpath("database")
    # rotor_file      = "H26F_scaled.csv"             # Rotor geometry
    data_path       = rotor_data_path 
    # data_path       = uns.def_data_path  # Path to rotor database
    xfoil           = true
    ncrit           = 9


    # Rotor Orientation
    # tilt            = 90.0
    yaw             = 0.0


    # Discretization
    n               = 20                 # Number of blade elements per blade
    r               = 1/5                # Geometric expansion of elements

    # Read radius of this rotor and number of blades
    R, B            = uns.read_rotor(rotor_file; data_path=data_path)[[1,3]]

    # Operating conditions
    rho             = 1.225                      # (kg/m^3) air density
    mu              = 1.81e-5                    # (kg/ms) air dynamic viscosity
    speedofsound    = 342.35                     # (m/s) speed of sound
    # magVinf         = J * RPM / 60 * (2*R)
    Vinf(X, t)      = magVinf*[cosd(AOA), 0, sind(AOA)] # (m/s) freestream velocity vector

    ReD             = 2*pi*RPM/60*R * rho/mu * 2*R      # Diameter-based Reynolds number
    Matip           = 2*pi*RPM/60*R / speedofsound      # Tip Mach number

    println("""
        RPM:    $(RPM)
        Vinf:   $(Vinf(zeros(3), 0)) m/s
        Matip:  $(round(Matip, digits=3))
        ReD:    $(round(ReD, digits=0))
    """)

    # Solver parameters
    VehicleType     = uns.UVLMVehicle      # Unsteady solver
    const_solution  = VehicleType==uns.QVLMVehicle  # Whether to assume that the
                                                # solution is constant or not
    
    nrevs           = 5                    # Number of revolutions in simulation
    nsteps_per_rev  = 36                   # Time steps per revolution
    nsteps          = nrevs * nsteps_per_rev # Number of time steps
    ttot            = nsteps/nsteps_per_rev / (RPM/60)  # (s) total simulation time
    
    # VPM particle shedding
    p_per_step      = 2                         # Sheds per time step
    shed_starting   = true                      # Whether to shed starting vortex
    shed_unsteady   = true                      # Whether to shed vorticity from unsteady loading
    max_particles   = ((2*n+1)*B)*nsteps*p_per_step + 1 # Maximum number of particles

    # Regularization
    sigma_rotor_surf= R/40                      # Rotor-on-VPM smoothing radius
    lambda_vpm      = 2.125                     # VPM core overlap
                                                # VPM smoothing radius
    sigma_vpm_overwrite = lambda_vpm * 2*pi*R/(nsteps_per_rev*p_per_step)

    # Rotor solver
    vlm_rlx         = 0.7                       # VLM relaxation <-- this also applied to rotors
    hubtiploss_correction = vlm.hubtiploss_nocorrection # Hub and tip loss correction

    # VPM solver
    vpm_viscous     = vpm.Inviscid()            # VPM viscous diffusion scheme

    if VehicleType == uns.QVLMVehicle
        # NOTE: If the quasi-steady solver is used, this mutes warnings regarding
        #       potential colinear vortex filaments. This is needed since the
        #       quasi-steady solver will probe induced velocities at the lifting
        #       line of the blade
        uns.vlm.VLMSolver._mute_warning(true)
    end

    # Generate rotor

    println("Generating geometry...")
    rotor = uns.generate_rotor(rotor_file; pitch=pitch,
                                            n=n, CW=CW, blade_r=r,
                                            altReD=[RPM, J, mu/rho],
                                            xfoil=xfoil,
                                            ncrit=ncrit,
                                            data_path=data_path,
                                            verbose=true,
                                            verbose_xfoil=false,
                                            plot_disc=false);
                            
    rotor_sys_pos = [0.0 ,0.0 ,0.0]
    rotor_sys_axis = uns.gt.rotation_matrix(0.0, tilt, yaw)
    vlm.setcoordsystem(rotor, rotor_sys_pos, rotor_sys_axis)
    
    # Generate vehicle
    println("Generating vehicle...")
    system = vlm.WingSystem()                   # System of all FLOWVLM objects
    vlm.addwing(system, "Rotor", rotor)
    
    rotors = [rotor]                            # Defining this rotor as its own system
    rotor_systems = (rotors, )                  # All systems of rotors
    
    wake_system = vlm.WingSystem()              # System that will shed a VPM wake
    
    if VehicleType != uns.QVLMVehicle
        vlm.addwing(wake_system, "Rotor", rotor)
    end

    vehicle = VehicleType(   system;
                            rotor_systems=rotor_systems,
                            wake_system=wake_system
                         );

    # Define maneuver
    Vvehicle(t) = zeros(3)
    
    anglevehicle(t) = zeros(3)
    
    RPMcontrol(t) = 1.0
    
    angles = ()                                 # Angle of each tilting system (none)
    RPMs = (RPMcontrol, )                       # RPM of each rotor system
    
    maneuver = uns.KinematicManeuver(angles, RPMs, Vvehicle, anglevehicle)

    # Simualtion definition
    Vref = 0.0                                  # Reference velocity to scale maneuver by
    RPMref = RPM                                # Reference RPM to scale maneuver by

    Vinit = Vref*Vvehicle(0)                    # Initial vehicle velocity
    Winit = pi/180*(anglevehicle(1e-6) - anglevehicle(0))/(1e-6*ttot)  # Initial angular velocity

    simulation = uns.Simulation(vehicle, maneuver, Vref, RPMref, ttot;
                                                    Vinit=Vinit, Winit=Winit);

    # Generate rotor monitor
    figs, figaxs = [], []
    monitor_rotor = uns.generate_monitor_rotors(rotors, J, rho, RPM, nsteps, AOA, pitch, tilt, yaw;
                                                t_scale=RPM/60,
                                                t_lbl="Revolutions",
                                                out_figs=figs,
                                                out_figaxs=figaxs,
                                                save_path=save_path,
                                                run_name=run_name,
                                                figname="rotor monitor")

    # Run the simulation
    println("Running simulation $run_name...")
    uns.run_simulation(simulation, nsteps;
                        # ----- SIMULATION OPTIONS -------------
                        Vinf=Vinf,
                        rho=rho, mu=mu, sound_spd=speedofsound,
                        # ----- SOLVERS OPTIONS ----------------
                        p_per_step=p_per_step,
                        max_particles=max_particles,
                        vpm_viscous=vpm_viscous,
                        sigma_vlm_surf=sigma_rotor_surf,
                        sigma_rotor_surf=sigma_rotor_surf,
                        sigma_vpm_overwrite=sigma_vpm_overwrite,
                        vlm_rlx=vlm_rlx,
                        hubtiploss_correction=hubtiploss_correction,
                        shed_unsteady=shed_unsteady,
                        shed_starting=shed_starting,
                        extra_runtime_function=monitor_rotor,
                        # ----- OUTPUT OPTIONS ------------------
                        save_path=save_path,
                        run_name=run_name,
                        );


    # Visualization with Paraview
    if paraview
        println("Calling Paraview for $run_name...")
        files = joinpath(save_path, run_name * "_pfield...xmf;")
        for bi in 1:B
            files *= run_name * "_Rotor_Blade$(bi)_loft...vtk;"
            files *= run_name * "_Rotor_Blade$(bi)_vlm...vtk;"
        end
        run(`paraview --data=$(files)`)
    end
end

# Function to schedule multiple simulations
function schedule_simulations()
    # List of simulation parameters for each job
    simulation_params = [

        # Training Cases:
        # ("rotor_dataset_H26F_rpm_1394.0_vinf_11.4_-0.5_92.7",    "H26F_scaled.csv", 1394.0, 11.4, 0.001, -1.6, -0.5, 92.7,  false),
        # ("rotor_dataset_H26F_rpm_1221.0_vinf_9.8_-1.7_86.7",     "H26F_scaled.csv", 1221.0, 9.8, 0.001, 7.3, -1.7, 86.7,    false),
        # ("rotor_dataset_H26F_rpm_1518.0_vinf_1.6_-1.0_86.0",     "H26F_scaled.csv", 1518.0, 1.6, 0.001, -1.3, -1.0, 86.0,   false),
        # ("rotor_dataset_H26F_rpm_1790.0_vinf_17.1_0.1_94.4",     "H26F_scaled.csv", 1790.0, 17.1, 0.001, 5.4, 0.1, 94.4,    false),
        # ("rotor_dataset_H26F_rpm_1883.0_vinf_4.5_-0.3_86.0",     "H26F_scaled.csv", 1883.0, 4.5, 0.001, 0.5, -0.3, 86.0,    false),
        # ("rotor_dataset_H26F_rpm_1873.0_vinf_16.7_0.0_85.9",     "H26F_scaled.csv", 1873.0, 16.7, 0.001, 3.7, 0.0, 85.9,    false),
        # ("rotor_dataset_H26F_rpm_1985.0_vinf_5.0_2.0_92.9",      "H26F_scaled.csv", 1985.0, 5.0, 0.001, 3.1, 2.0, 92.9,     false),
        # ("rotor_dataset_H26F_rpm_1246.0_vinf_7.3_-0.9_87.9",     "H26F_scaled.csv", 1246.0, 7.3, 0.001, 4.1, -0.9, 87.9,    false),
        # ("rotor_dataset_H26F_rpm_1095.0_vinf_4.3_0.8_91.1",      "H26F_scaled.csv", 1095.0, 4.3, 0.001, -0.8, 0.8, 91.1,    false),
        # ("rotor_dataset_H26F_rpm_1173.0_vinf_5.4_-1.6_90.8",     "H26F_scaled.csv", 1173.0, 5.4, 0.001, -1.5, -1.6, 90.8,   false),
        
        # ("rotor_dataset_H26F_rpm_1562.0_vinf_6.9_-1.4_86.1",     "H26F_scaled.csv", 1562.0, 6.9, 0.001, 6.0, -1.4, 86.1,    false),
        # ("rotor_dataset_H26F_rpm_1698.0_vinf_6.6_-1.9_91.7",     "H26F_scaled.csv", 1698.0, 6.6, 0.001, 0.7, -1.9, 91.7,    false),
        # ("rotor_dataset_H26F_rpm_1028.0_vinf_13.1_1.2_93.8",     "H26F_scaled.csv", 1028.0, 13.1, 0.001, 0.3, 1.2, 93.8,    false),
        # ("rotor_dataset_H26F_rpm_1206.0_vinf_19.3_-1.5_94.0",    "H26F_scaled.csv", 1206.0, 19.3, 0.001, 6.8, -1.5, 94.0,   false),
        # ("rotor_dataset_H26F_rpm_1121.0_vinf_17.4_-1.2_92.7",    "H26F_scaled.csv", 1121.0, 17.4, 0.001, -0.2, -1.2, 92.7,  false),
        # ("rotor_dataset_H26F_rpm_1026.0_vinf_8.6_-1.3_88.8",     "H26F_scaled.csv", 1026.0, 8.6, 0.001, 7.9, -1.3, 88.8,    false),
        # ("rotor_dataset_H26F_rpm_1368.0_vinf_12.5_1.5_85.1",     "H26F_scaled.csv", 1368.0, 12.5, 0.001, 7.1, 1.5, 85.1,    false),
        # ("rotor_dataset_H26F_rpm_1239.0_vinf_10.7_-0.8_89.3",    "H26F_scaled.csv", 1239.0, 10.7, 0.001, -1.3, -0.8, 89.3,  false),
        # ("rotor_dataset_H26F_rpm_1876.0_vinf_1.7_0.5_89.4",      "H26F_scaled.csv", 1876.0, 1.7, 0.001, 0.7, 0.5, 89.4,     false),
        # ("rotor_dataset_H26F_rpm_1614.0_vinf_9.7_0.7_92.9",      "H26F_scaled.csv", 1614.0, 9.7, 0.001, 2.8, 0.7, 92.9,     false),
        
        # ("rotor_dataset_H26F_rpm_1182.0_vinf_6.7_0.5_90.0",   "H26F_scaled.csv", 1182.0, 6.7, 0.001, 4.4, 0.5, 90.0,      false),
        # ("rotor_dataset_H26F_rpm_1210.0_vinf_7.2_-1.5_89.3",  "H26F_scaled.csv", 1210.0, 7.2, 0.001, 6.2, -1.5, 89.3,     false),
        # ("rotor_dataset_H26F_rpm_1094.0_vinf_6.1_0.5_88.2",   "H26F_scaled.csv", 1094.0, 6.1, 0.001, 3.5, 0.5, 88.2,      false),
        # ("rotor_dataset_H26F_rpm_1988.0_vinf_17.9_-1.9_91.4", "H26F_scaled.csv", 1988.0, 17.9, 0.001, 5.3, -1.9, 91.4,    false),
        # ("rotor_dataset_H26F_rpm_1160.0_vinf_12.1_0.1_92.1",  "H26F_scaled.csv", 1160.0, 12.1, 0.001, -1.3, 0.1, 92.1,    false),
        # ("rotor_dataset_H26F_rpm_1215.0_vinf_3.8_1.2_94.4",   "H26F_scaled.csv", 1215.0, 3.8, 0.001, 6.5, 1.2, 94.4,      false),
        # ("rotor_dataset_H26F_rpm_1014.0_vinf_10.9_-1.8_89.5", "H26F_scaled.csv", 1014.0, 10.9, 0.001, -1.9, -1.8, 89.5,   false),
        # ("rotor_dataset_H26F_rpm_1797.0_vinf_16.6_-0.9_91.6", "H26F_scaled.csv", 1797.0, 16.6, 0.001, 3.8, -0.9, 91.6,    false),
        # ("rotor_dataset_H26F_rpm_1835.0_vinf_8.7_1.6_91.5",   "H26F_scaled.csv", 1835.0, 8.7, 0.001, 4.0, 1.6, 91.5,      false),
        # ("rotor_dataset_H26F_rpm_1275.0_vinf_2.7_-0.5_92.8",  "H26F_scaled.csv", 1275.0, 2.7, 0.001, 7.2, -0.5, 92.8,     false),
        
        # ("rotor_dataset_H26F_rpm_1150.0_vinf_17.6_0.4_93.4",  "H26F_scaled.csv", 1150.0, 17.6, 0.001, 0.6, 0.4, 93.4,     false),
        # ("rotor_dataset_H26F_rpm_1319.0_vinf_9.2_-0.7_89.2",  "H26F_scaled.csv", 1319.0, 9.2, 0.001, 1.8, -0.7, 89.2,     false),
        # ("rotor_dataset_H26F_rpm_1072.0_vinf_5.9_1.1_91.5",   "H26F_scaled.csv", 1072.0, 5.9, 0.001, 7.9, 1.1, 91.5,      false),
        # ("rotor_dataset_H26F_rpm_1896.0_vinf_3.6_-0.8_87.9",  "H26F_scaled.csv", 1896.0, 3.6, 0.001, 2.7, -0.8, 87.9,     false),
        # ("rotor_dataset_H26F_rpm_1930.0_vinf_6.1_-1.5_87.3",  "H26F_scaled.csv", 1930.0, 6.1, 0.001, 7.8, -1.5, 87.3,     false),
        # ("rotor_dataset_H26F_rpm_1688.0_vinf_12.0_0.2_93.5",  "H26F_scaled.csv", 1688.0, 12.0, 0.001, 1.3, 0.2, 93.5,     false),
        # ("rotor_dataset_H26F_rpm_1241.0_vinf_12.1_-1.1_91.5", "H26F_scaled.csv", 1241.0, 12.1, 0.001, 0.6, -1.1, 91.5,    false),
        # ("rotor_dataset_H26F_rpm_1170.0_vinf_16.0_-0.7_91.0", "H26F_scaled.csv", 1170.0, 16.0, 0.001, 5.1, -0.7, 91.0,    false),
        # ("rotor_dataset_H26F_rpm_1456.0_vinf_10.0_0.7_88.7",  "H26F_scaled.csv", 1456.0, 10.0, 0.001, -1.5, 0.7, 88.7,    false),
        # ("rotor_dataset_H26F_rpm_1453.0_vinf_18.0_-0.2_91.8", "H26F_scaled.csv", 1453.0, 18.0, 0.001, 5.0, -0.2, 91.8,    false),

        # ("rotor_dataset_H26F_rpm_1854.0_vinf_18.4_1.6_89.0",    "H26F_scaled.csv", 1854.0, 18.4, 0.001, -0.7, 1.6, 89.0,    false),
        # ("rotor_dataset_H26F_rpm_1677.0_vinf_16.6_0.9_87.9",    "H26F_scaled.csv", 1677.0, 16.6, 0.001, 3.6, 0.9, 87.9,     false),
        # ("rotor_dataset_H26F_rpm_1138.0_vinf_6.3_0.8_87.6",     "H26F_scaled.csv", 1138.0, 6.3, 0.001, 5.8, 0.8, 87.6,      false),
        # ("rotor_dataset_H26F_rpm_1687.0_vinf_16.4_1.6_89.1",   "H26F_scaled.csv", 1687.0, 16.4, 0.001, 0.0, 1.6, 89.1,    false),
        # ("rotor_dataset_H26F_rpm_1774.0_vinf_8.7_1.2_87.8",     "H26F_scaled.csv", 1774.0, 8.7, 0.001, 7.8, 1.2, 87.8,      false),
        # ("rotor_dataset_H26F_rpm_1802.0_vinf_18.6_1.9_93.9",   "H26F_scaled.csv", 1802.0, 18.6, 0.001, 4.6, 1.9, 93.9,    false),
        # ("rotor_dataset_H26F_rpm_1025.0_vinf_19.4_1.7_92.2",    "H26F_scaled.csv", 1025.0, 19.4, 0.001, 0.9, 1.7, 92.2,     false),
        # ("rotor_dataset_H26F_rpm_1774.0_vinf_2.9_1.5_89.5",     "H26F_scaled.csv", 1774.0, 2.9, 0.001, 5.0, 1.5, 89.5,      false),
        # ("rotor_dataset_H26F_rpm_1029.0_vinf_6.3_-1.2_88.0",    "H26F_scaled.csv", 1029.0, 6.3, 0.001, 0.4, -1.2, 88.0,    false),
        # ("rotor_dataset_H26F_rpm_1407.0_vinf_6.8_1.3_92.6",     "H26F_scaled.csv", 1407.0, 6.8, 0.001, 1.4, 1.3, 92.6,      false),

        # ("rotor_dataset_H26F_rpm_1712.0_vinf_11.3_-0.2_88.3",   "H26F_scaled.csv", 1712.0, 11.3, 0.001, 6.3, -0.2, 88.3,    false),
        # ("rotor_dataset_H26F_rpm_1130.0_vinf_6.2_1.0_90.6",     "H26F_scaled.csv", 1130.0, 6.2, 0.001, 3.1, 1.0, 90.6,      false),
        # ("rotor_dataset_H26F_rpm_1765.0_vinf_8.6_0.0_85.1",     "H26F_scaled.csv", 1765.0, 8.6, 0.001, 4.4, 0.0, 85.1,      false),
        # ("rotor_dataset_H26F_rpm_1881.0_vinf_10.0_0.2_92.5",    "H26F_scaled.csv", 1881.0, 10.0, 0.001, 1.2, 0.2, 92.5,     false),
        # ("rotor_dataset_H26F_rpm_1333.0_vinf_12.7_-1.0_86.6",   "H26F_scaled.csv", 1333.0, 12.7, 0.001, 6.8, -1.0, 86.6,    false),
        # ("rotor_dataset_H26F_rpm_1731.0_vinf_11.0_-2.0_92.9",   "H26F_scaled.csv", 1731.0, 11.0, 0.001, 5.7, -2.0, 92.9,    false),
        # ("rotor_dataset_H26F_rpm_1393.0_vinf_15.0_0.4_90.8",    "H26F_scaled.csv", 1393.0, 15.0, 0.001, 0.8, 0.4, 90.8,     false),
        # ("rotor_dataset_H26F_rpm_1600.0_vinf_6.5_1.7_87.3",     "H26F_scaled.csv", 1600.0, 6.5, 0.001, 7.4, 1.7, 87.3,      false),
        # ("rotor_dataset_H26F_rpm_1761.0_vinf_19.1_-0.3_85.9",   "H26F_scaled.csv", 1761.0, 19.1, 0.001, 1.9, -0.3, 85.9,    false),
        # ("rotor_dataset_H26F_rpm_1958.0_vinf_9.0_-1.0_90.5",    "H26F_scaled.csv", 1958.0, 9.0, 0.001, 1.1, -1.0, 90.5,    false),

        # ("rotor_dataset_H26F_rpm_1939.0_vinf_10.8_-0.4_89.5",   "H26F_scaled.csv", 1939.0, 10.8, 0.001, -1.2, -0.4, 89.5,   false),
        # ("rotor_dataset_H26F_rpm_1288.0_vinf_6.0_0.6_91.2",     "H26F_scaled.csv", 1288.0, 6.0, 0.001, -1.4, 0.6, 91.2,     false),
        # ("rotor_dataset_H26F_rpm_1565.0_vinf_13.4_0.2_89.5",    "H26F_scaled.csv", 1565.0, 13.4, 0.001, 1.8, 0.2, 89.5,     false),
        # ("rotor_dataset_H26F_rpm_1789.0_vinf_1.3_1.8_92.8",     "H26F_scaled.csv", 1789.0, 1.3, 0.001, 6.7, 1.8, 92.8,      false),
        # ("rotor_dataset_H26F_rpm_1098.0_vinf_15.8_0.1_90.6",    "H26F_scaled.csv", 1098.0, 15.8, 0.001, 7.8, 0.1, 90.6,     false),
        # ("rotor_dataset_H26F_rpm_1704.0_vinf_19.4_1.0_90.4",    "H26F_scaled.csv", 1704.0, 19.4, 0.001, 0.1, 1.0, 90.4,     false),
        # ("rotor_dataset_H26F_rpm_1105.0_vinf_1.1_0.8_91.3",     "H26F_scaled.csv", 1105.0, 1.1, 0.001, 3.4, 0.8, 91.3,      false),
        # ("rotor_dataset_H26F_rpm_1636.0_vinf_19.0_0.5_90.6",    "H26F_scaled.csv", 1636.0, 19.0, 0.001, 7.0, 0.5, 90.6,     false),
        # ("rotor_dataset_H26F_rpm_1733.0_vinf_1.7_-0.2_91.1",    "H26F_scaled.csv", 1733.0, 1.7, 0.001, -1.7, -0.2, 91.1,    false),
        ("rotor_dataset_H26F_rpm_1683.0_vinf_13.7_-0.0_92.6",   "H26F_scaled.csv", 1683.0, 13.7, 0.001, 6.9, -0.0, 92.6,    false),
        
        # ("rotor_dataset_H26F_rpm_1280.0_vinf_5.4_-1.6_94.0",    "H26F_scaled.csv", 1280.0, 5.4, 0.001, -1.3, -1.6, 94.0,    false),
        # ("rotor_dataset_H26F_rpm_1280.0_vinf_5.7_0.8_92.9",     "H26F_scaled.csv", 1280.0, 5.7, 0.001, 6.2, 0.8, 92.9,      false),
        # ("rotor_dataset_H26F_rpm_1161.0_vinf_12.9_-1.4_92.3",   "H26F_scaled.csv", 1161.0, 12.9, 0.001, 1.4, -1.4, 92.3,    false),
        # ("rotor_dataset_H26F_rpm_1237.0_vinf_9.3_0.5_87.8",     "H26F_scaled.csv", 1237.0, 9.3, 0.001, -1.0, 0.5, 87.8,     false),
        # ("rotor_dataset_H26F_rpm_1788.0_vinf_2.6_-1.4_90.5",    "H26F_scaled.csv", 1788.0, 2.6, 0.001, 0.4, -1.4, 90.5,     false),
        # ("rotor_dataset_H26F_rpm_1199.0_vinf_19.8_1.0_86.3",    "H26F_scaled.csv", 1199.0, 19.8, 0.001, 7.3, 1.0, 86.3,     false),
        # ("rotor_dataset_H26F_rpm_1385.0_vinf_6.6_0.3_89.8",     "H26F_scaled.csv", 1385.0, 6.6, 0.001, -1.2, 0.3, 89.8,     false),
        # ("rotor_dataset_H26F_rpm_1335.0_vinf_4.3_1.3_93.4",     "H26F_scaled.csv", 1335.0, 4.3, 0.001, 0.1, 1.3, 93.4,      false),
        # ("rotor_dataset_H26F_rpm_1961.0_vinf_8.7_-1.0_93.9",    "H26F_scaled.csv", 1961.0, 8.7, 0.001, 3.8, -1.0, 93.9,     false),
        # ("rotor_dataset_H26F_rpm_1236.0_vinf_3.6_1.1_88.2",     "H26F_scaled.csv", 1236.0, 3.6, 0.001, 5.4, 1.1, 88.2,      false),

        
        # Testing Cases:
        # ("rotor_dataset_H26F_rpm_1113.0_vinf_8.7_-1.7_92.0",  "H26F_scaled.csv", 1113.0, 8.7, 0.001, -0.9, -1.7, 92.0,    false),
        # ("rotor_dataset_H26F_rpm_1809.0_vinf_9.7_-1.1_90.3",  "H26F_scaled.csv", 1809.0, 9.7, 0.001, 2.4, -1.1, 90.3,     false),
        # ("rotor_dataset_H26F_rpm_1961.0_vinf_19.0_-0.5_88.0", "H26F_scaled.csv", 1961.0, 19.0, 0.001, -0.2, -0.5, 88.0,   false),
        # ("rotor_dataset_H26F_rpm_1617.0_vinf_12.9_0.8_85.4",  "H26F_scaled.csv", 1617.0, 12.9, 0.001, 5.9, 0.8, 85.4,     false),
        # ("rotor_dataset_H26F_rpm_1638.0_vinf_15.7_-0.2_93.9", "H26F_scaled.csv", 1638.0, 15.7, 0.001, 3.5, -0.2, 93.9,    false),
        # ("rotor_dataset_H26F_rpm_1638.0_vinf_1.4_1.5_88.2",   "H26F_scaled.csv", 1638.0, 1.4, 0.001, 2.8, 1.5, 88.2,      false),
        # ("rotor_dataset_H26F_rpm_1700.0_vinf_9.6_-0.4_94.9",  "H26F_scaled.csv", 1700.0, 9.6, 0.001, 3.2, -0.4, 94.9,     false),
        # ("rotor_dataset_H26F_rpm_1143.0_vinf_11.8_-1.6_89.1", "H26F_scaled.csv", 1143.0, 11.8, 0.001, 7.5, -1.6, 89.1,    false),
        # ("rotor_dataset_H26F_rpm_1290.0_vinf_5.2_-0.9_94.9",  "H26F_scaled.csv", 1290.0, 5.2, 0.001, 1.7, -0.9, 94.9,     false),
        # ("rotor_dataset_H26F_rpm_1283.0_vinf_19.0_1.3_85.8",  "H26F_scaled.csv", 1283.0, 19.0, 0.001, 4.5, 1.3, 85.8,     false),
        #--------------------------------------------------------------------------------------

        # ("dataset_H26F_500_0.001_00_00_90", "H26F.csv", 500, 0.001, 0.0, 0.0, 90.0, false),
        # ("dataset_H26F_1250_0.001_00_00_90", "H26F.csv", 1250, 0.001, 0.0, 0.0, 90.0, false),
        # ("dataset_H26F_2000_0.001_00_00_90", "H26F.csv", 2000, 0.001, 0.0, 0.0, 90.0, false)

        # Add more jobs as needed...
    ]

    # Run each simulation sequentially
    for (run_name, rotor_file, RPM, magVinf, J, AOA, pitch, tilt, CW) in simulation_params
        run_single_simulation(run_name, rotor_file, RPM, magVinf, J, AOA, pitch, tilt, CW)
    end
end

# Call the scheduler to run all simulations
schedule_simulations()
