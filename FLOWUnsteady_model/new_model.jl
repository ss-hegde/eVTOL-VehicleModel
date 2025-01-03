import FLOWUnsteady as uns
import FLOWVLM as vlm
import FLOWVPM as vpm

run_name        = "eVTOL-model"      # Name of this simulation
save_path       = run_name                  # Where to save this simulation
paraview        = true                      # Whether to visualize with Paraview


add_rotors      = true
verbose         = true
v_lvl       = 0

# ----------------- GEOMETRY PARAMETERS ----------------------------------------
k                 = 0.1301                # Scaling factor

rotor_data_path = joinpath("database")
rotor_file      = "H26F_scaled.csv"             # Rotor geometry
data_path       = rotor_data_path         # Path to rotor database
pitch           = 0.0                       # (deg) collective pitch of blades
CWs             = [true, false]                     # Clock-wise rotation
xfoil           = true                      # Whether to run XFOIL
ncrit           = 9                         # Turbulence criterion for XFOIL
n_rotor         = 20            # Number of blade elements per blade
r_rotor         = 1/10                  # Geometric expansion between blade elements

tilt            = 90.0
yaw             = 0.0

# Read radius of this rotor and number of blades
R, B            = uns.read_rotor(rotor_file; data_path=data_path)[[1,3]]

# # Operating conditions
# RPM             = 1600                      # RPM
# J               = 0.15                   # Advance ratio Vinf/(nD)
# AOA             = 1.0                         # (deg) Angle of attack (incidence angle)

# rho             = 1.225                     # (kg/m^3) air density
# mu              = 1.81e-5                   # (kg/ms) air dynamic viscosity
# speedofsound    = 342.35                    # (m/s) speed of sound

# magVinf         = J*RPM/60*(2*R)
# Vinf(X, t)      = magVinf*[cosd(AOA), sind(AOA), 0] # (m/s) freestream velocity vector
# ReD             = 2*pi*RPM/60*R * rho/mu * 2*R      # Diameter-based Reynolds number
# Matip           = 2*pi*RPM/60*R / speedofsound      # Tip Mach number

# Rotor positions
# L1_rotor
x_pos_rL1       = 0.36 * k
y_pos_rL1       = -3.842 * k
z_pos_rL1       = 0.5 * k 

# L2_rotor
x_pos_rL2       = 2.28 * k
y_pos_rL2       = -3.842 * k
z_pos_rL2       = 0.5 * k 

# L3_rotor
x_pos_rL3       = 4.2 * k
y_pos_rL3       = -3.842 * k
z_pos_rL3       = 0.5 * k 

# L4_rotor
x_pos_rL4       = 6.54 * k
y_pos_rL4       = -3.842 * k
z_pos_rL4       = 0.5 * k 

# R1_rotor
x_pos_rR1       = 0.36 * k
y_pos_rR1       = 3.842 * k
z_pos_rR1       = 0.5 * k 

# R2_rotor
x_pos_rR2       = 2.28 * k
y_pos_rR2       = 3.842 * k
z_pos_rR2       = 0.5 * k

# R3_rotor
x_pos_rR3       = 4.2 * k
y_pos_rR3       = 3.842 * k
z_pos_rR3       = 0.5 * k 

# R4_rotor
x_pos_rR4       = 6.54 * k
y_pos_rR4       = 3.842 * k
z_pos_rR4       = 0.5 * k

# ----------------- SIMULATION PARAMETERS --------------------------------------

# Operating conditions
RPM             = 1600                      # RPM
J               = 0.0001                    # Advance ratio Vinf/(nD)
AOA             = 0                         # (deg) Angle of attack (incidence angle)

rho             = 1.071778                  # (kg/m^3) air density
mu              = 1.85508e-5                # (kg/ms) air dynamic viscosity
speedofsound    = 342.35                    # (m/s) speed of sound

# NOTE: For cases with zero freestream velocity, it is recommended that a
#       negligible small velocity is used instead of zero in order to avoid
#       potential numerical instabilities (hence, J here is negligible small
#       instead of zero)

magVinf         = J*RPM/60*(2*R)
Vinf(X, t)      = magVinf*[cos(AOA*pi/180), sin(AOA*pi/180), 0]  # (m/s) freestream velocity vector

ReD             = 2*pi*RPM/60*R * rho/mu * 2*R      # Diameter-based Reynolds number
Matip           = 2*pi*RPM/60 * R / speedofsound    # Tip Mach number

println("""
    RPM:    $(RPM)
    Vinf:   $(Vinf(zeros(3), 0)) m/s
    Matip:  $(round(Matip, digits=3))
    ReD:    $(round(ReD, digits=0))
""")

# ----------------- SOLVER PARAMETERS ------------------------------------------

# Aerodynamic solver
VehicleType     = uns.UVLMVehicle           # Unsteady solver
# VehicleType     = uns.QVLMVehicle         # Quasi-steady solver
const_solution  = VehicleType==uns.QVLMVehicle  # Whether to assume that the
                                                # solution is constant or not
# Time parameters
nrevs           = 10                        # Number of revolutions in simulation
nsteps_per_rev  = 36                        # Time steps per revolution
nsteps          = const_solution ? 2 : nrevs*nsteps_per_rev # Number of time steps
ttot            = nsteps/nsteps_per_rev / (RPM/60)       # (s) total simulation time

# VPM particle shedding
p_per_step      = 4                         # Sheds per time step
shed_starting   = false                     # Whether to shed starting vortex
shed_unsteady   = true                      # Whether to shed vorticity from unsteady loading
unsteady_shedcrit = 0.001                   # Shed unsteady loading whenever circulation
                                            #  fluctuates by more than this ratio
max_particles   = ((2*n_rotor+1)*B)*nsteps*p_per_step + 1 # Maximum number of particles

# Regularization
sigma_rotor_surf= R/10                      # Rotor-on-VPM smoothing radius
lambda_vpm      = 2.125                     # VPM core overlap
                                            # VPM smoothing radius
sigma_vpm_overwrite = lambda_vpm * 2*pi*R/(nsteps_per_rev*p_per_step)
sigmafactor_vpmonvlm= 1                     # Shrink particles by this factor when
                                            #  calculating VPM-on-VLM/Rotor induced velocities

# Rotor solver
vlm_rlx         = 0.5                       # VLM relaxation <-- this also applied to rotors
hubtiploss_correction = ((0.4, 5, 0.1, 0.05), (2, 1, 0.25, 0.05)) # Hub and tip correction

# VPM solver
vpm_integration = vpm.euler                 # VPM temporal integration scheme
# vpm_integration = vpm.rungekutta3

vpm_viscous     = vpm.Inviscid()            # VPM viscous diffusion scheme
# vpm_viscous   = vpm.CoreSpreading(-1, -1, vpm.zeta_fmm; beta=100.0, itmax=20, tol=1e-1)

vpm_SFS         = vpm.SFS_none              # VPM LES subfilter-scale model
# vpm_SFS       = vpm.SFS_Cd_twolevel_nobackscatter
# vpm_SFS       = vpm.SFS_Cd_threelevel_nobackscatter
# vpm_SFS       = vpm.DynamicSFS(vpm.Estr_fmm, vpm.pseudo3level_positive;
#                                   alpha=0.999, maxC=1.0,
#                                   clippings=[vpm.clipping_backscatter])
# vpm_SFS       = vpm.DynamicSFS(vpm.Estr_fmm, vpm.pseudo3level_positive;
#                                   alpha=0.999, rlxf=0.005, minC=0, maxC=1
#                                   clippings=[vpm.clipping_backscatter],
#                                   controls=[vpm.control_sigmasensor],
#                                   )

# NOTE: In most practical situations, open rotors operate at a Reynolds number
#       high enough that viscous diffusion in the wake is actually negligible.
#       Hence, it does not make much of a difference whether we run the
#       simulation with viscous diffusion enabled or not. On the other hand,
#       such high Reynolds numbers mean that the wake quickly becomes turbulent
#       and it is crucial to use a subfilter-scale (SFS) model to accurately
#       capture the turbulent decay of the wake (turbulent diffusion).

if VehicleType == uns.QVLMVehicle
    # Mute warnings regarding potential colinear vortex filaments. This is
    # needed since the quasi-steady solver will probe induced velocities at the
    # lifting line of the blade
    uns.vlm.VLMSolver._mute_warning(true)
end



# ----------------- WAKE TREATMENT ---------------------------------------------
# NOTE: It is known in the CFD community that rotor simulations with an
#       impulsive RPM start (*i.e.*, 0 to RPM in the first time step, as opposed
#       to gradually ramping up the RPM) leads to the hub "fountain effect",
#       with the root wake reversing the flow near the hub.
#       The fountain eventually goes away as the wake develops, but this happens
#       very slowly, which delays the convergence of the simulation to a steady
#       state. To accelerate convergence, here we define a wake treatment
#       procedure that suppresses the hub wake for the first three revolutions,
#       avoiding the fountain effect altogether.
#       This is especially helpful in low and mid-fidelity simulations.

suppress_fountain   = true                  # Toggle

# Supress wake shedding on blade elements inboard of this r/R radial station
no_shedding_Rthreshold = suppress_fountain ? 0.35 : 0.0

# Supress wake shedding for this many time steps
no_shedding_nstepsthreshold = 3*nsteps_per_rev

omit_shedding = []          # Index of blade elements to supress wake shedding

# Function to suppress or activate wake shedding
function wake_treatment_supress(sim, args...; optargs...)

    # Case: start of simulation -> suppress shedding
    if sim.nt == 1

        # Identify blade elements on which to suppress shedding
        for i in 1:vlm.get_m(rotor)
            HS = vlm.getHorseshoe(rotor, i)
            CP = HS[5]

            if uns.vlm.norm(CP - vlm._get_O(rotor)) <= no_shedding_Rthreshold*R
                push!(omit_shedding, i)
            end
        end
    end

    # Case: sufficient time steps -> enable shedding
    if sim.nt == no_shedding_nstepsthreshold

        # Flag to stop suppressing
        omit_shedding .= -1

    end

    return false
end


# ----------------- 1) VEHICLE DEFINITION --------------------------------------
println("Generating geometry...")


if add_rotors
        
    if verbose; println("\t"^(v_lvl+1)*"Generating Rotors"); end;        
    
    rotors_left = vlm.Rotor[]
    
    if verbose; println("\t"^(v_lvl+2)*"Generating Rotor L1"); end;
    rotor_L1 = uns.generate_rotor(rotor_file; pitch=pitch,
                                            n=n_rotor, CW=CWs[1], blade_r=r_rotor,
                                            altReD=[RPM, J, mu/rho],
                                            xfoil=xfoil,
                                            ncrit=ncrit,
                                            data_path=data_path,
                                            verbose=true,
                                            verbose_xfoil=false,
                                            plot_disc=false,
                                            save_polars=save_path
                                            );
    
    rotor_L1_pos = [x_pos_rL1, y_pos_rL1, z_pos_rL1]
    rotor_L1_axis = uns.gt.rotation_matrix(0.0, 93.0, 0.0)
    vlm.setcoordsystem(rotor_L1, rotor_L1_pos, rotor_L1_axis)
    
    if verbose; println("\t"^(v_lvl+2)*"Generating Rotor L2"); end;
    rotor_L2 = uns.generate_rotor(rotor_file; pitch=pitch,
                                            n=n_rotor, CW=CWs[2], blade_r=r_rotor,
                                            altReD=[RPM, J, mu/rho],
                                            xfoil=xfoil,
                                            ncrit=ncrit,
                                            data_path=data_path,
                                            verbose=true,
                                            verbose_xfoil=false,
                                            plot_disc=false,
                                            save_polars=save_path
                                            );
    
    rotor_L2_pos = [x_pos_rL2, y_pos_rL2, z_pos_rL2]
    rotor_L2_axis = uns.gt.rotation_matrix(0.0, 93.0, 0.0)
    vlm.setcoordsystem(rotor_L2, rotor_L2_pos, rotor_L2_axis)
    
    if verbose; println("\t"^(v_lvl+2)*"Generating Rotor L3"); end;
    rotor_L3 = uns.generate_rotor(rotor_file; pitch=pitch,
                                            n=n_rotor, CW=CWs[1], blade_r=r_rotor,
                                            altReD=[RPM, J, mu/rho],
                                            xfoil=xfoil,
                                            ncrit=ncrit,
                                            data_path=data_path,
                                            verbose=true,
                                            verbose_xfoil=false,
                                            plot_disc=false,
                                            save_polars=save_path
                                            );
    
    rotor_L3_pos = [x_pos_rL3, y_pos_rL3, z_pos_rL3]
    rotor_L3_axis = uns.gt.rotation_matrix(0.0, 93.0, 0.0)
    vlm.setcoordsystem(rotor_L3, rotor_L3_pos, rotor_L3_axis)
    
    if verbose; println("\t"^(v_lvl+2)*"Generating Rotor L4"); end;
    rotor_L4 = uns.generate_rotor(rotor_file; pitch=pitch,
                                            n=n_rotor, CW=CWs[2], blade_r=r_rotor,
                                            altReD=[RPM, J, mu/rho],
                                            xfoil=xfoil,
                                            ncrit=ncrit,
                                            data_path=data_path,
                                            verbose=true,
                                            verbose_xfoil=false,
                                            plot_disc=false,
                                            save_polars=save_path
                                            );
    
    rotor_L4_pos = [x_pos_rL4, y_pos_rL4, z_pos_rL4]
    rotor_L4_axis = uns.gt.rotation_matrix(0.0, 93.0, 0.0)
    vlm.setcoordsystem(rotor_L4, rotor_L4_pos, rotor_L4_axis)
    
    push!(rotors_left, rotor_L1)
    push!(rotors_left, rotor_L2)
    push!(rotors_left, rotor_L3)
    push!(rotors_left, rotor_L4)
    
    rotors_right = vlm.Rotor[]
    
    if verbose; println("\t"^(v_lvl+2)*"Generating Rotor R1"); end;
    rotor_R1 = uns.generate_rotor(rotor_file; pitch=pitch,
                                            n=n_rotor, CW=CWs[2], blade_r=r_rotor,
                                            altReD=[RPM, J, mu/rho],
                                            xfoil=xfoil,
                                            ncrit=ncrit,
                                            data_path=data_path,
                                            verbose=true,
                                            verbose_xfoil=false,
                                            plot_disc=false,
                                            save_polars=save_path
                                            );
    
    rotor_R1_pos = [x_pos_rR1, y_pos_rR1, z_pos_rR1]
    rotor_R1_axis = uns.gt.rotation_matrix(0.0, 93.0, 0.0)
    vlm.setcoordsystem(rotor_R1, rotor_R1_pos, rotor_R1_axis)
    
    if verbose; println("\t"^(v_lvl+2)*"Generating Rotor R2"); end;
    rotor_R2 = uns.generate_rotor(rotor_file; pitch=pitch,
                                            n=n_rotor, CW=CWs[1], blade_r=r_rotor,
                                            altReD=[RPM, J, mu/rho],
                                            xfoil=xfoil,
                                            ncrit=ncrit,
                                            data_path=data_path,
                                            verbose=true,
                                            verbose_xfoil=false,
                                            plot_disc=false,
                                            save_polars=save_path
                                            );
    
    rotor_R2_pos = [x_pos_rR2, y_pos_rR2, z_pos_rR2]
    rotor_R2_axis = uns.gt.rotation_matrix(0.0, 93.0, 0.0)
    vlm.setcoordsystem(rotor_R2, rotor_R2_pos, rotor_R2_axis)
    
    if verbose; println("\t"^(v_lvl+2)*"Generating Rotor R3"); end;
    rotor_R3 = uns.generate_rotor(rotor_file; pitch=pitch,
                                            n=n_rotor, CW=CWs[2], blade_r=r_rotor,
                                            altReD=[RPM, J, mu/rho],
                                            xfoil=xfoil,
                                            ncrit=ncrit,
                                            data_path=data_path,
                                            verbose=true,
                                            verbose_xfoil=false,
                                            plot_disc=false,
                                            save_polars=save_path
                                            );
    
    rotor_R3_pos = [x_pos_rR3, y_pos_rR3, z_pos_rR3]
    rotor_R3_axis = uns.gt.rotation_matrix(0.0, 93.0, 0.0)
    vlm.setcoordsystem(rotor_R3, rotor_R3_pos, rotor_R3_axis)
    
    if verbose; println("\t"^(v_lvl+2)*"Generating Rotor R4"); end;
    rotor_R4 = uns.generate_rotor(rotor_file; pitch=pitch,
                                            n=n_rotor, CW=CWs[1], blade_r=r_rotor,
                                            altReD=[RPM, J, mu/rho],
                                            xfoil=xfoil,
                                            ncrit=ncrit,
                                            data_path=data_path,
                                            verbose=true,
                                            verbose_xfoil=false,
                                            plot_disc=false,
                                            save_polars=save_path
                                            );
    
    rotor_R4_pos = [x_pos_rR4, y_pos_rR4, z_pos_rR4]
    rotor_R4_axis = uns.gt.rotation_matrix(0.0, 93.0, 0.0)
    vlm.setcoordsystem(rotor_R4, rotor_R4_pos, rotor_R4_axis)
    
    
    push!(rotors_right, rotor_R1)
    push!(rotors_right, rotor_R2)
    push!(rotors_right, rotor_R3)
    push!(rotors_right, rotor_R4)
end

# Vector of sigle rotor for monitors of each single rotor
if add_rotors    
    all_rotors = vlm.Rotor[]
    push!(all_rotors, rotor_L1)
    push!(all_rotors, rotor_L2)
    push!(all_rotors, rotor_L3)
    push!(all_rotors, rotor_L4)
    push!(all_rotors, rotor_R1)
    push!(all_rotors, rotor_R2)
    push!(all_rotors, rotor_R3)
    push!(all_rotors, rotor_R4)
end

# Generate rotor

println("Generating vehicle...")

# Generate vehicle
system = vlm.WingSystem()                   # System of all FLOWVLM objects

if add_rotors
    vlm.addwing(system, "Rotor_L1", rotor_L1)
    vlm.setRPM(rotor_L1, RPM)

    vlm.addwing(system, "Rotor_L2", rotor_L2)
    vlm.setRPM(rotor_L2, RPM)

    vlm.addwing(system, "Rotor_L3", rotor_L3)
    vlm.setRPM(rotor_L3, RPM)

    vlm.addwing(system, "Rotor_L4", rotor_L4)
    vlm.setRPM(rotor_L4, RPM)

    vlm.addwing(system, "Rotor_R1", rotor_R1)
    vlm.setRPM(rotor_R1, RPM)

    vlm.addwing(system, "Rotor_R2", rotor_R2)
    vlm.setRPM(rotor_R2, RPM)

    vlm.addwing(system, "Rotor_R3", rotor_R3)
    vlm.setRPM(rotor_R3, RPM)

    vlm.addwing(system, "Rotor_R4", rotor_R4)
    vlm.setRPM(rotor_R4, RPM)
end

if add_rotors
    rotor_systems = (rotors_left, rotors_right);
else
    rotor_systems = ();
end


wake_system = vlm.WingSystem()              # System that will shed a VPM wake
                                            # NOTE: Do NOT include rotor when using the quasi-steady solver
# if VehicleType != uns.QVLMVehicle
#     # vlm.addwing(wake_system, "Rotor", rotor)

#     vlm.addwing(wake_system, "Rotor_L1", rotor_L1)
#     vlm.addwing(wake_system, "Rotor_L2", rotor_L2)
#     vlm.addwing(wake_system, "Rotor_L3", rotor_L3)
#     vlm.addwing(wake_system, "Rotor_L4", rotor_L4)
#     vlm.addwing(wake_system, "Rotor_R1", rotor_R1)
#     vlm.addwing(wake_system, "Rotor_R2", rotor_R2)
#     vlm.addwing(wake_system, "Rotor_R3", rotor_R3)
#     vlm.addwing(wake_system, "Rotor_R4", rotor_R4)
# end

if add_rotors
    rotors = vcat(rotors_left, rotors_right);
end
if add_rotors
    if VehicleType==uns.VLMVehicle

        for (i, rotor) in enumerate(rotors)
            vlm.addwing(wake_system, "Rotor$i", rotor)
        end
        
    else
        uns.vlm.VLMSolver._mute_warning(true)
    end
end

vehicle = uns.VLMVehicle(   system;
                            rotor_systems=rotor_systems,
                            wake_system=wake_system
                         );


# ------------- 2) MANEUVER DEFINITION -----------------------------------------
# Non-dimensional translational velocity of vehicle over time
Vvehicle(t) = zeros(3)

# Angle of the vehicle over time
anglevehicle(t) = zeros(3)

# RPM control input over time (RPM over `RPMref`)
RPMcontrol(t) = 1.0

angles = ()                                 # Angle of each tilting system (none)
RPMs = (RPMcontrol, )                       # RPM of each rotor system

maneuver = uns.KinematicManeuver(angles, RPMs, Vvehicle, anglevehicle)


# ------------- 3) SIMULATION DEFINITION ---------------------------------------

Vref = 0.0                                  # Reference velocity to scale maneuver by
RPMref = RPM                                # Reference RPM to scale maneuver by
Vinit = Vref*Vvehicle(0)                    # Initial vehicle velocity
Winit = pi/180*(anglevehicle(1e-6) - anglevehicle(0))/(1e-6*ttot)  # Initial angular velocity

simulation = uns.Simulation(vehicle, maneuver, Vref, RPMref, ttot;
                                                    Vinit=Vinit, Winit=Winit);

# Restart simulation
restart_file = nothing

# NOTE: Uncomment the following line to restart a previous simulation.
#       Point it to a particle field file (with its full path) at a specific
#       time step, and `run_simulation` will start this simulation with the
#       particle field found in the restart simulation.

# restart_file = "/path/to/a/previous/simulation/rotorhover-example_pfield.360"


# ------------- 4) MONITORS DEFINITIONS ----------------------------------------
# Collect all monitors here
monitors = []
# Generate rotor monitor

if add_rotors

    for (si, rotor) in enumerate(vehicle.rotor_systems[1])

        monitor_name = "rotor_L$(si)"
        rotor_vec = [rotor, ]
        rotor_monitor_L = uns.generate_monitor_rotors(rotors, J, rho, RPM, nsteps, AOA, pitch, tilt, yaw;
                                                                t_scale=RPM/60,
                                                                t_lbl="Revolutions",
                                                                out_figs=figs,
                                                                out_figaxs=figaxs,
                                                                save_path=save_path,
                                                                run_name=run_name,
                                                                figname="rotor monitor - L")
        push!(monitors, rotor_monitor_L)

    end

    for (si, rotor) in enumerate(vehicle.rotor_systems[2])

        monitor_name = "rotor_R$(si)"
        rotor_vec = [rotor, ]
        rotor_monitor_R = uns.generate_monitor_rotors(rotors, J, rho, RPM, nsteps, AOA, pitch, tilt, yaw;
                                                                t_scale=RPM/60,
                                                                t_lbl="Revolutions",
                                                                out_figs=figs,
                                                                out_figaxs=figaxs,
                                                                save_path=save_path,
                                                                run_name=run_name,
                                                                figname="rotor monitor - R")
        push!(monitors, rotor_monitor_R)

    end

end

# Generate monitor of flow enstrophy (numerical stability)
monitor_enstrophy = uns.generate_monitor_enstrophy(;
                                            save_path=save_path,
                                            run_name=run_name,
                                            figname="enstrophy monitor"
                                            )

# Generate monitor of SFS model coefficient Cd
monitor_Cd = uns.generate_monitor_Cd(;
                                            save_path=save_path,
                                            run_name=run_name,
                                            figname="Cd monitor"
                                            )
# Concatenate monitors
monitors = uns.concatenate(monitor_enstrophy, monitor_Cd, monitors)


# ------------- 5) RUN SIMULATION ----------------------------------------------
println("Running simulation...")

# Concatenate monitors and wake treatment procedure into one runtime function
runtime_function = uns.concatenate(monitors, wake_treatment_supress)

# Run simulation
uns.run_simulation(simulation, nsteps;
                    # ----- SIMULATION OPTIONS -------------
                    Vinf=Vinf,
                    rho=rho, mu=mu, sound_spd=speedofsound,
                    # ----- SOLVERS OPTIONS ----------------
                    p_per_step=p_per_step,
                    max_particles=max_particles,
                    vpm_integration=vpm_integration,
                    vpm_viscous=vpm_viscous,
                    vpm_SFS=vpm_SFS,
                    sigma_vlm_surf=sigma_rotor_surf,
                    sigma_rotor_surf=sigma_rotor_surf,
                    sigma_vpm_overwrite=sigma_vpm_overwrite,
                    sigmafactor_vpmonvlm=sigmafactor_vpmonvlm,
                    vlm_rlx=vlm_rlx,
                    hubtiploss_correction=hubtiploss_correction,
                    shed_starting=shed_starting,
                    shed_unsteady=shed_unsteady,
                    unsteady_shedcrit=unsteady_shedcrit,
                    omit_shedding=omit_shedding,
                    extra_runtime_function=runtime_function,
                    # ----- RESTART OPTIONS -----------------
                    restart_vpmfile=restart_file,
                    # ----- OUTPUT OPTIONS ------------------
                    save_path=save_path,
                    run_name=run_name,
                    save_wopwopin=true,  # <--- Generates input files for PSU-WOPWOP noise analysis
                    );




# ----------------- 6) VISUALIZATION -------------------------------------------
if paraview
    println("Calling Paraview...")

    # Files to open in Paraview
    files = joinpath(save_path, run_name*"_pfield...xmf;")
    for bi in 1:B
        global files
        files *= run_name*"_Rotor_Blade$(bi)_loft...vtk;"
        files *= run_name*"_Rotor_Blade$(bi)_vlm...vtk;"
    end

    # Call Paraview
    run(`paraview --data=$(files)`)

end

