"""
    Generates the monitors of Vahana eVTOL simulation
"""
function generate_monitor_eMagicOne(vehicle, rho, RPMref, nsteps, save_path, Vinf;
                                                        add_wings=true,
                                                        wingmonitor_optargs=[])
    #include(joinpath("eMagicOne_vehicle.jl" ))


    airfoil_datapath = joinpath("database", "airfoils")

    # Collect all monitors here
    monitors = []

    # -------------------- WING MONITORS ---------------------------------------
    # Reference parameters for calculating coefficients
    # NOTE: make b, ar, and qinf equals to 1.0 to obtain dimensional force
    b_ref, ar_ref = 1.0, 1.0
    qinf = 1.0
    Jref = 1.0

    # Force axis labels
    CL_lbl = "Lift (N)"
    CD_lbl = "Drag (N)"

    # Directions of force components
    L_dir = [0, 0, 1]
    D_dir = [1, 0, 0]

    # Generate function that computes wing aerodynamic forces
    calc_aerodynamicforce_fun_wing = uns.generate_calc_aerodynamicforce(;
                                        add_parasiticdrag=true,
                                        add_skinfriction=true,
                                        airfoilpolar="xf-naca4412-il-1000000.csv",
					parasiticdrag_args = (read_path=airfoil_datapath,
                                                calc_cd_from_cl=false,
                                                add_skinfriction=true,
                                                Mach=nothing),
                                        )
    calc_aerodynamicforce_fun_canard = uns.generate_calc_aerodynamicforce(;
                                        add_parasiticdrag=true,
                                        add_skinfriction=true,
                                        airfoilpolar="xf-naca6412-il-1000000.csv",
					parasiticdrag_args = (read_path=airfoil_datapath,
                                                calc_cd_from_cl=false,
                                                add_skinfriction=true,
                                                Mach=nothing),
                                        )

    if add_wings

        # Main wing monitor
        monitor_name = "wing_main"
        mainwing_system = vlm.get_wing(vehicle.vlm_system, "Wing")
        mainwing_monitor = uns.generate_monitor_wing(mainwing_system, Vinf, b_ref, ar_ref,
                                                        rho, qinf, nsteps;
                                                        calc_aerodynamicforce_fun=calc_aerodynamicforce_fun_wing,
                                                        save_path=save_path,
                                                        run_name=monitor_name,
                                                        figname=monitor_name,
                                                        CL_lbl=CL_lbl,
                                                        CD_lbl=CD_lbl,
                                                        L_dir=L_dir,
                                                        D_dir=D_dir,
                                                        wingmonitor_optargs...)

        # Tandem wing monitor
        monitor_name = "Canarad"
        tandemwing_system = vlm.get_wing(vehicle.vlm_system, "Canard")
        tandemwing_monitor = uns.generate_monitor_wing(tandemwing_system, Vinf, b_ref, ar_ref,
                                                        rho, qinf, nsteps;
                                                        calc_aerodynamicforce_fun=calc_aerodynamicforce_fun_canard,
                                                        save_path=save_path,
                                                        run_name=monitor_name,
                                                        figname=monitor_name,
                                                        CL_lbl=CL_lbl,
                                                        CD_lbl=CD_lbl,
                                                        L_dir=L_dir,
                                                        D_dir=D_dir,
                                                        wingmonitor_optargs...)

        push!(monitors, mainwing_monitor)
        push!(monitors, tandemwing_monitor)
    end




    # -------------------- ROTOR MONITORS --------------------------------------
    for (si, rotors) in enumerate(vehicle.rotor_systems)

        monitor_name = "rotorsys$(si)"

        rotors_monitor = uns.generate_monitor_rotors(rotors, Jref, rho, RPMref,
                                                        nsteps;
                                                        save_path=save_path,
                                                        run_name=monitor_name,
                                                        figname=monitor_name,
                                                        save_init_plots=false)
        push!(monitors, rotors_monitor)
    end

    for (si, rotor) in enumerate(vehicle.rotor_systems[1])
 
	monitor_name = "rotor_L$(si)"
	rotor_vec = [rotor, ]
	rotor_monitor_L = uns.generate_monitor_rotors(rotor_vec, Jref, rho, RPMref,
	                                                nsteps;
	                                                save_path=save_path,
	                                                run_name=monitor_name,
	                                                figname=monitor_name,
	                                                save_init_plots=false)
	push!(monitors, rotor_monitor_L)

    end

    for (si, rotor) in enumerate(vehicle.rotor_systems[2])
 
	monitor_name = "rotor_R$(si)"
	rotor_vec = [rotor, ]
	rotor_monitor_R = uns.generate_monitor_rotors(rotor_vec, Jref, rho, RPMref,
	                                                nsteps;
	                                                save_path=save_path,
	                                                run_name=monitor_name,
	                                                figname=monitor_name,
	                                                save_init_plots=false)
	push!(monitors, rotor_monitor_R)

    end


    # -------------------- OTHER MONITORS --------------------------------------

    # State-variable monitor
    statevariable_monitor = uns.generate_monitor_statevariables(; save_path=save_path)

    # Global enstrophy monitor (numerical stability)
    monitor_enstrophy = uns.generate_monitor_enstrophy(; save_path=save_path)

    # Monitor of SFS model coefficient Cd
    monitor_Cd = uns.generate_monitor_Cd(; save_path=save_path)


    # -------------------- CONCATENATE MONITORS --------------------------------
    return uns.concatenate(statevariable_monitor, monitor_enstrophy, monitor_Cd, monitors...)
end