#=##############################################################################
# DESCRIPTION 
    1) Vehicle definition
=###############################################################################

"""
    Generates the geometry of eMagicOne aircraft
"""
function generate_eMagicOne_vehicle(;
                                    # VEHICLE OPTIONS
                                    xfoil       = true,             # Whether to run XFOIL
                                    n_factor::Int = 1,              # Discretization factor
                                    add_wings   = true,             # Whether to add the wings
                                    add_rotors  = true,             # Whether to add the rotors
                                    VehicleType = uns.VLMVehicle,   # Type of vehicle to generate (uns.QVLMVehicle for quasi-steady solver)
                                    data_path   = uns.def_data_path,# Database path
                                    # data_path_fuselage   = uns.examples_path,# Database path to the examples
                                    # OUTPUT OPTIONS
                                    run_name    = "eMagicOne",
                                    verbose     = true,
                                    v_lvl       = 0
                                    )

    ############################################################################
    # PARAMETERS
    ############################################################################

    if verbose; println("\t"^(v_lvl)*"Defining parameters..."); end;
    
    # MAIN Wing
    b_w               = 7.684                     # (m) span length
    ar_w              = 5.6919                    # Aspect ratio b/c_tip
    tr_w              = 0.8                       # Taper ratio c_tip/c_root
    twist_root_w      = 0.0                       # (deg) twist at root
    twist_tip_w       = 0.0                       # (deg) twist at tip
    lambda_w          = 0.0                      # (deg) sweep
    gamma_w           = 0.0                       # (deg) dihedral
    x_pos_w           = 4.505                     # (m) position of the wing in the x-direction
    y_pos_w           = 0.0                       # (m) position of the wing in the y-direction
    z_pos_w           = -0.2                       # (m) position of the wing in the z-direction
    n_w               = 24*n_factor           # Number of wing elements per side
    r_w               = 2.0
    
    central = false

    # Main-wing winglets
    b_wl            = b_w/6                 # (m) span of winglet from top to bottom
    ar_wl           = 3.0                   # Aspect ratio (b/c_tip)
    tr_wl           = (b_wl/ar_wl)/(b_w/ar_w/tr_w)  # Taper ratio (c_tip/c_root)
    twist_r_wl      = 0.0                   # (deg) twist at root
    twist_t_wl      = 0.0                   # (deg) twist at tip
    lambda_wl       = 60.0                  # (deg) sweep
    gamma_wl        = 0.0                  # (deg) dihedral
    n_wl            = 8*n_factor                   # Number of wing elements per side
    r_wl            = 2.0                   # Geometric expansion of wing elements
    x_pos_wl        = x_pos_w               # (m) position of the winglet in the x-direction
    y_pos_wl        = b_w/2                 # (m) position of the winglet in the y-direction
    z_pos_wl        = z_pos_w               # (m) position of the winglet in the z-direction

    # Canard
    b_c               = 7.684                     # (m) span length
    ar_c              = 15.368                    # Aspect ratio b/c_tip
    tr_c              = 1.0                       # Taper ratio c_tip/c_root
    twist_root_c      = 0.0                       # (deg) twist at root
    twist_tip_c       = 0.0                       # (deg) twist at tip
    lambda_c          = 0.0                      # (deg) sweep
    gamma_c           = 0.0                       # (deg) dihedral
    x_pos_c           = 1.14                     # (m) position of the wing in the x-direction
    y_pos_c           = 0.0                       # (m) position of the wing in the y-direction
    z_pos_c           = 0.35                       # (m) position of the wing in the z-direction
    n_c               = n_w                       # Number of wing elements per side  
    r_c               = r_w                         # Geometric expansion of wing elements

    # Fuselage
    # Path to DegenGeom file 
    # REMEMBER TO CHSNGE THE PATH
    # geom_path = joinpath(uns.examples_path, "aircraft-vsp", "fuselage_model_2_DegenGeom.csv")
    # boom_path = joinpath(uns.examples_path, "aircraft-vsp", "boom_model_DegenGeom.csv")

    geom_path = joinpath("aircraft_geometry", "fuselage_model_2_DegenGeom.csv")
    boom_path = joinpath("aircraft_geometry", "boom_model_DegenGeom.csv")
    rotor_data_path = joinpath("database")
    # Booms
    x_pos_booom          = 0.0                    # (m) position of the wing in the x-direction
    y_pos_boom           = 3.842                       # (m) position of the wing in the y-direction
    z_pos_boom           = 0.35                       # (m) position of the wing in the z-direction


    # Rotor
    # Rotor geometry
    rotor_file      = "H26F_mod.csv"             # Rotor geometry
    # rotor_file      = "DJI9443.csv"             # Rotor geometry
    # data_path       = uns.def_data_path         # Path to rotor database
    data_path       = rotor_data_path         # Path to rotor database
    pitch           = 0.0                       # (deg) collective pitch of blades
    CWs             = [true, false]                     # Clock-wise rotation
    xfoil           = xfoil                      # Whether to run XFOIL
    ncrit           = 9                         # Turbulence criterion for XFOIL
    n_rotor         = 7*n_factor            # Number of blade elements per blade
    r_rotor         = 1/20                  # Geometric expansion between blade elements
    
    # Read radius of this rotor and number of blades
    R, B            = uns.read_rotor(rotor_file; data_path=data_path)[[1,3]]

    # Operating conditions
    RPM             = 1600                      # RPM
    J               = 0.15                   # Advance ratio Vinf/(nD)
    AOA             = 1.0                         # (deg) Angle of attack (incidence angle)

    rho             = 1.225                     # (kg/m^3) air density
    mu              = 1.81e-5                   # (kg/ms) air dynamic viscosity
    speedofsound    = 342.35                    # (m/s) speed of sound

    magVinf         = J*RPM/60*(2*R)
    Vinf(X, t)      = magVinf*[cosd(AOA), sind(AOA), 0] # (m/s) freestream velocity vector
    ReD             = 2*pi*RPM/60*R * rho/mu * 2*R      # Diameter-based Reynolds number
    Matip           = 2*pi*RPM/60*R / speedofsound      # Tip Mach number

    # Rotor positions
    # L1_rotor
    x_pos_rL1       = 0.36
    y_pos_rL1       = -3.842
    z_pos_rL1       = 0.5 

    # L2_rotor
    x_pos_rL2       = 2.28
    y_pos_rL2       = -3.842
    z_pos_rL2       = 0.5 

    # L3_rotor
    x_pos_rL3       = 4.2
    y_pos_rL3       = -3.842
    z_pos_rL3       = 0.5 

    # L4_rotor
    x_pos_rL4       = 6.54
    y_pos_rL4       = -3.842
    z_pos_rL4       = 0.5 

    # R1_rotor
    x_pos_rR1       = 0.36
    y_pos_rR1       = 3.842
    z_pos_rR1       = 0.5 

    # R2_rotor
    x_pos_rR2       = 2.28
    y_pos_rR2       = 3.842
    z_pos_rR2       = 0.5

    # R3_rotor
    x_pos_rR3       = 4.2
    y_pos_rR3       = 3.842
    z_pos_rR3       = 0.5 

    # R4_rotor
    x_pos_rR4       = 6.54
    y_pos_rR4       = 3.842
    z_pos_rR4       = 0.5 

    ############################################################################
    # GENERATE COMPONENTS
    ############################################################################

    if verbose; println("\t"^(v_lvl)*"Generating components..."); end;

    # ------------ MAIN WING ---------------------------------------------
    # Generate wing
    if verbose; println("\t"^(v_lvl+1)*"Generating main wing assembly..."); end;

    wing = vlm.simpleWing(b_w, ar_w, tr_w, twist_root_w, lambda_w, gamma_w;
                            twist_tip=twist_tip_w, n=n_w, r=r_w, central=central);

    # Translate the wing into the position
    wing_pos = [x_pos_w, y_pos_w, z_pos_w]
    wing_axis = uns.gt.rotation_matrix(0.0, 0.0, 0.0)
    vlm.setcoordsystem(wing, wing_pos, wing_axis)

    # ------------ WINGLETS--- ---------------------------------------------
    # Generate Winglets
    winglet_R = vlm.simpleWing(b_wl, ar_wl, tr_wl, twist_r_wl, lambda_wl, gamma_wl;
    twist_tip=twist_t_wl, n=n_wl, r=r_wl)

    # Translate the winglet into the position
    winglet_R_pos = [x_pos_wl, y_pos_wl, z_pos_wl]
    winglet_R_axis = uns.gt.rotation_matrix(0.0, 0.0, 90.0)
    vlm.setcoordsystem(winglet_R, winglet_R_pos, winglet_R_axis)

    winglet_L = vlm.simpleWing(b_wl, ar_wl, tr_wl, twist_r_wl, lambda_wl, gamma_wl;
    twist_tip=twist_t_wl, n=n_wl, r=r_wl)

    # Translate the winglet into the position
    winglet_L_pos = [x_pos_wl, -y_pos_wl, z_pos_wl]
    winglet_L_axis = uns.gt.rotation_matrix(0.0, 0.0, -90.0)
    vlm.setcoordsystem(winglet_L, winglet_L_pos, winglet_L_axis)

    main_wing_fixed = vlm.WingSystem()
    main_wing_moving = vlm.WingSystem()
    main_wing = vlm.WingSystem()
    # Assemble the winglets to main wing
    vlm.addwing(main_wing_fixed, "Wing", wing)
    vlm.addwing(main_wing_moving, "Winglet_L", winglet_L)
    vlm.addwing(main_wing_moving, "Winglet_R", winglet_R)

    vlm.addwing(main_wing, "Fixed", main_wing_fixed)
    vlm.addwing(main_wing, "Moving", main_wing_moving)

    # ------------CANARD------ ---------------------------------------------
    if verbose; println("\t"^(v_lvl+1)*"Generating canard..."); end;

    # Generate canard
    canard = vlm.simpleWing(b_c, ar_c, tr_c, twist_root_c, lambda_c, gamma_c;
    twist_tip=twist_tip_c, n=n_c, r=r_c, central=central);

    # Translate the canard into the position
    canard_pos = [x_pos_c, y_pos_c, z_pos_c]
    canard_axis = uns.gt.rotation_matrix(0.0, 0.0, 0.0)
    vlm.setcoordsystem(canard, canard_pos, canard_axis)

    # ------------FUSELAGE---------------------------------------------------
    if verbose; println("\t"^(v_lvl+1)*"Generating Fuselage..."); end;
    # Genarate Fuselage
    # Import VSP Components from DegenGeom file
    comp_fus = uns.read_degengeom(geom_path);

    fuselage = uns.import_vsp(comp_fus[1])
    vertStabl = uns.import_vsp(comp_fus[2])

    comp_boom = uns.read_degengeom(boom_path);
    boom_R = uns.import_vsp(comp_boom[1])
    boom_L = uns.import_vsp(comp_boom[2])

    
    
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

    ############################################################################
    # DEFINE VEHICLE
    ############################################################################
    println("Generating vehicle...")
    # Generate vehicle
    system = vlm.WingSystem()                   # System of all FLOWVLM objects
    vlm.addwing(system, "WingSystem", main_wing)
    vlm.addwing(system, "Canard", canard)
    vlm.addwing(system, "vertStabl", vertStabl)
    # vlm.addwing(system, "Winglet_L", winglet_L)
    # vlm.addwing(system, "Winglet_R", winglet_R)

    tilting_systems = (main_wing_moving, )

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

    fuse_grid = uns.gt.MultiGrid(3)
    uns.gt.addgrid(fuse_grid, "Fuselage", fuselage)
    uns.gt.addgrid(fuse_grid, "Boom_L", boom_L)
    uns.gt.addgrid(fuse_grid, "Boom_R", boom_R)

    grids = [fuse_grid]

    if add_rotors
        rotor_systems = (rotors_left, rotors_right);
    else
        rotor_systems = ();
    end

    vlm_system = vlm.WingSystem();                        # System solved through VLM solver
    vlm.addwing(vlm_system, "Wing", wing)
    vlm.addwing(vlm_system, "Canard", canard)
    vlm.addwing(vlm_system, "vertStabl", vertStabl)
    vlm.addwing(vlm_system, "Winglet_L", winglet_L)
    vlm.addwing(vlm_system, "Winglet_R", winglet_R)

    # vlm.addwing(vlm_system, "MWing", main_wing)

    wake_system = vlm.WingSystem();                       # System that will shed a VPM wake
    vlm.addwing(wake_system, "Wing", wing)
    vlm.addwing(wake_system, "Canard", canard)
    vlm.addwing(wake_system, "vertStabl", vertStabl)
    vlm.addwing(wake_system, "Winglet_L", winglet_L)
    vlm.addwing(wake_system, "Winglet_R", winglet_R)

    # vlm.addwing(wake_system, "MWing", main_wing)

    if add_rotors
        rotors = vcat(rotors_left, rotors_right);
    end
    if add_rotors
        if VehicleType==uns.VLMVehicle
            # vlm.addwing(wake_system, "Rotor_L1", rotor_L1)
            # vlm.addwing(wake_system, "Rotor_L2", rotor_L2)
            # vlm.addwing(wake_system, "Rotor_L3", rotor_L3)
            # vlm.addwing(wake_system, "Rotor_L4", rotor_L4)
            # vlm.addwing(wake_system, "Rotor_R1", rotor_R1)
            # vlm.addwing(wake_system, "Rotor_R2", rotor_R2)
            # vlm.addwing(wake_system, "Rotor_R3", rotor_R3)
            # vlm.addwing(wake_system, "Rotor_R4", rotor_R4)
            for (i, rotor) in enumerate(rotors)
                vlm.addwing(wake_system, "Rotor$i", rotor)
            end
            
        else
            uns.vlm.VLMSolver._mute_warning(true)
        end
    end


    vehicle = uns.VLMVehicle(   system;
                                    tilting_systems = tilting_systems,
                                    rotor_systems=rotor_systems,
                                    vlm_system=vlm_system,
                                    wake_system=wake_system,
                                    grids=grids
                                );

    ############################################################################
    # VISUALIZE VEHICLE
    ############################################################################                        
    # ----------------- EXPORT GEOMETRY --------------------------------------------


    if isdir(save_path); rm(save_path, recursive=true, force=true); end
    mkdir(save_path)

    uns.vlm.setVinf(system, Vinf)
    str = uns.save_vtk(vehicle, run_name; path=save_path)

    # # Open up geometry in ParaView
    # str = joinpath(save_path, str)
    # run(`paraview -data=$(str)`)

    return vehicle
end



