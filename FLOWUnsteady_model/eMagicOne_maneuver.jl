"""
    Generates the eVTOL transition maneuver of eMagicOne aircraft
"""
function generate_maneuver_eMagicOne(; disp_plot=true,    # If true, it will display a plot of the maneuver
                                    add_rotors=true,    # Whether to add rotors contorl inputs to the maneuver
                                    V0=0.0001           # Initial non-dimensional velocity (slightly different than zero to avoid instabilities)
                                    )

    # NOTE: The following parameters are non-dimensional and scaled between 0
    #       and 1. Here, t=0 is the beginning of the maneuver and t=1 is the end.

    # Define end time of each stage
    #  Stage 1: [0,  t1] -> Take off
    #  Stage 2: [t1, t2] -> Transition
    #  Stage 3: [t2, t3] -> Cruise
    #  Stage 4: [t3, t4] -> Transition
    #  Stage 5: [t4, 1 ] -> Landing
    t1, t2, t3, t4 = 0.2, 0.3, 0.5, 0.6

    # Target velocity at each stage (ratio of cruise velocity)
    V1 = 0.001
    V2 = 0.10
    V3 = 1.00
    V4 = V2
    V5 = 0.5*V1

    # Target RPM at each stage (ratio of hover RPM)
    RPM1 = 1.05
    RPM2 = 1.50
    #RPM3 = 0.75
    RPM3 = 1.0
    RPM4 = 1.00
    RPM5 = 0.90

    # r_RPMh_stup = 3.50             # Ratio between stacked and main rotor RPM in hover
    # r_RPMh_tw = 0.75               # Ratio between tandem and main rotor RPM in hover

    # NOTE: -x is in the direction of flight and +z is climb with ground at z=0


    ############################################################################
    # AIRCRAFT VELOCITY
    ############################################################################
    """
        Receives a nondimensional time between 0 and 1, and returns the
        non-dimensional velocity vector of the aircraft at that instant.
    """
    function Vaircraft(t)

        # ------------ TAKE OFF ------------------------------------------------
        if t<t1

            # Weibull acceleration to target climb
            if t<t1/2
                val = t / (t1/2)
                Vz = V0 + V1*(1-exp(-(2*val)^5))

            # Weibull deceleration to hover
            else
                val = 1 - (t-t1/2) / (t1/2)
                Vz = V1*(1-exp(-(2*val)^5))
            end

            return [0, 0, Vz]

        # ------------ TRANSITION ----------------------------------------------
        elseif t<t2

            # Weibull acceleration to cruise
            val = 1.5 * (t-t1) / (t2-t1) + 0.25
            Vx = V3*(1-exp(-(val)^5))

            # Slight drop during transition
            val = (1/0.8) * (t-t1) / (t2-t1)
            Vz = 0.5 * (5 * val^4 * exp(-val^5))
            Vz = 0.2*V3 * Vz

            return [-Vx, 0, -Vz]

        # ------------ CRUISE --------------------------------------------------
        elseif t<t3

            return [-V3, 0, 0]

        # ------------ TRANSITION ----------------------------------------------
        elseif t<t4

            # Weibull deceleration to hover
            val = 1.5 * (1 - (t-t3) / (t4-t3))
            Vx = V3*(1-exp(-(val)^5))

            # Slight climb during transition
            val = 2.0 * (t-t3) / (t4-t3)
            Vz = 0.5 * (5 * val^4 * exp(-val^5))
            Vz = 0.2*V3 * Vz

            return [-Vx, 0, Vz]

        # ------------ LANDING -------------------------------------------------
        else

            # Weibull acceleration to target descend
            if (t-t4)<(1-t4)*0.35
                val = (t-t4) / ((1-t4)*0.35)
                Vz = V5*(1-exp(-(2*val)^5))
            # Weibull deceleration to hover
            else
                val = 1 - ((t-t4) - (1-t4)*0.35) / ((1-t4)*(1-0.35))
                Vz = V0 + V5*(1-exp(-(2*val)^10))
            end

            return [0, 0, -Vz]

        end
    end


    ############################################################################
    # AIRCRAFT ANGLES
    ############################################################################
    """
        Receives a nondimensional time between 0 and 1, and returns the attitude
        of the aircraft. The attitude is a vector indicating the angle of the
        vehicle frame with respect to each global axis.
    """
    function angleaircraft(t)

        # ------------ TAKE OFF ------------------------------------------------
        if t<t1

            return [0, 0, 0]

        # ------------ TRANSITION ----------------------------------------------
        elseif t<t2
            # Slightly pitching down during transition
            # val = 1.5 * (t-t1) / (t2-t1)
            # val3 = 1.5*val
            # val3 = (val3)^(val3 < 1 ? 3 : 1.5)
            # angle_aircraft = -15*(1.5/0.75 * val3^0.5 * exp(-(val3^1.5)))
            # return [0, angle_aircraft, 0]
            return [0, 0, 0]

        # ------------ CRUISE --------------------------------------------------
        elseif t<t3
            return [0, 0, 0]

        # ------------ TRANSITION ----------------------------------------------
        elseif t<t4
            # Weibull deceleration to hover
            # val = 1.5 * (1 - (t-t3) / (t4-t3))
            # # Slight climb during transition
            # val = 2.0 * (t-t3) / (t4-t3)
            # angle_aircraft = 15 * 0.5 * (5 * val^4 * exp(-val^5))
            # return [0, angle_aircraft, 0]
            return [0, 0, 0]

        # ------------ LANDING -------------------------------------------------
        else
            return [0, 0, 0]
        end
    end


    ############################################################################
    # CONTROL INPUT: TILTING SYSTEMS
    ############################################################################
    """
        Receives a nondimensional time between 0 and 1, and returns a vector
        with the angles of the main-wing system (in degrees) relative to the
        aircraft frame.

        This tilting system is made out of the tip-mounted tiltrotors on the
        main wing and winglets.
    """
    function angle_main(t)

        # ------------ TAKE OFF ------------------------------------------------
        if t<t1

            # Tiltrotors pointing up
            return [0, 0, 0]

        # ------------ TRANSITION ----------------------------------------------
        elseif t<t2

            # Tilt forward to cruise position
            # val = 1.5 * (t-t1) / (t2-t1)
            # val2 = val - 0.00
            # angle = 90 - 90*(1-exp(-(val2)^6))
            return [0, 0, 0]

        # ------------ CRUISE --------------------------------------------------
        elseif t<t3

            # Tiltrotors pointing forward in propeller mode
            return [0, 0, 0]

        # ------------ TRANSITION ----------------------------------------------
        elseif t<t4

            # Tilt back to hover position
            # val = 1.5 * (1 - (t-t3) / (t4-t3))
            # angle = 90 - 90*(1-exp(-(val)^5))
            return [0, 0, 0]

        # ------------ LANDING -------------------------------------------------
        else
            # Tiltrotors pointing up
            return [0, 0, 0]
        end
    end

    

    ###########################################################################
    # ROTORS RPM
    ############################################################################
    """
        Receives a nondimensional time between 0 and 1, and returns the RPM of
        main-wing tiltrotors normalized by the hover RPM.
    """
    function RPM_main_left(t)

        # ------------ TAKE OFF ------------------------------------------------
        if t<t1

            # Weibull acceleration to target climb
            if t<t1/2
                val = t / (t1/2)
                RPM_w = RPM1*(1.05-exp(-(5*val)^2))

            # Weibull deceleration to hover
            else
                val = 1 - (t-t1/2) / (t1/2)
                RPM_w = 1.0 + (RPM1-1.0)*(1-exp(-(2*val)^5))
            end

            return RPM_w

        # ------------ TRANSITION ----------------------------------------------
        elseif t<t2

            # Increase RPM to sustain forward flight and hover
            if (t-t1)<(t2-t1)*0.25
                val = (t-t1)/((t2-t1)*0.25)
                RPM_w = 1.0 + (RPM2-1.0)*(1-exp(-(2*val)^5))

            # Decrease RPM as tiltrotors transition to propeller mode
            else
                val = 1 - ((t2-t)/(t2-t1)-0.25) / (1-0.25)
                RPM_w = RPM2 + (RPM3-RPM2)*(1-exp(-(1.1*val)^5))
            end

            return RPM_w

        # ------------ CRUISE --------------------------------------------------
        elseif t<t3

            # Propeller mode
            RPM_w = RPM3
            return RPM_w

        # ------------ TRANSITION ----------------------------------------------
        elseif t<t4

            # Weibull acceleration to hover
            val = 1.5 * (1 - (t-t3) / (t4-t3))
            RPM_w = RPM4 + (RPM3-RPM4)*(1-exp(-(val)^3))
            return RPM_w

        # ------------ LANDING -------------------------------------------------
        else

            # Weibull deceleration to target descend
            if (t-t4)<(1-t4)*0.55
                val = (t-t4) / ((1-t4)*0.55)
                RPM_w = RPM4 + (RPM5-RPM4)*(1-exp(-(3*val)^5))

            # Weibull acceleration to hover
            else
                val = ((t-t4) - (1-t4)*0.55) / ((1-t4)*(1-0.55))
                val = val + 0.60
                RPM_w = RPM5*(3.0*(val)^(3.0-1)*exp(-(val)^5.0))/1.0
            end

            return RPM_w

        end
    end


    """
        Receives a nondimensional time between 0 and 1, and returns the RPM of
        main-wing tiltrotors normalized by the hover RPM.
    """
    function RPM_main_right(t)

        # ------------ TAKE OFF ------------------------------------------------
        if t<t1

            # Weibull acceleration to target climb
            if t<t1/2
                val = t / (t1/2)
                RPM_w = RPM1*(1.05-exp(-(5*val)^2))

            # Weibull deceleration to hover
            else
                val = 1 - (t-t1/2) / (t1/2)
                RPM_w = 1.0 + (RPM1-1.0)*(1-exp(-(2*val)^5))
            end

            return RPM_w

        # ------------ TRANSITION ----------------------------------------------
        elseif t<t2

            # Increase RPM to sustain forward flight and hover
            if (t-t1)<(t2-t1)*0.25
                val = (t-t1)/((t2-t1)*0.25)
                RPM_w = 1.0 + (RPM2-1.0)*(1-exp(-(2*val)^5))

            # Decrease RPM as tiltrotors transition to propeller mode
            else
                val = 1 - ((t2-t)/(t2-t1)-0.25) / (1-0.25)
                RPM_w = RPM2 + (RPM3-RPM2)*(1-exp(-(1.1*val)^5))
            end

            return RPM_w

        # ------------ CRUISE --------------------------------------------------
        elseif t<t3

            # Propeller mode
            RPM_w = RPM3
            return RPM_w

        # ------------ TRANSITION ----------------------------------------------
        elseif t<t4

            # Weibull acceleration to hover
            val = 1.5 * (1 - (t-t3) / (t4-t3))
            RPM_w = RPM4 + (RPM3-RPM4)*(1-exp(-(val)^3))
            return RPM_w

        # ------------ LANDING -------------------------------------------------
        else

            # Weibull deceleration to target descend
            if (t-t4)<(1-t4)*0.55
                val = (t-t4) / ((1-t4)*0.55)
                RPM_w = RPM4 + (RPM5-RPM4)*(1-exp(-(3*val)^5))

            # Weibull acceleration to hover
            else
                val = ((t-t4) - (1-t4)*0.55) / ((1-t4)*(1-0.55))
                val = val + 0.60
                RPM_w = RPM5*(3.0*(val)^(3.0-1)*exp(-(val)^5.0))/1.0
            end

            return RPM_w

        end
    end

    
    ############################################################################
    # MANEUVER OBJECT
    ############################################################################
    # Angle function of each tilting system
    angle = (angle_main, )

    # RPM function of each rotor system
    if add_rotors
        RPM = (RPM_main_left, RPM_main_right )
    else
        RPM = ()
    end

    # Aircraft velocity and angles (attitude)
    Vvehicle = Vaircraft
    anglevehicle = angleaircraft

    # Define non-dimensional maneuver
    maneuver = uns.KinematicManeuver(angle, RPM, Vvehicle, anglevehicle)

    # Plot maneuver for verification purposes
    if disp_plot
        uns.plot_maneuver(maneuver; tstages=[t1, t2, t3, t4])
    end

    return Vvehicle, anglevehicle, maneuver
end