import os                                   # type: ignore    
import numpy as np                          # type: ignore
import pandas as pd                         # type: ignore
from scipy.fft import fft                   # type: ignore
from scipy.interpolate import CubicSpline   # type: ignore
import torch                                # type: ignore
from torch.utils.data import Dataset        # type: ignore
import matplotlib.pyplot as plt             # type: ignore


from utility_functions import downsample_to_35, organize_data



# Condition function to filter subdirectories
def subdir_condition_wing(subdir_name):
    """
    Condition: Only process subdirectories whose names start with 'propeller-example_dji'.
    Modify this function to apply a specific filtering logic.
    """
    return subdir_name.startswith('eMO_hover')  # Change this condition as needed


class WingDataset(Dataset):
    def __init__(self, root_dir, af_model_ESCNN_Cl=None, af_model_ESCNN_Cd=None, airfoil_cl=None, airfoil_cd=None, device=None, wing_name=None,
                subdir_condition=None ,af_data_path='/mnt/e/Course_Materials/ROM/wing_model/wing_section/small_database_testing_csv'):
        """
        Args:
            root_dir (string): Root directory with subdirectories containing CSV files.
            af_data_path (string): Airfoil polar database directory contaiing airfoil polars of used airfoils.            
            subdir_condition (callable, optional): A function or condition to filter subdirectories by name.
        """
        self.root_dir = root_dir
        self.data = []
        self.targets = []
        self.time_data = []  # Store time data separately
        self.AOA_data = []
        self.v_inf_data = []
        self.cl_data = []
        self.cd_data = []

        self.cl_data = []
        self.cd_data = []
        self.fft_cl_r = []
        self.fft_cl_i = []
        self.fft_cd_r = []
        self.fft_cd_i = []
        
        self.subdir_condition = subdir_condition
        self.af_data_path = af_data_path
        self.wing_name = wing_name

        self.device = device
        self.af_model_ESCNN_Cl = af_model_ESCNN_Cl
        self.af_model_ESCNN_Cd = af_model_ESCNN_Cd
        self.airfoil_cl = airfoil_cl
        self.airfoil_cd = airfoil_cd


        # Traverse the root directory to gather data
        self._load_data()

    def _load_data(self):
        """
        Helper function to read CSV files from each subdirectory and extract relevant columns.
        """
        # Iterate through each subdirectory in the root directory
        for subdir, _, files in os.walk(self.root_dir):
            subdir_name = os.path.basename(subdir)
            
            # Apply subdirectory name condition
            if self.subdir_condition and not self.subdir_condition(subdir_name):
                continue

            for file in files:

                if file.endswith(self.wing_name+"_convergence.csv"):
                    # Load the CSV file
                    csv_path = os.path.join(subdir, file)
                    df = pd.read_csv(csv_path)
                    df = df[df['T'] <= 0.28]
                    
                    # Extract necessary columns for input features from wing_convergence.csv
                    time = df['T'].values  # Time
                    
                    v_inf = df['Vinf'].values
                    AOA = df['alpha(eff)'].values
                    theta = df['theta'].values
                    vz = df['vVehicle_z'].values
                                        
                    # Extract Cl and Cq for output variables
                    cl = df['CL'].values  # Lift coefficient (CL)
                    cd = df['CD'].values  # Drag coefficient (CD)

                    fft_cl = fft(cl)
                    fft_cl_real = np.real(fft_cl)
                    fft_cl_imag = np.imag(fft_cl)

                    fft_cd = fft(cd)
                    fft_cd_real = np.real(fft_cd)
                    fft_cd_imag = np.imag(fft_cd)

                    # Add the geometry parameters
                    geom_data_path = os.path.join(subdir, self.wing_name+'_geometry_data.csv')         # Read from wing_geometry_data.csv file - b, ar, tr, sweep, dihedral
                    geom_df = pd.read_csv(geom_data_path, nrows=1)

                    ones_empty = np.ones_like(time)                                         # Empty array of ones
                    
                    b_data = geom_df["b_wing"].values
                    b_data_array = geom_df["b_wing"].values * ones_empty                    # Making this an array simplifies things later
                    
                    ar_data = geom_df["ar_wing"].values
                    ar_data_array = geom_df["ar_wing"].values * ones_empty
                    
                    tr_data = geom_df["tr"].values
                    tr_data_array = geom_df["tr"].values * ones_empty

                    sweep_data = geom_df["lambda"].values
                    sweep_data_array = geom_df["lambda"].values * ones_empty
                    
                    gamma_data = geom_df["gamma"].values
                    gamma_data_array = geom_df["gamma"].values * ones_empty
                    
                    # Store the following in separate lists for easy access
                    self.time_data.append(time)
                    self.AOA_data.append(AOA)
                    self.v_inf_data.append(v_inf)
                    
                    self.cl_data.append(cl)
                    self.cd_data.append(cd)

                    self.fft_cl_r.append(fft_cl_real)
                    self.fft_cl_i.append(fft_cl_imag)
                    self.fft_cd_r.append(fft_cd_real)
                    self.fft_cd_i.append(fft_cd_imag)

                    # Extract the airfoil details
                    af_name = geom_df["airfoil "].values
                    extension = '.csv'
                    af_name = str(af_name[0]+extension)

                    af_name_new = af_name.split('-')
                    af_coordinate = str(af_name_new[1]+'_coordinates.dat')          # Name of airfoil coordinates. Adjust it according to the names being used
                    # print(af_coordinate)
                    
                    af_polar_data = pd.read_csv(os.path.join(self.af_data_path, af_name), skiprows=10)
                    # af_polar_data = af_polar_data[[(af_polar_data["Alpha"] >= -2) & (af_polar_data["Alpha"] <= 12)]]
                    # af_coordinate_data = pd.read_csv(os.path.join(self.af_data_path, af_coordinate), delim_whitespace=True)  # or use delimiter=','
                    af_coordinate_data = np.loadtxt(os.path.join(self.af_data_path, af_coordinate))

                    af_coordinate_x = af_coordinate_data[:,0]
                    af_coordinate_y = af_coordinate_data[:,1] 

                    AOA_af_polar = af_polar_data["Alpha"].values
                    cl_af_polar = af_polar_data["Cl"].values
                    cd_af_polar = af_polar_data["Cd"].values

                    # Fit a polynomial to Cl and Cd data
                    # Create cubic spline interpolation
                    spline_cl_airfoil = CubicSpline(AOA_af_polar, cl_af_polar)
                    spline_cd_airfoil = CubicSpline(AOA_af_polar, cd_af_polar)

                    cl_poly_coeff = np.polyfit(AOA_af_polar, cl_af_polar, deg=6)
                    cd_poly_coeff = np.polyfit(AOA_af_polar, cd_af_polar, deg=6)

                    cl_airfoil_calc = spline_cl_airfoil(AOA)
                    cd_airfoil_calc = spline_cd_airfoil(AOA)

                    #-------------------------------------------------------------------------------------------------------------------
                    # Note - To determine the airfoil aerodynamic coefficients, there are two approaaches.
                    # 1. Use the polar file from the database to fit a polynomial function and then estimate the Cl/Cd for any given AOA.
                    # 2. USe the pre-trained Neural Network to predict the Cl/Cd.  
                    #-------------------------------------------------------------------------------------------------------------------
                    
                    # Using approach 2 - Using RBF NN to predict the Cl and Cd
                    #-------------------------------------------------------------------------------------------------------------------
                
                    degree = 3
                    
                    input_sequence_cl_NN = []
                    input_sequence_cd_NN = []

                    # Prepare the data to input to ESCNN Model
                    x_escnn = np.array(downsample_to_35(af_coordinate_x)).reshape(1, -1)
                    y_escnn = np.array(downsample_to_35(af_coordinate_y)).reshape(1, -1)
                    aoa_escnn = np.array(AOA_af_polar).reshape(1, -1)

                    elements_ESCNN = organize_data(x_escnn, y_escnn, aoa_escnn)

                    if elements_ESCNN.shape != (0,):
                        input_escnn = elements_ESCNN

                        input_escnn = torch.tensor(input_escnn, dtype=torch.float32).to(self.device)

                        # Evaluate the model on test dataset
                        with torch.no_grad():
                            Cl_escnn_pred = self.af_model_ESCNN_Cl.forward(input_escnn)
                            Cd_escnn_pred = self.af_model_ESCNN_Cd.forward(input_escnn)
                            
                        Cl_escnn_pred = Cl_escnn_pred.cpu().detach().numpy()  # Convert tensor to numpy array
                        Cl_escnn_pred = Cl_escnn_pred.squeeze(1)

                        Cd_escnn_pred = Cd_escnn_pred.cpu().detach().numpy()  # Convert tensor to numpy array
                        Cd_escnn_pred = Cd_escnn_pred.squeeze(1)

                        plt_af_polar_comparison = False     # If needed for debugging
                        if plt_af_polar_comparison == True:
                            plt.figure()
                            # plt.plot(aoa_test[j], Cl_escnn_pred)
                            # # plt.plot(alphas_t[0], Cl_eval_org_scale)
                            # plt.plot(aoa_test[j], cl_test[j])

                            plt.legend(['NN Model - ESCNN', 'UIUC database'])
                            # plt.title(r'$C_d$ vs $\alpha$ for {} airfoil'.format(keyword))
                            plt.xlabel(r'AOA [$\alpha$]')
                            plt.ylabel(r'$C_l$')

                    else:
                        continue

                    af_coordinate_x_i = downsample_to_35(af_coordinate_x)
                    af_coordinate_y_i = downsample_to_35(af_coordinate_y)
                    AOA_af_polar_i = downsample_to_35(AOA_af_polar)
                    cl_af_polar_i = downsample_to_35(Cl_escnn_pred)
                    cd_af_polar_i = downsample_to_35(Cd_escnn_pred)

                    input_sequence_cl_NN = [
                                            af_coordinate_x_i, af_coordinate_y_i, AOA_af_polar_i, cl_af_polar_i
                                        ]
                    
                    input_sequence_cd_NN = [
                                            af_coordinate_x_i, af_coordinate_y_i, AOA_af_polar_i, cd_af_polar_i
                                        ]



                    input_sequence_cl_NN = np.array(input_sequence_cl_NN, dtype=float).reshape(1, -1)
                    input_sequence_cd_NN = np.array(input_sequence_cd_NN, dtype=float).reshape(1, -1)
                    
                    # print("Airfoil Cl NN - Input shape: ", input_sequence_cl_NN.shape)
                    # print("Airfoil Cd NN - Input shape: ", input_sequence_cd_NN.shape)

                    NN_cl_model_ip_data = torch.tensor(input_sequence_cl_NN, dtype=torch.float32).to(self.device) 
                    NN_cd_model_ip_data = torch.tensor(input_sequence_cd_NN, dtype=torch.float32).to(self.device) 

            
                    # Neural network model to predict the airfoil aerodynamic coeficients    
                    with torch.no_grad():  # Disable gradient computation for inference
                        predicted_coefficients_cl = self.airfoil_cl(NN_cl_model_ip_data)
                        predicted_coefficients_cd = self.airfoil_cd(NN_cd_model_ip_data)

                    predicted_af_coefficients_cl = predicted_coefficients_cl.cpu().detach().numpy()
                    predicted_af_coefficients_cd = predicted_coefficients_cd.cpu().detach().numpy()
                    
                    # print(predicted_af_coefficients_cl)

                    # polynomial_cl_new = np.poly1d(predicted_coefficients[0])
                    polynomial_cl_pred = np.poly1d(predicted_af_coefficients_cl[0])
                    polynomial_cd_pred = np.poly1d(predicted_af_coefficients_cd[0])

                    x_new = np.linspace(AOA_af_polar[0], AOA_af_polar[-1], 100)
                    # y_new_cl = polynomial_cl_new(x_new_cl)
                    y_new_cl = polynomial_cl_pred(x_new)
                    y_new_cd = polynomial_cd_pred(x_new)

                    

                    # # plt.figure()
                    plt_af_polar_NN = False
                    if plt_af_polar_NN == True:
                        plt.figure(figsize=(15, 5))
                        # plt.title(af_name)
                        
                        plt.subplot(1,2,1)
                        plt.plot(AOA_af_polar, cl_af_polar, color='red', label='UIUC Database')
                        plt.plot(x_new, y_new_cl, label=f'Polynomial Degree {degree} - NN')
                        plt.xlabel(r'$\alpha$')
                        plt.ylabel(r'$C_l$')
                        plt.title(af_name)
                        plt.legend()


                        plt.subplot(1,2,2)
                        plt.plot(AOA_af_polar, cd_af_polar, color='red', label='UIUC Database')
                        plt.plot(x_new, y_new_cd, label=f'Polynomial Degree {degree} - NN')
                        plt.xlabel(r'$\alpha$')
                        plt.ylabel(r'$C_d$')
                        plt.title(af_name)
                        plt.legend()
                    
                    
                    cl_af_NN = polynomial_cl_pred(AOA)
                    cd_af_NN = polynomial_cd_pred(AOA)


                    # print(cl_af_NN)
                    # print(cl_airfoil_calc)
                    plt_af_lvl_coeff = False
                    if plt_af_lvl_coeff==True:
                        plt.figure(figsize=(15, 5))
                        
                        
                        plt.subplot(1,2,1)
                        plt.plot(time, cl_airfoil_calc, color='red', label='UIUC Database')
                        plt.plot(time, cl_af_NN, label=f'Polynomial Degree {degree} - NN')
                        plt.xlabel(r'$\alpha$')
                        plt.ylabel(r'$C_l$')
                        plt.title(af_name)
                        plt.legend()

                        plt.subplot(1,2,2)
                        plt.plot(time, cd_airfoil_calc, color='red', label='UIUC Database')
                        plt.plot(time, cd_af_NN, label=f'Polynomial Degree {degree} - NN')
                        plt.xlabel(r'$\alpha$')
                        plt.ylabel(r'$C_d$')
                        plt.title(af_name)
                        plt.legend()
                    
                    #-------------------------------------------------------------------------------------------------------------------
                    # Calculate the Cl and Cd using the emperical equations using airfoil aerodynamic coefficients
                    # Calculate Wing Aerodynamic Coefficients

                    vel_comp_factor_cl = 0.01                # Assumed from trial & error. Depends on flow velocity and wing geometry
                    vel_comp_factor_cd = 0.00005              # Assumed from trial & error. Depends on flow velocity and wing geometry
                    e = 0.85                                # Oswald factor

                    # cl_wing_calc = cl_airfoil_calc / (1 + cl_airfoil_calc / (np.pi * ar_data * e)) * vel_comp_factor_cl

                    # cd_induced = (cl_wing_calc ** 2) / (np.pi * ar_data * e) * vel_comp_factor_cd

                    cl_wing_calc = cl_af_NN / (1 + cl_af_NN / (np.pi * ar_data * e)) * vel_comp_factor_cl

                    cd_induced = (cl_wing_calc ** 2) / (np.pi * ar_data * e) * vel_comp_factor_cd

                    # Step 5: Total drag coefficient for the wing
                    cd_wing_calc = cd_af_NN + cd_induced
                    # cd_wing_calc = cd_airfoil_calc + cd_induced

                    # For each simulation, the input sequence is structured as (n_timesteps, n_features)
                    sequence_inputs = []
                    sequence_outputs = []
                    for i in range(len(time)):
                        # Each time step has time, omega, and predefined variables: alpha, J, theta, yaw, tilt
                        input_data = [
                            time[i], AOA[i], v_inf[i], b_data_array[i], ar_data_array[i], tr_data_array[i], sweep_data_array[i], 
                            gamma_data_array[i], cl_wing_calc[i], cd_wing_calc[i]
                        ]
                        # output_data = [fft_cl_real[i], fft_cl_imag[i], fft_cd_real[i], fft_cd_imag[i]]
                        output_data = [cl[i], cd[i]]
                        
                        
                        sequence_inputs.append(input_data)
                        sequence_outputs.append(output_data)

                    sequence_inputs = np.array([sequence_inputs], dtype=float)
                    sequence_outputs = np.array([sequence_outputs], dtype=float)

                    # Append input sequence (n_timesteps, num_features) and output (Cl, Cd)
                    self.data.append(sequence_inputs)
                    self.targets.append(sequence_outputs)  # Append the whole Cl and Cd sequences


    def __len__(self):
        """
        Returns the total number of sequences in the dataset.
        """
        return len(self.data)

    def __getitem__(self, idx):
        """
        Returns a single sequence and its targets.
        """
        inputs = self.data[idx]  # Input sequence: (n_timesteps, n_features)
        targets = self.targets[idx]  # Output: (n_timesteps, 2)
        # return inputs, targets
        return torch.tensor(inputs, dtype=torch.float32), torch.tensor(targets, dtype=torch.float32)

    def get_variable(self, variable_name):
        """
        Returns a list of arrays for the specified variable.
        Args:
            'time' - timesteps
            
        """
        if variable_name == 'time':
            return self.time_data  # Return all time steps for each simulation
        elif variable_name == 'CL':
            return self.cl_data  # Return all omega (RPM) values for each simulation
        elif variable_name == 'CD':
            return self.cd_data
        elif variable_name == 'AOA':
            return self.AOA_data
        elif variable_name == 'Vinf':
            return self.v_inf_data  
        elif variable_name == 'fft_cl':
            return self.fft_cl_r, self.fft_cl_i 
        elif variable_name == 'fft_cd':
            return self.fft_cd_r, self.fft_cd_i 
        else:
            raise ValueError(f"Variable {variable_name} not supported.")
        

print("create_wing_dataset module loaded")