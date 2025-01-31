import os                                   # type: ignore    
import numpy as np                          # type: ignore
import pandas as pd                         # type: ignore
from scipy.fft import fft                   # type: ignore
import torch                                # type: ignore
from torch.utils.data import Dataset        # type: ignore

# Condition function to filter subdirectories
def subdir_condition_rotor(subdir_name):
    """
    Condition: Only process subdirectories whose names start with 'propeller-example_dji'.
    Modify this function to apply a specific filtering logic.
    """
    return subdir_name.startswith('eMO_hover')  # Change this condition as needed


class PropellerDataset(Dataset):
    # def __init__(self, root_dir, alpha, J, theta, yaw, tilt, subdir_condition=None):
    def __init__(self, root_dir, rotor_notation=None, subdir_condition=None):
        """
        Args:
            root_dir (string): Root directory with subdirectories containing CSV files.
            alpha (float): Angle of attack.
            J (float): Advance ratio.
            theta (float): Pitch.
            yaw (float): Yaw.
            tilt (float): Tilt.
            subdir_condition (callable, optional): A function or condition to filter subdirectories by name.
        """
        self.root_dir = root_dir
        self.data = []
        self.targets = []
        self.time_data = []  # Store time data separately
        self.omega_data = []  # Store omega (RPM) separately
        self.ct_data = []
        self.cq_data = []
        self.fft_ct_r = []
        self.fft_ct_i = []
        self.fft_cq_r = []
        self.fft_cq_i = []
        self.subdir_condition = subdir_condition

        self.rotor_notation = rotor_notation

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
                if file.endswith(self.rotor_notation+"_convergence.csv"):
                    # Load the CSV file
                    csv_path = os.path.join(subdir, file)
                    df = pd.read_csv(csv_path)
                    df = df[df['T'] < 0.28]
                    
                    # Extract necessary columns for input features
                    time = df['T'].values  # Time
                    omega = df['RPM_1'].values  # RPM
                    J = df['J'].values
                    AOA = df['AOA'].values
                    v_inf = df['v_inf'].values
                    pitch = df['Pitch (blade)'].values
                    tilt = df['Tilt'].values
                    yaw = df['Yaw'].values
                    
                    ref_angle = df['ref age (deg)'].values

                    # Store time and omega in separate lists for easy access
                    self.time_data.append(time)
                    self.omega_data.append(omega)
                    
                    # Extract CT and CQ as output variables
                    ct = df['CT_1'].values  # Thrust coefficient (CT)
                    cq = df['CQ_1'].values  # Torque coefficient (CQ)

                    fft_ct = fft(ct)
                    fft_ct_real = np.real(fft_ct)
                    fft_ct_imag = np.imag(fft_ct)

                    fft_cq = fft(cq)
                    fft_cq_real = np.real(fft_cq)
                    fft_cq_imag = np.imag(fft_cq)

                    
                    # Store ct and cq in separate lists for easy access
                    self.ct_data.append(ct)
                    self.cq_data.append(cq)
                    self.fft_ct_r.append(fft_ct_real)
                    self.fft_ct_i.append(fft_ct_imag)
                    self.fft_cq_r.append(fft_cq_real)
                    self.fft_cq_i.append(fft_cq_imag)


                    # For each simulation, the input sequence is structured as (n_timesteps, n_features)
                    sequence_inputs = []
                    sequence_outputs = []
                    for i in range(len(time)):
                        # Each time step has time, omega, and predefined variables: alpha, J, theta, yaw, tilt
                        input_data = [
                            # time[i], ref_angle[i], omega[i], (100*np.sin(omega[i]*time[i])), 
                            # (100*np.cos(omega[i]*time[i])), AOA[i], J[i], pitch[i], tilt[i], yaw[i]
                            time[i], omega[i], AOA[i], v_inf[i], (100*np.sin(omega[i]*time[i])), 
                            (100*np.cos(omega[i]*time[i])), J[i], pitch[i], tilt[i], yaw[i]
                            # omega[i], AOA[i], J[i], pitch[i], tilt[i], yaw[i]
                        ]
                        # output_data = [ct[i], cq[i], fft_ct_imag[i], fft_cq_imag[i]]
                        output_data = [ct[i], cq[i]]
                        # output_data = [fft_ct_real[i], fft_ct_imag[i], fft_cq_real[i], fft_cq_imag[i]]
                        
                        
                        sequence_inputs.append(input_data)
                        sequence_outputs.append(output_data)

                    sequence_inputs = np.array([sequence_inputs], dtype=float)
                    sequence_outputs = np.array([sequence_outputs], dtype=float)

                    # Append input sequence (n_timesteps, num_features) and output (CT, CQ)
                    self.data.append(sequence_inputs)
                    self.targets.append(sequence_outputs)  # Append the whole CT and CQ sequences


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
            'omega' - rotational velocity [RPM]
            'CT' - Thrust coefficient
            'CQ' - Torque coefficient
            'fft_ct_r' - Thrust ccoefficient FFT - real part
            'fft_ct_i' - Thrust ccoefficient FFT - imag part
            'fft_cq_r' - Torque ccoefficient FFT - real part
            'fft_cq_i' - Torque ccoefficient FFT - imag part

        """
        if variable_name == 'time':
            return self.time_data  # Return all time steps for each simulation
        elif variable_name == 'omega':
            return self.omega_data  # Return all omega (RPM) values for each simulation
        elif variable_name == 'CT':
            return self.ct_data 
        elif variable_name == 'CQ':
            return self.cq_data
        elif variable_name == 'fft_ct':
            return self.fft_ct_r, self.fft_ct_i 
        elif variable_name == 'fft_cq':
            return self.fft_cq_r, self.fft_cq_i  
        else:
            raise ValueError(f"Variable {variable_name} not supported.")