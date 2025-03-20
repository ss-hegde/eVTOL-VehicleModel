# %% [markdown]
# # eVTOL Vehicle Model Data Processing
# 
# This Jupyter Notebook is designed to process and analyze data from eVTOL (electric Vertical Take-Off and Landing) vehicle simulations. The notebook performs the following tasks:
# 
# 1. **Data Import and Initialization**:
#     - Imports necessary libraries such as `os`, `pandas`, and `numpy`.
#     - Sets the root directory for the data files.
#     - Initializes a counter to keep track of the number of folders processed.
# 
# 2. **Data Processing**:
#     - Iterates through each folder in the root directory.
#     - Extracts the velocity information from the folder name.
#     - Processes specific CSV files (`rotor_L1_convergence.csv`, `rotor_L2_convergence.csv`, etc.) within each folder.
#     - Adds a new column with the velocity value to each DataFrame.
#     - Saves the modified DataFrame, either overwriting the original file or saving it as a new file.
# 
# 3. **Output**:
#     - Prints the folder and file being processed.
#     - Displays the velocity value extracted from the folder name.
#     - Shows the number of folders processed.
# 
# This notebook helps in organizing and preparing the simulation data for further analysis and visualization.

# %%
import os
import pandas as pd
import numpy as np

root = "/mnt/e/eVTOL_model/eVTOL-VehicleModel/FLOWUnsteady_simulations/aircraft_data/training_data"

counter = 1
for folder in os.listdir(root):
    split_folder = folder.split("_")
    vel = split_folder[2]

    vel = vel.split("v")
    vel = float(vel[1])

    
    print("Processing the folder - ", folder)
    print("Velocity - ", vel)

    for file in os.listdir(os.path.join(root, folder)):
        if file.startswith("rotor_L1_convergence.csv"):
            print("Processing the file - ", file)   
            file_path = os.path.join(root, folder, file)
            df_L1 = pd.read_csv(file_path)

            # Adding a new column with velocity value (example)
            df_L1["v_inf"] = vel * np.ones_like(df_L1["T"])  # Replace with any logic you need

            # Save the modified DataFrame (overwrite or save as a new file)
            df_L1.to_csv(file_path, index=False)  # Overwrite the original file
            # df_L1.to_csv(file_path.replace(".csv", "_modified.csv"), index=False)  # Save as new file

        if file.startswith("rotor_L2_convergence.csv"):
            print("Processing the file - ", file)   
            file_path = os.path.join(root, folder, file)
            df_L2 = pd.read_csv(file_path)

            # Adding a new column with velocity value (example)
            df_L2["v_inf"] = vel * np.ones_like(df_L2["T"])  # Replace with any logic you need

            # Save the modified DataFrame (overwrite or save as a new file)
            df_L2.to_csv(file_path, index=False)  # Overwrite the original file
            # df_L1.to_csv(file_path.replace(".csv", "_modified.csv"), index=False)  # Save as new file

        if file.startswith("rotor_L3_convergence.csv"):
            print("Processing the file - ", file)   
            file_path = os.path.join(root, folder, file)
            df_L3 = pd.read_csv(file_path)

            # Adding a new column with velocity value (example)
            df_L3["v_inf"] = vel * np.ones_like(df_L3["T"])  # Replace with any logic you need

            # Save the modified DataFrame (overwrite or save as a new file)
            df_L3.to_csv(file_path, index=False)  # Overwrite the original file
            # df_L1.to_csv(file_path.replace(".csv", "_modified.csv"), index=False)  # Save as new file

        if file.startswith("rotor_L4_convergence.csv"):
            print("Processing the file - ", file)   
            file_path = os.path.join(root, folder, file)
            df_L4 = pd.read_csv(file_path)

            # Adding a new column with velocity value (example)
            df_L4["v_inf"] = vel * np.ones_like(df_L4["T"])  # Replace with any logic you need

            # Save the modified DataFrame (overwrite or save as a new file)
            df_L4.to_csv(file_path, index=False)  # Overwrite the original file
            # df_L1.to_csv(file_path.replace(".csv", "_modified.csv"), index=False)  # Save as new file

        if file.startswith("rotor_R1_convergence.csv"):
            print("Processing the file - ", file)   
            file_path = os.path.join(root, folder, file)
            df_R1 = pd.read_csv(file_path)

            # Adding a new column with velocity value (example)
            df_R1["v_inf"] = vel * np.ones_like(df_R1["T"])  # Replace with any logic you need

            # Save the modified DataFrame (overwrite or save as a new file)
            df_R1.to_csv(file_path, index=False)  # Overwrite the original file
            # df_L1.to_csv(file_path.replace(".csv", "_modified.csv"), index=False)  # Save as new file

        if file.startswith("rotor_R2_convergence.csv"):
            print("Processing the file - ", file)   
            file_path = os.path.join(root, folder, file)
            df_R2 = pd.read_csv(file_path)

            # Adding a new column with velocity value (example)
            df_R2["v_inf"] = vel * np.ones_like(df_R2["T"])  # Replace with any logic you need

            # Save the modified DataFrame (overwrite or save as a new file)
            df_R2.to_csv(file_path, index=False)  # Overwrite the original file
            # df_L1.to_csv(file_path.replace(".csv", "_modified.csv"), index=False)  # Save as new file

        if file.startswith("rotor_R3_convergence.csv"):
            print("Processing the file - ", file)   
            file_path = os.path.join(root, folder, file)
            df_R3 = pd.read_csv(file_path)

            # Adding a new column with velocity value (example)
            df_R3["v_inf"] = vel * np.ones_like(df_R3["T"])  # Replace with any logic you need

            # Save the modified DataFrame (overwrite or save as a new file)
            df_R3.to_csv(file_path, index=False)  # Overwrite the original file
            # df_L1.to_csv(file_path.replace(".csv", "_modified.csv"), index=False)  # Save as new file

        if file.startswith("rotor_R4_convergence.csv"):
            print("Processing the file - ", file)   
            file_path = os.path.join(root, folder, file)
            df_R4 = pd.read_csv(file_path)

            # Adding a new column with velocity value (example)
            df_R4["v_inf"] = vel * np.ones_like(df_R4["T"])  # Replace with any logic you need

            # Save the modified DataFrame (overwrite or save as a new file)
            df_R4.to_csv(file_path, index=False)  # Overwrite the original file
            # df_L1.to_csv(file_path.replace(".csv", "_modified.csv"), index=False)  # Save as new file

    print("Number of folders processed - ", counter)
    print("\n\n")
    counter += 1

            

    


