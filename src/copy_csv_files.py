# %%
import os
import shutil

def copy_csv_files(source_dir, destination_dir):
    """
    Copies all CSV files from the source directory to the destination directory
    while maintaining the folder structure.
    """
    for root, _, files in os.walk(source_dir):
        for file in files:
            if file.endswith(".csv"):
                # Get the relative path of the file from the source directory
                relative_path = os.path.relpath(root, source_dir)
                
                # Create the corresponding destination path
                dest_folder = os.path.join(destination_dir, relative_path)
                os.makedirs(dest_folder, exist_ok=True)

                # Copy the file
                src_file = os.path.join(root, file)
                dest_file = os.path.join(dest_folder, file)
                shutil.copy2(src_file, dest_file)

                print(f"Copied: {src_file} -> {dest_file}")




# %%
# Specify the source and destination directories
source_directory = "/mnt/e/eVTOL_model/eVTOL-VehicleModel/FLOWUnsteady_simulations/aircraft_data"  # Change this to your actual source path
destination_directory = "/mnt/e/eVTOL_model/eVTOL-VehicleModel/data/aircraft_data"  # Change this to your desired destination path

# Run the function
copy_csv_files(source_directory, destination_directory)


