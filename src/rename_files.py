## Rename files
# %%
import os

def rename_files_in_folders(root_dir):
    """
    Traverse all folders and rename specific files based on given rules.

    Args:
        root_dir (str): Root directory to start searching for files.
    """
    # Mapping of old filenames to new filenames
    rename_map = {
        "Canarad_convergence.csv": "eMO_Canard_convergence.csv",
        "canard_geometry_data.csv": "Canard_geometry_data.csv",
        "wing_main_convergence.csv": "eMO_wing_main_convergence.csv",
        "wing_geometry_data.csv": "wing_main_geometry_data.csv"
    }
    
    # Traverse through all subdirectories and files
    for dirpath, _, filenames in os.walk(root_dir):
        for filename in filenames:
            # Check if the file matches any of the keys in rename_map
            if filename in rename_map:
                old_file_path = os.path.join(dirpath, filename)
                new_file_path = os.path.join(dirpath, rename_map[filename])
                
                # Rename the file
                try:
                    os.rename(old_file_path, new_file_path)
                    print(f"Renamed: {old_file_path} -> {new_file_path}")
                except Exception as e:
                    print(f"Error renaming {old_file_path}: {e}")




