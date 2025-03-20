import pandas as pd
import matplotlib.pyplot as plt

num_data = 2
if num_data == 1:
    # Load the CSV file
    # file_path = "/mnt/e/eVTOL_model/eVTOL-VehicleModel/result/airfoil_model/Cd_evaluation_results_lambda_0.2.csv"  # Update the path if needed
    file_path = "/mnt/e/eVTOL_model/eVTOL-VehicleModel/result/rotor_model/rotorModel_evaluation_results.csv"
    df = pd.read_csv(file_path)

    # Convert the "Relative L2 Norm Cd (%)" column to numeric values
    df.iloc[:, 1] = pd.to_numeric(df.iloc[:, 1], errors='coerce')

    # Create a professional-looking histogram for the Relative L2 Error (%)
    plt.figure(figsize=(8, 5))

    # Use a refined style for professional appearance
    plt.style.use("seaborn-v0_8-whitegrid")

    # Plot histogram with increased bin size and refined aesthetics
    plt.hist(df.iloc[:, 1], bins=40, edgecolor='black', alpha=0.75, linewidth=1.2)

    # Labels and title with LaTeX formatting for an academic look
    plt.xlabel(r"Relative $L_2$ Error (%)", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.title(r"Distribution of Relative $L_2$ Error for $C_D$ Model", fontsize=16, fontweight="bold")

    # Customize ticks for better readability
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    # Remove the top and right borders for a cleaner academic look
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)

    # Save the figure in high-resolution format
    plt.savefig("Cd_L2_Error_Histogram.pdf", bbox_inches="tight", dpi=300)

    # Show the plot
    plt.show()

elif num_data == 2:
    # Load the CSV file
    # file_path = "/mnt/e/eVTOL_model/eVTOL-VehicleModel/result/wing_model/wingModel_evaluation_results.csv"  # Update the path if needed
    file_path = "/mnt/e/eVTOL_model/eVTOL-VehicleModel/result/rotor_model/rotorModel_evaluation_results.csv"
    df = pd.read_csv(file_path)

    # Convert the relevant columns to numeric values
    df.iloc[:, 1] = pd.to_numeric(df.iloc[:, 1], errors='coerce')  # Cl error
    df.iloc[:, 2] = pd.to_numeric(df.iloc[:, 2], errors='coerce')  # Cd error

    # Create the figure
    plt.figure(figsize=(8, 5))
    plt.style.use("seaborn-v0_8-whitegrid")

    # Plot overlapping histograms
    plt.hist(df.iloc[:, 1], bins=40, edgecolor='black', alpha=0.6, label=r"Relative $L_2$ Error for $C_T$", linewidth=1.2)
    plt.hist(df.iloc[:, 2], bins=40, edgecolor='black', alpha=0.6, color='tab:orange', label=r"Relative $L_2$ Error for $C_Q$", linewidth=1.2)

    # Labels and title
    plt.xlabel("Relative $L_2$ Error (%)", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.title(r"Histogram of Relative $L_2$ Error for $C_T$ and $C_Q$", fontsize=16, fontweight="bold")

    # Customize ticks and legend
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)

    # Remove top and right borders
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)

    # Save and show
    plt.savefig("Ct_Cq_L2_Error_Overlapping_Histogram.pdf", bbox_inches="tight", dpi=300)
    plt.show()