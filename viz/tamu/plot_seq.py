import subprocess
import sys
import os
import pandas as pd
import matplotlib.pyplot as plt

def create_and_save_plots(df, seqsimdir):
    # Set the figure size and adjust subplots
    plt.figure(figsize=(10, 8))

    # Titles and labels for plots
    metrics = {
        "Load_Shed": "Energy Imbalance (MWh)",
        "Line_Investments": "Total Line Investments (Increments)",
        "Storage_Investments": "Total Storage Investments (MW)",
        "Used_Renewables": "Total Used Renewables (MWh)"
    }

    # Create a plot for each metric
    for i, (metric, title) in enumerate(metrics.items(), start=1):
        plt.subplot(2, 2, i)  # Position each plot in a 2x2 grid
        plt.plot(df['Year'], df[metric], marker='o', linestyle='-')
        plt.title(title)
        plt.xlabel('Year')
        plt.ylabel(metric)
        plt.grid(True)
        plt.xticks(fontsize='small') 
        plt.xticks(rotation=90)  # Rotate x-axis ticks 90 degrees

    # Adjust layout to prevent overlap
    plt.tight_layout()

    # Save the figure to a file in the seqsimdir directory
    plt.savefig(os.path.join(seqsimdir, 'summary_plots.png'))

    # Optionally show the plot
    plt.show()

renewable_types = {"wind", "solar", "wind_offshore", "hydro"}

if __name__ == "__main__":
    inner_seqsimdir = sys.argv[1]

    seqsimdir = f"../../{inner_seqsimdir}"

    df = pd.DataFrame(columns=["Year",
                               "Load_Shed",
                               "Line_Investments",
                               "Storage_Investments",
                               "Used_Renewables"])

    for item in os.listdir(seqsimdir):
        # Construct the full path to the item
        simdir = os.path.join(seqsimdir, item)
        # Check if this item is a directory
        if os.path.isdir(simdir):
            # Plot investments
            subprocess.run(["python3", "plot_investments.py", os.path.join(inner_seqsimdir, item)], check=True)
            # Plot hourly generation
            subprocess.run(["python3", "plot_hourly_generation.py", os.path.join(inner_seqsimdir, item)], check=True)
            # Read CSVs
            energy = pd.read_csv(f"{simdir}/output/energy.csv")
            line_inv = pd.read_csv(f"{simdir}/output/line_investments.csv")
            storage_inv = pd.read_csv(f"{simdir}/output/storage_investments.csv")

            year = item
            load_shed = energy[energy['Energy_Imbalance'] < 0]['Energy_Imbalance'].sum()
            line_investments = line_inv['Upgrade_Lvl'].sum()
            storage_investments = storage_inv['Storage_Power'].sum()
            used_renewables = energy[[col for col in renewable_types if col in energy.columns]].sum(axis=1).sum()

            # add to the df
            new_row = pd.DataFrame({
                "Year": [year],
                "Load_Shed": [load_shed],
                "Line_Investments": [line_investments],
                "Storage_Investments": [storage_investments],
                "Used_Renewables": [used_renewables]
            })
            df = pd.concat([df, new_row], ignore_index=True)

    df.to_csv(f"{seqsimdir}/seq_summary_results.csv", index=False)
    # Call the function with your DataFrame and directory
    create_and_save_plots(df, seqsimdir)




        