import subprocess
import sys
import os
import pandas as pd
import matplotlib.pyplot as plt

def read_summary_to_row(df, year):
    # Convert the DataFrame to a dictionary and then to a single-row DataFrame
    row_dict = df.set_index('Variable').to_dict()['Value']
    row_dict['Year'] = year  # Add the year to the dictionary
    row_df = pd.DataFrame([row_dict])
    return row_df

def plot_stacked_cost(df, seqsimdir):
    df = df.sort_values(by='Year')
    years = df['Year']
    df.set_index('Year', inplace=True)
    # Select only the relevant cost columns
    cost_columns = ["over_generation_penalty", "under_served_penalty", "generation_costs", "storage_investment_costs", "line_investment_costs"]
    cost_columns_without_under_served = ["over_generation_penalty", "generation_costs", "storage_investment_costs", "line_investment_costs"]

    df_costs = df[cost_columns]
    df_costs_without_under_served = df[cost_columns_without_under_served]

    # Prepare the data for stackplot
    values = [df_costs[category] for category in cost_columns]
    values_without_under_served = [df_costs_without_under_served[category] for category in cost_columns_without_under_served]

    # Plotting
    fig, ax = plt.subplots(figsize=(10, 6))

    # Stacked area plot with baseline set to zero
    ax.stackplot(years, values, labels=cost_columns, baseline='zero')

    # Add labels and title
    ax.set_xlabel('Year')
    ax.set_ylabel('Costs')
    ax.set_title('Stacked Area Plot of Costs Over Years')
    ax.legend(loc='upper left')
    plt.xticks(fontsize='small') 
    plt.xticks(rotation=90)
    plt.savefig(os.path.join(seqsimdir, 'summary_costs.png'))

    # Plotting without under_served_penalty
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.stackplot(years, values_without_under_served, labels=cost_columns_without_under_served)
    ax.set_xlabel('Year')
    ax.set_ylabel('Costs')
    ax.set_title('Stacked Area Plot of Costs Over Years (Without Under Served Penalty)')
    ax.legend(loc='upper left')
    plt.xticks(fontsize='small')
    plt.xticks(rotation=90)
    plt.savefig(os.path.join(seqsimdir, 'summary_costs_without_under_served.png'))

    return 0

def create_and_save_plots(df, seqsimdir):
    df = df.sort_values(by='Year')
    # Set the figure size and adjust subplots
    plt.figure(figsize=(10, 8))

    # Titles and labels for plots
    metrics = {
        "Load_Shed": "Energy Imbalance (puh)",
        "Line_Investments": "Total Line Investments (Increments)",
        "Storage_Investments": "Total Storage Investments (pu)",
        "Used_Renewables": "Total Used Renewables (puh)",
    }

    # Create a plot for each metric
    for i, (metric, title) in enumerate(metrics.items(), start=1):
        plt.subplot(2, 2, i)  # Position each plot in a 2x2 grid
        plt.plot(df['Year'], df[metric], marker='o', linestyle='-', label=title)
        if metric == "Used_Renewables":
            plt.plot(df['Year'], df["Prod_Renewables"], marker='x', linestyle='--', label="Total Produced Renewables (puh)")
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
                               "Used_Renewables",
                               "Prod_Renewables"])
    
    summary_df = pd.DataFrame()

    for item in os.listdir(seqsimdir):
        # Construct the full path to the item
        simdir = os.path.join(seqsimdir, item)
        # Check if this item is a directory
        if os.path.isdir(simdir):
            # Plot investments
            # subprocess.run(["python3", "plot_investments.py", os.path.join(inner_seqsimdir, item)], check=True)
            # Plot hourly generation
            # subprocess.run(["python3", "plot_hourly_generation.py", os.path.join(inner_seqsimdir, item)], check=True)
            # Read CSVs
            energy = pd.read_csv(f"{simdir}/output/energy.csv")
            line_inv = pd.read_csv(f"{simdir}/output/line_investments.csv")
            storage_inv = pd.read_csv(f"{simdir}/output/storage_investments.csv")
            summary = pd.read_csv(f"{simdir}/output/summary_data.csv")

            year = item
            load_shed = energy[energy['Energy_Imbalance'] < 0]['Energy_Imbalance'].sum()
            line_investments = line_inv['Upgrade_Lvl'].sum()
            storage_investments = storage_inv['Storage_Power'].sum()
            used_renewables = energy[[col for col in renewable_types if col in energy.columns]].sum(axis=1).sum()
            prod_renewables = energy[[f"{col}_production" for col in renewable_types if f"{col}_production" in energy.columns]].sum(axis=1).sum()

            # add to the df
            new_row = pd.DataFrame({
                "Year": [year],
                "Load_Shed": [load_shed],
                "Line_Investments": [line_investments],
                "Storage_Investments": [storage_investments],
                "Used_Renewables": [used_renewables],
                "Prod_Renewables": [prod_renewables]
            })
            df = pd.concat([df, new_row], ignore_index=True)

            # add to the summary_df
            row_df = read_summary_to_row(summary, year)
            summary_df = pd.concat([summary_df, row_df], ignore_index=True)

    df.to_csv(f"{seqsimdir}/seq_summary_results.csv", index=False)
    # Call the function with your DataFrame and directory
    create_and_save_plots(df, seqsimdir)
    plot_stacked_cost(summary_df, seqsimdir)




        