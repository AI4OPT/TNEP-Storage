import matplotlib.pyplot as plt
import pandas as pd
import sys
import json
import numpy as np
import toml

def average_profile(profile_dict, num_hours):
    num_scenarios = len(profile_dict)
    scenario_sum = np.zeros(num_hours)
    # Sum up all the arrays across scenarios
    for scenario in profile_dict.values():
        scenario_sum += np.array(scenario)
    # Calculate the average of scenarios
    return scenario_sum / num_scenarios

simdir = sys.argv[1]
# Read data json
with open(f"../../{simdir}/data.json", 'r') as file:
    data = json.load(file)

with open(f'../../{simdir}/config.toml', 'r') as toml_file:
    config = toml.load(toml_file)
    
NUM_HOURS = config["num_hours"]
nonrenewable_types = set(config["nonrenewable_types"])
renewable_types = set(config["renewable_types"])

 # Compute the total nonrenewable capacity
total_pmax_nonrenewable = 0
for key in data["gen"]:
    # Check if the 'gen_type' of this generator is in the set of nonrenewable types
    if data["gen"][key]["gen_type"] in nonrenewable_types:
        # Add the 'pmax' value of this generator to the total
        total_pmax_nonrenewable += data["gen"][key].get("pmax", 0)

for rep_index in range(1, len(data["param"]["dates"]) + 1):
    date = data["param"]["dates"][rep_index - 1]
    # Read CSV file
    df_energy = pd.read_csv(f"../../{simdir}/output/{date}/energy.csv")

    # Initialize a NumPy array to store the sum of average loads for each hour (24 hours)
    total_load = np.zeros(NUM_HOURS)

    # Calculate the load per hour for each bus and sum them
    for key in data["bus"]:
        total_load += data["bus"][key]["load"][f"{rep_index}"]

    # Calculate renewable production
    renewable_production = np.zeros(NUM_HOURS)
    for key in data["gen"]:
        if data["gen"][key]["gen_type"] in renewable_types:
            renewable_production += data["gen"][key]["profile"][f"{rep_index}"]

    # Calculate wind and solar production
    wind_production = np.zeros(NUM_HOURS)
    solar_production = np.zeros(NUM_HOURS)
    for key in data["gen"]:
        if data["gen"][key]["gen_type"] == "wind":
            wind_production += data["gen"][key]["profile"][f"{rep_index}"]
        elif data["gen"][key]["gen_type"] == "solar":
            solar_production += data["gen"][key]["profile"][f"{rep_index}"]

    # Now compute the actual renewable output
    summed_renewable_df = df_energy.groupby('Hour')[list(renewable_types & set(df_energy.columns))].sum()
    renewable_outputs = summed_renewable_df.sum(axis=1).values
    # Now compute the actual nonrenewable output
    summed_nonrenewable_df = df_energy.groupby('Hour')[list(nonrenewable_types & set(df_energy.columns))].sum()
    nonrenewable_outputs = summed_nonrenewable_df.sum(axis=1).values
    # Now compute charge and discharge
    summed_charge_df = df_energy.groupby('Hour')['Charge'].sum()
    charge_amount = summed_charge_df.values
    summed_discharge_df = df_energy.groupby('Hour')['Discharge'].sum()
    discharge_amount = summed_discharge_df.values

    hours = np.arange(1, NUM_HOURS+1)
    # Plot: just hourly generation of every metric to give quick overview
    plt.figure(figsize=(10, 5))  # Set the figure size
    plt.plot(hours, total_load, label='Total Load', marker='o')  # Plot total load
    plt.plot(hours, renewable_production, label='Renewable Production', marker='o')  # Plot renewable production
    plt.plot(hours, renewable_outputs, label='Renewable Outputs', marker='o')  # Plot renewable outputs
    plt.plot(hours, np.array([total_pmax_nonrenewable] * NUM_HOURS), label='Nonrenewable Capacity', marker='o')  # Plot nonrenewable capacity
    plt.plot(hours, nonrenewable_outputs, label='Nonrenewable Outputs', marker='o')  # Plot renewable outputs
    plt.plot(hours, charge_amount, label='Charging Amount', marker='o')
    plt.plot(hours, discharge_amount, label='Discharge Amount', marker='o')

    plt.title(f'Energy Production and Load Over Time on {date}')
    plt.xlabel('Hour of the Day')
    plt.xlim(1,24)
    plt.ylabel('Energy (puh)')
    plt.ylim(bottom=0)
    plt.legend()
    plt.grid(True)
    plt.savefig(f'../../{simdir}/output/{date}/hourly_generation.png')
    plt.close()

    # Plot: resource adequacy - is there enough energy to provide for load for the day?
    plt.figure(figsize=(10, 5))  # Set the figure size
    plt.plot(hours, total_load, label='Total Load', marker='o')  # Plot total load
    plt.plot(hours, renewable_production + total_pmax_nonrenewable, label='Renewable Production + Nonrenewable Max', marker='o') 

    plt.title(f'Resource Adequacy on {date}, Extra Energy = {int(sum(renewable_production + total_pmax_nonrenewable) - sum(total_load))}')
    plt.xlabel('Hour of the Day')
    plt.xlim(1,24)
    plt.ylabel('Energy (puh)')
    plt.ylim(bottom=0)
    plt.legend()
    plt.grid(True)
    plt.savefig(f'../../{simdir}/output/{date}/resource_adequacy.png')
    plt.close()

    
    # Third plot: unserved energy
    ue = np.zeros(NUM_HOURS)
    df_filtered = df_energy[df_energy['Energy_Imbalance'] <= 0]
    df_grouped = df_filtered.groupby('Hour')['Energy_Imbalance'].sum().reset_index()
    for _, row in df_grouped.iterrows():
        hour = int(row['Hour'])
        ue[hour - 1] = -1 * row['Energy_Imbalance']

    plt.figure(figsize=(10,5))
    plt.plot(hours, ue, label='Unserved Energy', marker='o')
    plt.title('Unserved Energy over Time')
    plt.xlabel('Hour of the Day')
    plt.xlim(1,24)
    plt.ylabel('Energy (puh)')
    plt.ylim(bottom=0)
    plt.legend()
    plt.grid(True)
    plt.savefig(f'../../{simdir}/output/{date}/ue.png')
    plt.close()

    # Fourth plot: curtailed renewable energy
    plt.figure(figsize=(10,5))
    plt.plot(hours, renewable_production - renewable_outputs, label='Renewable Curtailment', marker='o')  # Plot renewable curtailment
    plt.title('Renewable Curtailment over Time')
    plt.xlabel('Hour of the Day')
    plt.xlim(1,24)
    plt.ylabel('Energy (puh)')
    plt.ylim(bottom=0)
    plt.legend()
    plt.grid(True)
    plt.savefig(f'../../{simdir}/output/{date}/curtailed_hourly.png')
    plt.close()

    # New plot: will show the cycling of storage throughout the day
    discharge_color = "#2C7FB8"  # Blue for discharge
    charge_color = "#E69F00"  # Orange for charge
    plt.figure(figsize=(10, 5))  # Set the figure size
    plt.plot(hours, charge_amount, label='Storage Charging', color=charge_color, marker='o') 
    plt.plot(hours, discharge_amount, label='Storage Discharging', color=discharge_color, marker='o') 
    plt.title(f'Total Charge/Discharge Over Time on {date}')
    plt.xlabel('Hour of the Day')
    plt.xlim(1,24)
    plt.ylabel('Energy (puh)')
    plt.ylim(bottom=0)
    plt.legend()
    plt.grid(True)
    plt.savefig(f'../../{simdir}/output/{date}/charge_discharge.png')
    plt.close()

    # New plot: will show charge and discharge
    # show load, wind (before curtailment), solar (before curtailment)
    load_avg = round(np.mean(total_load), 3)
    wind_avg = round(np.mean(wind_production), 3)
    solar_avg = round(np.mean(solar_production), 3)

    plt.figure(figsize=(10, 5))  # Set the figure size
    plt.plot(hours, total_load, label=f'Avg Load: {load_avg}', color='#377eb8', marker='o')  # Blue for load
    plt.plot(hours, wind_production, label=f'Avg Wind: {wind_avg}', color='#4daf4a', marker='o')  # Green for wind
    plt.plot(hours, solar_production, label=f'Avg Solar: {solar_avg}', color='#ff7f00', marker='o')  # Orange for solar
    plt.axhline(y=total_pmax_nonrenewable, color='#e41a1c', linestyle='--', label='Nonrenewable Cap')  # Dark red
    plt.ylim(0, 1000)
    plt.title(f'Load, Wind, Solar on {date}')
    plt.xlabel('Hour of the Day')
    plt.xlim(1,24)
    plt.ylabel('Energy (puh)')
    plt.ylim(bottom=0)
    plt.legend()
    plt.grid(True)
    plt.savefig(f'../../{simdir}/output/{date}/hourly_load_solar_wind.png')
    plt.close()

    # New plot: stacked plot that will show how load is being served
    # non-renewables, renewables (after curtailment), storage discharge, and load-shed in that order
    # charging amount should be subtracted
    total_load = np.array(total_load)
    renewable_outputs = np.array(renewable_outputs)  # After curtailment
    nonrenewable_outputs = np.array(nonrenewable_outputs)
    discharge_amount = np.array(discharge_amount)
    charge_amount = np.array(charge_amount)
    load_shed = np.maximum(0, total_load - (nonrenewable_outputs + renewable_outputs + discharge_amount))

    # Stack data for the layers: Non-Renewables, Renewables, Discharge
    y = np.row_stack([nonrenewable_outputs, renewable_outputs, discharge_amount, load_shed])
    y_stack = np.cumsum(y, axis=0)

    # Labels and colors for each layer
    labels = ["Non-Renewables", "Renewables (After Curtailment)", "Storage Discharge", "Load Shed"]
    colors = ["#8C510A", "#5AAE61", "#2C7FB8", "#E41A1C"]  # Dark brown, green, blue, bright red

    # Plotting the stackplot
    plt.figure(figsize=(10, 5))

    # Fill each layer and include in the legend
    plt.fill_between(hours, 0, y_stack[0, :], facecolor=colors[0], alpha=0.7, label=labels[0])
    plt.fill_between(hours, y_stack[0, :], y_stack[1, :], facecolor=colors[1], alpha=0.7, label=labels[1])
    plt.fill_between(hours, y_stack[1, :], y_stack[2, :], facecolor=colors[2], alpha=0.7, label=labels[2])
    plt.fill_between(hours, y_stack[2, :], y_stack[3, :], facecolor=colors[3], alpha=0.7, label=labels[3])

    # Plot total load as a line on top for reference
    plt.plot(hours, total_load, label='Total Load', color='black', linestyle='--', linewidth=1.0)

    # Set plot labels and title
    plt.title(f'Energy Sources Serving Load Over Time on {date}')
    plt.xlabel('Hour of the Day')
    plt.xlim(1,24)
    plt.ylabel('Energy (puh)')
    plt.ylim(bottom=0)
    plt.legend(loc='lower right', frameon=True)  # Add legend with colors
    plt.grid(True)

    # Save the plot
    plt.savefig(f'../../{simdir}/output/{date}/stacked_energy_serving_load.png')
    plt.close()