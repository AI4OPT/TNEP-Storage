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

    renewable_production = np.zeros(NUM_HOURS)
    for key in data["gen"]:
        if data["gen"][key]["gen_type"] in renewable_types:
            renewable_production += data["gen"][key]["profile"][f"{rep_index}"]

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
    plt.figure(figsize=(10, 5))  # Set the figure size
    plt.plot(hours, total_load, label='Total Load', marker='o')  # Plot total load
    plt.plot(hours, renewable_production, label='Renewable Production', marker='o')  # Plot renewable production
    plt.plot(hours, renewable_outputs, label='Renewable Outputs', marker='o')  # Plot renewable outputs
    plt.plot(hours, np.array([total_pmax_nonrenewable] * NUM_HOURS), label='Nonrenewable Capacity', marker='o')  # Plot nonrenewable capacity
    plt.plot(hours, nonrenewable_outputs, label='Nonrenewable Outputs', marker='o')  # Plot renewable outputs
    plt.plot(hours, charge_amount, label='Charging Amount', marker='o')
    plt.plot(hours, discharge_amount, label='Discharge Amount', marker='o')

    plt.title('Energy Production and Load Over Time')
    plt.xlabel('Hour of the Day')
    plt.ylabel('Energy (puh)')
    plt.legend()
    plt.grid(True)
    plt.savefig(f'../../{simdir}/output/{date}/hourly_generation.png')
    plt.close()

    # Second plot: stacked plot
    stacked_production = renewable_production + nonrenewable_outputs + discharge_amount
    stacked_outputs_only = renewable_outputs + nonrenewable_outputs
    stacked_outputs_with_discharge = renewable_outputs + nonrenewable_outputs + discharge_amount

    plt.figure(figsize=(10, 5))  # Set the figure size
    plt.plot(hours, total_load, label='Total Load', marker='o', linestyle='-', alpha=0.8)  # Plot total load
    plt.plot(hours, stacked_production, label='Production', marker='^', linestyle='--', alpha=0.6)  # Plot stacked production
    plt.plot(hours, stacked_outputs_only, label='After Curtailment Production', marker='x', linestyle=':', alpha=0.6)  # Plot stacked production
    plt.plot(hours, stacked_outputs_with_discharge, label='After Curtailment with Discharge', marker='o', linestyle='-.', alpha=0.6)  # Plot stacked production
    plt.title('Stacked Energy Production vs Total Load')
    plt.xlabel('Hour of the Day')
    plt.ylabel('Energy (puh)')
    plt.legend()
    plt.grid(True)
    plt.savefig(f'../../{simdir}/output/{date}/stacked_hourly_generation.png')
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
    plt.ylabel('Energy (puh)')
    plt.legend()
    plt.grid(True)
    plt.savefig(f'../../{simdir}/output/{date}/ue.png')
    plt.close()

    # Fourth plot: curtailed renewable energy
    plt.figure(figsize=(10,5))
    plt.plot(hours, renewable_production - renewable_outputs, label='Renewable Curtailment', marker='o')  # Plot renewable curtailment
    plt.title('Renewable Curtailment over Time')
    plt.xlabel('Hour of the Day')
    plt.ylabel('Energy (puh)')
    plt.legend()
    plt.grid(True)
    plt.savefig(f'../../{simdir}/output/{date}/curtailed_hourly.png')
    plt.close()

