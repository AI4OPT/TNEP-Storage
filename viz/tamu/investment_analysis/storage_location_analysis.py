import pandas as pd
import matplotlib.pyplot as plt
import sys
import numpy as np

simdir = "sim/r1/only_storage/2044"
controldir = "sim/r1/no_upgrades/run_seq/2044"
# simdir = sys.argv[1]
# controldir = sys.argv[2]

def print_storage_insights(simdir, controldir):
    df_storage = pd.read_csv(f"../../../{simdir}/output/storage_investments.csv")
    df_energy = pd.read_csv(f"../../../{controldir}/output/energy.csv")

    merged_data = pd.merge(df_energy, df_storage, on=['Node_Index', 'Node_Name', 'Lat', 'Lon'], how='left')

    nodes_with_storage = merged_data[(merged_data['Storage_Power'] > 0) & (merged_data['Storage_Energy'] > 0)]
    nodes_without_storage = merged_data[(merged_data['Storage_Power'] == 0) & (merged_data['Storage_Energy'] == 0)]

    def sum_cols(df):
        exclude_columns = ['Node_Name', 'Node_Index', 'Lat', 'Lon', 'Hour', 'Storage_Power', 'Storage_Energy']
        # Columns to sum
        columns_to_sum = [col for col in df.columns if col not in exclude_columns]
        # Group by 'Node_Index' and sum the relevant columns
        grouped_df = df.groupby('Node_Index')[columns_to_sum].sum().reset_index()
        # Add the non-summed columns back (taking the first occurrence)
        non_summed_columns = df.groupby('Node_Index')[['Node_Name', 'Lat', 'Lon', 'Storage_Power', 'Storage_Energy']].first().reset_index()
        # Merge the summed and non-summed dataframes
        result_df = pd.merge(non_summed_columns, grouped_df, on='Node_Index', suffixes=('', '_sum'))
        return result_df

    nodes_ws = sum_cols(nodes_with_storage)
    nodes_wos = sum_cols(nodes_without_storage)

    count_imbalance_ws = np.round((nodes_ws['Energy_Imbalance'] != 0).sum(), 2)
    count_imbalance_wos = np.round((nodes_wos['Energy_Imbalance'] != 0).sum(), 2)

    print(f"Percent Energy Imbalance in Nodes with Storage: {count_imbalance_ws} out of {len(nodes_ws)}")
    print(f"Percent Energy Imbalance in Nodes without Storage: {count_imbalance_wos} out of {len(nodes_wos)}")

    def calc_curtailment(df):
        columns_to_process = ['wind', 'solar', 'wind_offshore', 'hydro']
        # Subtract production columns from the corresponding columns without suffix
        for col in columns_to_process:
            production_col = f'{col}_production'
            df[f'{col}_difference'] = df[production_col] - df[col]
        df['total_difference'] = df[[f'{col}_difference' for col in columns_to_process]].sum(axis=1)
        df['curtailment_or_imbalance'] = (df['total_difference'] != 0) | (df['Energy_Imbalance'] != 0)
        return df

    nodes_ws = calc_curtailment(nodes_ws)
    nodes_wos = calc_curtailment(nodes_wos)

    count_curtail_ws = np.round((nodes_ws['total_difference'] != 0).sum(), 2)
    count_curtail_wos = np.round((nodes_wos['total_difference'] != 0).sum(), 2)

    print(f"Percent Count Curtailment in Nodes with Storage: {count_curtail_ws} out of {len(nodes_ws)}")
    print(f"Percent Count Curtailment in Nodes without Storage: {count_curtail_wos} out of {len(nodes_wos)}")

    count_both_ws = np.round((nodes_ws['curtailment_or_imbalance'] != 0).sum(), 2)
    count_both_wos = np.round((nodes_wos['curtailment_or_imbalance'] != 0).sum(), 2)

    print(f"Percent Count Curtailment/Imbalance in Nodes with Storage: {count_both_ws} out of {len(nodes_ws)}")
    print(f"Percent Count Curtailment/Imbalance in Nodes without Storage: {count_both_wos} out of {len(nodes_wos)}")

    return nodes_ws, nodes_wos