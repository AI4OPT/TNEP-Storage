import pandas as pd
import json
import sys
import numpy as np

CONGESTION_THRESHOLD = 0.99

simdir = "sim/r1/only_transmission/2048"
controldir = "sim/r1/no_upgrades/run_seq/2048"
# simdir = sys.argv[1]
# controldir = sys.argv[2]

def print_line_insights(simdir, controldir):
    df_line_investments = pd.read_csv(f"../../../{simdir}/output/line_investments.csv")
    df_flows = pd.read_csv(f"../../../{controldir}/output/flow.csv")
    df_merged = pd.merge(df_flows, df_line_investments[['Branch_Index', 'Upgrade_Lvl']], on='Branch_Index', how='left')
    df_merged['Congested'] = (df_merged['Power_Flow'].abs() / df_merged['Rate_A']) > CONGESTION_THRESHOLD
    congested_hours = df_merged.groupby('Branch_Index').agg({
        'Congested': 'sum',
        'Rate_A': 'first', 
        'Upgrade_Lvl': 'first'
    }).reset_index()
    congested_hours.columns = ['Branch_Index', 'Congested_Hours', 'Rate_A', 'Upgrade_Lvl']
    filtered_congested_upgraded = congested_hours[congested_hours['Upgrade_Lvl'] > 0]
    filtered_congested_nonupgraded = congested_hours[(congested_hours['Upgrade_Lvl'] == 0) & (congested_hours['Congested_Hours'] > 0)]
    
    return filtered_congested_upgraded, filtered_congested_nonupgraded




