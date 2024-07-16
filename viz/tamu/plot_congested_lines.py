import plotly.graph_objects as go
import pandas as pd
import json
import sys

CONGESTION_THRESHOLD = 0.99

simdir = sys.argv[1]
# Read CSV and JSON
df_flows = pd.read_csv("../../{}/output/flow.csv".format(simdir))
df_line_investments = pd.read_csv(f"../../{simdir}/output/line_investments.csv")
df_merged = pd.merge(df_flows, df_line_investments[['Branch_Index', 'Upgrade_Lvl']], on='Branch_Index', how='left')

with open(f"../../{simdir}/data.json", 'r') as file:
    data = json.load(file)

increment = data["param"]["cap_upgrade_increment"]
df_merged['New_Rate_A'] = df_merged['Rate_A'] + df_merged['Upgrade_Lvl'] * increment

df_old_congestion = df_merged[abs(df_merged['Power_Flow'] / df_merged['Rate_A']) > CONGESTION_THRESHOLD]
df_still_congested = df_merged[abs(df_merged['Power_Flow'] / df_merged['New_Rate_A']) > CONGESTION_THRESHOLD]

# Create a Plotly figure
fig = go.Figure()

# Add traces for each line that meets the threshold
for _, row in df_old_congestion.iterrows():
    fig.add_trace(go.Scattergeo(
        lon=[row['Lon1'], row['Lon2']],
        lat=[row['Lat1'], row['Lat2']],
        mode='lines',
        line=dict(width=1, color='blue'),  # Adjust line width as needed
        showlegend=False
    ))

for _, row in df_still_congested.iterrows():
    fig.add_trace(go.Scattergeo(
        lon=[row['Lon1'], row['Lon2']],
        lat=[row['Lat1'], row['Lat2']],
        mode='lines',
        line=dict(width=1, color='red'),  # Adjust line width as needed
        showlegend=False
    ))

# Update the layout of the figure
fig.update_layout(
    title=f'Congested Lines at Threshold {CONGESTION_THRESHOLD}',
    geo_scope='usa',
    showlegend=False
)
fig.write_image("../../{}/visual/congested_lines.png".format(simdir))