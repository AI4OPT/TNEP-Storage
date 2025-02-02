import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import sys
import json

CONGESTION_THRESHOLD = 0.99

simdir = sys.argv[1]
date = sys.argv[2]

# Read CSV and JSON
df_flows = pd.read_csv(f"../../../{simdir}/output/{date}/flow.csv")
df_lines = pd.read_csv(f"../../../{simdir}/output/line_investments.csv")
df_merged = pd.merge(df_flows, df_lines[['Branch_Index', 'Upgrade_Lvl']], on='Branch_Index', how='left')

with open(f"../../../{simdir}/data.json", 'r') as file:
    data = json.load(file)

increment = data["param"]["cap_upgrade_increment"]
if data["param"].get("cap_percent", False):
    df_merged['New_Rate_A'] = df_merged['Rate_A'] * (1.0 + df_merged['Upgrade_Lvl'] * increment)
else:
    df_merged['New_Rate_A'] = df_merged['Rate_A'] + df_merged['Upgrade_Lvl'] * increment

# Filter the datasets
df_old_congestion = df_merged[(abs(df_merged['Power_Flow'] / df_merged['Rate_A']) > CONGESTION_THRESHOLD) & 
                              (abs(df_merged['Power_Flow'] / df_merged['New_Rate_A']) <= CONGESTION_THRESHOLD)]
df_still_congested = df_merged[abs(df_merged['Power_Flow'] / df_merged['New_Rate_A']) > CONGESTION_THRESHOLD]
df_other = df_merged[~df_merged['Branch_Index'].isin(df_old_congestion['Branch_Index']) & 
                     ~df_merged['Branch_Index'].isin(df_still_congested['Branch_Index'])]

df_old_congestion = df_old_congestion.loc[
    df_old_congestion.groupby('Branch_Index')['Power_Flow'].apply(lambda x: x.abs().idxmax())
]
df_still_congested = df_still_congested.loc[
    df_still_congested.groupby('Branch_Index')['Power_Flow'].apply(lambda x: x.abs().idxmax())
]
df_other = df_other.loc[
    df_other.groupby('Branch_Index')['Power_Flow'].apply(lambda x: x.abs().idxmax())
]

# Create a Plotly figure
fig = go.Figure()

# Function to add lines to the map
def add_lines(df, color="black", width=0.2):
    traces = []
    for _, row in df.iterrows():
        traces.append(go.Scattermapbox(
            lat=[row['Lat1'], row['Lat2']],
            lon=[row['Lon1'], row['Lon2']],
            mode='lines',
            showlegend=False,
            line=dict(color=color, width=width)
        ))
    return traces

fig.add_traces(add_lines(df_still_congested, color="red", width=1.5)) 
fig.add_traces(add_lines(df_old_congestion, color="blue", width=1.2))   
fig.add_traces(add_lines(df_other, color="black", width=0.2)) 

fig.update_layout(
    showlegend=False,
    mapbox=dict(
        style="carto-positron",  # Ensure a valid style is used
        center=dict(lat=31.5, lon=-99.5),  # Adjust to desired center
        zoom=4,
    )
)

fig.write_image(f"../../../{simdir}/visual/congestion_{date}.png")