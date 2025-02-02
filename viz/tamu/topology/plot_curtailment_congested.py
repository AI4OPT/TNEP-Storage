import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import sys
import json

simdir = sys.argv[1]
date = sys.argv[2]
CONGESTION_THRESHOLD = 0.99
hour_start = int(sys.argv[3]) if len(sys.argv) > 3 else None
hour_end = int(sys.argv[4]) if len(sys.argv) > 4 else None

df_storage = pd.read_csv(f"../../../{simdir}/output/storage_investments.csv")
df_lines = pd.read_csv(f"../../../{simdir}/output/line_investments.csv")
df_energy = pd.read_csv(f"../../../{simdir}/output/{date}/energy.csv")
production_cols = [col for col in df_energy.columns if col.endswith('_production')]

# Filter by hour range if specified
if hour_start is not None and hour_end is not None:
    df_energy = df_energy[(df_energy['Hour'] >= hour_start) & (df_energy['Hour'] <= hour_end)]

# Calculate curtailment for each production type
for prod_col in production_cols:
    original_col = prod_col.replace('_production', '')
    curtailment_col = f'Curtailment_{original_col}'
    df_energy[curtailment_col] = df_energy[prod_col] - df_energy[original_col]

# Calculate total curtailment across all production types for each node and hour
curtailment_cols = [f'Curtailment_{col.replace("_production", "")}' for col in production_cols]
df_energy['Total_Curtailment'] = df_energy[curtailment_cols].sum(axis=1)

# Read CSV and JSON
df_flows = pd.read_csv(f"../../../{simdir}/output/{date}/flow.csv")
# Filter by hour range if specified
if hour_start is not None and hour_end is not None:
    df_flows = df_flows[(df_flows['Hour'] >= hour_start) & (df_flows['Hour'] <= hour_end)]

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

# Plot for Total Curtailment
fig = px.scatter_mapbox(
    df_energy,
    lat='Lat',
    lon='Lon',
    size='Total_Curtailment',               # Size based on total curtailment
    color='Total_Curtailment',              # Color based on total curtailment
    color_continuous_scale='Burg_r',     # Use a reversed viridis color scale (smaller curtailments are darker)
    center=dict(lat=31.5, lon=-99.5),       # Center the map on Texas
    zoom=4,                                 # Set zoom level
    labels={'Total_Curtailment': 'Total Curtailment (100 MWh)'},  # Label for the color legend
    mapbox_style="carto-positron",          # Map style
)

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

fig.add_traces(add_lines(df_old_congestion, color="blue", width=1.2))   
fig.add_traces(add_lines(df_still_congested, color="red", width=1.5)) 
fig.add_traces(add_lines(df_other, color="black", width=0.2)) 

hour_range_suffix = f"_{hour_start}-{hour_end}" if hour_start is not None and hour_end is not None else ""
fig.write_image(f"../../../{simdir}/visual/curtailment_congestion_{date}{hour_range_suffix}.png")