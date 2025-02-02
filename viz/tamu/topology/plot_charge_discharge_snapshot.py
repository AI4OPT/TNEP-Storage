import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import sys
import json

simdir = sys.argv[1]
date = sys.argv[2]
hour = sys.argv[3]

df_storage = pd.read_csv(f"../../../{simdir}/output/storage_investments.csv")
df_lines = pd.read_csv(f"../../../{simdir}/output/line_investments.csv")
df_energy = pd.read_csv(f"../../../{simdir}/output/{date}/energy.csv")

df_energy = df_energy[(df_energy['Hour'] == hour) & ((df_energy['Charge'] != 0) | (df_energy['Discharge'] != 0))]

# Add a dummy row for consistent scaling
dummy_row = {'Lat': 0, 'Lon': 0, 'Charge': 30, 'Discharge': 30, 'Hour': hour}
df_energy = pd.concat([df_energy, pd.DataFrame([dummy_row])], ignore_index=True)

# Create scatter mapbox
fig = go.Figure()

# Plot Charge points as red circles with size proportional to the Charge amount
fig.add_trace(go.Scattermapbox(
    lat=df_energy[df_energy['Charge'] > 0]['Lat'],
    lon=df_energy[df_energy['Charge'] > 0]['Lon'],
    mode='markers',
    marker=dict(
        size=df_energy[df_energy['Charge'] > 0]['Charge'],  # Set size based on Charge
        sizemode='area',  # Size circles by area for better proportionality
        color='red'
    ),
    name='Charge'
))

# Plot Discharge points as blue circles with size proportional to the Discharge amount
fig.add_trace(go.Scattermapbox(
    lat=df_energy[df_energy['Discharge'] > 0]['Lat'],
    lon=df_energy[df_energy['Discharge'] > 0]['Lon'],
    mode='markers',
    marker=dict(
        size=df_energy[df_energy['Discharge'] > 0]['Discharge'],  # Set size based on Discharge
        sizemode='area',  # Size circles by area for better proportionality
        color='blue'
    ),
    name='Discharge'
))

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

# fig.add_traces(add_lines(df_lines, color="black", width=0.2)) 

# Map settings
fig.update_layout(
    mapbox=dict(
        style="carto-positron",
        center=dict(lat=31.5, lon=-99.5),
        zoom=4
    ),
    margin={"r":0,"t":0,"l":0,"b":0},
    showlegend=False
)


fig.write_image(f"../../../{simdir}/visual/charge_discharge_{date}_{hour}.png")