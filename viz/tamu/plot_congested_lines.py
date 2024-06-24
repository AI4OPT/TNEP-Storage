import plotly.graph_objects as go
import pandas as pd
import sys

CONGESTION_THRESHOLD = 0.99

simdir = sys.argv[1]
# Read CSV file
df_flows = pd.read_csv("../../{}/output/flow.csv".format(simdir))

df_lines = df_flows[abs(df_flows['Power_Flow'] / df_flows['Rate_A']) > CONGESTION_THRESHOLD]

# Create a Plotly figure
fig = go.Figure()

# Add traces for each line that meets the threshold
for _, row in df_lines.iterrows():
    fig.add_trace(go.Scattergeo(
        lon=[row['Lon1'], row['Lon2']],
        lat=[row['Lat1'], row['Lat2']],
        mode='lines',
        line=dict(width=1, color='blue'),  # Adjust line width as needed
        showlegend=False
    ))

# Update the layout of the figure
fig.update_layout(
    title=f'Congested Lines at Threshold {CONGESTION_THRESHOLD}',
    geo_scope='usa',
    showlegend=False
)
fig.write_image("../../{}/visual/congested_lines.png".format(simdir))