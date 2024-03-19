import plotly.graph_objects as go
import pandas as pd
import sys

# Assuming simdir is provided, for example, by command line arguments
simdir = sys.argv[1]

# Read CSV file
df_flow = pd.read_csv(f"../../{simdir}/output/flow.csv")

# Normalize Power Flow for line width
max_flow = df_flow['Power_Flow'].abs().max()
df_flow['Normalized_Flow'] = df_flow['Power_Flow'].abs() / max_flow * 10  # Adjust scale for visibility

# Define color function based on directionality and power flow
def red_or_blue(lon1, lon2, power_flow):
    if lon1 <= lon2 and power_flow > 0:
        return 'red'
    elif lon1 > lon2 and power_flow <= 0:
        return 'red'
    else:
        return 'blue'

# Initialize figure with frames for each hour
frames = []
for hour in sorted(df_flow['Hour'].unique()):
    frame_data = []
    df_hour = df_flow[df_flow['Hour'] == hour]
    for index, row in df_hour.iterrows():
        color = red_or_blue(row['Lon1'], row['Lon2'], row['Power_Flow'])
        frame_data.append(go.Scattergeo(
            lon=[row['Lon1'], row['Lon2']],
            lat=[row['Lat1'], row['Lat2']],
            mode='lines',
            line=dict(width=row['Normalized_Flow'], color=color),
            hoverinfo='text',
            text=f"Hour: {row['Hour']}, Power Flow: {row['Power_Flow']}",
        ))
    frames.append(go.Frame(data=frame_data, name=str(hour)))

# Set the initial frame
fig = go.Figure(data=frames[0].data)

# Add frames to the figure
fig.frames = frames

# Slider setup
sliders = [{
    'steps': [{'args': [[frame.name], {'frame': {'duration': 300, 'redraw': True}, 'mode': 'immediate', 'transition': {'duration': 300}}],
               'label': frame.name, 'method': 'animate'} for frame in frames],
    'transition': {'duration': 300},
    'x': 0.1, 'xanchor': 'left', 'y': 0, 'yanchor': 'top'
}]

# Update layout with slider and other configurations
fig.update_layout(
    sliders=sliders,
    title_text='Power Flow Visualization with Directionality, Magnitude, and Time Animation',
    showlegend=False,
    geo=dict(scope='usa', projection_type='albers usa', showland=True, landcolor='rgb(217, 217, 217)'),
    updatemenus=[{
        'type': 'buttons',
        'buttons': [
            {'label': 'Play', 'method': 'animate', 'args': [None, {'frame': {'duration': 500, 'redraw': True}, 'fromcurrent': True, 'transition': {'duration': 300}}]},
            {'label': 'Pause', 'method': 'animate', 'args': [[None], {'frame': {'duration': 0, 'redraw': False}, 'mode': 'immediate', 'transition': {'duration': 0}}]}
        ],
        'direction': 'left',
        'pad': {'r': 10, 't': 87},
        'showactive': False,
        'x': 0.1, 'xanchor': 'right', 'y': 0, 'yanchor': 'top'
    }]
)

# Save the figure as an HTML file
fig.write_html(f"../../{simdir}/visual/flow_animation.html")
