import plotly.graph_objects as go
import pandas as pd
import sys

# Assuming simdir is provided, for example, by command line arguments
simdir = sys.argv[1]
additional_file = sys.argv[2] if len(sys.argv) > 2 else None



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

def create_figure(df, title_text):
    frames = []
    for hour in sorted(df['Hour'].unique()):
        frame_data = []
        df_hour = df[df['Hour'] == hour]
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

    fig = go.Figure(data=frames[0].data)
    fig.frames = frames

    # Slider setup
    sliders = [{
        'steps': [{'args': [[frame.name], 
                            {'frame': {'duration': 300, 'redraw': True},
                             'mode': 'immediate',
                             'transition': {'duration': 300}}],
                   'label': frame.name, 'method': 'animate'} for frame in frames],
        'transition': {'duration': 300},
        'x': 0.1, 'xanchor': 'left', 'y': 0, 'yanchor': 'top'
    }]

    # Update layout with slider and other configurations
    fig.update_layout(
        sliders=sliders,
        title_text=title_text,
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
    
    return fig


# Create the figure for the primary data
fig_flow = create_figure(df_flow, "TAMU transport flows")
fig_flow.write_html(f"../../{simdir}/visual/TAMU_flow_animation.html")

if additional_file:
    df_additional = pd.read_csv(f"../../{additional_file}")
    df_additional['Normalized_Flow'] = df_additional['Power_Flow'].abs() / max_flow * 10  # Use the same max_flow for normalization
    fig_additional = create_figure(df_additional, "EIA transport flows")
    fig_additional.write_html(f"../../{simdir}/visual/EIA_flow_animation.html")

