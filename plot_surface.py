"""Plots 3D dotted and triangulated surface
"""

import surface
import logging
import util
import plotly.graph_objects as go
import sys
import dash
from dash import html
from dash import dcc
import numpy as np

if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', level=logging.INFO)

    conf = util.parse_argv.parse(sys.argv)

    n_dots, dots = surface.dots.read('/wrk3/tests/dotted-surface.dat')
    n_edges, edges = surface.edges.read('/wrk3/tests/edges-triangulated-surface.dat')

    dotted_srf = go.Scatter3d(x=dots[:, 0],
                              y=dots[:, 1],
                              z=dots[:, 2],
                              marker=go.scatter3d.Marker(size=1),
                              opacity=0.7,
                              mode='markers')

    x = np.zeros(shape=0)
    y = np.zeros(shape=0)
    z = np.zeros(shape=0)
    counter = 0
    for edge in edges:
        start = edge[0:3]
        end = edge[3:]
        if counter < 30000:
            x = np.append(x, start[0])
            y = np.append(y, start[1])
            z = np.append(z, start[2])
            x = np.append(x, end[0])
            y = np.append(y, end[1])
            z = np.append(z, end[2])
            # Keep these:
            x = np.append(x, None)
            y = np.append(y, None)
            z = np.append(z, None)
        counter += 1
    edges_srf = go.Scatter3d(x=x, y=y, z=z, mode='lines', marker=go.scatter3d.Marker(size=2))

    figure = go.Figure([edges_srf])

    app = dash.Dash()
    app.layout = html.Div(
        style={'height': '1300px', 'margin': '0', 'border': '0', 'padding': '0'},
        children=[dcc.Graph(style={'height': '100%'}, figure=figure)]
    )

    app.run_server(debug=True)
