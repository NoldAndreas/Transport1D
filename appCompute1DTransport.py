import json

import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.express as px
import pandas as pd
import numpy as np

from SolveTransport import Transport

#************************
#************************
#Parameters
#************************
supply_vs_demand = 1.0; #define rate of over- or undersupply
L                = 1000;#length of domain (in mu_m)
N      = 200;
TR = Transport('params_standard.json',L,supply_vs_demand);
TR.params_input['name'] = 'somatic, passive protein transport';

TransportOptions = ['somatic, passive protein transport',\
                    'somatic, active protein transport',\
                    'dendritic, passive mRNA transport',\
                    'dendritic, active mRNA transport'];
#************************
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)



#fig.update_traces(marker_size=10)

app.layout = html.Div([
    dcc.Graph(
        id='graphOutput'
    ),
    html.Div([
        dcc.Dropdown(
            id='transport-option',
            options=[{'label': i, 'value': i} for i in TransportOptions],
            value=TransportOptions[1]
        )
    ]),
    html.Br(),
    html.Div([
        dcc.Slider(
            id='L-slider',
            min=100,
            max=1000,
            value=200,
            marks={str(l): str(l) for l in [100,200,300,400,500,600,700,800,900,1000]},
            step=None
            )]
    ),
    html.Br(),
    html.Button('Compute', id='button'),
    html.Div(id='my-output')
])


@app.callback(
    dash.dependencies.Output('graphOutput', 'figure'),
    dash.dependencies.Output('my-output', 'children'),
    dash.dependencies.Input('button', 'n_clicks'),
    dash.dependencies.State('transport-option', 'value'),
    dash.dependencies.State('L-slider', 'value'))
def update_output(n_clicks, value,valueL):
    #fig.update_layout(clickmode='event+select')
    if((n_clicks != None) and  (n_clicks>0)):
        TR.params_input['name'] = value;
        TR.params_input['L'] = valueL;
        TR.SolveTransportEqs(N);
        df = TR.GetSolution();
    else:
        df = TR.GetSolution();
    df = df.melt(id_vars=['x'],var_name="Quantity",value_name="Value")
    fig = px.line(df, x="x", y="Value",facet_row="Quantity",color="Quantity");
    fig.update_yaxes(matches=None);

    return fig,str(valueL)


if __name__ == '__main__':
    app.run_server(debug=True)
