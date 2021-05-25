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
    html.H1(children='mRNA and Protein transport in neuronal dendrites'),
    dcc.Graph(
        id='graphOutput'
    ),

    html.H6(children="General properties:"),
    html.Div([
        dcc.Dropdown(
            id='transport-option',
            options=[{'label': i, 'value': i} for i in TransportOptions],
            value=TransportOptions[1]
        )
    ]),
    html.Div([
        html.P('Length of domain', style={
            'textAlign': 'center'
        }),
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
    html.Div([

    html.Div([
        html.H6(children="mRNA properties:"),

        html.P(children="Halflife (seconds):"),dcc.Input(id="input-halflife-m", type="number", value=25200),

        html.P(children="Resting state diffusion constant (mu m^2/seconds):"),dcc.Input(id="input-D-m", type="number", value=0.1,step=0.1),

        html.P(children="Active transport velocity (mu m/seconds):"),dcc.Input(id="input-v-m", type="number", value=1),

        html.P(children="Alpha_m (per second):"),dcc.Input(id="alpha-m", type="number", value=0.3,step=0.1),

        html.P(children="Beta plus (per second):"),dcc.Input(id="beta-plus-m", type="number", value=0.3,step=0.1),

        html.P(children="Beta minus (per second):"),dcc.Input(id="beta-minus-m", type="number", value=0.4,step=0.1),
    ],style={'width': '49%', 'display': 'inline-block'}),

    #************************************************
    html.Div([
        html.H6(children="Protein properties:"),

        html.P(children="Halflife (seconds):"),dcc.Input(id="input-halflife-p", type="number", value=432000),

        html.P(children="Resting state diffusion constant (mu m^2/seconds):"),dcc.Input(id="input-D-p", type="number", value=0.1,step=0.1),

        html.P(children="Active transport velocity (mu m/seconds):"),dcc.Input(id="input-v-p", type="number", value=1),

        html.P(children="Alpha_m (per second):"),dcc.Input(id="alpha-p", type="number", value=0.3,step=0.1),

        html.P(children="Beta plus (per second):"),dcc.Input(id="beta-plus-p", type="number", value=0.3,step=0.1),

        html.P(children="Beta minus (per second):"),dcc.Input(id="beta-minus-p", type="number", value=0.4,step=0.1),
    ],style={'width': '49%', 'display': 'inline-block'}),
    ]),

    html.Br(),
    html.Button('Compute', id='button'),
    html.Div(id='my-output')
])


@app.callback(
    dash.dependencies.Output('graphOutput', 'figure'),
    dash.dependencies.Output('my-output', 'children'),
    dash.dependencies.Input('button', 'n_clicks'),
    dash.dependencies.State('transport-option', 'value'),
    dash.dependencies.State('L-slider', 'value'),
    dash.dependencies.State('input-halflife-m', 'value'),
    dash.dependencies.State('input-D-m', 'value'),
    dash.dependencies.State('input-v-m', 'value'),
    dash.dependencies.State('alpha-m', 'value'),
    dash.dependencies.State('beta-plus-m', 'value'),
    dash.dependencies.State('beta-minus-m', 'value'),
    dash.dependencies.State('input-halflife-p', 'value'),
    dash.dependencies.State('input-D-p', 'value'),
    dash.dependencies.State('input-v-p', 'value'),
    dash.dependencies.State('alpha-p', 'value'),
    dash.dependencies.State('beta-plus-p', 'value'),
    dash.dependencies.State('beta-minus-p', 'value'),
    )
def update_output(n_clicks, value,valueL,halflife_m,Dm,vm,alpham,betaPlusm,betaMinusm,\
                                         halflife_p,Dp,vp,alphap,betaPlusp,betaMinusp):
    #fig.update_layout(clickmode='event+select')
    if((n_clicks != None) and  (n_clicks>0)):
        TR.params_input['name'] = value;
        TR.params_input['L'] = valueL;
        TR.params_input['halflife_m'] = float(halflife_m);#float(halflife_m),
        TR.params_input['D_m'] = float(Dm);
        TR.params_input['v_m'] = float(vm);
        TR.params_input['alpha_m'] = float(alpham);
        TR.params_input['beta_plus_m'] = float(betaPlusm);
        TR.params_input['beta_minus_m'] = float(betaMinusm);
        TR.params_input['halflife_p'] = float(halflife_p);
        TR.params_input['D_p'] = float(Dp);
        TR.params_input['v_p'] = float(vp);
        TR.params_input['alpha_p'] = float(alphap);
        TR.params_input['beta_plus_p'] = float(betaPlusp);
        TR.params_input['beta_minus_p'] = float(betaMinusp);
        TR.SolveTransportEqs(N);
        df = TR.GetSolution();
    else:
        df = TR.GetSolution();
    df = df.melt(id_vars=['x'],var_name="Quantity",value_name="Value")
    fig = px.line(df, x="x", y="Value",facet_row="Quantity",color="Quantity");
    fig.update_yaxes(matches=None);

    return fig,str(TR.params_input)#TR.params_input)#fig


if __name__ == '__main__':
    app.run_server(debug=True)
