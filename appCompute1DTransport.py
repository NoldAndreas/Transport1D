import json

import dash
import dash_table
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import pandas as pd
import plotly.express as px
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

resultsToDisplay = ['Particles in active transport',\
                    'Excess proteins translated',\
                    'Excess mRNA transcribed'];

paramsToDisplay = ['L','halflife_m','uptake_mRNA','D_m','v_m',\
                    'uptake_proteins','halflife_p','D_p','v_p',\
                    'boundaryCondition'];

colsToDisplay = [{"name":"idx","id":'idx'}];
for p in paramsToDisplay:
    colsToDisplay.append({"name":p,"id":p});


rows_results = [];
rows_params  = [];
#************************
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)



#fig.update_traces(marker_size=10)

app.layout = html.Div([
    html.H1(children='mRNA and Protein transport in neuronal dendrites'),

    html.H3(children="Parameters:"),
    html.H6(children="General properties:"),

    dcc.RadioItems(
        id='boundaryCondition',
        options=[
            {'label': 'Fix influx', 'value': 'fix_influx'},
            {'label': 'Fix end concentration to zero', 'value': 'fix_EndConcentration'},
        ],
        value='fix_influx',
        labelStyle={'display': 'inline-block'}
    ) ,

    html.Div([
        dcc.Dropdown(
            id='transport-option',
            options=[{'label': i, 'value': i} for i in TransportOptions],
            value=TransportOptions[0]
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

        html.P('Synaptic uptake / localization:'),
        dcc.RadioItems(
            id='uptake_mRNA',
            options=[
                {'label': 'linear', 'value': 'linear'},
                {'label': 'tanh', 'value': 'tanh'},
                {'label': 'constant', 'value': 'const'},
            ],
            value='const',
            labelStyle={'display': 'inline-block'}
        ),

        html.P(children="Halflife (seconds):"),dcc.Input(id="input-halflife-m", type="number", value=25200),

        html.P(children="Resting state diffusion constant (mu m^2/seconds):"),dcc.Input(id="input-D-m", type="number", value=0.1,step=0.1),

        html.P(children="Active transport velocity (mu m/seconds):"),dcc.Input(id="input-v-m", type="number", value=1),
    ],style={'width': '49%', 'display': 'inline-block'}),

    #************************************************
    html.Div([
        html.H6(children="Protein properties:"),

        html.P('Synaptic uptake / localization:'),
        dcc.RadioItems(
            id='uptake_proteins',
            options=[
                {'label': 'linear', 'value': 'linear'},
                {'label': 'tanh', 'value': 'tanh'},
                {'label': 'constant', 'value': 'const'},
            ],
            value='const',
            labelStyle={'display': 'inline-block'}
        ),

        html.P(children="Halflife (seconds):"),dcc.Input(id="input-halflife-p", type="number", value=432000),

        html.P(children="Resting state diffusion constant (mu m^2/seconds):"),dcc.Input(id="input-D-p", type="number", value=0.1,step=0.1),

        html.P(children="Active transport velocity (mu m/seconds):"),dcc.Input(id="input-v-p", type="number", value=1),
    ],style={'width': '49%', 'display': 'inline-block'}),
    ]),

    html.Br(),
    html.Button('Compute', id='button'),

    #html.Div(id='my-output-parameters'),

    html.H3(children="Results:"),

    html.Div([
    html.Div(
    dcc.Graph(
        id='graphOutput'
    ),style={'width': '49%', 'display': 'inline-block'}),
    html.Div(
    dcc.Graph(
        id='graph-costs-1'
    ),style={'width': '49%', 'display': 'inline-block'}),
    ]),

    html.Label('Parameters of all computations'),
    dash_table.DataTable(
        id='tableParameters',
        columns=colsToDisplay,
        #[{'name': 'Column 1', 'id': 'column1'},{'name': 'Column 2', 'id': 'column2'}]
    ),

    #html.Div(id='my-output-results',style={'whiteSpace': 'pre-line'}),
    html.Div(
        [html.P(children='Site by Andreas Nold'),
        dcc.Link('Github', href='https://github.com/NoldAndreas/Transport1D',refresh=True)],
        style={'textAlign':'center'}),

])


@app.callback(
    dash.dependencies.Output('graphOutput', 'figure'),
    dash.dependencies.Output('graph-costs-1', 'figure'),
    #dash.dependencies.Output('my-output-parameters', 'children'),
    dash.dependencies.Output('tableParameters', component_property='data'),
    #dash.dependencies.Output('my-output-results', 'children'),
    dash.dependencies.Input('button', 'n_clicks'),
    dash.dependencies.State('boundaryCondition', 'value'),
    dash.dependencies.State('uptake_mRNA', 'value'),
    dash.dependencies.State('uptake_proteins', 'value'),
    dash.dependencies.State('transport-option', 'value'),
    dash.dependencies.State('L-slider', 'value'),
    dash.dependencies.State('input-halflife-m', 'value'),
    dash.dependencies.State('input-D-m', 'value'),
    dash.dependencies.State('input-v-m', 'value'),
    dash.dependencies.State('input-halflife-p', 'value'),
    dash.dependencies.State('input-D-p', 'value'),
    dash.dependencies.State('input-v-p', 'value'),
    )
def update_output(n_clicks, valueboundaryCondition,value_uptakemRNA,value_uptakeProtein,value,valueL,halflife_m,Dm,vm,\
                                         halflife_p,Dp,vp):
    #fig.update_layout(clickmode='event+select')
    if((n_clicks != None) and  (n_clicks>0)):
        TR.params_input['name'] = value;
        TR.params_input['L'] = valueL;
        TR.params_input['halflife_m'] = float(halflife_m);#float(halflife_m),
        TR.params_input['D_m'] = float(Dm);
        TR.params_input['v_m'] = float(vm);
        TR.params_input['halflife_p'] = float(halflife_p);
        TR.params_input['D_p'] = float(Dp);
        TR.params_input['v_p'] = float(vp);
        TR.params_input['boundaryCondition'] = valueboundaryCondition;
        TR.params_input['uptake_mRNA'] = value_uptakemRNA;
        TR.params_input['uptake_proteins'] = value_uptakeProtein;

        TR.SolveTransportEqs(N);
        df = TR.GetSolution();
        #rows_results.append((TR.GetParametersAndResults(combineParameter=True)).copy());
        rows_results.append(TR.GetResults(resultsToDisplay));
        rows_params.append(TR.GetParameters(paramsToDisplay));
    else:
        df = TR.GetSolution();

    df = df.melt(id_vars=['x'],var_name="Quantity",value_name="Value")

    fig = px.line(df, x="x", y="Value",facet_row="Quantity",color="Quantity",title="Protein and mRNA concentrations (last computation):");
    fig.update_yaxes(matches=None);

    if(len(rows_results) > 0):
        df_res        = pd.DataFrame(rows_results);
        df_res["idx"] = np.asarray(df_res.index);

        for p in df_res:
            if(p != "idx"):
                df_res[p] = np.maximum(1e-6,df_res[p]);
        #df_res["Particles in active transport (per synaptic protein)"] = np.maximum(1e-6,np.asarray(df_res["Particles in active transport (per synaptic protein)"]));
        df_p   = df_res.melt(id_vars=["idx"],var_name="Measure",value_name="Value");
        fig2   = px.strip(df_p,x="Value",color="idx",y="Measure",log_x=True,title="Energies of all computations:");
        fig2.layout.update(showlegend=False,xaxis = dict(
            tickmode = 'array',
            tickvals = [1e-6,1e-4,1e-2,1],
            ticktext = ['<1e-6','1e-4','1e-2','1']
        ))

    else:
        df_res = pd.DataFrame([]);
        fig2   = px.strip(df_res,title="Energies of all computations (per synaptic protein):");

    df_params        = pd.DataFrame(rows_params);
    df_params['idx'] = df_params.index;
#    dictParams = {};#
    #dictParams['column1'] = 1.0;
    #dictParams['column2'] = 1.0;
    return fig,fig2,df_params.to_dict("records")#str(TR.params_input),format(df_res.to_string())#TR.params_input)#fig


if __name__ == '__main__':
    app.run_server(debug=True)
