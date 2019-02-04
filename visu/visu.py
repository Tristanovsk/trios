from plotly.graph_objs import *
import matplotlib as mpl
import numpy as np

import os
import pandas as pd
import numpy as np
import flask
import json
import base64
import datetime
import io
import plotly
import plotly.plotly as py
import plotly.graph_objs as go
import matplotlib as mpl
from textwrap import dedent as d

import re

import dash
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
import dash_table_experiments as dte
from datetime import datetime

app = dash.Dash()
# app.css.append_css('data/aeronet_layout.css')
styles = {
    'pre': {
        'border': 'thin lightgrey solid',
        'overflowX': 'scroll'
    }
}

# ------------------------------------------------------
# layout section
# ------------------------------------------------------

app.layout = html.Div([
    html.Div([
        html.H1(
            'Above-water data visualization',
            className='eight columns',
            style={'display': 'inline-block'}),

        #   file selection box
        dcc.Upload(
            id='upload-data',
            children=html.Div([
                'Drag and Drop or ',
                html.A('Select Files')
            ]),
            style={
                'width': '40%',
                'height': '60px',
                'lineHeight': '60px',
                'borderWidth': '1px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'center',
                'margin': '10px',
                'float': 'right'
            },
            # Allow multiple files to be uploaded
            multiple=False)],
        style={'margin-top': '0',
               'width': '100%', 'display': 'inline-block',
               }),
    html.Div([
        html.H4('File...', id='filename',style={'float': 'left','width': '60%'}),
        html.Div([
            html.H4('Color variable:',style={'margin-bottom': '0','width':'50%'}),
            html.Div([
                dcc.Dropdown(
                    id='color-column',
                    value='sza',

                ),
                # dcc.RadioItems(
                #     id='color-type',
                #     options=[{'label': i, 'value': i} for i in ['Linear', 'Log']],
                #     value='Linear',
                #     labelStyle={'width': '25%','display': 'inline-block'})],
            ],
            style={'width': '50%','display': 'inline-block'})],
            style={'width': '40%',
                   'margin-top': '0',
                   'display': 'inline-block',
                   'margin-left': '0%',
                   'float': 'right'})],
        style={'margin-block-end': '7%'}),

    # Spectrum graphs
    html.Div([
        html.Div([

            html.H4('Spectral parameter 1'),
            dcc.Dropdown(
                id='spectrum1',
                value='Ed'),
        ],
            style={'width': '48.9%',
                   'float': 'left', }),

        html.Div([

            html.H4('Spectral parameter 2'),
            dcc.Dropdown(
                id='spectrum2',
                value='Lsky'),
        ],
            style={'width': '48.9%',
                   'float': 'right', }),],
        style={'width':'100%','margin-block-start': '1%'}),


        html.Div([
            dcc.Graph(id='graph1')],
            # className='eight columns',
            style={'width': '59.9%',
                   'margin-top': '0',
                   'display': 'inline-block',

                   }),
        html.Div([
            dcc.Graph(id='graph2')],
            # className='eight columns',
            style={'width': '59.9%',
                   'float': 'left',
                   'margin-top': '1%',
                   'display': 'inline-block',

                   }),
         # className='row'


    # Spectrum graphs
    html.Div([
        html.Div([
            dcc.Graph(id='graph3')],
            # className='eight columns',
            style={'width': '59.9%',
                   'margin-top': '1%',
                   'display': 'inline-block',

                   }),

    ],

        # style={'display': 'inline-block'},
        # className='row'
    ),

    html.Div([

        html.H4('Spectral parameter'),
        dcc.Dropdown(
            id='spectrum3',
            value='Lt'),
    ],
        style={'width': '48.9%',
               'float': 'left', }),


    html.Div([
        dcc.Markdown(d("""
                **Selection Data**

                Choose the lasso or rectangle tool in the graph's menu
                bar and then select points in the graph.
            """)),
        # html.Pre(id='selected-data', style=styles['pre']),
    ], className='three columns'),

    # hidden signal value
    html.Div(id='dataset', style={'display': 'none'}),
],

    style={
        'width': '90%',
        'fontFamily': 'Sans-Serif',
        'margin-left': 'auto',
        'margin-right': 'auto'})




def figure_spectrum_v1(df, column_name, color_column_name):
    dff = df
    parameters = df.loc[:, (color_column_name)].values
    dff = dff.loc[:, ([column_name])]
    wl = dff.columns.droplevel().values
    dff = dff.stack(level=['wl'])
    norm = mpl.colors.Normalize(vmin=np.nanmin(parameters), vmax=np.nanmax(parameters))
    # create a ScalarMappable and initialize a data structure
    cmap = mpl.cm.Spectral
    s_m = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
    s_m.set_array([])
    i = 0
    trace = []
    for date, x in dff.groupby(level=0):
        trace.append(Scattergl(
            x=wl,  # spectrum,
            y=x[column_name].values,
            text="depth="+str(float(parameters[i])/100)+" m",#x.index.get_level_values(0),
            mode='lines',
            marker={
                'size': 7,
                'opacity': 0.5,
                # 'color': 'rgba({}, {}, {}, {})'.format(*s_m.to_rgba(parameters[i]).flatten()),
                # x.unique(),#color': df.index.get_level_values(0),
                'line': {'width': 0.5, 'color': 'white'},
            },
            line=Line(color='rgba({}, {}, {}, {})'.format(*s_m.to_rgba(parameters[i])), width=2),
            showlegend=False
        ))
        i = i + 1

    # spectrum = df[label['aod']].stack()
    return {
        'data': trace,
        'layout': Layout(
            xaxis={
                'title': 'Wavelength (nm)',

            },
            yaxis={
                'title': column_name,

            },
            margin={'l': 50, 'b': 40, 't': 20, 'r': 50},
            hovermode='closest',

            height=400,
            font=dict(color='#CCCCCC'),
            titlefont=dict(color='#CCCCCC', size='14'),

            plot_bgcolor="#191A1A",
            paper_bgcolor="#020202",
        )
    }


def figure_spectrum(df, column_name, color_column_name='depth'):


    dff = df
    if color_column_name=='depth':
        parameters = df.loc[:, (color_column_name)].values
        dff = dff.loc[:, ([column_name,color_column_name])]
        dff.set_index('depth', append=True, inplace=True)

    wl = dff.columns.droplevel().values
    dff = dff.stack(level=['wl'])
    # norm = mpl.colors.Normalize(vmin=np.nanmin(parameters), vmax=np.nanmax(parameters))
    # # create a ScalarMappable and initialize a data structure
    # cmap = mpl.cm.Spectral
    # s_m = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
    # s_m.set_array([])
    i = 0
    trace = []
    colors = ['blue',  'green','grey', 'yellow', 'orange', 'red', 'purple']


    group = np.array(dff.index.get_level_values(color_column_name),dtype=float)/100

    opts = []
    for i in range(0, len(colors)):
        opt = {'target': np.unique(group)[i], 'value': dict(marker=dict(color=colors[i]))}
        opts.append(opt)

    aggs = ["all","avg","median","mode","rms","stddev","min","max"]
    agg_func = []
    for i in range(0, len(aggs)):
        if i == 0:
            agg = dict(
                args=['transforms[0].aggregations[0].func',0],
                label=aggs[i],
                method='restyle'
            )
        else:
            agg = dict(
                args=[ 'transforms[0].aggregations[0].func',aggs[i]],
                label=aggs[i],
                method='restyle'
            )

        agg_func.append(agg)

    trace=[dict(
        type = 'scatter',
        x=dff.index.get_level_values('wl'),  # spectrum,
        y=dff[column_name].values,
        text=dff.index.get_level_values(0),#group,#x.index.get_level_values(0),
        hoverinfo = 'text',
        name=dff.index.get_level_values(0),
        mode='lines',
        marker={
            'size': 7,
            'opacity': 0.5,
            # 'color': 'rgba({}, {}, {}, {})'.format(*s_m.to_rgba(parameters[i]).flatten()),
            # x.unique(),#color': df.index.get_level_values(0),
            'line': {'width': 0.5, 'color': 'white'},
        },
        #line=Line(color='rgba({}, {}, {}, {})'.format(*s_m.to_rgba(parameters[i])), width=2),
        #showlegend=False,

        transforms = [


#            {'type': 'groupby', 'groups': dff.index.get_level_values(0),'showlegend': False},

            {'type': 'aggregate', 'groups': dff.index.get_level_values('wl'), 'aggregations': [
             dict(target='y', func='avg', enabled = True),]},
            {'type': 'groupby', 'groups': group, 'styles': opts},
         ]
    )]

    # spectrum = df[label['aod']].stack()
    return {
        'data': trace,
        'layout':  dict(
            xaxis={
                'title': 'Wavelength (nm)',

            },
            yaxis={
                'title': column_name,

            },
            margin={'l': 50, 'b': 40, 't': 20, 'r': 50},
            hovermode='closest',

            height=400,
            font=dict(color='#CCCCCC'),
            titlefont=dict(color='#CCCCCC', size='14'),

            plot_bgcolor="#191A1A",
            paper_bgcolor="#020202",
            updatemenus = [dict(
                x = 0.85,
                y = 1.15,
                xref = 'paper',
                yref = 'paper',
                yanchor = 'top',
                active = 0,
                showactive = True,
                buttons = agg_func
  )]
        )
    }


def figure_profile(df, var_col='Luz', depth_col='depth_Luz'):


    dff = df

    dff = dff.loc[:, ([var_col,depth_col])]
    dff.set_index(depth_col, append=True, inplace=True)

    wl = dff.columns.droplevel().values
    dff = dff.stack(level=['wl'])
    norm = mpl.colors.Normalize(vmin=np.nanmin(wl), vmax=np.nanmax(wl))
    # create a ScalarMappable and initialize a data structure
    cmap = mpl.cm.Spectral
    s_m = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
    s_m.set_array([])
    i = 0
    trace = []
    colors = ['blue',  'green','grey', 'yellow', 'orange', 'red', 'purple']


    group = np.array(dff.index.get_level_values('wl'),dtype=int)

    opts = []
    for i in range(0, len(colors)):
        opt = {'target': np.unique(group)[i], 'value': dict(marker=dict(color=colors[i]))}
        opts.append(opt)

    aggs = ["all","avg","median","mode","rms","stddev","min","max"]
    agg_func = []
    for i in range(0, len(aggs)):
        if i == 0:
            agg = dict(
                args=['transforms[0].aggregations[0].func',0],
                label=aggs[i],
                method='restyle'
            )
        else:
            agg = dict(
                args=[ 'transforms[0].aggregations[0].func',aggs[i]],
                label=aggs[i],
                method='restyle'
            )

        agg_func.append(agg)

    trace=[dict(
        type = 'scattergl',
        x=dff[var_col].values,  # spectrum,
        y=dff.index.get_level_values(1),
        text=dff.index.get_level_values(2),#group,#x.index.get_level_values(0),
        hoverinfo = 'text',
        name=dff.index.get_level_values(0),
        mode='markers',
        marker=dict(color=group,
                    showscale=True,size= 7,
                    colorscale='Jet',
                    opacity= 0.5),
        # marker={
        #     'size': 7,
        #     'opacity': 0.5,
        #     'color': 'rgba({}, {}, {}, {})'.format(*s_m.to_rgba(wl).flatten()),
        #     # x.unique(),#color': df.index.get_level_values(0),
        #     'line': {'width': 0.5, 'color': 'white'},
        # },
        #line=Line(color='rgba({}, {}, {}, {})'.format(*s_m.to_rgba(wl)), width=2),
        #showlegend=False,

#         transforms = [
#
#
# #            {'type': 'groupby', 'groups': dff.index.get_level_values(0),'showlegend': False},
#
#             # {'type': 'aggregate', 'groups': dff.index.get_level_values('wl'), 'aggregations': [
#             #  dict(target='y', func='avg', enabled = True),]},
#             {'type': 'groupby', 'groups': group, 'styles': opts},
#          ]
    )]

    # spectrum = df[label['aod']].stack()
    return {
        'data': trace,
        'layout':  dict(
            xaxis={
                'title': 'Wavelength (nm)',

            },
            yaxis={
                'title': column_name,

            },
            margin={'l': 50, 'b': 40, 't': 20, 'r': 50},
            hovermode='closest',

            height=400,
            font=dict(color='#CCCCCC'),
            titlefont=dict(color='#CCCCCC', size='14'),

            plot_bgcolor="#191A1A",
            paper_bgcolor="#020202",
            updatemenus = [dict(
                x = 0.85,
                y = 1.15,
                xref = 'paper',
                yref = 'paper',
                yanchor = 'top',
                active = 0,
                showactive = True,
                buttons = agg_func
  )]
        )
    }

@app.callback(Output('graph1', 'figure'),
              [Input('spectrum1', 'value'),
               Input('color-column', 'value'),
               Input('dataset', 'children')])
def spectrum_figure(column_name, color_column_name,void):

    return figure_spectrum(df, 'Edz','depth')

@app.callback(Output('graph2', 'figure'),
              [Input('spectrum2', 'value'),
               Input('color-column', 'value'),
               Input('dataset', 'children')])
def spectrum_figure(column_name, color_column_name,void):

    return figure_spectrum(df, 'Luz','depth')

@app.callback(Output('graph3', 'figure'),
              [Input('spectrum3', 'value'),
               Input('color-column', 'value'),
               Input('dataset', 'children')])
def spectrum_figure(column_name, color_column_name,void):

    return figure_spectrum(df, 'Ed','depth')

app.run_server()
