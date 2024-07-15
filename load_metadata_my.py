from dash import Dash, html, dcc, Input, Output, callback
import plotly.express as px
import pandas as pd
import numpy as np
import base64
from skimage import io
import os

protoDUNElogo = '/exp/dune/app/users/gvittist/offline_dqm/work/protoDUNE_logo.png'
blank = '/nashome/x/xning/blank.png'
path='/nashome/x/xning/Pictures/debug/'

# external_stylesheets = ['https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css']
# app = Dash(external_stylesheets=external_stylesheets)


def b64_image(image_filename):
    with open(image_filename, 'rb') as f:
        image = f.read()
    return 'data:image/png;base64,' + base64.b64encode(image).decode('utf-8')


def list_metadata():
    metadata_array=[]
    #metadata_array = np.load("/nashome/x/xning/runinfo.npz")
    # data = np.load("/nashome/x/xning/runinfo2.npz")
    # metadata_array.append(data)
    #for x in os.listdir("./data_to_server/"):
    for x in os.listdir("/nashome/x/xning/runinfo/"):
       if x.endswith(".npz"):
           data = np.load("/nashome/x/xning/runinfo/"+x)
           metadata_array.append(data)

    return metadata_array


app = Dash(__name__)

colors = {
    'background': '#111111',
    'text': '#7FDBFF'
}

app.layout = html.Div([

    html.Div(style={'backgroundColor': colors['background']}, children=[
        html.Img(src=b64_image(protoDUNElogo),
        style={'width': '30%', 'float': 'center'})
    ]),
    
    html.Div(style={'backgroundColor': colors['background']}, children=[

        html.H1(children='Number of Ionization Charge at Anode Planes', style={
        'textAlign': 'center',
        'color': colors['text']
        })
    ]),

    html.Div(style={'backgroundColor': colors['background']}, children=[
    
        html.Div([
            'Run',
            dcc.Dropdown(
                id='run-number',
                style={'color' : colors['background'],'width':'100%'}
            ),
            'Event',
            dcc.Dropdown(
                id='event-number',
                style={'color' : colors['background'],'width':'100%'}
            ),


            dcc.Graph(
                id='waveform_ori',
            ),
            dcc.Graph(
                id='waveform_ANf',
            ),

            dcc.Graph(
                id='baseline',
            ),

            dcc.Graph(
                id='rms',
            ),

            html.Div(className='row',children='Which Plane'),
            dcc.RadioItems(options=['u', 'v', 'w'], value='u', id='plane', inline=True),
            html.Div(className='row',children='Which APA'),
            dcc.RadioItems(options=['1', '2', '3','4'], value='1', id='APA',inline=True),

            dcc.Graph(
                id='wf_sep',
            ),
            dcc.Graph(
                id='wf_sep_ANf',
            ),            
            dcc.Graph(
                id='fft',
            ),
            dcc.Graph(
                id='cov_ori',
            ),
            dcc.Graph(
                id='cov_ANf',
            ),

            dcc.Interval(
            id='interval-component',
            interval=1*60000, # in milliseconds
            n_intervals=0
            )     
             
        #], style={'color': colors['text'], 'width': '100%', 'display': 'inline-block', 'marginTop': '20px'}),
        ], style={'color': colors['text'],'width':'100%', 'marginTop': '20px','fontSize': '20px'}),

     ])
])

@callback(
    Output('run-number', 'options'),
    Input('interval-component', 'n_intervals')
    )
def update_run_bar(n):
    
    all_data = list_metadata()

    all_run = [rn['RunNumber'] for rn in all_data]

  #  for j in range(0, len(all_data)):
   #     all_run.append(all_data[j]["RunNumber"])

    return all_run

@callback(
    Output('event-number', 'options'),
    Input('run-number', 'value'),
    Input('interval-component', 'n_intervals')
    )
def update_event_bar(which_run,n):

    all_data = list_metadata()

    event_list=[]

    for ev in all_data:
        if ev['RunNumber'] == which_run:
            event_list = ev['EventNumbers']

   # event_list=np.array(event_list).astype(str)

    return event_list

@callback(
    Output('waveform_ori', 'figure'),
    Input('run-number', 'value'),
    Input('event-number', 'value')
    )
def update_store_data(which_run, which_event):

    runname = str(which_run)
    eventname = str(which_event)

    if (runname == "None") or (eventname == "None"):
        
        pathtoimage = protoDUNElogo
    
        img=io.imread(pathtoimage)

        fig = px.imshow(img)
        fig.update_layout(coloraxis_showscale=False, 
        width=1700, height=600, margin=dict(
            l=0,
            r=0,
            b=0,
            t=0,
            pad=0
        ))
        fig.update_xaxes(showticklabels=False)
        fig.update_yaxes(showticklabels=False)
        
        return fig
    
    else:
        #pathtoimage = "/exp/dune/app/users/gvittist/offline_dqm/work/storage/notcompressed/run_"+runname+"_event_"+eventname+".png"
        pathtoimage = path+"raw_"+runname+"_"+eventname+".png"
    
        img=io.imread(pathtoimage)

        fig = px.imshow(img)
        fig.update_layout(coloraxis_showscale=False, 
        width=1650, height=630, margin=dict(
            l=0, r=0, b=0, t=0, pad=0))
        fig.update_xaxes(showticklabels=False)
        fig.update_yaxes(showticklabels=False)

        return fig

@callback(
    Output('waveform_ANf', 'figure'),
    Input('run-number', 'value'),
    Input('event-number', 'value')
    )
def update_store_data1(which_run, which_event):

    runname = str(which_run)
    eventname = str(which_event)

    if (runname == "None") or (eventname == "None"):
        
        pathtoimage = blank
    
        img=io.imread(pathtoimage)

        fig = px.imshow(img)
        fig.update_layout(coloraxis_showscale=False, 
        width=1700, height=600, margin=dict(
            l=0,
            r=0,
            b=0,
            t=0,
            pad=0
        ))
        fig.update_xaxes(showticklabels=False)
        fig.update_yaxes(showticklabels=False)
        
        return fig
    
    else:
        #pathtoimage = "/exp/dune/app/users/gvittist/offline_dqm/work/storage/notcompressed/run_"+runname+"_event_"+eventname+".png"
        pathtoimage = path+"raw_ANf_"+runname+"_"+eventname+".png"
    
        img=io.imread(pathtoimage)

        fig = px.imshow(img)
        fig.update_layout(coloraxis_showscale=False, 
        width=1650, height=630, margin=dict(
            l=0, r=0, b=0, t=0, pad=0))
        fig.update_xaxes(showticklabels=False)
        fig.update_yaxes(showticklabels=False)

        return fig

@callback(
    Output('baseline', 'figure'),
    Input('run-number', 'value'),
    Input('event-number', 'value')
    )
def update_store_data(which_run, which_event):

    runname = str(which_run)
    eventname = str(which_event)

    if (runname == "None") or (eventname == "None"):
        
        pathtoimage = blank
    
        img=io.imread(pathtoimage)

        fig = px.imshow(img)
        fig.update_layout(coloraxis_showscale=False, 
        width=1700, height=600, margin=dict(
            l=0,
            r=0,
            b=0,
            t=0,
            pad=0
        ))
        fig.update_xaxes(showticklabels=False)
        fig.update_yaxes(showticklabels=False)
        
        return fig
    
    else:
        #pathtoimage = "/exp/dune/app/users/gvittist/offline_dqm/work/storage/notcompressed/run_"+runname+"_event_"+eventname+".png"
        pathtoimage = path+"baseline_"+runname+"_"+eventname+".png"
    
        img=io.imread(pathtoimage)

        fig = px.imshow(img)
        fig.update_layout(coloraxis_showscale=False, 
        width=1650, height=630, margin=dict(
            l=0, r=0, b=0, t=0, pad=0))
        fig.update_xaxes(showticklabels=False)
        fig.update_yaxes(showticklabels=False)

        return fig

@callback(
    Output('rms', 'figure'),
    Input('run-number', 'value'),
    Input('event-number', 'value')
    )
def update_store_data2(which_run, which_event):

    runname = str(which_run)
    eventname = str(which_event)

    if (runname == "None") or (eventname == "None"):
        
        pathtoimage = blank
    
        img=io.imread(pathtoimage)

        fig = px.imshow(img)
        fig.update_layout(coloraxis_showscale=False, 
        width=1700, height=600, margin=dict(
            l=0,
            r=0,
            b=0,
            t=0,
            pad=0
        ))
        fig.update_xaxes(showticklabels=False)
        fig.update_yaxes(showticklabels=False)
        
        return fig
    
    else:
        #pathtoimage = "/exp/dune/app/users/gvittist/offline_dqm/work/storage/notcompressed/run_"+runname+"_event_"+eventname+".png"
        pathtoimage = path+"rms_"+runname+"_"+eventname+".png"
    
        img=io.imread(pathtoimage)

        fig = px.imshow(img)
        fig.update_layout(coloraxis_showscale=False, 
        width=1650, height=630, margin=dict(
            l=0, r=0, b=0, t=0, pad=0))
        fig.update_xaxes(showticklabels=False)
        fig.update_yaxes(showticklabels=False)

        return fig

@callback(
    Output('fft', 'figure'),
    Input('run-number', 'value'),
    Input('event-number', 'value'),
    Input('plane', 'value'),
    Input('APA', 'value')
    )
def update_store_data3(which_run, which_event,which_plane,which_APA):

    runname = str(which_run)
    eventname = str(which_event)
    APA = str(which_APA)
    Plane = str(which_plane)

    if (runname == "None") or (eventname == "None"):
        
        pathtoimage = blank
    
        img=io.imread(pathtoimage)

        fig = px.imshow(img)
        fig.update_layout(coloraxis_showscale=False, 
        width=1700, height=600, margin=dict(
            l=0,
            r=0,
            b=0,
            t=0,
            pad=0
        ))
        fig.update_xaxes(showticklabels=False)
        fig.update_yaxes(showticklabels=False)
        
        return fig
    
    else:
        #pathtoimage = "/exp/dune/app/users/gvittist/offline_dqm/work/storage/notcompressed/run_"+runname+"_event_"+eventname+".png"
        pathtoimage = path+"fft"+APA+Plane+"_"+runname+"_"+eventname+".png"
    
        img=io.imread(pathtoimage)

        fig = px.imshow(img)
        fig.update_layout(coloraxis_showscale=False, 
        width=1650, height=630, margin=dict(
            l=0, r=0, b=0, t=0, pad=0))
        fig.update_xaxes(showticklabels=False)
        fig.update_yaxes(showticklabels=False)

        return fig

@callback(
    Output('cov_ori', 'figure'),
    Input('run-number', 'value'),
    Input('event-number', 'value'),
    Input('plane', 'value'),
    Input('APA', 'value')
    )
def update_store_data4(which_run, which_event,which_plane,which_APA):

    runname = str(which_run)
    eventname = str(which_event)
    APA = str(which_APA)
    Plane = str(which_plane)

    if (runname == "None") or (eventname == "None"):
        
        pathtoimage = blank
    
        img=io.imread(pathtoimage)

        fig = px.imshow(img)
        fig.update_layout(coloraxis_showscale=False, 
        width=1700, height=600, margin=dict(
            l=0,
            r=0,
            b=0,
            t=0,
            pad=0
        ))
        fig.update_xaxes(showticklabels=False)
        fig.update_yaxes(showticklabels=False)
        
        return fig
    
    else:
        #pathtoimage = "/exp/dune/app/users/gvittist/offline_dqm/work/storage/notcompressed/run_"+runname+"_event_"+eventname+".png"
        pathtoimage = path+"CovMatx_"+APA+Plane+"_"+runname+"_"+eventname+".png"
    
        img=io.imread(pathtoimage)

        fig = px.imshow(img)
        fig.update_layout(coloraxis_showscale=False, 
        width=1650, height=630, margin=dict(
            l=0, r=0, b=0, t=0, pad=0))
        fig.update_xaxes(showticklabels=False)
        fig.update_yaxes(showticklabels=False)

        return fig

@callback(
    Output('cov_ANf', 'figure'),
    Input('run-number', 'value'),
    Input('event-number', 'value'),
    Input('plane', 'value'),
    Input('APA', 'value')
    )
def update_store_data5(which_run, which_event,which_plane,which_APA):

    runname = str(which_run)
    eventname = str(which_event)
    APA = str(which_APA)
    Plane = str(which_plane)

    if (runname == "None") or (eventname == "None"):
        
        pathtoimage = blank
    
        img=io.imread(pathtoimage)

        fig = px.imshow(img)
        fig.update_layout(coloraxis_showscale=False, 
        width=1700, height=600, margin=dict(
            l=0,
            r=0,
            b=0,
            t=0,
            pad=0
        ))
        fig.update_xaxes(showticklabels=False)
        fig.update_yaxes(showticklabels=False)
        
        return fig
    
    else:
        #pathtoimage = "/exp/dune/app/users/gvittist/offline_dqm/work/storage/notcompressed/run_"+runname+"_event_"+eventname+".png"
        pathtoimage = path+"cov_ANfMatx_"+APA+Plane+"_"+runname+"_"+eventname+".png"
    
        img=io.imread(pathtoimage)

        fig = px.imshow(img)
        fig.update_layout(coloraxis_showscale=False, 
        width=1650, height=630, margin=dict(
            l=0, r=0, b=0, t=0, pad=0))
        fig.update_xaxes(showticklabels=False)
        fig.update_yaxes(showticklabels=False)

        return fig

@callback(
    Output('wf_sep', 'figure'),
    Input('run-number', 'value'),
    Input('event-number', 'value'),
    Input('plane', 'value'),
    Input('APA', 'value')
    )
def update_store_data3(which_run, which_event,which_plane,which_APA):

    runname = str(which_run)
    eventname = str(which_event)
    APA = str(which_APA)
    Plane = str(which_plane)

    if (runname == "None") or (eventname == "None"):
        
        pathtoimage = blank
    
        img=io.imread(pathtoimage)

        fig = px.imshow(img)
        fig.update_layout(coloraxis_showscale=False, 
        width=1700, height=600, margin=dict(
            l=0,
            r=0,
            b=0,
            t=0,
            pad=0
        ))
        fig.update_xaxes(showticklabels=False)
        fig.update_yaxes(showticklabels=False)
        
        return fig
    
    else:
        #pathtoimage = "/exp/dune/app/users/gvittist/offline_dqm/work/storage/notcompressed/run_"+runname+"_event_"+eventname+".png"
        pathtoimage = path+"raw_apa"+APA+"_"+Plane+"_"+runname+"_"+eventname+".png"
    
        img=io.imread(pathtoimage)

        fig = px.imshow(img)
        fig.update_layout(coloraxis_showscale=False, 
        width=1650, height=630, margin=dict(
            l=0, r=0, b=0, t=0, pad=0))
        fig.update_xaxes(showticklabels=False)
        fig.update_yaxes(showticklabels=False)

        return fig

@callback(
    Output('wf_sep_ANf', 'figure'),
    Input('run-number', 'value'),
    Input('event-number', 'value'),
    Input('plane', 'value'),
    Input('APA', 'value')
    )
def update_store_data3(which_run, which_event,which_plane,which_APA):

    runname = str(which_run)
    eventname = str(which_event)
    APA = str(which_APA)
    Plane = str(which_plane)

    if (runname == "None") or (eventname == "None"):
        
        pathtoimage = blank
    
        img=io.imread(pathtoimage)

        fig = px.imshow(img)
        fig.update_layout(coloraxis_showscale=False, 
        width=1700, height=600, margin=dict(
            l=0,
            r=0,
            b=0,
            t=0,
            pad=0
        ))
        fig.update_xaxes(showticklabels=False)
        fig.update_yaxes(showticklabels=False)
        
        return fig
    
    else:
        #pathtoimage = "/exp/dune/app/users/gvittist/offline_dqm/work/storage/notcompressed/run_"+runname+"_event_"+eventname+".png"
        pathtoimage = path+"raw_ANf_apa"+APA+"_"+Plane+"_"+runname+"_"+eventname+".png"
    
        img=io.imread(pathtoimage)

        fig = px.imshow(img)
        fig.update_layout(coloraxis_showscale=False, 
        width=1650, height=630, margin=dict(
            l=0, r=0, b=0, t=0, pad=0))
        fig.update_xaxes(showticklabels=False)
        fig.update_yaxes(showticklabels=False)

        return fig

if __name__ == '__main__':
    #app.run(debug=True,port=8051)
    app.run(debug=True)
