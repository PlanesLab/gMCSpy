import os
import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc, dash_table
import dash_html_components as html
from dash.dependencies import Input, Output
import pandas as pd
import plotly.graph_objects as go
import networkx as nx
import base64

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP], suppress_callback_exceptions=True)

# Define the SVG content for the home page
# Specify the path to your SVG file
# print working directory

os.getcwd()

svg_file_path = './dashboard/assets/logoGMCS.PNG'
with open(svg_file_path, 'rb') as f:
    # Read the raw bytes of the SVG file
    raw_bytes = f.read()

    # Base64 encode the SVG file bytes for HTML embedding   
    logo_base64 = base64.b64encode(raw_bytes).decode('ascii')


# Define the folder path (replace this with your desired folder path)
folder_path = "./logs"

# Define the files we want to extract information from
details = []
ganttDataframes = {}
solutionsDataframes = {}
checksDataframes = {}

# Function to get the list of elements from a specified folder path (first level only)
def get_elements_from_path(folder_path):
    elements = []
    for entry in os.scandir(folder_path):
        if entry.is_dir():
            try:
                processFileAdress = f"{entry.path}/process/process.log"
                processFile = open(processFileAdress, "r")
                # Read the first line of the file
                details.append(processFile.readline().strip())
                ganntId = '/' + entry.name + '/'
                ganttDataframes[ganntId] = read_and_format_processfile(processFileAdress)
                solutionFileAdress = f"{entry.path}/solutions/solutions.log"
                solutionsDataframes[ganntId] = read_and_format_solutionsfile(solutionFileAdress)
                try:
                    checksFileAdress = f"{entry.path}/checks/checkedSolutions.csv"
                    checksDataframes[ganntId] = read_and_format_checksfile(checksFileAdress)
                except Exception as e:
                    checksDataframes[ganntId] = None
                    print(entry.name + 'No checks files')
            except FileNotFoundError:
                print(entry.name + 'No process file')
                continue
            elements.append(entry.name)
            
    return elements

# Function to parse details
def parse_key_element_pairs(text):
    pairs = {}
    elements = text.split(', ')
    for element in elements:
        key, value = element.split(': ')
        if value.isdigit():
            value = int(value)
        elif value.replace('.', '', 1).isdigit():
            value = float(value)
        elif value.lower() == 'true':
            value = True
        elif value.lower() == 'false':
            value = False
        pairs[key.strip()] = value
    return pairs

# read and format process file for gantt chart
def read_and_format_processfile(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    df_list = []
    prev_finish_time = 0
    # Skip the first line and read the rest
    for line in lines[1:]:
        try:
            name = line.split(':')[0]
            time = prev_finish_time + float(line.split(':')[1].strip())
        except:
            continue
        
        if name=='totalTime':
            time = time - prev_finish_time
            prev_finish_time = 0
        # Create a dictionary for each task
        task_dict = {
            'Task': name,
            'Start': prev_finish_time,
            'Finish': time
        }

        df_list.append(task_dict)
        prev_finish_time=time

    df = pd.DataFrame(df_list)
    return df

# read and format solutions file
def read_and_format_solutionsfile(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

        df_list = []
       
        for line in lines[1:]:
            try:
                geneSet = line.split(':')[0]
                reactionSet = line.split(':')[1].strip()
                geneSetkey = geneSet.strip('[]').split(',')
                geneSetkey = frozenset(geneSetkey)
            except:
                continue
            
            df_dict = {
                'GeneSetKey': geneSetkey,
                'GeneSet': geneSet,
                'Reactions': reactionSet,
                'GeneSetLength': len(geneSet.split(',')),
                'ReactionsLength': len(reactionSet.split(','))
            }

            df_list.append(df_dict)


        df = pd.DataFrame(df_list)
        return df

# read and format checks fiile 
def read_and_format_checksfile(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

        df_list = []
       
        for line in lines[1:]:
            try:
                isGMCS = line.split(',')[1]
                indexC =   line.index('{')
                checkDict = line[indexC:]
                checkDict = checkDict.strip()
                checkDict = eval(checkDict)
                geneSet = checkDict['solution']
                isCutSet = checkDict['isCutSet']
                cutSetValue = checkDict['cutSetValue']
                isMinimal = checkDict['isMinimal']
                setGeneTested = checkDict['setGeneTested']
                objVals = checkDict['objVals']
                
                df_dict = {
                    'GeneSetKey': geneSet,
                    'isGMCS': isGMCS,
                    'isCutSet': isCutSet,
                    'cutSetValue': cutSetValue,
                    'isMinimal': isMinimal,
                    'setGeneTested': setGeneTested,
                    'objVals': objVals
                }
                df_list.append(df_dict)


            except:
                continue
            
        df = pd.DataFrame(df_list)
            
            
    return df

# Create a card layout for each element with clickable links
def create_card_layout(elements, details):
    cards = []
    for element in elements:
        card_content = [
            dbc.CardHeader(element),
            dbc.CardBody([
                dbc.FormText(details[elements.index(element)]),
            ]),
        ]
        # Make each card clickable and link to the corresponding detail page
        card = dbc.Card([dbc.CardLink(card_content[0], href=f"{element}/"), dbc.CardBody(card_content[1])], className="my-3")
        cards.append(card)
    return cards



# Get the list of elements from the specified folder
elements = get_elements_from_path(folder_path)


data = {
    'GeneSet': [],
    'Reactions': [],
    'GeneSetLength': [],
    'ReactionsLength': [],
    'Checked': [],
}

df = pd.DataFrame(data)



# Create the overall layout of the app

app.layout = html.Div(
    className="container",
    children=[
        dcc.Store(id='solDf-store'),
        dcc.Location(id="url", refresh=False),
        html.H1(children=[
        html.Center(
            html.Img(src='data:image/png;base64,{}'.format(logo_base64), style={'height':'50%', 'width':'50%'})
        )
    ]),
        html.Div(id="page-content"),
        html.Div(dcc.Checklist(
                id='checkbox-checks'), style={'display': 'none'})
    ],
)

def updateRadioCheck(value):
    
    return dcc.Checklist(
                id='checkbox-checks',
                options=[
                    {'label': '  Include Checks Process', 'value': 'includeChecks'}
                ],
                value=value
            )


# Callback to render content based on the URL
@app.callback(Output("page-content", "children"), 
              [Input("url", "pathname"),
               Input('checkbox-checks', 'value')])
def render_page_content(pathname, includeChecks):
    if pathname is None:
        return "404 - Page not found"

    if pathname.startswith("/gMCSpy_"):
        # Extract the selected directory name from the URL
        modelName = pathname.split("gMCSpy")[1]
        modelName = modelName.split("gurobi")[0]
        # Replace this part with code to generate the scatter plots or detailed information
        ganttDf = ganttDataframes[pathname]
        print(includeChecks)
        if includeChecks is None:
            includeChecks = []
        radioChecks = updateRadioCheck(includeChecks)
        ganttFig = createGantt(ganttDf, includeChecks)
        solDf = solutionsDataframes[pathname]
        checkDf = checksDataframes[pathname]
        if checkDf is not None:
            mergedDf = pd.merge(solDf, checkDf, how='left', on=['GeneSetKey'])
        else:
            mergedDf = solDf
            mergedDf['isGMCS'] = 'Not Checked'
        # select only the columns we want to show in the table
        solDf = mergedDf[['GeneSet', 'GeneSetKey', 'Reactions', 'GeneSetLength', 'ReactionsLength', 'isGMCS']]
        tableLayout = update_layout(solDf)
        
        
        
        return html.Div([
            dcc.Store(id='solDf-store', data=solDf.to_dict('records')),
            dbc.Row([
                dbc.Col(dbc.Button("Back", id="back-button", color="primary", href="/"), width="auto"),
                dbc.Col(html.H3(f"Model: {modelName}"), width=True),
            ]),
            # Add your scatter plots or detailed information here
            dbc.Row([
                dbc.Col(dcc.Graph(
                    id='timeline-plot',
                    figure=ganttFig
                ), width=10),
                dbc.Col(radioChecks, width=2, align='center'),
            ]),
        ]), tableLayout

    # By default, return the home page content (list of cards)
    return *create_card_layout(elements, details), html.Div(dcc.Checklist(id='checkbox-checks'), style={'display': 'none'})



def createGantt(ganttDf, includeChecks=False):
        # Groupped by task and sum the Finishtime, keep the min start time
        # Group by 'Task', sum the 'Finish' times, and keep the minimum 'Start' time
        ganttDf = ganttDf.groupby('Task').agg({'Start': 'min', 'Finish': 'max'}).reset_index()
        # Define the specific tasks that should be the first and last
        first_task = 'GMatrix'
        last_task = 'checks'
        # Filter out the rows corresponding to the first and last tasks
        first_task_row = ganttDf[ganttDf['Task'] == first_task]
        last_task_row = ganttDf[ganttDf['Task'] == last_task]
        # Remove the first and last tasks from the original dataframe
        ganttDf = ganttDf[ganttDf['Task'] != first_task]
        ganttDf = ganttDf[ganttDf['Task'] != last_task]
        # Concatenate the first task row, the main dataframe, and the last task row in the desired order
        if 'includeChecks' in includeChecks:
            ganttDf = pd.concat([first_task_row, ganttDf, last_task_row])
        else:
            ganttDf = pd.concat([first_task_row, ganttDf])
        # Reset the index to have continuous integer index values
        ganttDf = ganttDf.reset_index(drop=True)
        fig = go.Figure()
        # Iterate through each row in the DataFrame and add a horizontal bar trace
        for index, row in ganttDf[::-1].iterrows():
            fig.add_trace(go.Bar(
                y=[row['Task']],
                x=[row['Finish'] - row['Start']],
                base=[row['Start']],
                name=row['Task'],
                orientation='h'
            ))
        # Update layout
        fig.update_layout(
            title='Task Durations',
            xaxis_title='Duration',
            yaxis_title='Task',
            xaxis_tickformat='.3f',  # Format X-axis to show up to 3 decimal places
            showlegend=False
        )

        # Update layout
        fig.update_layout(
            title='Calculation process',
        )

        # Adjust figure margins and height for better alignment
        height_per_task = 100  # You can adjust this value as needed
        fig.update_layout(
            margin=dict(l=10, r=20, t=100),
            height=len(ganttDf) * height_per_task  # Adjust the height based on the number of tasks
        )
        return fig


def createSolTable(df):
   # Convert lists in the 'GeneSet' and 'Reactions' columns to comma-separated strings
    df['GeneSet'] = df['GeneSet'].apply(lambda x: ', '.join(x))
    df['Reactions'] = df['Reactions'].apply(lambda x: ', '.join(x))

    maxGeneSetLength = df['GeneSetLength'].max()
    maxReactionsLength = df['ReactionsLength'].max() 
    

def update_layout(dataframe):
    # Create a new layout with updated components
    genesDropdownValues = list(set().union(*dataframe.GeneSetKey))
    dataframe.drop(['GeneSetKey'], axis=1, inplace=True)
    reactionStep = 1
    if dataframe['ReactionsLength'].max() > 20:
       reactionStep = 20 
    updated_layout = dbc.Container([
        dbc.Row([
            dbc.Col([
                dcc.Dropdown(
                    id='geneSet_drop',
                    options=[{'label': x, 'value': x} for x in sorted(genesDropdownValues)],
                    multi=True
                )
            ], width=2),
            dbc.Col([
                dcc.Dropdown(
                    id='reactions_drop',
                    options=[{'label': x, 'value': x} for x in sorted(dataframe.Reactions.unique())],
                    multi=True
                )
            ], width=2),
            dbc.Col([
                dcc.Slider(
                    id='gene_slider',
                    min=1,
                    max=dataframe['GeneSetLength'].max(),
                    value=1,
                    step=1,
                    tooltip={"placement": "bottom", "always_visible": True}
                )
            ], width=2),
            dbc.Col([
                dcc.Slider(
                    id='react_slider',
                    min=1,
                    max=dataframe['ReactionsLength'].max(),
                    value=1,
                    step=reactionStep,
                    tooltip={"placement": "bottom", "always_visible": True}
                )
            ], width=2),
            dbc.Col([], width=2),
        ], justify="between", className='mt-3 mb-4'),

        dash_table.DataTable(
            id='gene-table',
            columns=[
                {'name': 'GeneSet', 'id': 'GeneSet'},
                {'name': 'Reactions', 'id': 'Reactions'},
                {'name': 'GeneSetLength', 'id': 'GeneSetLength'},  
                {'name': 'ReactionsLength', 'id': 'ReactionsLength'},
                {'name': 'isGMCS', 'id': 'isGMCS'},
            ],
            data=dataframe.to_dict('records'),
            page_size=20,

            style_data={
                'width': '150px', 'minWidth': '150px', 'maxWidth': '150px',
                'overflow': 'hidden',
                'textOverflow': 'ellipsis',
            }
        ),
    ], id='table-layout', style={'display': 'block', 'paddingBottom': '20px'})
    
    return updated_layout


@app.callback(
    Output('gene-table', 'data'),
    Input('gene_slider', 'value'),
    Input('react_slider', 'value'),
    Input('solDf-store', 'data'),
    Input('geneSet_drop', 'value'),
    Input('reactions_drop', 'value')
)
def update_table_data(gene_value, react_value, dfStored, geneSet_value, reactions_value):
    df = pd.DataFrame(dfStored)
    filtered_df = df[(df['GeneSetLength'] >= gene_value) & (df['ReactionsLength'] >= react_value)]
        
    if geneSet_value:
        if len(geneSet_value) == 1:
            gene_value = str(geneSet_value).strip('[]').strip("'")
            # Check the column and check if the value is in the list
            filtered_df = filtered_df[filtered_df['GeneSet'].str.contains(gene_value)]
        elif len(geneSet_value) > 1:
            emptyDf = pd.DataFrame()
            for gene in geneSet_value:
                interDf = filtered_df[filtered_df['GeneSet'].str.contains(str(gene))]
                if emptyDf.empty:
                    emptyDf = interDf
                else:
                    emptyDf = emptyDf.merge(interDf, how='inner')
            filtered_df = emptyDf
    if reactions_value:
        filtered_df = filtered_df[filtered_df['Reactions'].isin(reactions_value)]
        
    return filtered_df.to_dict('records')



if __name__ == "__main__":
    app.run_server(debug=True)
