import pandas
from bioservices.kegg import KEGG
import networkx as nx
from visJS2jupyter import visJS_module
import math
import xml.etree.ElementTree as ET
import spectra
from copy import deepcopy
import mygene
mg = mygene.MyGeneInfo()

def pathwayVisualization(KEGG_id, path_to_csv, redirect=True, compound=False):
    """
    The pathwayVisualization function returns a graph visualization based on user input
    
    Args:
        KEGG_id     (str): string specifying KEGG pathway ID to visualize
        path_to_csv (str): string specifying data to overlay on graph
        redirect    (bool): True to split nodes into their components. Defaults to True
        compound    (bool): True to display compounds (such as Ca2+). Defaults to False
        
    Returns:
        A graph visualization using the visjs_network function from visjs_2_jupyter
    """
    
    s = KEGG()
    result = s.parse_kgml_pathway(KEGG_id)
    
    ETroot = parsingXML(KEGG_id, s)
    
    G=nx.DiGraph()
    
    max_id, compound_array = addNodes(G, result)
    setCoord(G, ETroot)
    
    if redirect is False:
        getNodeSymbols(G, s, compound)
    else:
        parent_list, parent_dict = splitNodes(G, s, max_id)
    
    complex_array, component_array, node_dict, comp_dict = undefNodes(G, ETroot)
    
    if redirect is False:
        addEdges(G, result, component_array, node_dict)
    else:
        addAndRedirectEdges(G, result, complex_array, component_array, parent_list, parent_dict, node_dict, comp_dict)
    
    #add reactions to graph
    addReaction(G, ETroot)
    
    edge_to_name = dict()
    for edge in G.edges():
        if G.edge[edge[0]][edge[1]]['name'] == 'phosphorylation':
            edge_to_name[edge] = G.edge[edge[0]][edge[1]]['value']
        elif G.edge[edge[0]][edge[1]]['name'] == 'dephosphorylation':
            edge_to_name[edge] = G.edge[edge[0]][edge[1]]['value']
        elif 'dephosphorylation' in G.edge[edge[0]][edge[1]]['name']:
            edge_to_name[edge] = G.edge[edge[0]][edge[1]]['name'].replace('dephosphorylation', '-p')
        elif 'phosphorylation' in G.edge[edge[0]][edge[1]]['name']:
            edge_to_name[edge] = G.edge[edge[0]][edge[1]]['name'].replace('phosphorylation', '+p')
        else:
            edge_to_name[edge] = G.edge[edge[0]][edge[1]]['name']
            
        edge_to_name[edge] = edge_to_name[edge].replace('activation, ', '')
        edge_to_name[edge] = edge_to_name[edge].replace('inhibition, ', '')
        edge_to_name[edge] = edge_to_name[edge].replace('activation', '')
        edge_to_name[edge] = edge_to_name[edge].replace('inhibition', '')

    #edges are transparent
    edge_to_color = dict()
    for edge in G.edges():
        if 'activation' in G.edge[edge[0]][edge[1]]['name']:
            edge_to_color[edge] = 'rgba(26, 148, 49, 0.3)' #green
        elif 'inhibition' in G.edge[edge[0]][edge[1]]['name']:
            edge_to_color[edge] = 'rgba(255, 0, 0, 0.3)' #red
        else:
            edge_to_color[edge] = 'rgba(0, 0, 255, 0.3)' #blue
    
    #for graph with split nodes
    if redirect is True:
        #remove undefined nodes from graph
        G.remove_nodes_from(complex_array)

        #remove nodes with more than one gene
        G.remove_nodes_from(parent_list)

    if compound is False:
        #remove compound nodes
        G.remove_nodes_from(compound_array)
        
    node_to_symbol = dict()
    for node in G.node:
        if G.node[node]['type'] == 'map':
            node_to_symbol[node] = G.node[node]['gene_names']
        else:
            if 'symbol' in G.node[node]:
                node_to_symbol[node] = G.node[node]['symbol']
            elif 'gene_names'in G.node[node]:
                node_to_symbol[node] = G.node[node]['gene_names']
            else: 
                node_to_symbol[node] = G.node[node]['name']
            
    # getting name of nodes
    node_to_gene = dict()
    for node in G.node:
        node_to_gene[node] = G.node[node]['gene_names']
            
    # getting x coord of nodes
    node_to_x = dict()
    for node in G.node:
        node_to_x[node] = G.node[node]['x']
    
    # getting y coord of nodes
    node_to_y = dict()
    for node in G.node:
        node_to_y[node] = G.node[node]['y']
    
    id_to_log2fold = log2FoldChange(G, path_to_csv)
    
    # Create color scale with negative as green and positive as red
    my_scale = spectra.scale([ "green", "#CCC", "red" ]).domain([ -4, 0, 4 ])
    
    # color nodes based on log2fold data
    node_to_color = dict()
    
    for node in G.nodes():

        if node in id_to_log2fold:
            node_to_color[node] = my_scale(id_to_log2fold[node][0]).hexcode

        else:
            node_to_color[node] = '#f1f1f1'

    # getting nodes in graph
    nodes = G.nodes()
    numnodes = len(nodes)
    node_map = dict(zip(nodes,range(numnodes)))  # map to indices for source/target in edges
    
    # getting edges in graph
    edges = G.edges()
    numedges = len(edges)

    # dictionaries that hold per node and per edge attributes
    nodes_dict = [{"id":node_to_gene[n],"degree":G.degree(n),"color":node_to_color[n], "node_shape":"box",
                 "node_size":10,'border_width':1, "id_num":node_to_symbol[n], "x":node_to_x[n], "y":node_to_y[n]} for n in nodes]

    edges_dict = [{"source":node_map[edges[i][0]], "target":node_map[edges[i][1]], 
                  "color":edge_to_color[edges[i]], "id":edge_to_name[edges[i]], "edge_label":'',
                 "hidden":'false', "physics":'true'} for i in range(numedges)]        

    # html file label for first graph (must manually increment later)
    time = 1700

    # create graph here
    #return G
    return visJS_module.visjs_network(nodes_dict, edges_dict, time_stamp = time, node_label_field = "id_num", 
                               edge_width = 3, border_color = "black", edge_arrow_to = True, edge_font_size = 15,
                               edge_font_align= "top", physics_enabled = False, graph_width = 1000, graph_height = 1000)
    
def parsingXML(KEGG_id, s):
    """
    The parsingXML function parses the KEGG metabolic pathway specified by the KEGG ID, and returns an element tree object
    
    Args:
        KEGG_id (str): string specifying KEGG pathway ID to visualize
        s       (KEGG): KEGG object with built in functions
        
    Returns:
        An element tree object containing the XML version of the KEGG metabolic pathway
    """
    kgml_data = s.get(KEGG_id, "kgml")

    #convert from unicode to string
    kgml_str = kgml_data.encode('utf-8')

    # convert to XML
    root = ET.fromstring(kgml_str)
    return root
    
def addNodes(G, res):
    max_id = 0
    compound_array = []
    # add nodes to networkx graph and keep track of the compound nodes
    for entry in res['entries']:
        G.add_node(entry['id'], entry)
        if entry['type'] == "compound":
            compound_array.append(entry['id'])
        if int(entry['id']) > int(max_id):
            max_id = int(entry['id'])

    return max_id, compound_array

def setCoord(G, root):
    # get x and y coordinates for each node
    seen_coord = set()
    for entry in root.findall('entry'):
        node_id = entry.attrib['id']
        graphics = entry.find('graphics')
        if (graphics.attrib['x'], graphics.attrib['y']) in seen_coord:
            G.node[node_id]['x'] = (int(graphics.attrib['x']) + .1) * 2.5
            G.node[node_id]['y'] = (int(graphics.attrib['y']) + .1) * 2.5
            seen_coord.add((G.node[node_id]['x'], G.node[node_id]['y']))

        else:
            seen_coord.add((graphics.attrib['x'], graphics.attrib['y']))
            G.node[node_id]['x'] = int(graphics.attrib['x']) * 2.5
            G.node[node_id]['y'] = int(graphics.attrib['y']) * 2.5

#get node symbols
def getNodeSymbols(G, s, compound):
    # get symbol of each node
    for node, data in G.nodes(data=True):

        if data['type'] == 'gene':
            # one gene only
            if ' ' not in data['name']:
                gene_symbol = data['gene_names'].split(',', 1)[0]
                G.node[node]['symbol'] = gene_symbol    # this value will be displayed upon hover
                G.node[node]['label'] = gene_symbol     # this value will be displayed on node
            #multiple genes
            else: 
                result = data['name'].split("hsa:")
                result = ''.join(result)
                result = result.split()
                
                # generate label and symbol of genes
                for index, gene in enumerate(result):
                    if index == 0:
                        gene_symbol = str(entrez_to_symbol(gene))
                        gene_label = gene_symbol

                        # determine family of the genes
                        family = ""
                        replace = True
                        for i in gene_symbol:
                            if i.isalpha():
                                family = "".join([family,i])

                    # check if we can use family as label
                    else:
                        gene_symbol = gene_symbol + ', ' + str(entrez_to_symbol(gene))
                        if replace == True:
                            if not str(entrez_to_symbol(gene)).startswith(family):
                                replace = False

                    # generate label to be used if genes do not have a family
                    if index < 3:
                        gene_label = gene_symbol
                    if index == 3:
                        gene_label = gene_label + '...'  
                        
                # assign value to be displayed upon hover
                G.node[node]['symbol'] = gene_symbol
                
                # assign value to be displayed on node
                if replace == True:
                    G.node[node]['label'] = family
                else:
                    G.node[node]['label'] = gene_label
                    
        if compound is True:
            if data['type'] == 'compound':
                if isinstance(s.parse(s.get(data['name'])), unicode):
                    compound_array = s.parse(s.get(data['name'])).split()
                    name_index = compound_array.index('NAME')
                    gene_symbol = [compound_array[name_index+1]]
                else:
                    gene_symbol = s.parse(s.get(data['name']))['NAME']
                G.node[node]['gene_names'] = ' '.join(gene_symbol)
                G.node[node]['symbol'] = gene_symbol[0].replace(';', '')
            
#get node symbols and split nodes with multiple genes
def splitNodes(G, s, max_id):
    # get symbol of each node
    parent_to_list = dict()
    parent_list = []
    temp_node_id_array = []
    for node, data in G.nodes(data=True):

        if data['type'] == 'gene':
            # one gene only
            if ' ' not in data['name']:
                G.node[node]['symbol'] = data['gene_names'].split(',', 1)[0]
            # multiple genes
            else: 
                parent_list.append(node)
                gene_list = []
                result = data['name'].split("hsa:")
                result = ''.join(result)
                result = result.split()   # array of genes
                for index, gene in enumerate(result):
                    max_id = max_id + 1
                    gene_list.append(str(max_id))
                    temp_data = deepcopy(data)
                    temp_data['name'] = deepcopy(gene)
                    temp_data['symbol'] = deepcopy(str(entrez_to_symbol(gene)))
                    temp_data['parent_id'] = deepcopy(data['id'])
                    temp_data['id'] = str(max_id)

                    # we can change spacing based on # of characters in name
                    if len(result) == 2:
                        temp_data['x'] = deepcopy(temp_data['x']) + index*70 - (35 * (len(result)-1))
                    elif len(result) == 3:
                        temp_data['x'] = deepcopy(temp_data['x']) + index*70 - (35 * (len(result)-1))
                    #figure out the odd # case
                    else:
                        middle = math.ceil(len(result)/2)
                        if index < middle:
                            temp_data['x'] = deepcopy(temp_data['x']) + index*70 - (35 * (middle-1))
                            temp_data['y'] = deepcopy(temp_data['y']) - 13
                        else:
                            temp_data['x'] = deepcopy(temp_data['x']) + (index - middle)*70 - (35 * (len(result)-middle-1))
                            temp_data['y'] = deepcopy(temp_data['y']) + 13
                    G.add_node(str(max_id), temp_data)
                    
                    if index == 0:
                        gene_symbol = str(entrez_to_symbol(gene))
                    else:
                        gene_symbol = gene_symbol + ', ' + str(entrez_to_symbol(gene))
                G.node[node]['symbol'] = gene_symbol
                #G.node[node]['gene_id_list'] = gene_list
                parent_to_list[node] = gene_list
        if data['type'] == 'compound':
            if isinstance(s.parse(s.get(data['name'])), unicode):
                compound_array = s.parse(s.get(data['name'])).split()
                name_index = compound_array.index('NAME')
                gene_symbol = [compound_array[name_index+1]]
            else:
                gene_symbol = s.parse(s.get(data['name']))['NAME']
            G.node[node]['gene_names'] = ' '.join(gene_symbol)
            G.node[node]['symbol'] = gene_symbol[0].replace(';', '')
    return (parent_list, parent_to_list)

def undefNodes(G, root):
    #handling undefined nodes
    comp_dict = dict()       # key: complex; element: list of components
    node_to_comp = dict()    # key: component; element: complex
    comp_array_total = []    # array containing all component nodes
    complex_array = []       # array containing all undefined complexes
    for entry in root.findall('entry'):
        #array to store components of undefined nodes
        component_array = []
        if entry.attrib['name'] == 'undefined':
            node_id = entry.attrib['id']
            complex_array.append(node_id)

            #find components
            for index, component in enumerate(entry.iter('component')):
                component_array.append(component.get('id'))   
                #check to see which elements are components
                comp_array_total.append(component.get('id'))
                node_to_comp[component.get('id')] = node_id

            #store into node dictionary
            G.node[node_id]['component'] = component_array
            comp_dict[node_id] = component_array

            # store gene names
            gene_name_array = []
            for index, component_id in enumerate(component_array):
                if index == 0:
                    gene_name_array.append(G.node[component_id]['gene_names'])
                else:
                    gene_name_array.append('\n' + G.node[component_id]['gene_names'])

            G.node[node_id]['gene_names'] = gene_name_array

            # store gene symbols
            gene_symbol_array = []
            for index, component_id in enumerate(component_array):
                if index == 0:
                    gene_symbol_array.append(G.node[component_id]['symbol'])
                else:
                    gene_symbol_array.append('\n' + G.node[component_id]['symbol'])

            G.node[node_id]['symbol'] = gene_symbol_array
    return (complex_array, comp_array_total, node_to_comp, comp_dict)

#add edges without redirection
def addEdges(G, res, comp_array_total, node_to_comp):
    edge_list = []
    edge_pairs = []
    # add edges to networkx graph
    # redirect edges to point to undefined nodes containing components in order to connect graph
    for edge in res['relations']:
        source = edge['entry1']
        dest = edge['entry2']
        if (edge['entry1'] in comp_array_total) == True: 
            source = node_to_comp[edge['entry1']]
        if (edge['entry2'] in comp_array_total) == True:
            dest = node_to_comp[edge['entry2']] 
        edge_list.append((source, dest, edge))
        edge_pairs.append((source,dest))
        #check for duplicates
        if (source, dest) in G.edges():
            name = []
            value = []
            link = []
            name.append(G.edge[source][dest]['name'])
            value.append(G.edge[source][dest]['value'])
            link.append(G.edge[source][dest]['link'])
            name.append(edge['name'])
            value.append(edge['value'])
            link.append(edge['link'])
            G.edge[source][dest]['name'] = '\n'.join(name)
            G.edge[source][dest]['value'] = '\n'.join(value)
            G.edge[source][dest]['link'] = '\n'.join(link)
        #check for self loops
        elif source == dest:
            continue
        else:
            G.add_edge(source, dest, edge)

#Add and redirect edges
def addAndRedirectEdges(G, res, complex_array, comp_array_total, parent_list, parent_to_list, node_to_comp, comp_dict):
    # add edges to networkx graph
    for edge in res['relations']:
        edge_pairs = []
        source_redirected = False
        #redirect edges
        #find edges that point to complex
        source_component_list = [edge['entry1']]
        dest_component_list = [edge['entry2']]
        if (edge['entry1'] in complex_array) == True:
            source_component_list = comp_dict[edge['entry1']]
        if (edge['entry2'] in complex_array) == True:
            dest_component_list = comp_dict[edge['entry2']]

        #find edges that point to component
        if (edge['entry1'] in comp_array_total) == True: 
            source_undef_complex = node_to_comp[edge['entry1']]
            #look for other components in that complex
            source_component_list = comp_dict[source_undef_complex]
            source_redirected = True
        if (edge['entry2'] in comp_array_total) == True:
            dest_undef_complex = node_to_comp[edge['entry2']]
            #look for other components in that complex
            dest_component_list = comp_dict[dest_undef_complex]

            if source_redirected == True:   
                #look for self loops
                if (source_undef_complex == dest_undef_complex):
                    comp_nodes = comp_dict[source_undef_complex]
                    for node in comp_nodes:
                        # stop redirecting edges to this undefined node
                        comp_array_total.append(node)
                    # undo redirection of edge
                    source_component_list = [edge['entry1']]
                    dest_component_list = [edge['entry2']]

        #splitting nodes
        temp_source_array = deepcopy(source_component_list)
        for source in source_component_list:
            if source in parent_list:
                gene_list = parent_to_list[source]
                temp_source_array.remove(source)
                for gene in gene_list:
                    temp_source_array.append(gene)
        source_component_list = deepcopy(temp_source_array)

        temp_dest_array = deepcopy(dest_component_list)
        for dest in dest_component_list:
            if dest in parent_list:
                gene_list = parent_to_list[dest]
                temp_dest_array.remove(dest)
                for gene in gene_list:
                    temp_dest_array.append(gene)
        dest_component_list = deepcopy(temp_dest_array)

        # store edge pairs in array
        for source in source_component_list:
            for dest in dest_component_list:
                edge_pairs.append((source,dest))

        #check for duplicates
        for source, dest in edge_pairs:
            if (source, dest) in G.edges():

                # undo edge redirection
                if G.edge[source][dest]['name'] == edge['name']:
                    temp_source = G.edge[source][dest]['entry1']
                    temp_dest = G.edge[source][dest]['entry2']
                    temp_edge = G.edge[source][dest]
                    G.add_edge(temp_source, temp_dest, temp_edge)
                    G.remove_edge(source, dest)
                    G.add_edge(edge['entry1'], edge['entry2'], edge)
                    break
                else:
                    name = []
                    value = []
                    link = []
                    name.append(G.edge[source][dest]['name'])
                    value.append(G.edge[source][dest]['value'])
                    link.append(G.edge[source][dest]['link'])
                    name.append(edge['name'])
                    value.append(edge['value'])
                    link.append(edge['link'])
                    G.edge[source][dest]['name'] = '\n'.join(name)
                    G.edge[source][dest]['value'] = '\n'.join(value)
                    G.edge[source][dest]['link'] = '\n'.join(link)
            else:
                G.add_edge(source, dest, edge)

def addReaction(G, ETroot):
    for entry in ETroot.findall('reaction'):
        entryDict = dict()
        entryDict['reaction'] =  entry.attrib['name']
        entryDict['name'] =  entry.attrib['type']
        #source
        substrate = entry.find('substrate')
        source = substrate.attrib["id"]
        entryDict['entry1'] = substrate.attrib["name"]
        #dest
        product = entry.find('product')
        dest = product.attrib["id"]
        entryDict['entry2'] = product.attrib["name"]
        G.add_edge(source, dest, entryDict)
        
def log2FoldChange(G, path_to_csv):
    # log2FoldChange 
    DE_genes_df = pandas.read_csv(path_to_csv)

    short_df = DE_genes_df[['_id', 'Ensembl', 'log2FoldChange']]
    
    short_df.to_dict('split')
    
    gene_to_log2fold = dict()
    for entry in short_df.to_dict('split')['data']:
        if type(entry[0]) is float:
            if math.isnan(entry[0]):
                gene_to_log2fold[entry[1]] = entry[2]
            else:
                gene_to_log2fold[entry[0]] = entry[2]
        else:
            gene_to_log2fold[entry[0]] = entry[2]

    id_to_log2fold = dict()

    for id_, node in G.nodes(data=True):

        log2fold_array = []
        if node['name'] == 'undefined':
            continue
        elif node['type'] == 'map':
            continue
        else:
            result = node['name'].split("hsa:")
            result = ''.join(result)
            result = result.split()
            for item in result:
                if item in gene_to_log2fold.keys():
                    log2fold_array.append(gene_to_log2fold[item])
            if len(log2fold_array) > 0:
                id_to_log2fold[node['id']] = log2fold_array
    return id_to_log2fold


def entrez_to_symbol(query_gene):
    result = mg.getgene(query_gene, fields='symbol')
    return result['symbol']
