from ete3 import Tree, NodeStyle, TreeStyle, TextFace, CircleFace, ImgFace, RectFace, faces
import sys, copy
from collections import defaultdict
from operator import itemgetter

class MyTextFace(TextFace):
        def __init__(self, *args, **kwargs):
                TextFace.__init__(self, *args, **kwargs)
                self.margin_right = 10

#728 clusters (MCC)
clusters = {'GN1': ['EBOV|EMLab-RT|KG167||GIN|Boke|Boke-Dabis_Kasongony|MinION_LQ10|2015-06-29', 'EBOV|EMLab-RT|EM_COY_2015_017021||GIN|Fria|Banguine-|MinION_LQ05|2015-05-25'], 'SL3': ['EBOV|EMLab-RT|EM_FORE_2015_2878||GIN|Forecariah|Kaleah-Kindoyah|MinION_LQ05|2015-10-24', 'EBOV|EMLab-RT|EM_GUI_2015_004674||GIN|Siguiri|Lero|MinION_LQ05|2015-03-26'], 'Boke': ['EBOV|EMLab-RT|EM_COY_2015_013731||GIN|Coyah|Lansanaya-|MinION_LQ10|2015-03-14', 'EBOV|EMLab-RT|EM_COY_2015_015802||GIN|Coyah|Maneah-Kalokhoyah|MinION_LQ05|2015-04-14']}

#729 clusters (ML)
#clusters = {'GN1': ['EBOV|EMLab-RT|EM_COY_2015_017021||GIN|Fria|?|MinION_LQ05|2015-05-25', 'EBOV|EMLab-RT|KG35||GIN|Boke|?|MinION_LQ10|2015-06-08'], 'SL3': ['EBOV|EMLab-RT|EM-FORE_2015_1160||GIN|Forecariah|?|MinION_LQ10|2015-07-13', 'EBOV|EMLab-RT|EM_COY_2015_017057||GIN|Dubreka|?|MinION_LQ05|2015-05-27'], 'Boke': ['EBOV|EMLab-RT|EM_COY_2015_014098||GIN|Conakry|?|MinION_LQ05|2015-03-26', 'EBOV|EMLab-RT|KG35||GIN|Boke|?|MinION_LQ10|2015-06-08']}

#923 clusters
#clusters = {'GN1': ['EBOV|EMLab-RT|GUI_CTS_2015_0051||GIN|Boke|?|MinION_LQ10|2015-06-21', 'EBOV|EMLab-RT|EM_COY_2015_013857||GIN|Forecariah|?|MinION_LQ05|2015-03-18'], 'SL3': ['EBOV|EMLab-RT|EM_COY_2015_017135||GIN|Dubreka||MinION_LQ10|2015-05-29', 'EBOV|EMLab-RT|EM-FORE_2015_1160||GIN|Forecariah|?|MinION_LQ10|2015-07-13'], 'Boke': ['EBOV|EMLab-RT|EM_COY_2015_014102||GIN|Conakry|?|MinION_LQ05|2015-03-26', 'EBOV|EMLab-RT|GUI_CTS_2015_0051||GIN|Boke|?|MinION_LQ10|2015-06-21']}

#set size of nodes
size = {'big': 2.0, 'small': 10.0}

#prefec colours
colours = ['#E31A1C', '#33A02C', '#1F78B4', '#B15928', '#6A3D9A', '#FF7F00', '#FB9A99', '#B2DF8A', '#A6CEE3', '#FFFF99', '#CAB2D6', '#FDBF6F', '#D3D3D3']

#read metadata
def get_meta(metadata, big_tree):
    leaf_names = [n.get_leaf_names()[0].strip("'") for n in big_tree]
    for each in leaf_names:
        cols = each.split('|')
        if cols[2] == 'MinION':
            loc_strings = cols[4].split('-')
            metadata[each] = {'country': cols[3], 'date': cols[5], 'short_id': cols[1], 'prefec': loc_strings[0], 'subpre': loc_strings[1], 'instrument': 'MinION'}
        else:
            if cols[4].startswith('G'):
                cols[4] = 'GUI'
            metadata[each] = {'country': cols[4], 'date': cols[8], 'short_id': cols[2], 'prefec': cols[5], 'subpre': ' ', 'instrument': 'other'}
    return metadata

def get_meta_new(metadata, big_tree):
    leaf_names = [n.get_leaf_names()[0].strip("'") for n in big_tree]
    for each in leaf_names:
        cols = each.split('|')
        if cols[1] == 'EMLab-RT':
            metadata[each] = {'country': cols[4], 'date': cols[8], 'short_id': cols[2], 'prefec': cols[5], 'subpre': '', 'instrument': 'MinION', 'group': cols[1]}
        else:
            metadata[each] = {'country': cols[4], 'date': cols[8], 'short_id': cols[2], 'prefec': cols[5], 'subpre': '', 'instrument': 'other', 'group': cols[1]}
    return metadata        

def get_colours(clusters, tree, colours):    
    #get a list of prefectures for both clusters
    both_leaves = []
    for c in [key for key in list(clusters.keys()) if key in ['SL3', 'GN1']]:
        b = ["'" + clusters[c][0] + "'", "'" + clusters[c][1] + "'"]
        for a in tree.get_common_ancestor(b).get_leaves():
            both_leaves.append(a.name[1:-1])
    
    #count the number of samples of each prefecture and assign colours
    counts = defaultdict(int)
    colourDict= {}
    for each in both_leaves:
        #print each, metadata[each]['instrument'], metadata[each]['prefec']
        if metadata[each]['instrument'] == 'MinION':
            counts[metadata[each]['prefec']] += 1    
    for n, (key, value) in enumerate(sorted(list(counts.items()), key=itemgetter(1), reverse=True)):
        colourDict[key] = colours[n]
    for each in list(counts.keys()):
        print('%s\t%s\t%s' %(each, counts[each], colourDict[each]))
    return colourDict

#render tree function
def render_tree(tree, mode, cluster, colourDict, width, position):
    duplicates = ["'EBOV|EMLab-RT|KG12||GIN|Boke|?|MinION_LQ05|2015-05-27'", "'EBOV|EMLab-RT|KG45||GIN|Boke|?|MinION_LQ10|2015-06-09'", "'EBOV|EMLab-RT|KG90||GIN|Boke|?|MinION_LQ05|2015-06-19'", "'EBOV|EMLab-RT|KG91||GIN|Boke|?|MinION_LQ05|2015-06-20'"]
    print('Running %s: %s cluster' %(mode, cluster))
    if mode == 'small':
        #delete unwanted leaves
        keep_leaves = []
        b = ["'" + clusters[cluster][0] + "'", "'" + clusters[cluster][1] + "'"]
        for a in tree.get_common_ancestor(b).get_leaves():
            keep_leaves.append(a.name)
        delete_leaves = [leaf for leaf in tree.get_leaf_names() if leaf not in keep_leaves]
        #if cluster == 'Boke':
        #    delete_leaves.extend(duplicates)
        print('Keeping %s leaves' %len(keep_leaves))
        for leaf in delete_leaves:
            if tree.search_nodes(name=leaf)[0]:
                n = tree.search_nodes(name=leaf)[0]
                #if leaf in duplicates:
                #    n.delete(preserve_branch_length=True)
                #else:
                n.delete()
                #print 'Removed %s' %(n.get_leaf_names()[0])
        tree.ladderize()
        tree.write(outfile=sys.argv[1] + '_' + cluster + '.nwk')

    for n in tree.get_leaves():
        display_name = n.get_leaf_names()[0][1:-1]
        country = metadata[display_name]['country']
        instrument = metadata[display_name]['instrument']
        #if cluster == 'Boke':
        if mode == 'big':
            n.add_face(MyTextFace(metadata[display_name]['short_id'], ftype="Helvetica", fsize=size[mode]), column=0, position="aligned")
            n.add_face(MyTextFace(metadata[display_name]['prefec'], ftype="Helvetica", fsize=size[mode]), column=1, position="aligned")
            n.add_face(MyTextFace(metadata[display_name]['date'], ftype="Helvetica", fsize=size[mode]), column=2, position="aligned")
        if country == 'GUI' or country == 'GIN':
            if instrument == 'MinION':
                C = CircleFace(radius=size[mode]/2, color=colourDict[metadata[display_name]['prefec']])
            else:
                C = CircleFace(radius=size[mode]/2, color='#F1F1F1')
        elif country == 'SLE':
            if instrument == 'MinION':
                C = RectFace(width=size[mode], height=size[mode], bgcolor=colourDict[metadata[display_name]['prefec']], fgcolor='#FFFFFF')
            else:
                C = RectFace(width=size[mode], height=size[mode], bgcolor='#D3D3D3', fgcolor='#FFFFFF')
        elif country == 'LIB' or country == 'LBR':
            if instrument == 'MinION':
                C = RectFace(triangle=True, width=size[mode], height=size[mode], bgcolor=colourDict[metadata[display_name]['prefec']], fgcolor='#FFFFFF')
            else:
                C = RectFace(triangle=True, width=size[mode], height=size[mode], bgcolor='#939393', fgcolor='#FFFFFF')
        n.add_face(C, column=1, position=position)
    tree.render(file_name=sys.argv[1] + '_' + cluster + '.pdf', tree_style=ts, w=width)

big_tree = Tree(sys.argv[1])
mode = sys.argv[2]
metadata = {}
metadata = get_meta_new(metadata, big_tree)
colourDict = get_colours(clusters, big_tree, colours)

#remove dodgy sample
big_tree.search_nodes(name="'EBOV|EMLab-RT|IPDPFHGINSP_GUI_2015_5339||GIN|Conakry|?|MinION_LQ05|2015-04-08'")[0].delete(preserve_branch_length=True)
#root the same as the MCC tree
ancestor = big_tree.get_common_ancestor("'EBOV|EMLab|EM_079422|KR817187|GIN|Macenta|?||2014-03-27'","'EBOV|EMLab|Gueckedou-C05|KJ660348|GIN|Gueckedou|?||2014-03-19'")
big_tree.set_outgroup(ancestor)
big_tree.ladderize()

ts = TreeStyle()
ts.show_leaf_name = False
#ts.show_branch_support = True
ts.scale = 100000
if mode == 'small':
    ts.scale = 750000

#add legend
for each in list(colourDict.keys()):
    ts.legend.add_face(CircleFace(radius=size[mode]/2, color=colourDict[each]), column=0)
    ts.legend.add_face(TextFace(each, ftype="Helvetica", fsize=size[mode]), column=1)
ts.legend.add_face(CircleFace(radius=size[mode]/2, color='#F1F1F1'), column=0)
ts.legend.add_face(TextFace('Guinea', ftype="Helvetica", fsize=size[mode]), column=1)
ts.legend.add_face(RectFace(width=size[mode], height=size[mode], bgcolor='#D3D3D3', fgcolor='#FFFFFF'), column=0)
ts.legend.add_face(TextFace('Sierra Leone', ftype="Helvetica", fsize=size[mode]), column=1)
ts.legend.add_face(RectFace(triangle=True, width=size[mode], height=size[mode], bgcolor='#939393', fgcolor='#FFFFFF'), column=0)
ts.legend.add_face(TextFace('Liberia', ftype="Helvetica", fsize=size[mode]), column=1)

#reset nodes
ns = NodeStyle()
ns['size'] = 0
ns['hz_line_width'] = 2
ns['vt_line_width'] = 2
if mode =='big':
    ns['hz_line_width'] = 0
    ns['vt_line_width'] = 0
for n in big_tree.traverse():
        n.set_style(ns)

#render tree
if mode == 'big':
    cluster = 'big'
    render_tree(big_tree, mode, cluster, colourDict, width=2000, position='float')
elif mode == 'small':
    for cluster in list(clusters.keys()):
        tree = copy.deepcopy(big_tree)
        render_tree(tree, mode, cluster, colourDict, width=4000, position='branch-right')
else:
    print('Mode not recognised: %s' %mode)
