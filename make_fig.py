import numpy as np 
from unet import utils
from unet.sim_measures import jaccard_pixelwise, jaccard_roiwise
from helper_fxns import ijroi
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from shapely.geometry import Polygon, MultiPolygon
from shapely.ops import cascaded_union
from keras.preprocessing.image import img_to_array, load_img
from scipy.misc import imsave

import pandas as pd
import itertools
import json
import sys
import os

from bokeh.plotting import figure, show, output_file, ColumnDataSource
from bokeh.palettes import Category20_20 as palette
from bokeh.palettes import brewer








def read_json(data_file):
    #INPUT:
    #   data_file: path to json file
    #OUTPUT:
    #   return: data from json file

    with open(data_file, "r") as read_file:
        data = json.load(read_file)
        return(data)


def convert_to_single_polygon(multi_polygon_object, gene_name):
    # INPUT:
    #     multi_polygon_object: list of polygon objects
    #     gene_name: names of the genes being smushed together
    # OUTPUT:
    #     return: single polygon object

    single_poly = MultiPolygon(Polygon(p.exterior) for p in multi_polygon_object)
    single_poly = cascaded_union([Polygon(component.exterior).buffer(4).buffer(-3) for component in single_poly])
    if single_poly.type == "Polygon":
        return(single_poly)
    else:
        print("!!!!!!YOU NEED TO FIX THE IMAGEJ ROI FOR: ", gene_name, " AT: ", single_poly.bounds)
        return(None)




def shape_objs_from_roi_file(roi_file):
    roi = ijroi.read_roi_zip(roi_file)
    ply_list = []
    for r in roi:
        coords = r[1]
        coord_list = []
        for p in coords:
            pt = (p[0], p[1])
            coord_list.append(pt)
        ply = Polygon(coord_list)
        ply = ply.buffer(0)
        if ply.type == "MultiPolygon":
            ply = convert_to_single_polygon(ply, roi_file)
        
        if ply.area > 0:
            ply_list.append(ply)

    return(ply_list)


def shape_obj_from_df_coords(coords):
    # INPUT: coordinate list for a polygon
    # OUTPU: a shapely.Polygon object

    ys = coords[0]
    xs = coords[1]
    coord_list = []
    for y, x in zip(ys, xs):
        coord_list.append((y,x))
    ply = Polygon(coord_list)
    ply.buffer(0)
    
    return(ply)


         

def get_mean_pix_in_bb(img, bb_bounds):
    # print(bb_bounds)
    y1 = int(bb_bounds[0])
    y2 = int(bb_bounds[2])
    x1 = int(bb_bounds[1])
    x2 = int(bb_bounds[3])
    bb_img = img[y1:y2,x1:x2,:]
    bb_mean = np.average(bb_img)
    # print(bb_mean)
    # print(bb_img.shape)
    return(bb_mean)


def redraw_masks(roi_df, ch_dict, img_shape, mask_save_dir):
    from PIL import Image, ImageDraw
    # arr = np.zeros(img_shape)
    img_shape = (img_shape[1], img_shape[0])
    for ch_num, ch_name in ch_dict.items():
        img = Image.new('L', img_shape, 0)
        ch_df = roi_df.loc[roi_df['gene'] == ch_name]
        img = Image.new('L', img_shape, 0)
        for index, row in ch_df.iterrows():
            ply = shape_obj_from_df_coords(row['polygon_coords'])
            ys = np.array(ply.exterior.xy)[0].tolist()
            xs = np.array(ply.exterior.xy)[1].tolist()
            ply_coords = list(zip(xs, ys))
            # print(ply_coords)
            ImageDraw.Draw(img).polygon(ply_coords, outline=255, fill=255)
        mask = np.array(img)
        outfile = os.path.join(mask_save_dir, ch_name + '_test_mask.tif')
        imsave(outfile, mask)
        # print(mask.shape)
        # print(np.max(mask))



def df_from_roi_file(roi_zip_files, root_file_path, ch_dict, area_threshold):
    #INPUT:
    #   roi_zip_files: names of the roi files in the root folder
    #   root_file_path: path to the enclosing folder which contains roi/img/csv files
    #   ch_dict: key=index of channel in multichannel tiff img.  value=channel name.  I omitted dapi.
    #   area_threshold: minimum roi size to pass
    #OUTPUT:
    #   return: pd.dataframe of all rois of selected channels that meet minimum size requirement



    # roi_zip_files.sort()
    # gene_names = [s.split('_')[0] for s in roi_zip_files]
    # gene_names.sort()
    column_names = ['gene', 'channel_number', 'bb_coords', 'polygon_coords', 'polygon_area', 'mean_intensity_in_bb', 'centroid', 'index_num_coexpressed', 'channel_nums_coexpressed', 'channel_names_coexpressed']
    img_df = pd.DataFrame(columns=column_names)
    # for name, roi_file in zip(gene_names, roi_zip_files):
        # full_roi_file_path = os.path.join(root_file_path, roi_file)
        # ply_list = shape_objs_from_roi_file(full_roi_file_path)
    for ch, name in ch_dict.items():
        file_name = name + '_RoiSet.zip'
        full_roi_file_path = os.path.join(root_file_path, file_name)
        ply_list = shape_objs_from_roi_file(full_roi_file_path)
        for p in ply_list:
            if p.area > area_threshold:
                load_dict = {}
                load_dict['gene'] = name
                load_dict['bb_coords'] = list(p.bounds)
                load_dict['channel_number'] = ch
                load_dict['polygon_coords'] = list(p.exterior.xy)
                load_dict['polygon_area'] = p.area
                load_dict['mean_intensity_in_bb'] = get_mean_pix_in_bb(dummy_img, list(p.bounds))
                load_dict['centroid'] = list(p.centroid.coords)
                img_df = img_df.append(load_dict, ignore_index=True)
                # print(np.array(p.exterior.xy))
    # print(img_df)
    return(img_df)        



def find_coexpression_indexs(roi_df, channel_dict, iou_threshold):
    #INPUT:
    #   roi_df: pd.DataFrame of all rois
    #   channel_dict: key=index of channel in multichannel tiff img.  value=channel name.  I omitted dapi.
    #   iou_threshold: minimum IntersectionOverUnion value to count as coexpression
    #OUTPUT:
    #   return: pd.DataFrame same as 'roi_df' with extra columns for 'index_num_coexpressed' 'channel_nums_coexpressed' 'channel_names_coexpressed'

    # add new columns to df
    roi_df = roi_df
    roi_df['index_num_coexpressed'] = roi_df['index_num_coexpressed'].astype(object)
    roi_df['channel_nums_coexpressed'] = roi_df['channel_nums_coexpressed'].astype(object)
    roi_df['channel_names_coexpressed'] = roi_df['channel_names_coexpressed'].astype(object)
    
    # generate list of all possible combinations of coexpression n=2
    gene_list = list(channel_dict.values())
    comb_gene_list = list(itertools.combinations(gene_list, 2))
    

    for comb in comb_gene_list:
        g1_df = roi_df.loc[roi_df['gene'] == comb[0]]
        g2_df = roi_df.loc[roi_df['gene'] == comb[1]]
        for index1, row1 in g1_df.iterrows():
            poly1 = shape_obj_from_df_coords(row1['polygon_coords'])
            for index2, row2 in g2_df.iterrows():
                poly2 = shape_obj_from_df_coords(row2['polygon_coords'])

                # get IoU value IF the plygons intersect.  saves allot of time.
                if poly1.intersects(poly2):
                    intersection = poly1.intersection(poly2).area
                    union = (poly1.area+poly2.area) - intersection
                    iou = intersection/union

                    # if IoU > iou_threshold update index_num_coexpressed list in df
                    if iou > iou_threshold:
                        coex_prev_index_1 = roi_df.at[index1, 'index_num_coexpressed']
                        if np.isnan(coex_prev_index_1).all():
                            index_list_1 = [index2]
                            roi_df.at[index1, 'index_num_coexpressed'] = index_list_1

                        else:
                            coex_prev_index_1.append(index2)
                            roi_df.at[index1, 'index_num_coexpressed'] = coex_prev_index_1


                        coex_prev_index_2 = roi_df.at[index2, 'index_num_coexpressed']
                        if np.isnan(coex_prev_index_2).all():
                            index_list_2 = [index1]
                            roi_df.at[index2, 'index_num_coexpressed'] = index_list_2
                        else:
                            coex_prev_index_2.append(index1)
                            roi_df.at[index2, 'index_num_coexpressed'] = coex_prev_index_2

    return(roi_df)



def fill_in_df_holes(roi_df, ch_dict):
    #INPUT:
    #   roi_df: coexpression roi pd.DataFrame
    #   ch_dict: key=index of channel in multichannel tiff img.  value=channel name.  I omitted dapi.
    #OUTPUT:
    #   roi_df: df same as roi_df, with 'channel_number' and 'gene' columns filled in.
    #   count_df: df of counts of coexpression

    roi_df = roi_df

    # create an empty array to count coexpression into
    n_ch = int(max(ch_dict.keys()))+1
    count_arr = np.zeros((n_ch, n_ch))
    ch_keys = list(ch_dict.keys())

    for index, row in roi_df.iterrows():
        coex_indxs = row['index_num_coexpressed']

        # update count_arr
        count_coord_primary = int(row['channel_number'])
        count_arr[(count_coord_primary, count_coord_primary)] += 1

        if isinstance(coex_indxs, list):
            ch_num_list = []
            ch_name_list = []
            for inx in coex_indxs:
                ch_num = roi_df.at[inx, 'channel_number']
                ch_name = roi_df.at[inx, 'gene']
                ch_num_list.append(int(ch_num))
                ch_name_list.append(ch_name)

                # update count_arr
                count_arr[(count_coord_primary, int(ch_num))] += 1

            roi_df.at[index, 'channel_nums_coexpressed'] = ch_num_list
            roi_df.at[index, 'channel_names_coexpressed'] = ch_name_list
    

    # convert count_arr to a named dataframe of columns and indexes of ch_name
    count_df = pd.DataFrame(columns=ch_keys, index=ch_keys)
    for k1 in ch_keys:
        for k2 in ch_keys:
            count_df.at[k1,k2] = count_arr[(int(k1), int(k2))]
    count_df = count_df.rename(columns=ch_dict)
    count_df = count_df.rename(index=ch_dict)
    
    return(roi_df, count_df)





def format_df_for_plot(df, ch_dict):
    #INPUT:
    #   df: pd.DataFrame of coexpression
    #   ch_dict: key=index of channel in multichannel tiff img.  value=channel name.  I omitted dapi.
    #OUTPUT:
    #   return: 

    orig_df = df
    
    # list of 
    completed_indexes = []

    column_names = ['gene_list', 'ch_num_list', 'index_nums', 'merged_polygon_coords']
    reform_df = pd.DataFrame(columns=column_names)


    for inx, row in orig_df.iterrows():
        loader_dict = {}
        if inx not in completed_indexes:
            # cells w/o coexpression get loaded with data about itself
            if np.isnan(row['index_num_coexpressed']).all():
                loader_dict['gene_list'] = [row['gene']]
                loader_dict['ch_num_list'] = [row['channel_number']]
                loader_dict['index_nums'] = [inx]
                loader_dict['merged_polygon_coords'] = row['polygon_coords']
                reform_df = reform_df.append(loader_dict, ignore_index=True)
                completed_indexes.append(inx)

            # cells w/ coexpression get their df.indexs updated all at the same time
            else:
                poly_obj_list = []
                co_indexes = row['index_num_coexpressed']
                co_indexes.append(inx)
                gene_list = []
                ch_num_list = []
                
                for co_inx in co_indexes:
                    gene_list.append(orig_df.iloc[co_inx]['gene'])
                    ch_num_list.append(orig_df.iloc[co_inx]['channel_number'])
                    poly_coords = orig_df.iloc[co_inx]['polygon_coords']
                    poly_obj_list.append(shape_obj_from_df_coords(poly_coords))

                # converty multiple overlapping polygons into a single polygon.  this can cause some issues!!
                merged_poly = convert_to_single_polygon(poly_obj_list, gene_name=co_indexes)
                if merged_poly != None:
                    loader_dict['merged_polygon_coords'] = list(merged_poly.exterior.xy)
                    loader_dict['gene_list'] = gene_list
                    loader_dict['ch_num_list'] = ch_num_list
                    loader_dict['index_nums'] = co_indexes
                    reform_df = reform_df.append(loader_dict, ignore_index=True)
                    completed_indexes = completed_indexes + co_indexes
                else:
                    print("ERRRORRRRRRRR")
                    print(poly_obj_list)
                    print(co_indexes)
                    print(gene_list)

    return(reform_df)






def plot_roi_df_sankey(roi_df, ch_dict, primary_node, secondary_node):
    # REF: http://holoviews.org/reference/elements/bokeh/Sankey.html
    # INPUT:
    #   roi_df: the same df used for plot_roi_df()
    #   primary_name: name of the gene you want of the left of the sankey plot
    # OUTPUT:

    import holoviews as hv
    from holoviews import opts, dim
    hv.extension('bokeh')
    # sankey = hv.Sankey([
    # ['A', 'X', 5],
    # ['A', 'Y', 7],
    # ['A', 'Z', 6],
    # ['B', 'X', 2],
    # ['B', 'Y', 9],
    # ['B', 'Z', 4],
    # ['B', 'P', 0]]
    # )
    # sankey.opts(width=600, height=400)


    # pairs = [(k, v) for k, v in ch_dict.items()]
    # # print(pairs)
    # nodes = hv.Dataset(pairs, 'index', 'label')
    
    edges_d = {}
    edges = []

    for inx, row in roi_df.iterrows():
        if primary_node in row['gene_list']:
            if secondary_node in row['gene_list']:
                if (primary_node, secondary_node) not in edges_d:
                    edges_d[(primary_node, secondary_node)] = 1
                else:
                    edges_d[(primary_node, secondary_node)] = edges_d[(primary_node, secondary_node)] + 1
            else:    
                if (primary_node) not in edges_d:
                    edges_d[(primary_node)] = 1
                else:
                    edges_d[(primary_node)] = edges_d[(primary_node)] + 1
    print(edges_d)
  


    
    # hv.save(sankey, 'sankey.html')




def plot_roi_df(roi_df):
    color_dict = {}
    comb_dict = {}
    ys = []
    xs = []
    genes = []
    color_list = []
    for inx, row in roi_df.iterrows():
        # coords = row['centroid'][0]
        # ys.append(int(coords[0]))
        # xs.append(int(coords[1]))
        # print(type(row['merged_polygon_coords'][0].tolist()))
        ys.append(np.negative(row['merged_polygon_coords'][0]).tolist())
        xs.append(row['merged_polygon_coords'][1].tolist())
        genes.append(row['gene_list'])

    unique_comb = [list(x) for x in set(tuple(x) for x in genes)]
    colors = itertools.cycle(palette)

    for comb_n, comb in zip(range(len(unique_comb)), unique_comb):
        comb_dict[comb_n] = comb

    for comb_n, color in zip(range(len(unique_comb)), colors):
        color_dict[comb_n] = color
    
    for gene in genes:
        c = None
        for comb_n, comb in comb_dict.items():
            if comb == gene:
                c = color_dict[comb_n]
        color_list.append(c)
    # print(color_list)

    # colors = brewer['Spectral'][len(roi_df.gene_list.unique())]
    # colormap = {i: colors[i] for i in roi_df.gene_list.unique()}
    # colors = [colormap[x] for x in roi_df.gene_list]

    source = ColumnDataSource(data=dict(
        x=xs,
        y=ys,
        genes=genes,
        colors=color_list,
    ))
    TOOLTIPS = [
        ("Genes", "@genes"),
    ]

    TOOLS="pan,wheel_zoom,reset,save,box_select,hover"

    p = figure(tools=TOOLS, tooltips=TOOLTIPS)

    p.patches('x','y',color='colors',fill_alpha=0.6,source=source)
    # p.legend.click_policy='mute'
    show(p)




iou_threshold = 0.4
size_threshold = 100

DUMMY_IMG_PATH = './processed_images/test2/antiGFP_max-z.png'
dummy_img = load_img(DUMMY_IMG_PATH, color_mode='grayscale')
dummy_img = img_to_array(dummy_img)
print(dummy_img.shape)

DUMMY_JSON_PATH = './processed_images/test2/image_data.json'  
j = read_json(DUMMY_JSON_PATH)
ch_dictionary = j['channel_dictionary']

### removes dapi from calculations ###
del ch_dictionary['0']



file_path = None
if len(sys.argv) > 1:
    file_path = str(sys.argv[1])

files = os.listdir(file_path)
roi_files = [j for j in files if 'RoiSet.zip' in j]


df = df_from_roi_file(roi_files, file_path, ch_dictionary, size_threshold)
redraw_masks(df, ch_dictionary, dummy_img.shape[:2], './processed_images/test2/masks/')
df = find_coexpression_indexs(df, ch_dictionary, iou_threshold)
df, count_df = fill_in_df_holes(df, ch_dictionary)
plot_df = format_df_for_plot(df, ch_dictionary)
print(count_df)

### NEEDS ALLOT OF WORK ###
# plot_roi_df_sankey(plot_df, ch_dictionary, 'virSYFP', 'gad1')

### KEEP ####
# plot_roi_df(plot_df)

df.to_csv('./processed_images/test2/index.csv')
count_df.to_csv('./processed_images/test2/coexpression_counts.csv')
