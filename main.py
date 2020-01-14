import argparse
import os

import numpy as np
import pandas as pd
from skimage import io
from skimage.transform import resize
from skimage.filters import gaussian
import tensorflow as tf

from df_resources import unet, preproc, utils, pixelshift37, rundmc


def process(config):
    TILE_SHAPE=(540,540)
    PADDING=(184,184)
    SEED=0
    EL_SIZE=[600, 600] #micrometers
    BATCH_NORM=True
    LAMBDA=50 #50
    V_BAL=0.1 #0.1
    SIGMA_BAL=10 #10 
    SIGMA_SEP=6 #6



    image_path_list = []
    if config.multi_img_dir is not None:
        image_path_list = [os.path.normpath(os.path.join(config.multi_img_dir, f)) for f in os.listdir(config.multi_img_dir)]
    else:
        image_path_list.append(os.path.normpath(config.image_path))
    

    for f in image_path_list:
        if config.output_path == None:
            just_file_name = os.path.basename(f).split('.')[0]
            output_path = os.path.normpath(os.path.join('./output', just_file_name))
        else:
            output_path = config.output_path
        
        if not os.path.exists(output_path):
            os.makedirs(output_path)

        img = io.imread(f)
        if config.allign_channels:
            img = rundmc.pixel_shift_3d(img, dichroic_dictionary=config.dichroic_dict, json_file=config.pixelshift_json_path)
            img = rundmc.correct_shift_upsampling(img)
        img = rundmc.convert_16_to_8(img)

        io.imsave(os.path.join(output_path, os.path.basename(f)), img)
        # io.imshow(img[1,:,:])
        # io.show()
        if len(img.shape) == 3:
            img_list = [img[i] for i in range(img.shape[0])]
        else:
            img_list = [img]
        
        img_sizes = [i.shape for i in img_list]
        
        for n in range(len(img_list)):
            if img_list[n].shape != (config.resize_y, config.resize_x):
                reshaped_img = resize(img_list[n], (config.resize_y, config.resize_x), anti_aliasing=True)
            if config.gauss_level is not None:
                reshaped_img = gaussian(reshaped_img, sigma=config.gauss_level)
            img_list[n] = reshaped_img

        
        img_list = [np.expand_dims(img, axis=2) for img in img_list]
        X_test = np.empty(((0,) + img_list[0].shape))
        X_test = np.append(X_test, np.array(img_list), axis=0)

        data_test = [{'rawdata': i, 'element_size_um': EL_SIZE} for i in X_test]

        test_generator = preproc.TileGenerator(data = data_test,
                                       instancelabels=None,
                                       tile_shape=TILE_SHAPE,
                                       padding=PADDING,
                                       n_classes=2,
                                       border_weight_sigma_px=SIGMA_SEP,
                                       border_weight_factor=LAMBDA,
                                       foreground_background_ratio=V_BAL)

        model = unet.Unet2D(snapshot=config.model_path, 
                n_channels=1, 
                n_classes=2, 
                n_levels=4,
                batch_norm = BATCH_NORM,
                upsample=False,
                relu_alpha=0.1,
                n_features=64, name='U-Net')

        prediction = model.predict(test_generator)

        masks = []
        for ch in range(len(img_sizes)):
            mask = prediction[0][ch][...,1]
            if img_sizes[ch] != mask.shape:
                mask=resize(mask, img_sizes[ch], anti_aliasing=True)
            mask = utils.label_mask(mask, threshold=config.prediction_threshold, min_pixel=config.min_pix_per_cell)
            masks.append(mask)
            mask_save_path = os.path.join(output_path, (config.channel_list[ch]+'_mask.png'))
            io.imsave(mask_save_path, mask)

        # iou, max_a, max_b = utils.iou_mapping(masks[2], masks[1], min_roi_size=200)
        # print(iou, max_a, max_b)
        # print(iou.shape)
        all_poly_df = pd.DataFrame(columns=['cell_n', 'gene_name', 'row_pixels', 'col_pixels'])
        for n, gene_name in enumerate(config.channel_list):
            print('Making ROIs for: ', gene_name)
            polygons = rundmc.polygons_from_labels(masks[n], gene_name)
            all_poly_df = all_poly_df.append(polygons, ignore_index=True)
        # print(all_poly_df)
        csv_save_path = os.path.join(output_path, 'cell_coordinates.csv')
        all_poly_df.to_csv(csv_save_path)




        






if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    ### GET PATH DATA
    parser.add_argument('--image_path', type=str, default='./images/Ms_dapi_CN1500_gad594_slc546_pv647_30x_[1]_V1_BG.tif')
    parser.add_argument('--multi_img_dir', type=str)
    parser.add_argument('--output_path', type=str)
    parser.add_argument('--model_path', type=str, default='./df_resources/full_model.h5')

    ### GET VARIOUS ARGUMENTS
    parser.add_argument('--gauss_level', type=float, default=0.1)
    parser.add_argument('--resize_x', type=int, default=2040)
    parser.add_argument('--resize_y', type=int, default=2040)
    parser.add_argument('--channel_list', nargs='+', default=['dapi',
                                                                'syfp',
                                                                'gad1',
                                                                'vip',
                                                                'pvalb']) #enter into command line like (no commas): dapi syfp ...
    
    parser.add_argument('--allign_channels', type=bool, default=False)
    parser.add_argument('--pixelshift_json_path', type=str, default='./df_resources/shiftset.json')
    parser.add_argument('--dichroic_dict', action = type('', (argparse.Action, ), dict(__call__ = lambda a, p, n, v, o: getattr(n, a.dest).update(dict([v.split('=')])))),
                                                        default={0:'dmqb',
                                                                1:'dm4',
                                                                2:'dm4',
                                                                3:'dmqb',
                                                                4:'dmqb'}) # $ python argparse_dict_argument.py --arg a=b --arg aa=bb
                                                                            # Namespace(env={'a': 'b', 'aa': 'bb'})

    parser.add_argument('--prediction_threshold', type=float, default=0.5)
    parser.add_argument('--min_pix_per_cell', type=int, default=400)

    



    config = parser.parse_args()
    process(config)