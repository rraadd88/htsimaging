## test examples using papermill
import logging
from pathlib import Path
import os
from os.path import exists

from testbook import testbook
import papermill as pm

print(">>>>>>>>>>>>>>>>> current dir:"+os.getcwd())
if not os.getcwd().endswith("/examples"):
    os.chdir('./examples/')
    print(">>>>>>>>>>>>>>>>> current dir:"+os.getcwd())

@testbook("_make_demo_data.ipynb", execute=True)
def test_make_demo_data(tb):
    pass # execute only because tests are present in the notebook itself
    return

# parameters
input_dir_path='./inputs/'

def test_protein_abundance_and_normalization(
    input_path='./protein_abundance_and_normalization.ipynb',
    ):
    output_dir_path= './outputs/' + Path(input_path).stem + '/'
    os.makedirs(output_dir_path,exist_ok=True)

    parameters=dict(
    ## parameters
        input_path=f'{input_dir_path}/image_intensity_cells.npy',
        segmented_image_path=f'{input_dir_path}/image_segmented_cells.npy',
        output_path=f'{output_dir_path}/01_gfpby_cell.tsv',
        misaligned_fraction_max=0.9,
        edge_off=20,
    )
    
    logging.info(parameters)
    Path(output_dir_path).mkdir(exist_ok=True)

    _=pm.execute_notebook(
        input_path=input_path,
        output_path='./outputs/' + Path(input_path).name,
        parameters=parameters,
        # report_mode=True,
        kernel_name='test',
        execution_timeout=360,
    )
    assert exists(parameters['output_path'])
    
def test_protein_abundance_by_marker_location(
    input_path='./protein_abundance_by_marker_location.ipynb',
    ):
    output_dir_path= './outputs/' + Path(input_path).stem + '/'
    os.makedirs(output_dir_path,exist_ok=True)

    parameters=dict(
    ## parameters
        ## parameters
        input_path=f'{input_dir_path}/image_marker_cells.npy', ## path to the image from channel marking a subcellular localization
        output_path=f'{output_dir_path}/01_gfpby_cell_and_marker.tsv',

        image_intensity_path=f'{input_dir_path}/image_intensity_cells.npy', ## path to the image with channel which is to be used to caluculate the abundance
        regions_path=f'{input_dir_path}/image_segmented_cells.npy', ## segmented regions (dtype: bool)
        marker_intensity_min_quantile=0.975,
        pixels_per_cell_min=100,
        non_marker_intensity_quantile_off=0.02,

        force=False,
        test=True,
    )

    logging.info(parameters)
    Path(output_dir_path).mkdir(exist_ok=True)

    _=pm.execute_notebook(
        input_path=input_path,
        output_path='./outputs/' + Path(input_path).name,
        parameters=parameters,
        # report_mode=True,
        kernel_name='test',
        execution_timeout=360,
    )
    assert exists(parameters['output_path'])
