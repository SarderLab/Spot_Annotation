"""


Generating 10x Visium Spatial Transcriptomics spot annotations

Inputs:
    - Spot coordinates = csv file with barcodes, image_row, and image_col columns
        - Actual image rows and columns values are flipped (image_row = image x direction, image_col = image y direction) in previous examples
    - Aligned Spot spatial transcriptomics information = csv file(s) with barcodes as column names and rows as cell type/gene acronyms
    - Cell types definitions = file containing main cell type, subtypes included, and their states for cell types
        - csv file with subtypes and states columns delimited by periods for each subtype and it's state within each main cell type
    - Output type = str value (case insensitive)
        - xml = Aperio ImageScope format
        - geojson = GeoJSON format
        - json = Histomics large-image annotation format

"""

import os
import sys
sys.path.append('..')
import json
import argparse
import numpy as np
import pandas as pd
import girder_client
from shapely.geometry import Polygon, Point
#from geojson import Feature, dump
from wsi_annotations_kit import wsi_annotations_kit as wak
import rpy2.robjects as robjects
from tqdm import tqdm
import os
import shutil


class SpotAnnotation:
    def __init__(self,
                 rds_file,
                 definitions_file: str,
                 image_id,
                 gc):
        
        self.rds_file = rds_file
        self.definitions_file = definitions_file

        self.image_id = image_id
        self.gc = gc

        # Reading in csv files from paths
        self.definitions = pd.read_csv(self.definitions_file)

        self.process_rds()

        # Determining MPP from coordinates
        self.mpp = self.calculate_mpp()

        self.barcodes = list(self.coordinates.index)
        print(f'number of barcodes: {len(self.barcodes)}')

        self.omics_data = self.process_omics()
        self.spot_annotations = self.process_spots()

        self.save()

    def calculate_mpp(self):

        # Finding minimum distance spots from first spot and using that to determine MPP
        # spot centroids are 100um apart and spot diameters are 55um
        spot_x_coords = self.coordinates['imagecol'].tolist()
        spot_y_coords = self.coordinates['imagerow'].tolist()

        # Test spot is the first one
        test_spot = [spot_x_coords[0],spot_y_coords[0]]
        #print(f'test_spot: {test_spot}')

        # Spot distances
        spot_dists = np.array([self.distance(test_spot, [x, y]) for x, y in zip(spot_x_coords, spot_y_coords)])
        #print(f'Number of spot dists: {np.shape(spot_dists)}')
        spot_dists = np.unique(spot_dists[spot_dists > 0])
        #print(spot_dists)
        #print(f'Number of unique spot_dists: {np.shape(spot_dists)}')
        min_spot_dist = np.min(spot_dists)
        #print(f'Minimum spot distance: {min_spot_dist}')

        # Minimum distance between the test spot and another spot = 100um (same as doing 100/min_spot_dist)
        mpp = 1/(min_spot_dist/100)
        #print(f'calculated mpp = {mpp}')

        return mpp

    def distance(self, point1, point2):
        # Distance between 2 points
        return (((point1[0]-point2[0])**2)+((point1[1]-point2[1])**2))**0.5

    def process_rds(self):

        # Processing RDS file containing spot coordinates and omics data
        robjects.r('library(Seurat)')
        robjects.r('library(stringr)')

        robjects.r('''
                    # Function to extract normalized cell type fractions from RDS file
                    get_cell_norms <- function(rds_file){
                        read_rds_file <- readRDS(rds_file)

                        if ("predsubclassl2" %in% names(read_rds_file@assays)){
                            cell_type_fract <- GetAssayData(read_rds_file@assays[["predsubclassl2"]])
                        } else if ("pred_subclass_l2" %in% names(read_rds_file@assays)){
                            cell_type_fract <- GetAssayData(read_rds_file@assays[["pred_subclass_l2"]])
                        }

                        # Normalizing so that the columns sum to 1
                        cell_type_fract <- cell_type_fract[1:nrow(cell_type_fract)-1,]
                        cell_type_norm <- sweep(cell_type_fract,2,colSums(cell_type_fract),'/')
                        cell_type_norm[is.na(cell_type_norm)] = 0
                   
                        spot_coords <- read_rds_file@images[["slice1"]]@coordinates
                   
                        # Getting UMI counts
                        umi_count <- data.frame(read_rds_file@meta.data[["nCount_Spatial"]],row.names=colnames(read_rds_file@assays[["Spatial"]]))

                        return(c(cell_type_norm,spot_coords,umi_count))
                    }
                    
                    ''')

        get_cell_norms = robjects.globalenv['get_cell_norms']
        cell_norm_output = get_cell_norms(self.rds_file)

        # Converting R dataframes to pandas dataframes
        with (robjects.default_converter + robjects.pandas2ri.converter).context():
            cell_type_norms = robjects.conversion.get_conversion().rpy2py(cell_norm_output[0])
            spot_coordinates = robjects.conversion.get_conversion().rpy2py(cell_norm_output[1])
            umi_counts = robjects.conversion.get_conversion().rpy2py(cell_norm_output[2])

        self.omics = cell_type_norms
        self.coordinates = spot_coordinates
        self.umi_counts = umi_counts

    def process_omics(self):

        # Aggregating counts information based on counts definition file
        sub_types_list = self.definitions['Sub_Types'].tolist()
        cell_states_list = self.definitions['Cell_States'].tolist()
        main_types_list = self.definitions['Main_Types'].tolist()

        # Initializing dictionary to hold compressed counts data
        # Keys will be:
        #   - main_cell_types = percentage of each main cell type for each spot
        #    - {main_cell_type}
        #       - pct_subtypes = percentage of each cell subtype present for each main cell type for each spot
        #       - pct_states = percentage of each cell state present for each main cell type for each spot

        slide_compressed_counts = {}
        slide_compressed_counts['main_cell_types'] = pd.DataFrame()

        for m_idx, main_name in enumerate(main_types_list):
            
            # Have to split the sub-types and states for each main cell type row
            sub_list = sub_types_list[m_idx].split('.')
            state_list = cell_states_list[m_idx].split('.')

            slide_compressed_counts[main_name] = {}
            
            # Narrowing down omics df to only rows containing a sub-type for this main cell type
            sub_df = self.omics.loc[self.omics.index.isin(sub_list)]
            
            # Percentage of this main cell type is the sum of all the normalized percentages of it's sub-types (total will still sum to 1
            sub_sum_df = sub_df.sum(axis=0)
            if slide_compressed_counts['main_cell_types'].empty:
                slide_compressed_counts['main_cell_types'] = pd.DataFrame(data=sub_sum_df.values, columns = [main_name],index=list(sub_sum_df.index))
            else:
                new_df = pd.DataFrame(data = sub_sum_df.values, columns=[main_name], index=list(sub_sum_df.index))
                slide_compressed_counts['main_cell_types'] = pd.concat([slide_compressed_counts['main_cell_types'],new_df],axis=1)

            # Normalizing the sub-types sum dataframe gets percentage of each subtype
            pct_count_df = (sub_df/sub_sum_df).fillna(0)
            slide_compressed_counts[main_name]['pct_subtypes'] = pct_count_df

            # Initializing state percentage dataframe
            state_pct_df = pct_count_df.copy()

            # Setting the index of the dataframe only if the sub-type that has that state is present
            state_pct_df.index = [state_list[i] for i in range(len(state_list)) if sub_list[i] in list(state_pct_df.index)]

            # Combining states with the same name
            state_pct_df = state_pct_df.groupby(level=0).sum()

            # Re-normalizing if certain cell states were removed, replacing inf/nan with zero
            state_pct_df = (state_pct_df/state_pct_df.sum(axis=0)).fillna(0)

            # Sorting index for consistency
            state_pct_df = state_pct_df.sort_index()

            slide_compressed_counts[main_name]['pct_states'] = state_pct_df
    
        # Adding in rows for missing main cell types
        if slide_compressed_counts['main_cell_types'].shape[1]<len(main_types_list):
            add_column_names = [i for i in main_types_list if i not in slide_compressed_counts['main_cell_types'].columns]
            for add in add_column_names:
                added_df = pd.DataFrame({add: [0]*slide_compressed_counts['main_cell_types'].shape[0]}, index=list(slide_compressed_counts['main_cell_types'].index))
                slide_compressed_counts['main_cell_types'] = pd.concat([slide_compressed_counts['main_cell_types'], added_df], axis=1)

                # Add these into pct_states and pct_subtypes?

        slide_compressed_counts['main_cell_types'] = slide_compressed_counts['main_cell_types'].sort_index(axis=1)

        return slide_compressed_counts

    def process_spots(self):

        # Initializing annotations using wak
        spot_annotations = wak.Annotation(mpp=self.mpp)
        spot_annotations.add_names(['Spots'])

        # Iterating through barcodes and creating equal sized spots centered on coordinates
        spot_pixel_diameter = int((1/self.mpp)*55)
        spot_pixel_radius = int(spot_pixel_diameter/2)
        for b in tqdm(self.barcodes):

            # Pulling out that row from coordinates
            b_coords = self.coordinates.loc[b]
            b_x = b_coords['imagecol']
            b_y = b_coords['imagerow']
            
            # Creating shapely annotation for this spot
            spot_poly = Point(b_x, b_y).buffer(spot_pixel_radius)

            # Extracting -omics information from previous processing step
            spot_main_cell_types = self.omics_data['main_cell_types'].loc[b].to_frame()
            spot_main_cell_types.columns = ['Main_Cell_Types']
            spot_main_cell_types_dict = spot_main_cell_types.to_dict()

            # Cell State info
            spot_cell_state_dict = {}
            spot_cell_subtype_dict = {}
            for m in list(spot_main_cell_types_dict['Main_Cell_Types'].keys()):

                spot_cell_state = self.omics_data[m]['pct_states'][b].to_frame()
                spot_cell_state.columns = [m]
                spot_cell_state_dict[m] = spot_cell_state.to_dict()[m]

                spot_cell_subs = self.omics_data[m]['pct_subtypes'][b].to_frame()
                spot_cell_subs.columns = [m]
                spot_cell_subs_dict = spot_cell_subs.to_dict()
                spot_cell_subtype_dict[m] = spot_cell_subs_dict[m]

            # Properties written to annotations = 
            # Main_Cell_Types = aggregated cell subtypes
            # Cell_States = cell state distribution for each main cell type
            # Cell_Subtypes = cell subtype distribution for each main cell type
            # All_Subtypes = raw cell subtype distribution
            # UMI_Count = QC metric for number of reads recorded for each spot
            spot_properties = {
                'Main_Cell_Types': spot_main_cell_types_dict['Main_Cell_Types'],
                'Cell_States': spot_cell_state_dict,
                'Cell_Subtypes': spot_cell_subtype_dict,
                'All_Subtypes': self.omics[b].to_dict(),
                'UMI_Count': self.umi_count.loc[b].values
                }
            # Adding spot to annotations, crs is just the origin since we are using non-scaled centroid points
            # name is set to be the barcode
            spot_annotations.add_shape(
                poly=spot_poly,
                box_crs=[0, 0],
                structure='Spots',
                name=b,
                properties=spot_properties)

        # Print str property of Annotation object
        print(spot_annotations)
        return spot_annotations

    def save(self):

        annot = wak.Histomics(self.spot_annotations)
        _ = self.gc.post(f'annotation/item/{self.image_id}',
                        data=json.dumps(annot.json),
                        headers={
                            'X-HTTP-Method': 'POST',
                            'Content-Type':'application/json'
                            }
                        )
        print('uploating layers')
        print('annotation uploaded...\n')

