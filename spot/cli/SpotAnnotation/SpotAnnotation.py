import os
import sys
from ctk_cli import CLIArgumentParser
import girder_client

sys.path.append("..")

from spot_code.Generate_Spot_Annotations import SpotAnnotation

def main(args):  

    # Defining main inputs into SpotAnnotation from arguments
    basedir = args.basedir
    rds_file = args.rds_file
    definitions_file = args.definitions_file
    girderApiUrl = args.girderApiUrl
    girderToken = args.girderToken
    input_files = args.input_files

    # simple processing of input arguments
    image_name = input_files.split(os.sep)[-1]
    gc = girder_client.GirderClient(apiUrl=girderApiUrl)
    gc.setToken(girderToken)

    # Getting image id
    print(f'rds_file: {rds_file}')
    print(f'basedir: {basedir}')
    print(f'definitions_file: {definitions_file}')
    print(f'input_files: {input_files}')
    folder_items = gc.get(f'/resource/{basedir}/items?token={girderToken}',parameters={'type':'folder','limit':10000})
    item_names = [i['name'] for i in folder_items]

    image_id = folder_items[item_names.index(image_name)]['_id']
    # instantiating SpotAnnotation object, automatically outputs annotations to image_name
    SpotAnnotation(rds_file,definitions_file,image_id,gc)


if __name__ == "__main__":
    main(CLIArgumentParser().parse_args())