import os
import sys
from ctk_cli import CLIArgumentParser
import girder_client

sys.path.append("..")

from spot_code.Generate_Spot_Annotations import SpotAnnotation

def main(args):  

    # Defining main inputs into SpotAnnotation from arguments

    counts_file = args.counts_file
    organ = args.organ
    definitions_file = args.definitions_file
    girderApiUrl = args.girderApiUrl
    girderToken = args.girderToken
    input_files = args.input_files
    gene_selection_method = {
        'method': args.gene_selection_method,
        'n': args.n,
        'list': args.gene_list
    }

    # simple processing of input arguments
    gc = girder_client.GirderClient(apiUrl=girderApiUrl)
    gc.setToken(girderToken)

    # Getting image id
    for a in vars(args):
        print(f'{a}: {getattr(args,a)}')

        
    # instantiating SpotAnnotation object, automatically outputs annotations to image_name
    SpotAnnotation(
        counts_file,
        definitions_file,
        input_files,
        gc,
        organ,
        gene_selection_method
    )


if __name__ == "__main__":
    main(CLIArgumentParser().parse_args())