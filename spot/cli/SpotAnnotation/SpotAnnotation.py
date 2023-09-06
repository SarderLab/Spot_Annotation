import os
import sys
from ctk_cli import CLIArgumentParser


sys.path.append("..")


def main(args):  
    
    cmd = "python3 ../spot_code/Generate_Spot_Annotations.py   --basedir {}  --rds_file {} --definitions_file {}  --girderApiUrl {} --girderToken {} \
             --input_files {}".format(args.basedir, args.rds_file, args.definitions_file, args.girderApiUrl, args.girderToken, args.input_files)
    print(cmd)
    sys.stdout.flush()
    os.system(cmd)  


if __name__ == "__main__":
    main(CLIArgumentParser().parse_args())