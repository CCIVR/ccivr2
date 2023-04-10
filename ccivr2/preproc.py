import argparse
import os
from ccivr2.__init__ import __version__

def getargs():
    
    parser = argparse.ArgumentParser(
        description='Identification of various kinds of antisense-oriented transcripts including total cis-NATs, nearby head-to-head, and nearby tail-to-tail based on genome annotation'
        )

    parser.add_argument('input',
        help='the path of the CSV file to read'
        )

    parser.add_argument('-o', '--output',
        help=
            'the path of the directory to save results, '
            'if not specified, files would be stored in the same location as input file.'
        )

    parser.add_argument('-v', '--version', action='version',
        version=__version__,
        help='show the version number and exit'
        )

    args = parser.parse_args()

    return args


def getmode():

    print(
        "Mode 1: CisNats (embedded, fully-overlapped, head-to-head, tail-to-tail)" + '\n'
        "Mode 2: TSS comparison" + '\n'
        "Mode 3: TES comparison"
    )

    # initial value (dummy)
    tss_range = ["0","0"]
    tes_range = ["0","0"]

    loop1 = 0
    loop2 = 0

    def test_bigsmall(r):
        if int(r[0]) < int(r[1]):
            return 1
        else:
            print("Incorrect answer: the following expression must be satisfied; min < max")
            return 0

    while loop1 == 0:
        ans = input("Select identification mode [1/2/3]: ")

        if ans == "1":
            loop1 = 1
            
        elif ans == "2":
            while loop2 == 0:
                tss_range = input("Set the range of AS_TSS position relative to S_TSS [min,max]: ").split(",")
                loop2 = test_bigsmall(tss_range)

            loop1 = 1

        elif ans == "3":
            while loop2 == 0:
                tes_range = input("Set the range of AS_TES position relative to S_TES [min,max]: ").split(",")
                loop2 = test_bigsmall(tes_range)

            loop1 = 1

        else:
            print("Incorrect answer")

    class ModeInformation:
        def __init__(self,mode,range_tss_distance:list,range_tes_distance:list):
            self.mode = mode
            self.tss = range_tss_distance
            self.tes = range_tes_distance

    setting = ModeInformation(
        mode = ans,
        range_tss_distance = tss_range,
        range_tes_distance = tes_range
        )

    return setting


def get_prior_inf():

    class Paths:
        def __init__(self,input,output:list):
            self.input = input
            self.output = output

    class SettingInformation:
        def __init__(self,paths:Paths,mode,tss_range,tes_range):
            self.paths = paths
            self.mode = mode
            self.tss = tss_range
            self.tes = tes_range

    args = getargs()
    modeset = getmode()

    # form path
    # input
    input_file = os.path.abspath(args.input)

    # output:location
    if args.output:
        output_dir = os.path.abspath(args.output)
    else:
        output_dir = os.path.dirname(input_file)
    # ourput:name
    if modeset.mode == "1":
        terms = modeset.mode
    elif modeset.mode == "2":
        terms = modeset.mode + '_({r})'.format(r=','.join(modeset.tss))
    elif modeset.mode == "3":
        terms = modeset.mode + '_({r})'.format(r=','.join(modeset.tes))

    dirname = 'ccivr_output_mode' + terms

    ccivr_output = os.path.join(output_dir, dirname)
    os.makedirs(ccivr_output, exist_ok=True)

    output_file_table = os.path.join(ccivr_output, 'Table.csv')
    output_file_summary = os.path.join(ccivr_output, 'Summary.csv')

    path_inf = Paths(
        input = input_file, 
        output = [ccivr_output, output_file_table, output_file_summary]
        )

    inf = SettingInformation(
        paths = path_inf,
        mode = modeset.mode,
        tss_range = modeset.tss,
        tes_range = modeset.tes
        )

    return inf