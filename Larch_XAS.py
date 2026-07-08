"""
File: Larch_XAS
Name: Cheng-Chu Chung
----------------------------------------
TODO: XAS data processing
Github source: https://github.com/xraypy/xraylarch
Color palettes for Python: https://jiffyclub.github.io/palettable/#palette-interface
"""
print('Import packages...')
import time
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
print('Import xraydb package...')
from xraydb import guess_edge, xray_edge
print('Import larch package...')
from larch.xafs import pre_edge, find_e0
from larch.io import read_ascii, write_ascii, merge_groups, read_athena, create_athena, write_group
import athena_project
from athena_project import read_athena
import palettable as pltt
from collections import defaultdict
from pprint import pprint
import configparser
import matplotlib
import re
from colorama import Fore, Back, Style, init
print("Before, Backend used by matplotlib is: ", matplotlib.get_backend())
matplotlib.rcParams['backend'] = 'wxAgg'
print("After, Backend used by matplotlib is: ", matplotlib.get_backend())

# Constant
"""
FILE_TYPE instructions:
'.prj' for merging fluorescence scans
'' or '.dat' for merging single scan
'.txt' for plotting scans and you will need to have a .ini file to set up the plotting parameters
---------------------------------------------
TRANSMISSION_MODE instructions:
'Auto' is default to process the transmission data based on the notation in the raw file
Otherwise True or False to decide what type of data you want to export (transmission or fluorescence scans)
True for transmission scans; False for fluorescence scans
"""
FILE_TYPE = '.txt'   # <------------------------------------------------------------------------------------- data type
TRANSMISSION_MODE = 'Auto'
INPUT_PATH = r"G:\Other computers\我的 筆記型電腦 MSI\Research data\SSID\202411\20241104 BMM AE\Ni-b59-03-CrCuNiCr-AE\Output_files"    # <----------------------- Data folder input
OUTPUT_PATH = Path(f'{INPUT_PATH}\Output_files')

# Merged Constant
SKIP_SCANS = ['Mn_b58_04_ScVMnSc_AE_19p858000000000004_300_002']  # [] if scans are good or just add scans you want to exclude
IF_NOR = False                                     # Do normalization for fluorescence scans
ADD_DEV = False                                    # Add plus and minus standard deviation lines for fluorescence scans
SHOW_DATA_INFORMATION = False                      # List athena parameters, such as atomic symbol, edge, label, etc.

# Plot Constant for .txt
"""
You could set FILE_INDEX = 0, SAMPLE_LIST = [], STANDARD_LIST = [], 
SAMPLE_LABEL = [], ENERGY_RANGE = () as a default for your first try.
"""
CONFIG_FILE = r"G:\Other computers\我的 筆記型電腦 MSI\Research data\SSID\202411\20241104 BMM AE\Ni-b59-03-CrCuNiCr-AE\Output_files\Ni-b59-03-CrCuNiCr-AE.ini"   # <-------------------- .ini setting for plotting or leave it blank for data preprocessing
IF_SAVE = False  # Save the plot or not, so you can set IF_SAVE=False if you don't want to save the plot

config = configparser.ConfigParser()
if Path(CONFIG_FILE).is_file():
    config.read(CONFIG_FILE, encoding="utf8")
    is_ini = True
    print("=================")
    print("Use .ini input")
    print("-----------------")
else:
    is_ini = False
    print("=================")
    print("Manually input so please remove all eval commands if error occurs or file is not found")
    print("-----------------")

FILE_INDEX = config.getint('samples', 'file_index') if is_ini else 0                         # Which file in the file list you want to plot
SAMPLE_LIST = eval(config['samples']['sample_list']) if is_ini else []                                    # [] for default or [1, 7, 5, 3] for a index list you want to plot
STANDARD_LIST = eval(config['samples']['standard_list']) if is_ini else []                                # [] if none or [5, 3] in the SAMPLE_LIST become dash lines
SAMPLE_LABEL = eval(config['legends']['sample_label']) if is_ini else []                                  # [] for default or add a specific name list
FIGURE_SIZE = eval(config['format']['figure_size']) if is_ini else (6.4, 4.8)                             # Cheng-Hung uses (6, 7.5), but the default is (6.4, 4.8)
PALETTE = eval(config['format']['palette']) if is_ini else pltt.colorbrewer.diverging.Spectral_4_r        # pld.Spectral_4_r  # _r if you want to reverse the color sequence
CMAP = PALETTE.mpl_colormap                                                                               # .mpl_colormap attribute is a continuous, interpolated map
LINEWIDTH = eval(config['format']['linewidth']) if is_ini else 2                                          # Line width
FRAMELINEWIDTH = eval(config['format']['framelinewidth']) if is_ini else 2                                # Frame line width
COLOR_INCREMENT = eval(config['format']['color_increment']) if is_ini else 0
OFFSET = eval(config['format']['offset']) if is_ini else 0                                                # Value you want to add to an y offset for each curve.
ENERGY_RANGE = eval(config['format']['energy_range']) if is_ini else ()                                   # () for default, (18900, 19150) for Nb, (4425, 4625) for Sc
Y_RANGE = eval(config['format']['y_range']) if is_ini else ()                                             # () for default
ENERGY_INTERVAL = eval(config['format']['energy_interval']) if is_ini else 0                              # This parameter works only when you set a ENERGY_RANGE
# IF_SAVE = eval(config['format']['if_save']) if is_ini else True                                           # Save the plot or not
OUTPUT_FILENAME = eval(config['format']['output_filename']) if is_ini else "Default"
NUM_COLUMN = 1
DETECTOR_INDEX_HEAD, DETECTOR_INDEX_TAIL = 7, 14    # 2024 cycle 3 updated the number of detectors from 4 to 7
ADD_REFERENCE = True # Add a reference prj file for each sample, so you can set add_reference=False if you don't want to have it
IF_DELETE_SINGLE_PRJ = True  # Delete single prj files after merging, so you can set IF_DELETE_SIGNLE_PRJ=False if you want to keep them

def main():
    files = Path(INPUT_PATH).glob(f'*{FILE_TYPE}')
    # new_merge_project = create_athena(f'default.prj')     # Call the athena_project.py in larch.io
    new_merge_project = athena_project.create_athena(f'default.prj')  # Call athena_project.py in current folder

    if FILE_TYPE != '.txt':
        if not OUTPUT_PATH.exists():
            OUTPUT_PATH.mkdir()    # Create an output folder to save all generated data/files

    if FILE_TYPE == '.prj':
        process_prj_file(files, new_merge_project)

    elif FILE_TYPE == '.txt':
        plot_xas(files)

    else:
        # new_merge_project = athena_project.create_athena(f'default.prj')  # Call athena_project.py in current folder

        # Create each prj
        create_transmission_prj(files, add_reference=ADD_REFERENCE)  # Sometimes the reference has, so you can set add_reference=False

        # Create a group prj containing all sample data
        create_transmission_prj_group(new_merge_project)

        if IF_DELETE_SINGLE_PRJ:
            print("\n==============================")
            print('Delete single prj files')
            print("------------------------------")
            delete_multiple_file_types(OUTPUT_PATH, extensions=["*merged*", "*reference*"])


def plot_xas(files):
    """
    :param files: txt, a txt file with organized dataset output from Athena
    :return: None
    """
    f_list = []
    print("==============================")
    print('Files')
    print("------------------------------")

    for index, file in enumerate(files):
        f_list.append(file)
        print(index, file)

    print("\n==============================")
    print(f'Data column in file number {FILE_INDEX} ---> {f_list[FILE_INDEX].name}')
    print("------------------------------")
    file = read_ascii(f_list[FILE_INDEX])
    file_keys = file.__dir__()  # Make a list
    for index, key in enumerate(file_keys):
        print(f'{index:>2} {key}')

    energy = getattr(file, file_keys[6])

    # Do the plotting
    f1, ax1 = plt.subplots(1, 1, figsize=FIGURE_SIZE)
    print("\n==============================")
    print('Index   Filename')
    print("------------------------------")
    increment = 0   # Increment for offset
    if len(SAMPLE_LIST) == 0:
        color_idx = np.linspace(0, 1, len(file_keys)-7+COLOR_INCREMENT)   # All plots have their own color
        for i in range(len(file_keys)-7):   # Start from the first sample name because file key 0-6 are data information
            sample_index = i + 7
            sample_name = file_keys[sample_index]
            sample_label = sample_name
            # if len(SAMPLE_LABEL) > i:
            #     sample_label = SAMPLE_LABEL[i]
            # else:
            #     sample_label = sample_name
            mu = getattr(file, file_keys[sample_index])
            ax1.plot(energy, mu + OFFSET * increment, linewidth=LINEWIDTH, color=CMAP(color_idx[increment+COLOR_INCREMENT]), label=sample_label)
            increment += 1
            print('{:>2} {}'.format(sample_index, sample_name))
    else:
        color_idx = np.linspace(0, 1, len(SAMPLE_LIST)+COLOR_INCREMENT)   # Only the plots you want their own color
        for sample_index in SAMPLE_LIST:
            sample_name = file_keys[sample_index]
            if len(SAMPLE_LABEL) > SAMPLE_LIST.index(sample_index):
                sample_label = SAMPLE_LABEL[SAMPLE_LIST.index(sample_index)]
            else:
                sample_label = sample_name
            mu = getattr(file, file_keys[sample_index])
            if sample_index in STANDARD_LIST:
                ax1.plot(energy, mu + OFFSET * increment, '--', linewidth=LINEWIDTH, color=CMAP(color_idx[increment+COLOR_INCREMENT]), label=sample_label)
            else:
                ax1.plot(energy, mu + OFFSET * increment, linewidth=LINEWIDTH, color=CMAP(color_idx[increment+COLOR_INCREMENT]), label=sample_label)
            increment += 1
            print('{:>2} {}'.format(sample_index, sample_name))

    # Frame linewidth
    spineline = ['left', 'right', 'top', 'bottom']
    for direction in spineline:
        ax1.spines[direction].set_linewidth(FRAMELINEWIDTH)

    # Plot format
    if ENERGY_RANGE == ():
        ax1.set_xlim(energy.min() // 1 + 1, energy.max() // 1 - 1)
        plt.xticks(fontsize=14)
    else:
        ax1.set_xlim(ENERGY_RANGE)
        plt.xticks(np.arange(ENERGY_RANGE[0], ENERGY_RANGE[1], step=ENERGY_INTERVAL), fontsize=14)
    if Y_RANGE != ():
        ax1.set_ylim(Y_RANGE)
    plt.title(OUTPUT_FILENAME, fontsize=20, pad=15) \
        if OUTPUT_FILENAME != "" \
        else plt.title(Path(CONFIG_FILE).stem, fontsize=20, pad=15)
    x_label = r'$\mathregular{Energy\ (eV)}$'
    y_label = r'$\mathregular{Normalized\ \mu(E)}$'     # XANES normalization doesn't have x
    # plt.yticks([])  # Disable ticks
    ax1.tick_params(width=FRAMELINEWIDTH)
    ax1.set_xlabel(x_label, fontsize=18)
    ax1.set_ylabel(y_label, fontsize=18, labelpad=10)
    # plt.rcParams["axes.linewidth"] = 5
    plt.legend(loc='lower right', framealpha=1, frameon=False, fontsize=14, ncol=NUM_COLUMN, reverse=True)
    plt.tight_layout()
    if IF_SAVE:
        plt.yticks([])  # Disable ticks
        config_file_location = Path(CONFIG_FILE).parent
        output_filename = check_filename_repetition(OUTPUT_FILENAME, config_file_location) \
            if OUTPUT_FILENAME != "" \
            else check_filename_repetition(Path(CONFIG_FILE).stem, config_file_location)
        plt.savefig("{}/{}.png".format(Path(config_file_location), output_filename), dpi=300, transparent=False)
    plt.show()


def process_prj_file(files, new_merge_project):
    # new_merge_project = create_athena(f'default.prj')     # Call the athena_project.py in larch.io
    # new_merge_project = athena_project.create_athena(f'default.prj')  # Call athena_project.py in current folder

    for index, file_prj in enumerate(files):
        if 'Created' not in file_prj.name:
            merge_scan(file_prj, new_merge_project)

    new_merge_project.save(f'{OUTPUT_PATH}/Created_group.prj')
    print('\n=================================================================================')
    print(f'Save merge project into ---> Created_group.prj')
    print('=================================================================================')

    # Call Path again to grab the created athena project
    files = OUTPUT_PATH.glob(f'*.prj')
    calibrate_energy(files)


def merge_scan(file_prj, new_merge_project):
    """
    :param file_prj: prj file, a prj file from the current folder
    :param new_merge_project: prj file, an empty prj file to store merged data
    :return: None
    """
    filename = file_prj.name
    print(f'\nWhether the file "{filename}" exists?', file_prj.exists())

    # Read Athena project
    scans = read_athena(f'{file_prj}')
    scans_grouplist = []    # Each scan is a group

    # Print each scan name in the prj file
    print("\n==============================")
    print('Scan name')
    print("------------------------------")
    for name, group in scans._athena_groups.items():
        if name not in SKIP_SCANS:
            scans_grouplist.append(group)
            print(name, group)

    # Print scan information
    first_scan_information = scans_grouplist[0]
    if SHOW_DATA_INFORMATION:
        print("\n==============================")
        print(f'Athena parameters in {first_scan_information.label}')
        print("------------------------------")

        show_data_information(first_scan_information)

    # Plot scans
    fig, ax = plt.subplots()
    for index, scan in enumerate(scans_grouplist):
        if IF_NOR:  # Do normalization
            pre_edge(scan.energy, scan.mu, group=scan)
            plt.plot(scan.energy, scan.flat, label=scan.label)
        else:
            plt.plot(scan.energy, scan.mu, label=scan.label)

    # Do the merge and plot
    merges = merge_groups(scans_grouplist)
    if IF_NOR:  # Do normalization
        pre_edge(merges.energy, merges.mu, group=merges)    # <---------- Automatically define white region
        plt.plot(merges.energy, merges.flat, label=f'{first_scan_information.label[:-4]}_merged')
        if ADD_DEV:
            plt.plot(merges.energy, merges.flat + merges.mu_std * merges.flat / merges.mu, '-',
                     label=f'{first_scan_information.label[:-4]}_merged+std')
            plt.plot(merges.energy, merges.flat - merges.mu_std * merges.flat / merges.mu, '-',
                     label=f'{first_scan_information.label[:-4]}_merged-std')
    else:
        # Add e0, atom, and edge information
        e0 = find_e0(merges.energy, mu=merges.mu, group=merges)
        atsym, edge = guess_edge(merges.e0)
        merges.atsym = atsym
        merges.edge = edge

        plt.plot(merges.energy, merges.mu, label=f'{first_scan_information.label[:-4]}_merged')
        if ADD_DEV:
            plt.plot(merges.energy, merges.mu + merges.mu_std, '-',
                     label=f'{first_scan_information.label[:-4]}_merged+std')
            plt.plot(merges.energy, merges.mu - merges.mu_std, '-',
                     label=f'{first_scan_information.label[:-4]}_merged-std')

    if SHOW_DATA_INFORMATION:
        print("\n==============================")
        print(f'Merged parameters in {first_scan_information.label[:-4]}_merged')
        print("------------------------------")

        show_data_information(merges)

        print('Atomic system:', merges.atsym)
        print('E0:', merges.e0)
        print('edge:', merges.edge)
        print('Xray edge:', xray_edge(merges.atsym, merges.edge)[0])

    # Replace '-' with '_' because '-' will cause error
    scan_name = f'{first_scan_information.label[:-4]}_merged'.replace('-', '_').replace(' ', '_').replace('.', 'p')
    new_merge_project.add_group(merges, scan_name)  # Add the merge into the new prj

    # Plotting format
    ax.set_xlim((merges.energy.min() // 1 + 1, merges.energy.max() // 1 - 1))
    plt.xlabel('$\mathregular{Energy\ (eV)}$', fontsize=12)
    plt.ylabel('$\mathregular{\chi\mu(E)}$', fontsize=12)
    plt.title(f'{first_scan_information.atsym} {first_scan_information.edge}-edge')
    plt.legend(loc='lower right', framealpha=1, frameon=False)
    if IF_SAVE:
        plt.savefig("{}/{}.png".format(Path(INPUT_PATH), first_scan_information.label[:-4]), dpi=300, transparent=False)
        print('\n=================================================================================')
        print(f'Save merge image into ---> {first_scan_information.label[:-4]}.png')
        print('=================================================================================')
    plt.close()


def create_transmission_prj(files, add_reference=True):
    """
    :param files: txt file, a txt file from the current folder
    :return: None
    """
    scan_dictionary = defaultdict(list)
    print("==============================")
    print('Index Files')
    print("------------------------------")

    # Create scan dictionary and each item contains energy, reference, scan1, scan2, scan3, etc...
    for index, scan in enumerate(files):
        scan = scan.resolve()  # Make the path absolute, resolving any symlinks
        scanname = scan.name

        if scanname[-3:].isnumeric() or scanname[-3:] == 'dat':   # <---- file type .001, .002, .003, or 0001.dat, etc.
            print(index, scanname)
            scan = read_ascii(scan)
            scan_header = str(scan.header)
            scan_plot_hint_index = scan_header.find('# Scan.plot_hint')
            plot_hint = scan_header[scan_plot_hint_index: scan_header.find("'", scan_plot_hint_index)]

            if SHOW_DATA_INFORMATION:
                print("\n==============================")
                print(f'Scan attributes in {scanname}')
                print("------------------------------")
                print('Attributes:', scan.__dir__())
                print('Header:')
                pprint(scan._members()['header'])
                print('Data columns:', scan.array_labels)
                show_data_information(scan)

            # Append energy, mu!!!
            if scanname[-3:] == 'dat':
                space_index = scanname.find(' ', -11)
                sample_name = scanname[:space_index].replace('-', '_').replace('(', '').replace(')', '').replace(' ',
                                                                                                                 '_').replace('.', 'p')

                if f'{sample_name}_energy_mu' not in scan_dictionary:
                    scan_dictionary[f'{sample_name}_energy_mu'] = []
                    scan_dictionary[f'{sample_name}_energy_mu'].append(scan.energy)                   # <--- Energy
                    scan_dictionary[f'{sample_name}_energy_mu'].append(np.log(scan.it / scan.ir))     # <--- Reference
                tens_digit = int(scanname[space_index + 1:-4]) // 10 * 10
                units_digit = int(scanname[space_index + 1:-4]) % 10
                if f'{sample_name}_00{tens_digit + units_digit}' not in SKIP_SCANS:
                    scan_dictionary[f'{sample_name}_energy_mu'].append(np.log(scan.i0 / scan.it))   # <--- Transmission

            else:
                sample_name = scanname[:-4].replace('-', '_').replace('(', '').replace(')', '').replace(' ', '_').replace('.', 'p')
                if f'{sample_name}_energy_mu' not in scan_dictionary:
                    scan_dictionary[f'{sample_name}_energy_mu'] = []
                    scan_dictionary[f'{sample_name}_energy_mu'].append(scan.energy)                   # <--- Energy
                    scan_dictionary[f'{sample_name}_energy_mu'].append(np.log(scan.it / scan.ir))     # <--- Reference
                tens_digit = int(scanname[-2:]) // 10 * 10
                units_digit = int(scanname[-2:]) % 10
                if f'{sample_name}_00{tens_digit + units_digit}' not in SKIP_SCANS:
                    if TRANSMISSION_MODE == 'Auto':
                        if 'ln' in plot_hint:
                            scan_dictionary[f'{sample_name}_energy_mu'].append(
                                np.log(scan.i0 / scan.it))  # <------------------------------------------- Transmission
                        else:
                            print('Detectors:')
                            fluorescence_total_counts = 0
                            for if_index in range(DETECTOR_INDEX_HEAD, DETECTOR_INDEX_TAIL):   # 2024 cycle 3 updated number of detectors from 4 to 7
                                fluorescence_detector = scan.array_labels[if_index]
                                print(fluorescence_detector)
                                fluorescence_counts = getattr(scan, fluorescence_detector)
                                fluorescence_total_counts += fluorescence_counts
                            scan_dictionary[f'{sample_name}_energy_mu'].append(
                                fluorescence_total_counts / scan.i0)  # <--------------------------------- Fluorescence
                    else:       # User may want to decide which data type they want to export
                        if TRANSMISSION_MODE:
                            scan_dictionary[f'{sample_name}_energy_mu'].append(
                                np.log(scan.i0 / scan.it))  # <------------------------------------------- Transmission
                        else:
                            print('Detectors:')
                            fluorescence_total_counts = 0
                            for if_index in range(DETECTOR_INDEX_HEAD, DETECTOR_INDEX_TAIL):   # 2024 cycle 3 updated number of detectors from 4 to 7
                                fluorescence_detector = scan.array_labels[if_index]
                                print(fluorescence_detector)
                                fluorescence_counts = getattr(scan, fluorescence_detector)
                                fluorescence_total_counts += fluorescence_counts
                            scan_dictionary[f'{sample_name}_energy_mu'].append(
                                fluorescence_total_counts / scan.i0)  # <--------------------------------- Fluorescence
                    print('')

    # Append merged data, so each item will contain energy, reference, scan1, scan2, scan3, etc... and a merged scan.
    print("\n==============================")
    print('Scan plots and you may copy the scan name below into the SKIP_SCANS to exclude the scan')
    print("------------------------------")

    for sample_data in scan_dictionary:
        fig, ax = plt.subplots(1, 1, figsize=FIGURE_SIZE)   # Each figure should have its own format
        merge = np.mean(scan_dictionary[sample_data][2:], axis=0)   # <------------------------------- Merged array
        scan_dictionary[sample_data].append(merge)
        energy = scan_dictionary[sample_data][0]        # <------------------------------------------- Energy array
        reference = scan_dictionary[sample_data][1]     # <------------------------------------------- Reference array

        sample_name = sample_data[:-10]

        # Find the float-like part using regex
        match = re.search(r'\d+p\d+', sample_name)
        if match:
            float_str = match.group().replace('p', '.')
            rounded = round(float(float_str), 4)  # round to 3 decimal places
            rounded_str = str(rounded).replace('.', 'p')  # format back to original style
            sample_name = sample_name.replace(match.group(), rounded_str)
            print(sample_name)
        else:
            print("No float-like part found.")

        write_ascii("{}/{}_merged.prj".format(Path(OUTPUT_PATH), f'{sample_name}'), energy, merge,
                    label='energy mu',
                    header=['energy', 'mu'])
        if add_reference:
            write_ascii("{}/{}_reference.prj".format(Path(OUTPUT_PATH), f'{sample_name}'), energy, reference,
                        label='energy mu',
                        header=['energy', 'mu'])

        # Plot
        skip_times = 0
        for mu_index in range(1, len(scan_dictionary[sample_data])):
            mu = scan_dictionary[sample_data][mu_index]
            # Label
            if mu_index == len(scan_dictionary[sample_data])-1:
                label = f'{sample_name}_merged'
            else:
                if f'{sample_name}_00{mu_index - 1 + skip_times}' in SKIP_SCANS:
                    skip_times += 1
                label = f'{sample_name}_00{mu_index - 1 + skip_times}'

            if mu_index > 1:
                print(label)
                plt.plot(energy, mu, label=label)

        # Plotting format
        ax.set_xlim(energy.min() // 1 + 1, energy.max() // 1 - 1)
        plt.xticks(fontsize=14)
        plt.title(sample_name, fontsize=20)
        x_label = r'$\mathregular{Energy\ (eV)}$'
        y_label = r'$\mathregular{x\mu(E)}$'
        plt.yticks(fontsize=14)
        ax.set_xlabel(x_label, fontsize=18)
        ax.set_ylabel(y_label, fontsize=18)
        plt.legend(loc='lower right', framealpha=1, frameon=False)
        plt.tight_layout()
        if IF_SAVE:
            plt.savefig("{}/{}.png".format(Path(INPUT_PATH), f'{sample_name}'), dpi=300, transparent=False)
            print('\n=================================================================================')
            print(f'Save figures into ---> {sample_name}.png')
            print('=================================================================================')
            print('')
        plt.close('all')


def create_transmission_prj_group(new_merge_project):
    files = OUTPUT_PATH.glob(f'*.prj')  # Import an input path to read prj files which have been processed
    for index, file_prj in enumerate(files):
        if '.prj' in file_prj.name and 'Created' not in file_prj.name:
            group = read_ascii(f'{file_prj}')  # <------------------------- take care, not read_athena
            group.filename = group.filename[:-4]  # <------------------------- take care, rename group filename

            if SHOW_DATA_INFORMATION:
                print("\n==============================")
                print(f'Scan attributes in {group.filename}')
                print("------------------------------")
                show_data_information(group)

            # Replace special characters because they might cause an error
            sample_name = f'{group.filename}'.replace('-', '_').replace('(', '').replace(')', '').replace('.', 'p')
            print(f'Sample name: {sample_name}')
            # new_merge_project.add_group(group, sample_name)
            try:
                new_merge_project.add_group(group, sample_name)
            except Exception as e:
                # Handle the exception and print an error message with color, Style.RESET_ALL to reset the color
                print("\n==============================")
                print(f"{Back.BLUE}Failed to add group {group} \n with sample {sample_name}: {e}")
                print("------------------------------")
                # Continue with the rest of your code

    new_merge_project.save(f'{OUTPUT_PATH}/Created_transmission_group.prj')
    print('=================================================================================')
    print(f'Save merge project into ---> Created_transmission_group.prj')
    print('=================================================================================')

    # Call Path again to grab the created athena project
    files = OUTPUT_PATH.glob(f'*.prj')
    calibrate_energy(files)


def calibrate_energy(files):
    """
    :param files: files, prj, or no-type files from the current folder
    :return: None
    """
    data_list = []

    # Add e0, atsym, edge, xrayedge information into the created group
    for file in files:
        if 'Created' in file.name and 'calibration' not in file.name:
            print("\n==============================")
            print(f'Calculate energy shift for {file.name}')
            print("------------------------------")
            group = read_athena(file)
            for name, data in group._athena_groups.items():
                e0 = find_e0(data.energy, mu=data.mu, group=data)
                atsym, edge = guess_edge(data.e0)
                data.atsym = atsym
                data.edge = edge
                data.xrayedge = xray_edge(data.atsym, data.edge)[0]
                print(name, data)
                print('e0:', e0)
                print(f'{data.atsym} energy edge:', data.xrayedge)
                if 'foil' in data.label or 'reference' in data.label:    # <------------------------- reference keyword
                    data.energyshift = data.xrayedge - data.e0
                    print('Energy shift:', data.energyshift)
                data_list.append(data)
                print('')

    first_scan_information = data_list[0]

    if SHOW_DATA_INFORMATION:
        print("==============================")
        print(f'Athena parameters in {first_scan_information.label}')
        print("------------------------------")

        show_data_information(first_scan_information)

        print('')

    reference_energy_shift_dictionary = defaultdict(list)

    # Calculate calibrated energy
    print("\n==============================")
    print(f'Add "energy shifts" attributes to the reference group')
    print("------------------------------")
    for index, data in enumerate(data_list):
        sample_name = data.label
        if 'foil' in sample_name:
            reference_energy_shift_dictionary[sample_name] = data.energyshift
            print(index, sample_name)
        elif 'reference' in sample_name:
            reference_energy_shift_dictionary[sample_name] = data.energyshift
            print(index, sample_name)

    new_merge_project = athena_project.create_athena(f'default.prj')  # Call athena_project.py in current folder

    print("\n==============================")
    print(f'Energy shift of the reference')
    print("------------------------------")
    pprint(reference_energy_shift_dictionary)

    print("\n==============================")
    print(f'Energy calibration')
    print("------------------------------")
    reference_checklist = []
    for index, data in enumerate(data_list):
        for reference_name in reference_energy_shift_dictionary:
            if data.label[:-9] in reference_name and data.label[:-6] in reference_name:
                print(index, data.label)
                energy_shift = reference_energy_shift_dictionary[reference_name]
                # print('Energy before:', data.energy[0])
                print('Energy E0 before:', data.e0)
                print('Energy shift:', energy_shift)
                data.energy = data.energy + energy_shift    # <------------------------------------- Energy calibration
                # print('Energy after:', data.energy[0])
                print('Energy E0 after:', find_e0(data.energy, mu=data.mu, group=data))
                reference_checklist.append(reference_name)
                if np.abs(energy_shift) > 3:
                    print('WARNING: Energy shift is too large, please check the reference. '
                                          'Otherwise, please use Created_group.prj to do the calibration manually.')
                print('')

        # Replace special characters because they might cause error
        filename = f'{data.label}'.replace('-', '_').replace('(', '').replace(')', '').replace('.', 'p')
        new_merge_project.add_group(data, filename)

    new_merge_project.save(f'{OUTPUT_PATH}/Created_group_with_calibration.prj')
    print('=================================================================================')
    print(f'Save merge project into ---> Created_group_with_calibration.prj')
    print('=================================================================================')


def show_data_information(group):
    """
    :param group: group, a group contains x-ray information, such as energy, mu, norm, etc...
    :return: None
    """
    for scan_attribute in dir(group):
        print(scan_attribute, type(getattr(group, scan_attribute)))


def check_filename_repetition(output_filename, directory):
    """
    :param output_filename: string, output filename
    :return: string, new output filename
    """
    print("\n==============================")
    print('Check filename repetition')
    print("------------------------------")
    files = Path(directory).glob(f'*.png')
    png_list = []
    for index, file in enumerate(files):
        png_list.append(file.name[:-4])

    print(output_filename)
    while output_filename in png_list:
        output_filename = output_filename + '_1'
        print(output_filename)
    return output_filename

def delete_multiple_file_types(directory, extensions=["*.prj", "*.png"]):
    """
    Delete multiple file types from a directory.

    :param directory: str or Path, the directory containing generated files
    :param extensions: list of str, file extensions to filter files for deletion
    """
    output_path = Path(directory)
    for ext in extensions:
        files_to_delete = output_path.glob(ext)
        for file in files_to_delete:
            try:
                file.unlink()
                print(f"Deleted file: {file}")
            except Exception as e:
                print(f"Failed to delete {file}: {e}")




if __name__ == '__main__':
    main()