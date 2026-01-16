from importlib import import_module
import subprocess
from pathlib import Path
import sys
import numpy as np


def output_terminal(cmd, print_output=True, skip_error=False, print_cmd=False,
                    **kwargs):
    """
    Runs a command in a terminal and save the output in a list
    of strings

    Parameters
    ==========
    cmd: str
        bash command to be executed in the terminal.
    print_output: bool. Default=True
        True for printing the output besides of returning it. Default False.
    skip_error: bool. Default=False
        True for continuing running although the command fails. Default False.
    **kwargs:
        additional options for subprocess.Popen
    print_cmd: Default=False
        print the command that was just executed.

    Return
    ======
    (list) [#linesStr] output of the executed command, line by line.
    """
    if print_cmd:
        print(cmd)
    p = subprocess.Popen(cmd,
                         shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         text=True,
                         **kwargs)

    out = ""
    output = ""
    while not (output == '' and p.poll() is not None):
        output = p.stdout.readline()
        if output:
            out += output
            if print_output:
                print(output, end="")
    return_code = p.wait()

    if not skip_error:
        assert not return_code, f"ERROR executing the command \"{cmd}\"  " + \
            "with output_terminal with the next message:\n" + \
            p.stderr.read().strip()

    return out


pymodules = {
    'F_max_stretch': 'sith.utils.tools',
    'distance': 'sith.utils.tools',
    'all_xyz2pdb': 'sith.utils.tools',
    'conf2pdb': 'sith.utils.tools',
    'diff_bonds': 'sith.utils.tools',
    'extract_bonds': 'sith.utils.tools',
    'change_distance': 'sith.utils.tools',
    'shake_except': 'sith.utils.tools',
    'proline_state': 'sith.g09_stretching.sith_tools',
    'gen_randpep': 'sith.g09_stretching.sith_tools',
    'protonate': 'sith.g09_stretching.protonate',
    'log2xyz2': 'sith.g09_stretching.g09_xyz',
    'log2xyz': 'sith.g09_stretching.g09_xyz',
    'reduce_structs': 'sith.g09_stretching.from_extreme.info_from_opt',
    'continuous_e2e': 'sith.g09_stretching.from_extreme.info_from_opt',
    'info_from_opt': 'sith.g09_stretching.from_extreme.info_from_opt',
    'change_prolines_dofs': 'sith.g09_stretching.change_dofs',
}

sh_executers = {
    'find_blocks': './utils/find_blocks.sh',
    'clean_ds': './utils/clean_ds.sh',
    'basics': './utils/basics.sh',
    'workflow': './g09_stretching/workflow.sh',
    'swap_atoms_in_com': './g09_stretching/swap_atoms_in_com.sh',
    'stretching': './g09_stretching/stretching.sh',
    'simplify_path': './g09_stretching/simplify_path.sh',
    'restart_g09': './g09_stretching/restart_g09.sh',
    'proline_mod': './g09_stretching/proline_mod.sh',
    'workflow_from_extreme': './g09_stretching/from_extreme/workflow_from_extreme.sh',
    'opt_and_forces': './g09_stretching/from_extreme/opt_and_forces.sh',
    'extr_dofs': './g09_stretching/from_extreme/extr_dofs.sh',
    'continuous_path': './g09_stretching/from_extreme/continuous_path.sh',
    'find_forces': './g09_stretching/find_forces.sh',
    'extract_forces': './g09_stretching/extract_forces.sh',
    'create_peptide': './g09_stretching/create_peptide.sh',
}

other_files = {
}


def _read_arguments():
    """
    Function that reads args and kwargs.

    Return
    ======
    (tuple) args
    """
    if len(sys.argv) == 1:
        return ((), {})
    argument = '_reader_args'
    values = ''
    args_dict = {}
    for entry in np.array(sys.argv[1:]):
        if '--' == entry[:2]:
            if argument == '_reader_args':
                args_dict[argument] = values.split(" ")[2:]
                argument = entry.replace('--', '')
                values = ''
            else:
                args_dict[argument] = eval(values)
                argument = entry.replace('--', '')
                values = ''
        else:
            values += ' ' + entry
    if argument == '_reader_args':
        args_dict[argument] = values.split(" ")[2:]
    else:
        args_dict[argument] = eval(values)
    args = args_dict['_reader_args']
    del args_dict['_reader_args']
    
    if '--' == sys.argv[1][:2]:
        args = ()

    return (args, args_dict)


def main():
    """
    This function run each time sith is called from the terminal.

    Return
    ======
    (None)
    """
    # Help menu of this code
    if sys.argv[1] == '-h' or sys.argv[1] == '--help' or sys.argv[1] == 'help':
        functions = list(pymodules.keys()) + list(sh_executers.keys())
        functions.append('tests')
        functions.sort()

        print("\n"
              "'sith' contains the next set of tools:\n\n")
        for function in functions:
            print("    -   " + function)

        print("\nTo use one of the functions from the terminal, use:\n\n"
              "  'sith <function> <arg1> <arg2>...'\n\n"
              "For detailed information of any function, use \"-h\" as first"
              " argument (<arg1>).")

    elif sys.argv[1] == 'tests':
        testdir = Path(__file__).parent
        cmd = f"cd {str(testdir)}/../tests ; pytest -v --color=yes" + \
            ' '.join(sys.argv[2:])
        output_terminal(cmd)

    # python module from terminal
    elif sys.argv[1] in pymodules.keys():
        module = import_module(pymodules[sys.argv[1]])
        method = getattr(module, sys.argv[1])

        if '-h' in sys.argv:
            print(method.__doc__)
        elif '-path' in sys.argv:
            print(pymodules[sys.argv[1]])
        else:
            arguments = _read_arguments()
            output = method(*arguments[0], **arguments[1])
            if output is not None:
                print(output)

    # bash codes
    elif sys.argv[1] in sh_executers.keys():
        if '-path' in sys.argv[2:]:
            path = (str(Path(__file__).parent)[:-3]
                    + sh_executers[sys.argv[1]][2:])
            print(path)
        else:
            command = str(Path(__file__).parent)[:-3] + \
                sh_executers[sys.argv[1]][2:] + ' ' + \
                ' '.join(sys.argv[2:])

            output_terminal(command, print_output=True)

    # other files
    elif sys.argv[1] in other_files.keys():
        print(str(Path(__file__).parent)[:-3] + other_files[sys.argv[1]][2:])

    # own path
    elif sys.argv[1] == 'path':
        print(str(Path(__file__).parent)[:-3])
    
    # open documentation
    elif sys.argv[1] == 'doc':
        command = "xdg-open " + str(Path(__file__).parent)[:-3] + \
            "../../docs/html/index.html"
        output_terminal(command)

    # Not recognized keyword
    else:
        print(f"ERROR: keyword {sys.argv[1]} not recognized as part of"
              " sith. Use 'sith -h' to see the options you can use.")


if __name__ == "__main__":
    main()
