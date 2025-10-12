import subprocess


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
            if print_output and len(output.strip()) != 0:
                print(output.strip())
    return_code = p.wait()

    if not skip_error:
        assert not return_code, f"ERROR executing the command \"{cmd}\"  " + \
            "with output_terminal with the next message:\n" + \
            p.stderr.read().strip()

    return out
