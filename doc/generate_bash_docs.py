"""
Generate the rst pages that document the bash scripts of sith.

For every ``*.sh`` file under ``src/sith`` two pages are written:

- ``modules/bash_rsts_doc/<name>.rst``: the help text of the script, i.e. the
  text printed by its ``print_help`` function (``sith <name> -h``).
- ``modules/bash_rsts_scripts/<name>.rst``: the full source of the script.

The help text is extracted statically from the source instead of executing
``<script> -h``, because running it lets the shell expand the backticks of
``:bashscript:`...``` roles as command substitutions.

This module is executed automatically by ``conf.py`` on every Sphinx build,
and it can also be run by hand::

    python doc/generate_bash_docs.py
"""
from pathlib import Path
import re
import sys

DOC_DIR = Path(__file__).resolve().parent
SRC_DIR = DOC_DIR.parent / 'src' / 'sith'
MODULES_DIR = DOC_DIR / 'modules'
HELP_DIR = MODULES_DIR / 'bash_rsts_doc'
SCRIPTS_DIR = MODULES_DIR / 'bash_rsts_scripts'

# keeps :bashscript:`...` roles, escapes any other backtick.
ROLE = re.compile(r':bashscript:`[^`]*`')


def extract_help(script: Path) -> list:
    """
    Extract the help text of a bash script from its ``print_help`` function.

    Parameters
    ==========
    script: Path
        path to the bash script.

    Return
    ======
    (list) lines of the help text, without the leading and trailing empty
    lines.
    """
    lines = script.read_text().splitlines()
    try:
        start = next(i for i, line in enumerate(lines)
                     if re.match(r'\s*print_help\s*\(\)\s*\{', line))
    except StopIteration:
        raise ValueError(f"{script} has no print_help function")

    opener = lines[start + 1].strip()
    body = []
    if opener.startswith('cat'):
        # cat << EOF ... EOF
        delimiter = opener.split('<<')[1].strip().strip('\'"-')
        for line in lines[start + 2:]:
            if line.strip() == delimiter:
                break
            body.append(line)
    elif opener.startswith('echo'):
        # echo [-e] " ... "
        first = opener.split('"', 1)[1]
        if first:
            body.append(first)
        for line in lines[start + 2:]:
            if line.rstrip().endswith('"') and not line.rstrip()\
                    .endswith('\\"'):
                last = line.rstrip()[:-1]
                if last:
                    body.append(last)
                break
            body.append(line)
    else:
        raise ValueError(f"Unrecognized print_help format in {script}")

    # undo the shell escapes inside double quotes (\\ \" \$ \`), as echo does.
    body = [re.sub(r'\\([\\"$`])', r'\1', line.rstrip()) for line in body]
    while body and not body[0]:
        body.pop(0)
    while body and not body[-1]:
        body.pop()
    return body


def escape_rst(line: str) -> str:
    """
    Escape the characters of a help line that rst would interpret as inline
    markup, except for :bashscript: roles.

    Parameters
    ==========
    line: str
        line of the help text.

    Return
    ======
    (str) escaped line.
    """
    def escape(text):
        text = text.replace('\\', '\\\\')
        text = text.replace('*', '\\*').replace('`', '\\`')
        text = text.replace('|', '\\|')
        # word_ would be a reference in rst
        text = re.sub(r'(\w)_(?=\W|$)', r'\1\\_', text)
        return text

    out = ''
    last = 0
    for match in ROLE.finditer(line):
        out += escape(line[last:match.start()]) + match.group(0)
        last = match.end()
    return out + escape(line[last:])


def help_page(name: str, relpath: str, help_lines: list) -> str:
    """
    rst page with the help text of a script.
    """
    bar = '=' * len(name)
    body = '\n'.join(('      ' + escape_rst(line)).rstrip() or ''
                     for line in help_lines)
    return (f"{bar}\n{name}\n{bar}\n\n"
            ".. container:: bash-script-title\n\n"
            f"   :ref:`[script] <{name}>` **sith/{relpath}**\n\n"
            ".. container:: bash-script-doc\n\n"
            "   .. line-block::\n\n"
            f"{body}\n")


def script_page(name: str, relpath: str) -> str:
    """
    rst page with the source of a script.
    """
    title = f"Script of sith {name}"
    bar = '=' * len(title)
    return (f".. _{name}:\n\n"
            f"{bar}\n{title}\n{bar}\n\n"
            f".. literalinclude:: ../../../src/sith/{relpath}\n"
            "   :language: bash\n")


def write_if_changed(path: Path, content: str) -> bool:
    """
    Write content only if it differs, so Sphinx does not rebuild unchanged
    pages.
    """
    if path.exists() and path.read_text() == content:
        return False
    path.write_text(content)
    return True


def included_scripts() -> set:
    """
    Names of the scripts included in some modules/sith*.rst page.
    """
    names = set()
    for rst in MODULES_DIR.glob('sith*.rst'):
        names.update(re.findall(r'include:: bash_rsts_doc/(\w+)\.rst',
                                rst.read_text()))
    return names


def main() -> int:
    HELP_DIR.mkdir(parents=True, exist_ok=True)
    SCRIPTS_DIR.mkdir(parents=True, exist_ok=True)

    scripts = sorted(SRC_DIR.rglob('*.sh'))
    names = [s.stem for s in scripts]
    duplicated = {n for n in names if names.count(n) > 1}
    assert not duplicated, f"duplicated script names: {duplicated}"

    changed = []
    for script in scripts:
        relpath = script.relative_to(SRC_DIR).as_posix()
        help_lines = extract_help(script)
        if write_if_changed(HELP_DIR / f'{script.stem}.rst',
                            help_page(script.stem, relpath, help_lines)):
            changed.append(f'bash_rsts_doc/{script.stem}.rst')
        if write_if_changed(SCRIPTS_DIR / f'{script.stem}.rst',
                            script_page(script.stem, relpath)):
            changed.append(f'bash_rsts_scripts/{script.stem}.rst')

    # pages of scripts that do not exist anymore
    for folder in (HELP_DIR, SCRIPTS_DIR):
        for rst in folder.glob('*.rst'):
            if rst.stem not in names:
                print(f"[generate_bash_docs] WARNING: {rst} documents a "
                      "script that does not exist.", file=sys.stderr)

    missing = set(names) - included_scripts()
    for name in sorted(missing):
        print(f"[generate_bash_docs] WARNING: {name}.sh is not included in "
              "any doc/modules/sith*.rst page.", file=sys.stderr)

    for page in changed:
        print(f"[generate_bash_docs] updated {page}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
