"""
A util function to clean the gromacs files.
"""

import os
import sys

def set_opls(file_name: str):
    """
    For the given gromacs file change to the opls combination rules.
    """

    with open(file_name) as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if "comb-rule" in line:
                index = i
                break
        # now sub
        lines[i + 1] = "  1           3               yes            0.5         0.5\n"
    # now put this in a new file
    with open(file_name, "w")as out:
        out.writelines(lines)


def main():
    """
    Find all gromacs top files and run set opls
    """
    for top in os.listdir("."):
        if top.endswith("top"):
            set_opls(top)


if __name__ == "__main__":
    main()