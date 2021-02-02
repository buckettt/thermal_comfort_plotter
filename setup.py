# -*- coding: utf-8 -*-

# A very simple setup script to create a single executable
#
# hello.py is a very simple 'Hello, world' type script which also displays the
# environment in which the script runs
#
# Run the build process by running the command 'python setup.py build'
#
# If everything works well you should find a subdirectory in the build
# subdirectory that contains the files needed to run the script without Python

import sys, os
from cx_Freeze import setup, Executable

base = None
#if sys.platform == 'win32': # Commentedout to include console.
#    base = 'Win32GUI'

executables = [
    Executable('thermalcomfort_plotter.py',
	base=base,
	icon=r'Y:\git_projects\mf_toolbox\dev\mf_modules\res\MF_O_trans.ico')
]

additional_mods = ["numpy.core._methods", "numpy.lib.format"]
exclude_mods = ["babel", "mf_modules.mf_pyIES", "scipy", "PyQt5", "tornado", "zmq", "sphinx", "sphinx_rtd_theme", "psutil", "notebook", "nbconvert", "lxml", "cryptography", "bottleneck", "pandas", "markupsafe", "prompt_toolkit", "jedi", "chardet"]

build_exe_options = {"excludes": exclude_mods, "includes": additional_mods}

os.environ['TCL_LIBRARY'] = r'C:\ProgramData\Anaconda3\tcl\tcl8.6'
os.environ['TK_LIBRARY'] = r'C:\ProgramData\Anaconda3\tcl\tk8.6'

setup(name='Thermal Comfort Plotter',
      version='0.1',
	  includes = ['os'],
	  options = {"build_exe": build_exe_options},
      description='Plot different comfort metrics based on various parameters.',
      executables=executables
      )