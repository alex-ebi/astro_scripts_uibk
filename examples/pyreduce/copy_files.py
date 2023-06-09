"""
This script copies the pyreduce config files to the pyreduce directory.
Enter the path to your pyreduce directory which contains the instruments, settings, masks and wavecal directories.
e.g. python3 copy_files.py ~/PyReduce/pyreduce
"""

from pathlib import Path
import shutil

pyreduce_dir = Path(input('Enter the path to your pyreduce directory which contains the instruments, settings, '
                          'masks and wavecal directories.'))

copy_dict = {'echelle_uibk.json': 'instruments',
             'echelle_uibk_example.py': 'instruments',
             'settings_ECHELLE_UIBK.json': 'settings',
             'mask_echelle_uibk.fits': 'masks',
             'echelle_uibk_2D.npz': 'wavecal'}

for file, target in copy_dict.items():
    shutil.copy(file, pyreduce_dir / target)
