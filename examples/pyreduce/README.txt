You have to copy those files into the pyreduce package to make the reduction of ECHELLE_UIBK work.
For convenience you can use copy_files.py

echelle_uibk.json -> pyreduce/instruments/
echelle_uibk.py -> pyreduce/instruments/
settings_ECHELLE_UIBK.json -> pyreduce/settings/
mask_echelle_uibk.fits -> pyreduce/masks/
echelle_uibk_2D.npz -> pyreduce/wavecal/