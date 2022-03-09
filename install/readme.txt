Instructions for set-up the python environment to work on Windows
1) Install miniconda.
1) Open a cmd prompt and navigate to the install folder on the github repo.
2) Open the yml file  with a text editor and modify last line (prefix) with the path of your miniconda folder.
2) Execute the command: conda env create -f kalpana.yml
5) Install the 4 wheel files available on the install folder. Follow the next order: gdal, fiona, shapely and rtree. For installing, execute: pip install name_of_the_wheel_file (use tab for auto-complete).
6) Install geopandas using pip: pip install geopandas