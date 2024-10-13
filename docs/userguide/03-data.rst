Data Preparation
================

Download VISp data
------------------

In this tutorial we work with a filtered subset of the spot data of
the murine primary visual cortex (VISp) profiled using multiplexed smFISH.

To access the full dataset, please refer to this Zenodo record: https://zenodo.org/record/3478502

First, download the data:

::

   curl "https://s3.amazonaws.com/starfish.data.spacetx/spacetx-website/data/smFISH_Allen/s3_spot_table.csv" -o s3_spot_table.csv

Load data into python
---------------------

Let’s start with loading our python packages:

::

   import numpy as np
   import pandas as pd
   import matplotlib.pyplot as plt
   import ssam

Now we can load the mRNA spot table. Each row describes one mRNA spot
and the columns contain its coordinates and target gene. We load the
required columns into a dataframe:

::

   df = pd.read_csv("s3_spot_table.csv", usecols=['rotated_x', 'rotated_y', 'gene'])

We are only interested in the coordinates of the spots and the gene names,
so we filter the dataframe accordingly using the ``usecols`` argument.

Now we rename the columns to match the format SSAM expects:

::

   df = df.rename(columns={'rotated_x': 'x', 'rotated_y': 'y'})

SSAM expects the columns to be named ``x``, ``y``, ``z`` and ``gene`` (``z`` is optional).

If your dataset is organized differently, you will have to reshape it
before continuing with the next steps.

Transform Data
--------------

As a next step we normalize the coordinates of the spots. We want the
origin to be at the top left corner of the tissue, so we subtract the
minimum x and y coordinates from all spots:

::

   df.x -= df.x.min()
   df.y -= df.y.min()

SSAM expects the coordinates to be in micrometers, so we need to
multiply the coordinates by the pixel size of the image.

In this dataset, the coordinates are already in micrometers, so we don't
need to do this. If your data is in a different unit, you will have to
convert it into micrometers by multiplying with the conversion factor.

::

   # Below is not needed for this dataset, but to show how it could be done
   conversion_factor = 1.0 # change this if your data is in a different unit
   df.x *= um_per_pixel
   df.y *= um_per_pixel


Create the ``SSAMDataset`` and ``SSAMAnalysis`` objects
-------------------------------------------------------

For the SSAM analysis we need to create a ``SSAMDataset`` and ``SSAMAnalysis`` objects.

The ``SSAMDataset`` object will contain all analysis results and intermediate data,
while the ``SSAMAnalysis`` object will perform the analysis steps.

We create a new ``SSAMDataset`` object by providing a name for the
dataset. This will initialize a new directory (as a ``Zarr`` store)
where all the data will be stored.

::

   # Below will create a new directory named "ssam_msmfish" in the current working directory
   ds = ssam.SSAMDataset("ssam_msmfish")

Next, we create a ``SSAMAnalysis`` object. This object will perform the
analysis steps on the dataset.

::

   analysis = ssam.SSAMAnalysis(ds, verbose=True)

The ``verbose=True`` argument let us see the progress of the analysis.

It is also possible to specify the number of cores to use for the
analysis steps. This is useful for speeding up the analysis on
multi-core machines:

::

   analysis = ssam.SSAMAnalysis(ds, ncores=10, verbose=True) # use 10 cores

Now we can start the analysis with the `kernel density
estimation <kde.md>`__ step.