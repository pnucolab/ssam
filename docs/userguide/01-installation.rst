Installation
============

A step-by-step guide
--------------------

The easiest way to prepare a python environment for SSAM is using
`conda <https://docs.conda.io/projects/conda/en/latest/user-guide/install/>`__.
Keeping python projects in isolated environments prevents dependency
version conflicts or conflicts with your OS installation of python which
usually depends on older versions incompatible with current scientific
packages.

Create your environment:

::

   conda create -n ssam python=3.11 numpy=2.0.2 gxx_linux-64

Remember to activate before using it:

::

   conda activate ssam


Finally we switch to pip:

::

   pip install ssam

(Optional) Using ``sctransform`` for normalization
--------------------------------------------------

As of SSAM v1.1.0, SSAM uses log normalization by default. Read this article for hints on choosing an appropriate normalization method for your data: https://www.nature.com/articles/s41592-023-01814-1

If you want to use ``sctransform`` to normalize your vectors, you can additionally install ``R`` and the R packages ``sctransform`` and ``feather``.

First, install R:

::

   conda install -c r r-base

Then, install the required R packages. Open R and type:

::

   install.packages("sctransform")
   install.packages("feather")

And install ``pyarrow`` package for interoperation between R and Python:

::

   pip install pyarrow

Next we can download and prepare our :doc:`data <02-data>`__.

SSAM’s source code
------------------

In case you want to work with `SSAM’s source
code <https://github.com/pnucolab/ssam>`__, it is also hosted on github.
