import io
import setuptools
try:
    import numpy as np
except ImportError:
    print("Please install Numpy first. e.g. pip install numpy")
    exit(1)

# Every SIMD kernel in c/utils.cpp (AVX-512F on x86, SVE on aarch64) carries a function-level target attribute,
# so no arch flags are required.
module_utils = setuptools.extension.Extension('ssam.utils', sources=["c/utils.cpp"], extra_compile_args=["-fopenmp"], extra_link_args=["-fopenmp"], include_dirs=[np.get_include()])

with io.open("README.rst", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ssam",
    version="1.1.2",
    author="Jeongbin Park",
    author_email="jeongbin.park@pusan.ac.kr",
    description="SSAM",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pnucolab/ssam",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)",
        "Operating System :: POSIX",
    ],
    ext_modules = [module_utils],
    install_requires=[
        "numpy>=1.26.4,<=2.0.2",
        "zarr<3",
        "networkx",
        "numba",
        "scipy",
        "pandas",
        "matplotlib",
        "seaborn",
        "scikit-learn",
        "umap-learn",
        "leidenalg",
        "sparse",
        "scikit-image",
        "dask",
        "dask[array]",
        "opencv-python",
        "packaging",
        "tqdm",
    ]
)
