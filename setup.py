from distutils.core import setup

setup(
    name = "jupyterPathway",
    packages = ["jupyterPathway"],
    version = "0.0.4",
    description= "Canonical pathway visualization and analysis",
    url = "https://github.com/ucsd-ccbb/jupyter-pathway",
    author="Brin Rosenthal (sbrosenthal@ucsd.edu), Lilith Huang (lihuang@ucsd.edu)",
    author_email="sbrosenthal@ucsd.edu",
    keywords = ['Jupyter notebook', 'interactive', 'pathway', 'KEGG', 'metabolic'],
    license = 'MIT',
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
    ]
)
