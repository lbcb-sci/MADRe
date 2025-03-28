from setuptools import setup, find_packages

setup(
    name="madre",
    version="1.0.0",
    py_modules=["MADRe"],
    packages=["src"],
    package_dir={"src": "src"},
    entry_points={
        "console_scripts": [
            "madre = MADRe:main",
            "read-classification = src.ReadClassification:main",
            "database-reduction = src.DatabaseReduction:main",
            "calculate-abundances = src.CalculateAbundances:main"
        ]
    },
    install_requires=[
        "scikit-learn"
    ],
    author="Josipa Lipovac",
    author_email="josipa.lipovac@fer.unizg.hr",
    description="Strain-level metagenomic classification with Metagenome Assembly driven Database Reduction approach.",
    url="https://github.com/lbcb-sci/MADRe",
)
