import os
from setuptools import setup

version_py = os.path.join(os.path.dirname(__file__), "artic", "version.py")
version = open(version_py).read().strip().split("=")[-1].replace('"', "").strip()
long_description = """
``artic`` is a pipeline for working with virus sequencing data sequenced with nanopore
"""

HERE = os.path.dirname(__file__)

setup(
    name="artic",
    version=version,
    requires=["python (>=3.5)"],
    packages=["artic"],
    author="Nick Loman",
    description="A toolset for working with nanopore sequencing data",
    long_description=long_description,
    url="https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html",
    package_dir={"artic": "artic"},
    package_data={"artic": []},
    zip_safe=False,
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "artic=artic.pipeline:main",
            "align_trim=artic.align_trim:main",
            "align_trim_n=artic.align_trim_n:main",
            "margin_cons=artic.margin_cons:main",
            "margin_cons_medaka=artic.margin_cons_medaka:main",
            "vcfextract=artic.vcfextract:main",
            "artic_vcf_merge=artic.vcf_merge:main",
            "artic_vcf_filter=artic.vcf_filter:main",
            "artic_make_depth_mask=artic.make_depth_mask:main",
            "artic_fasta_header=artic.fasta_header:main",
            "artic_mask=artic.mask:main",
            "artic_get_stats=artic.artic_mqc:main",
        ],
    },
    author_email="n.j.loman@bham.ac.uk",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
