from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="pm7_enhanced_charged",
    version="0.2.0",
    author="Debjyoti Bhattacharya",
    author_email="bhattadeb34@gmail.com",  # Update with your email
    description="PM7 semi-empirical quantum chemistry calculator with charge support and proton affinity calculations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bhattadeb34/pm7_enhanced_charged",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    python_requires=">=3.7",
    install_requires=[
        "numpy>=1.19.0",
        "pandas>=1.1.0",
        "rdkit>=2020.09.1",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
        ],
    },
)
