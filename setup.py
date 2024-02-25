from setuptools import find_packages, setup

__docs__ = "NovAna (Novelty Analysis) is a cheminformatics tool that " \
           "that extends the decomposition of molecules into " \
           "Bemis-Murcko frameworks (J. Med. Chem. 1996, 39, 15, 2887-2893)."

# read the contents of your README file
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="novana",
    packages=find_packages("src"),
    version="0.2.0",
    description=__docs__,
    long_description=long_description,
    long_description_content_type='text/markdown',
    author="Gian Marco Ghiandoni",
    author_email="ghiandoni.g@gmail.com",
    license="MIT",
    install_requires=["rdkit"],
    setup_requires=["pytest-runner", "flake8"],
    tests_require=["pytest==7.1.2"],
    test_suite="test",
    package_dir={"novana": "src/novana"},
)
