from setuptools import setup, find_packages


def get_version():
    return '0.0.1'

with open("README.md") as handle:
    readme = handle.read()

setup(
    name='varanno',
    version=get_version(),
    author='Alyss Flynn',
    author_email='alyss.flynn@me.com',
    description='Varanno is a prototype variant annotation tool written in python. This was developed for a take home interview project with Tempus AI.',
    long_description=readme,
    license='LICENSE',
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    package_data={
        "varanno": ["*.txt"],
        "varanno.data": ["*"],
    },
    python_requires=">=3.10",
    install_requires=[
        "biopython==1.84",
        "numpy==2.1.2",
        "pandas==2.2.3",
    ]
)
