# Tempus-Project / varanno

Varanno is a prototype variant annotation tool written in python. This was developed for a take home interview project with Tempus AI. 

## Getting started

Clone the project and cd into the project directory.

Create and activate a virtual environment:
```
$ python3 -m venv .venv
$ source .venv/bin/activate
```

Install project requirements:
```
$ python -m pip install --upgrade pip
$ python -m pip install -e .'[dev]'
```

### Run tests

Run tests with 
`$ pytest`

## Usage

### Run from console

You can run varanno from the command line, specifying an input VCF file and output directory for the annotation data and metadata files:

`$ python -m varanno -f "{vcf_input_file}" -o "{output_directory}"`

Or run from an interactive console to inspect the data dynamically. 

