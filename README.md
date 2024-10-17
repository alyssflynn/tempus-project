# Tempus-Project (varanno)

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

The script will create the `output_directory` if it doesn't already exist, and generates 3-4 files when processing the VCF file:
1. `annotations.csv` contains the annotations for each variant in the VCF file.
2. `metadata.json` contains a JSON of the VCF file headers.
3. `tmp.log` contains a log of the actions taken during script execution. 
4. `errors.log` contains a list of errors encountered while running the script (only present if any occured).


### Python package
You can also use the varanno package to process and interact with the VCF data in other formats like pandas dataframes. 

```python3
import pandas as pd
from varanno import Reader, Record, VariantAnnotation

vcf_file = "path/to/vcf.txt"

reader = Reader(vcf_file)
reader.load_records()

metadata_df = pd.DataFrame(reader.meta_structs())
records_df = pd.DataFrame(reader.records)

annotations = list(reader.annotation_generator())
annotations_df = pd.DataFrame(annotations)
```

An important thing to note is that the annotation process involves a lot of calls to the VEP API, so you should avoid running annotations more than once. A good strategy would be to run the script in the console, and load the results into a dataframe (or your format of choice) one the process is complete.

```bash
python -m varanno -f "path/to/vcf.txt" -o "output"
...
```

```python3
import pandas as pd
from varanno import Reader, Record, VariantAnnotation

vcf_file = "path/to/vcf.txt"
annotations_csv = "output/annotations.csv"

reader = Reader(vcf_file)
reader.load_records()

metadata_df = pd.DataFrame(reader.meta_structs())
records_df = pd.DataFrame(reader.records)

annotations_df = pd.read_csv(annotations_csv)
```
