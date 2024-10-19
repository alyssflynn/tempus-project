# Tempus-Project (varanno)

Varanno is a prototype variant annotation tool written in python. This was developed for a take home interview project with Tempus AI. It uses a streaming pattern to ingest and process the data in order to minimize memory usage.

A CSV containing variant annotations for the provided VCF file is located [here](/annotation_results/annotations.csv).

[See project notes below!](#notes)

## Getting started

**This package requires python version 3.10.0 or greater!** Please verify your system meets this requirement with: `python3 --version`

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

### Running tests

This project uses pytest. Use this command to run all varanno tests:
```
$ pytest
```

### Other commands

run linter: `$ ruff format src/ --diff`

run mypy: `$ mypy src/`

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
$ python -m varanno -f "path/to/vcf.txt" -o "output"
```
...then after waiting for it to complete...
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

## Notes

### Dependency avoidance
There are a lot of existing tools for reading VCF files out there, but I avoided using these because that seemed more in the spirit of the project. In reality I suspect rolling your own file reader isn't the best approach.

### VEP API - GRCh37 vs GRCh38
The prompt specifies using the VEP API available at this endpoint: https://rest.ensembl.org/#VEP. This API was yielding pretty unreliable results (or reliably poor results), as more than 95% of requests returned 400s. 

I did some digging and realized that this may be due to the genome assembly version--the VEP API currently uses GRCh38, but the provided VCF file appears to have been based on the GRCh37 assembly (based on the metadata `'refFile': '/data/ref_genome/human_g1k_v37.fasta'`). I then switched to using a different VEP API version ("https://grch37.rest.ensembl.org"), which dramatically improved request success rate. 

The package currently hard-codes the URL (and therefore the API version used in vep requests), but In a future version I'd probably want to make the package more extensible so you could toggle between assembly versions, and/or find a way to interpret assembly version from the file metadata. 

### VEP query speed
The current process runs fairly slowly. I was originally using the single-HGVS endpoint and switched to the bulk endpoint to save time (and coupled this to making the reader behave like a generator). 

Despite this, I'm still unsatisfied with the runtime. I opted to leave it as is for the sake of time, but for a longer-term projects I'd want to try a few different approaches to speed things up. One idea would be to switch to using `asynio/aiohttp`, as I've had a lot of success with it in the past for long chains of API calls like this, or alternatively take an offline approach and download the data needed to run this locally.    

## Batch size
The current batch size is 50 by default, but this can be modified in code. The VEP API specs say that this value can be increased up to 300, although higher batch size may not necessarily mean faster runtime as VEP API request latency seems to scale noticeably with batch size. An interesting next step would be to benchmark various batch sizes and optimize for the best runtime.  

## VEP hgvs API bulk endpoint issue
In my first pass, I assumed the VEP hgvs bulk endpoint (which accepts a list of `hgvs_notations` strings and returns a JSON array of objects) was returning an output list equal in length to the input list. It turns out the results for query strings it can't process are filtered out of the output entirely, so it is not guaranteed that the output list is equal in length to the input list provided. This resulted in dropped rows and annotation outputs that were misaligned with the input records. 

To fix this I added a condition to realign inputs with outputs if the lists were of different lengths. For a longer term project I'd want to add more robust error handling and more useful annotations for input rows which failed certain steps.

## Multi-allelic variant rows
The input VCF file contained some rows which contained commas in the ALT column, like this:
```
1	91859795	.	TATGTGA	CATGTGA,CATGTGG	2962	PASS	BRF=0.23;FR=0.5000,0.5000;HP=3;HapScore=2;MGOF=1;MMLQ=37;MQ=59.76;NF=115,51;NR=93,44;PP=2962,2892;QD=20;SC=TTGCCAGCAATATGTGATAAG;SbPval=0.54;Source=Platypus;TC=209;TCF=115;TCR=94;TR=208,95;WE=91859809;WS=91859785	GT:GL:GOF:GQ:NR:NV	1/2:-1,-1,-1:1:99:209,209:208,95
```
The VCF processor currently still attempts to generate an HGVS notation string for these rows mostly as usual, which doesn't appear to be a valid format. The VEP API therefore generally returns no data for these rows (see "VEP hgvs API bulk endpoint issue" above).

I wasn't sure how to best handle these records, so rather than filtering them out I decided to handle processing failures gracefully and fill unknown output annotations with empty or "unknown" values. This way the input record would still be represented in the output file with whatever annotations we were able to generate successfully (e.g. the coverage value).

In a real-life project situation, I'd ask for clarification on how to handle rows like these during the requirements-gathering process of the project, e.g.:
- How should these rows be handled?
- Should these rows be split into separate variants and annotated independently?
- [What is the/Is there a] correct way to format the HGVS notation string to query the VEP API for notation data on these multi-allelic variants?
