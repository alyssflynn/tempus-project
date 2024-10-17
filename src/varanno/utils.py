import re


VCF_META_KEYVAL = r"^##(?P<key>\w+)=(?P<value>.+)$"


VCF_META_STRUCT = (
    r'^##(?P<key>\w+)=<'                                      # start
    r'ID=(?P<id>\w+)'                                         # required (string)
    r'(,Number=(?P<number>[.0-9A-Z]+))?'                      # required (int, ., or char)
    r'(,Type=(?P<type>Integer|Float|Flag|Character|String))?' # required (string enum)
    r',Description="(?P<description>.+)"'                     # required (quoted string)
    r'(,Source="(?P<source>\w+)")?'                           # optional (quoted string)
    r'(,Version="(?P<version>\d+)")?'                         # optional (quoted int) NOTE: always int?
    r'>$'                                                     # end
)


def cast_float(value: str):
    """Attempts to cast a value as a float, otherwise returns it unchanged."""
    try:
        return float(value)
    except ValueError:
        return value


def parse_record_info(info: str):
    """Parses a VCF record's INFO value into a dictionary.
    
    E.g.
    ```
    BRF=0.16;FR=1.0000;HP=1;HapScore=1;MGOF=3;MMLQ=33;MQ=59.75;NF=89;NR=67;PP=2965;QD=20;SC=CACTTTCCTCATCCACTTTGA;SbPval=0.58;Source=Platypus;TC=160;TCF=90;TCR=70;TR=156;WE=1158639;WS=1158621
    ```
    becomes
    ```
    {
        'BRF': '0.16', 'FR': '1.0000', 'HP': '1', 'HapScore': '1', 'MGOF': '3', 
        'MMLQ': '33', 'MQ': '59.75', 'NF': '89', 'NR': '67', 'PP': '2965', 'QD': '20', 
        'SC': 'CACTTTCCTCATCCACTTTGA', 'SbPval': '0.58', 'Source': 'Platypus', 'TC': '160', 
        'TCF': '90', 'TCR': '70', 'TR': '156', 'WE': '1158639', 'WS': '1158621'
    }    
    ```
    """
    return dict(re.findall(r"\b(\w+)=([^;]+);?", info))


def parse_format_sample(formatstr:str, samplestr:str):
    """Parses a VCF record's FORMAT and SAMPLE values into a dictionary. 

    E.g.
    ```
    FORMAT              sample
    GT:GL:GOF:GQ:NR:NV  1/1:-300.0,-57.5,0.0:3:99:194:193
    ```
    becomes
    ```
    {'GT': '1/1', 'GL': '-300.0,-43.88,0.0', 'GOF': '3', 'GQ': '99', 'NR': '160', 'NV': '156'}
    ```
    """
    return dict(zip(
        formatstr.split(":"),
        samplestr.split(":")
    ))
