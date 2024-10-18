import pytest 
from varanno.parse import parse_record_info, parse_format_sample


@pytest.fixture
def sample_record_info():
    return "BRF=0.16;FR=1.0000;HP=1;HapScore=1;MGOF=3;MMLQ=33;MQ=59.75;NF=89;NR=67;PP=2965;QD=20;SC=CACTTTCCTCATCCACTTTGA;SbPval=0.58;Source=Platypus;TC=160;TCF=90;TCR=70;TR=156;WE=1158639;WS=1158621"


def test_parse_record_info(sample_record_info):
    assert parse_record_info(sample_record_info) == {
        'BRF': '0.16', 'FR': '1.0000', 'HP': '1', 'HapScore': '1', 
        'MGOF': '3', 'MMLQ': '33', 'MQ': '59.75', 'NF': '89', 'NR': '67', 
        'PP': '2965', 'QD': '20', 'SC': 'CACTTTCCTCATCCACTTTGA', 'SbPval': '0.58', 
        'Source': 'Platypus', 'TC': '160', 'TCF': '90', 'TCR': '70', 'TR': '156', 
        'WE': '1158639', 'WS': '1158621'
    }


def test_parse_format_sample():
    formatstr = "GT:GL:GOF:GQ:NR:NV"
    samplestr = "1/1:-300.0,-43.88,0.0:3:99:160:156"
    assert parse_format_sample(formatstr, samplestr) == {
        'GT': '1/1', 'GL': '-300.0,-43.88,0.0', 'GOF': '3', 'GQ': '99', 'NR': '160', 'NV': '156'
    }
