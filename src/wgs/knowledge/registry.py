"""All knowledge sources, in refresh order."""

from .clingen import ClinGen
from .clinvar import ClinVar
from .ensembl import Ensembl
from .population import GnomAD, ThousandGenomes

SOURCES = [Ensembl(), ClinVar(), ClinGen(), ThousandGenomes(), GnomAD()]
