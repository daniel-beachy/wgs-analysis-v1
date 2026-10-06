"""All knowledge sources, in refresh order."""

from .clingen import ClinGen
from .clinvar import ClinVar
from .curated import HPO, AcmgSF, Mondo, Orphanet
from .ensembl import Ensembl
from .population import GnomAD, ThousandGenomes
from .predictors import AlphaMissense, GeneConstraint, Revel

SOURCES = [Ensembl(), ClinVar(), ClinGen(), ThousandGenomes(), GnomAD(), HPO(), Mondo(), Orphanet(), AcmgSF(),
           GeneConstraint(), AlphaMissense(), Revel()]
