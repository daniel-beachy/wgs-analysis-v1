from wgs.hgvs import protein


def test_missense_and_nonsense():
    assert protein("314N>314D", "missense") == "p.Asn314Asp"
    assert protein("73W>73*", "stop_gained") == "p.Trp73Ter"
    assert protein("741L", "synonymous") == "p.Leu741="
    assert protein("1M>1L", "start_lost") == "p.Met1?"
    assert protein("204*>204R", "stop_lost") == "p.Ter204Argext*?"


def test_frameshift_names_first_changed_residue():
    assert protein("302SRNG*>302SQEREE*", "frameshift") == "p.Arg303GlnfsTer6"
    assert protein("635PPPRT*>635PPSH*", "frameshift") == "p.Pro637SerfsTer3"


def test_inframe():
    assert protein("128GKK>128G", "inframe_deletion") == "p.Lys129_Lys130del"
    assert protein("92SSGS>92S", "inframe_deletion") == "p.Ser93_Ser95del"
    assert protein("1M>1MT", "inframe_insertion") == "p.Met1delinsMetThr"


def test_unparseable_is_none():
    assert protein(None) is None
    assert protein("weird") is None
