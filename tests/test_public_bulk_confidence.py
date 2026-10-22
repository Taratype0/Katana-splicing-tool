import pytest
from scripts.audit_public_bulk_confidence import suppa_se_structure


def test_suppa_coordinates_and_strand():
    result = suppa_se_structure("G;SE:chr1:10-20:30-40:-")
    assert result["upstreamEE"] == 10
    assert result["exonStart_0base"] == 19
    assert result["exonEnd"] == 30
    assert result["downstreamES"] == 39
    assert result["strand"] == "-"
    with pytest.raises(ValueError):
        suppa_se_structure("G;SE:chr1:10-20:19-40:+")
