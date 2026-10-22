import numpy as np
import pandas as pd

from scripts.benchmark_original_marvel import audit_se


def test_marvel_gates_and_length_factor(tmp_path):
    def write(name, values):
        pd.DataFrame({"tran_id": ["e1", "e2", "e3", "e4"], "cell": values}).to_csv(tmp_path / name, sep="\t", index=False)
    write("SE_sj.included.1.tsv", [20, 1, 100, 0])
    write("SE_sj.included.2.tsv", [20, 1, 1, 0])
    write("SE_sj.excluded.tsv", [20, 1, 20, 20])
    write("psi_SE.tsv", [.5, np.nan, np.nan, 0])
    result = audit_se(tmp_path)
    assert result["measurable_event_cells"] == 2
    assert result["missingness_identical"]
