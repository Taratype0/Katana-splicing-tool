import duckdb
import numpy as np
import pandas as pd

from katana_next.science.motiftrace_exact import full_family_adjust


def test_disk_grouped_corrections_match_full_vectors():
    # Ties need LAST ranks for BH, FIRST ranks for Holm.
    p = np.r_[np.repeat([1e-20, .0002, .01, .2, 1.], [5, 37, 111, 32, 400]),
              np.random.default_rng(71).uniform(size=1000)]
    frame = pd.DataFrame({'p': p})
    m = len(p)
    con = duckdb.connect()
    con.execute('SET threads=1')
    con.register('frame', frame)
    con.execute('CREATE TABLE probabilities AS SELECT p, count(*) AS n FROM frame GROUP BY p')
    con.execute('CREATE TABLE ranked AS SELECT *, sum(n) OVER (ORDER BY p ROWS UNBOUNDED PRECEDING) AS r FROM probabilities')
    harmonic = np.sum(1. / np.arange(1, m + 1))
    actual = con.execute(f'''SELECT p,
        least(1., min(p*{m}/r) OVER (ORDER BY p DESC ROWS UNBOUNDED PRECEDING)) AS bh,
        least(1., {harmonic}*min(p*{m}/r) OVER (ORDER BY p DESC ROWS UNBOUNDED PRECEDING)) AS by_q,
        least(1., max(p*({m}-r+n)) OVER (ORDER BY p ROWS UNBOUNDED PRECEDING)) AS holm
        FROM ranked''').df().set_index('p')
    for column, expected in zip(['bh', 'by_q', 'holm'], full_family_adjust(p)):
        np.testing.assert_allclose(actual.loc[p, column].to_numpy(), expected, rtol=2e-14, atol=1e-15)
    con.close()
