import numpy as np
from katana_next.science import rbp_short_extension as original
from katana_next.science import rbp_short_extension_cached as cached


def test_cached_sequence_and_junction_opportunities_are_exact():
    rng=np.random.default_rng(20261015)
    sequences=["".join(rng.choice(list("ACGTN"),size=i)) for i in range(1,95)]
    for pattern in ("AAA","[AC]G[CT]","ACTTACTT","A","AAAAAAAAAAAAAAAAAAAAAAAAA"):
        for junctions in (None,[(len(s)//3,2*len(s)//3) for s in sequences]):
            a=original.consensus_scores(sequences,pattern,junctions)
            b=cached.consensus_scores(sequences,pattern,junctions)
            np.testing.assert_equal(a[0],b[0])
            assert a[1]==b[1]
            for x,y in zip(a[2],b[2]):
                np.testing.assert_equal(x,y)


def test_cached_rank_test_is_exact_with_ties_missing_and_constant_columns():
    rng=np.random.default_rng(20261015)
    a=rng.integers(0,10,size=(61,30)).astype(float)
    b=rng.integers(0,10,size=(53,30)).astype(float)
    a[:,0]=b[:,0]=0
    a[:5,1]=np.nan
    b[:40,2]=np.nan
    a[:,3]=b[:,3]=np.nan
    for x,y in zip(original.rank_test_matrix(a,b),cached.rank_test_matrix(a,b)):
        np.testing.assert_equal(x,y)
