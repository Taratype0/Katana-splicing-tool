import numpy as np
from katana_next.science.heatmap_style import row_scaled, pastel_diverging


def test_scale_and_missing():
    raw=np.array([[.1,.2,.3],[.5,.5,np.nan],[np.nan,np.nan,np.nan]])
    result=row_scaled(raw)
    np.testing.assert_allclose(result[0],[-1,0,1],atol=1e-12)
    np.testing.assert_allclose(result[1,:2],[0,0])
    assert np.isnan(result[1,2]) and np.isnan(result[2]).all()
    assert raw[0,0]==.1


def test_pastel_endpoints():
    cmap=pastel_diverging()
    assert cmap(0.)[2]>cmap(0.)[0]
    assert cmap(1.)[0]>cmap(1.)[2]
    assert min(cmap(.5)[:3])>.99
