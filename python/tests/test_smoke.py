import pytest

import miradac


def test_version():
    assert miradac.__version__ == "0.1.0"


@pytest.mark.skipif(not miradac.HAS_SYMBOLIC, reason="numeric-only build")
def test_has_symbolic():
    assert miradac.HAS_SYMBOLIC is True
