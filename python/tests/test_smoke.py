import pytest

import miradac


def test_version():
    from importlib.metadata import version
    assert miradac.__version__ == version("miradac")


@pytest.mark.skipif(not miradac.HAS_SYMBOLIC, reason="numeric-only build")
def test_has_symbolic():
    assert miradac.HAS_SYMBOLIC is True
