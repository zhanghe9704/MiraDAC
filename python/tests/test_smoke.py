import miradac


def test_version():
    assert miradac.__version__ == "0.1.0"


def test_has_symbolic():
    assert miradac.HAS_SYMBOLIC is True
