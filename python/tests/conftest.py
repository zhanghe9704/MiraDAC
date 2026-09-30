import pytest

import miradac as da


@pytest.fixture
def env():
    """Default env at order 4 with 3 variables; cleared afterwards."""
    da.init(4, 3, 400)
    yield
    da.clear()
