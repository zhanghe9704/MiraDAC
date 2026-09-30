import pytest

import miradac as da


def test_init_clear_cycles():
    da.init(4, 3, 100)
    assert (da.max_order(), da.current_order(), da.nvars()) == (4, 4, 3)
    assert da.full_length() == 35
    assert da.poolsize() == 100
    assert da.count() + da.remain() == 100
    x = da.var(0) * da.var(1)
    assert x.order == 4
    da.clear()

    da.init(2, 2, 50)
    assert (da.max_order(), da.nvars(), da.full_length(), da.poolsize()) == (2, 2, 6, 50)
    y = da.var(1) * da.var(1)
    assert y.element([0, 2]) == 1.0
    da.clear()

    with pytest.raises(da.EnvError):
        da.count()


def test_count_tracks_vectors(env):
    n = da.count()
    v = da.NDA(1.0)
    assert da.count() == n + 1
    del v
    assert da.count() == n


def test_eps(env):
    old = da.get_eps()
    try:
        da.set_eps(1e-10)
        assert da.get_eps() == 1e-10
        with pytest.raises(ValueError):
            da.set_eps(-1.0)
        assert da.get_eps() == 1e-10
    finally:
        da.set_eps(old)


def test_change_and_restore_order(env):
    assert da.change_order(2) is True
    assert (da.current_order(), da.max_order(), da.full_length()) == (2, 4, 10)
    assert da.change_order(5) is False
    assert da.current_order() == 2
    da.restore_order()
    assert da.current_order() == 4


def test_order_context_nests(env):
    with da.order(3):
        assert da.current_order() == 3
        with da.order(2):
            assert da.current_order() == 2
            assert (da.var(0) * da.var(1) * da.var(2)).iszero()
        assert da.current_order() == 3
        assert (da.var(0) * da.var(1) * da.var(2)).element([1, 1, 1]) == 1.0
    assert da.current_order() == 4


def test_order_context_same_order_nested(env):
    with da.order(3):
        with da.order(3):
            pass
        assert da.current_order() == 3
    assert da.current_order() == 4


def test_order_context_restores_on_error(env):
    with pytest.raises(KeyError):
        with da.order(3):
            with da.order(1):
                raise KeyError
    assert da.current_order() == 4
    with pytest.raises(KeyError):
        with da.order(2):
            raise KeyError
    assert da.current_order() == 4


def test_order_context_rejects_higher_order(env):
    with pytest.raises(ValueError):
        with da.order(5):
            pass
    assert da.current_order() == 4
