import subprocess
import sys
import textwrap

import pytest

import miradac as da


def test_env_error_is_runtime_error():
    assert issubclass(da.EnvError, RuntimeError)


def test_use_after_clear_raises():
    da.init(4, 3, 100)
    x = da.var(0) + 1.0
    da.clear()
    with pytest.raises(da.EnvError):
        x + x
    with pytest.raises(da.EnvError):
        x.con
    with pytest.raises(da.EnvError):
        da.exp(x)
    with pytest.raises(da.EnvError):
        str(x)
    with pytest.raises(da.EnvError):
        da.NDA(1.0)
    with pytest.raises(da.EnvError):
        da.var(0)

    # A live vector meeting one of a cleared env: check_env's logic_error.
    da.init(4, 3, 100)
    try:
        y = da.var(1)
        with pytest.raises(da.EnvError):
            y + x
    finally:
        da.clear()
    del x  # destructor of a retired-env vector must not crash


def test_pool_exhaustion_raises():
    da.init(2, 2, 10)
    try:
        keep = []
        with pytest.raises(RuntimeError, match="Run out of vectors"):
            for _ in range(11):
                keep.append(da.NDA(1.0))
        assert da.remain() == 0
        keep.pop()                                   # a slot freed at exhaustion
        assert da.remain() == 1                      # goes back on the free list
        keep.append(da.NDA(2.0))
        assert da.remain() == 0
    finally:
        da.clear()


def test_pool_exhaustion_on_returned_result_does_not_abort():
    code = textwrap.dedent("""
        import miradac as da
        da.init(2, 2, 10)
        keep = [da.NDA(1.0) for _ in range(7)]
        try:
            a = da.var(0)
        except RuntimeError:
            pass
    """)
    r = subprocess.run([sys.executable, "-c", code], capture_output=True)
    assert r.returncode == 0


def test_var_out_of_range(env):
    with pytest.raises(IndexError):
        da.var(99)
    with pytest.raises(IndexError):
        da.var(-1)
    with pytest.raises(IndexError):
        da.base[3]


def test_value_errors(env):
    with pytest.raises(ValueError):
        da.asin(da.NDA(2.0))          # std::domain_error
    with pytest.raises(ValueError):
        da.var(0) / 0.0               # std::invalid_argument
    with pytest.raises(ValueError):
        da.var(0).element([-1, 0, 0])


def test_no_env():
    with pytest.raises(da.EnvError):
        da.NDA()
    with pytest.raises(da.EnvError):
        da.current_order()
