from circhic import utils
import numpy as np
import pytest


def test_convert_xy_to_thetar():
    lengths = np.array([42])
    random_state = np.random.RandomState(seed=42)
    n = 100
    x = random_state.randint(0, lengths.sum(), n)
    y = random_state.randint(0, lengths.sum(), n)

    # Flush test to check the code runs
    theta, r = utils.convert_xy_to_thetar((x, y), lengths)
