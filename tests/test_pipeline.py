import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / 'IS_proximity_publication.py'
spec = importlib.util.spec_from_file_location('is_proximity', SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_interval_distance_linear_and_circular():
    assert int(module.interval_distances(10, 20, np.array([30]), np.array([40]))[0]) == 10
    assert int(module.interval_distances(10, 20, np.array([15]), np.array([25]))[0]) == 0
    assert int(module.interval_distances(2, 5, np.array([95]), np.array([98]), 100)[0]) == 4


def test_canonical_reverse_order():
    assert module.canonical('IS26', ['a', 'b'], 'IS903') == module.canonical('IS903', ['b', 'a'], 'IS26')


def test_example_loads():
    root = Path(__file__).resolve().parents[1] / 'examples'
    t = module.load_targets(root / 'amr.xlsx', 'AMR')
    i = module.load_is(root / 'is.csv', r'^IS')
    lengths = module.read_lengths(root / 'lengths.csv')
    module.validate_lengths(t, i, lengths, True)
    nearest, nearby = module.analyze_proximity(t, i, lengths, True, [1000, 5000, 10000], 10000)
    assert len(nearest) == len(t)
    assert nearest['Distance_bp'].isna().sum() >= 1
    assert len(nearby) >= 1
    mods, fams = module.call_modules(t, i, 30000)
    assert isinstance(mods, pd.DataFrame)
    assert isinstance(fams, pd.DataFrame)


def test_permutation_reproducibility():
    root = Path(__file__).resolve().parents[1] / 'examples'
    t = module.load_targets(root / 'amr.xlsx', 'AMR')
    i = module.load_is(root / 'is.csv', r'^IS')
    lengths = module.read_lengths(root / 'lengths.csv')
    a = module.permutation_test(t, i, lengths, True, [1000], 20, 42)
    b = module.permutation_test(t, i, lengths, True, [1000], 20, 42)
    pd.testing.assert_frame_equal(a, b)
