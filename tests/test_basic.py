import ast
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.Utils.calculations import mylog


@pytest.mark.parametrize(
    ('first', 'second', 'expected_first', 'expected_second'),
    [
        ([1.0, 10.0, 100.0], [1.0, -1.0, 100.0], [0.0, 2.0], [0.0, 2.0]),
        ([1.0, 10.0, 100.0, -1.0, 10000.0], [-1.0, 10.0, 100.0, 1000.0, 100000.0], [1.0, 2.0, 4.0], [1.0, 2.0, 5.0]),
        ([0.0, np.nan, np.inf, 10.0, -1.0], [1.0, 2.0, 3.0, 100.0, 4.0], [1.0], [2.0]),
        ([0.0, -1.0, np.nan], [1.0, 2.0, 3.0], [], []),
        ([], [], [], []),
    ],
)
def test_mylog_filters_pairs_jointly_without_domain_warnings(first, second, expected_first, expected_second):
    with np.errstate(divide='raise', invalid='raise'):
        actual_first, actual_second = mylog(first, second)
    np.testing.assert_array_equal(actual_first, expected_first)
    np.testing.assert_array_equal(actual_second, expected_second)


@pytest.mark.parametrize('dtype', [np.float32, np.float64])
def test_mylog_preserves_all_positive_values_and_input_arrays(dtype):
    first = np.array([0.1, 1.0, 2.0, 10.0, 100.0], dtype=dtype)
    second = np.array([0.01, 0.2, 3.0, 100.0, 1000.0], dtype=dtype)
    original_first, original_second = first.copy(), second.copy()

    actual_first, actual_second = mylog(first, second)

    np.testing.assert_array_equal(actual_first, np.log10(first))
    np.testing.assert_array_equal(actual_second, np.log10(second))
    assert actual_first.dtype == first.dtype
    assert actual_second.dtype == second.dtype
    np.testing.assert_array_equal(first, original_first)
    np.testing.assert_array_equal(second, original_second)


def test_mylog_respects_masks_in_both_inputs():
    first = np.ma.array([1.0, 10.0, 100.0], mask=[False, True, False])
    second = np.ma.array([10.0, 100.0, 1000.0], mask=[True, False, False])

    actual_first, actual_second = mylog(first, second)

    np.testing.assert_array_equal(actual_first, [2.0])
    np.testing.assert_array_equal(actual_second, [3.0])


def test_mylog_handles_matching_dataframe_column_positions():
    frame = pd.DataFrame({'first': [1.0, 10.0, 100.0], 'second': [10.0, -1.0, 1000.0]}, index=[8, 3, 5])
    actual_first, actual_second = mylog(frame['first'], frame['second'])
    np.testing.assert_array_equal(actual_first, [0.0, 2.0])
    np.testing.assert_array_equal(actual_second, [1.0, 3.0])


def test_mylog_rejects_shapes_that_would_broadcast_unpaired_observations():
    with pytest.raises(ValueError, match='same shape'):
        mylog(np.ones((3, 1)), np.ones(3))


@pytest.mark.parametrize('filename', ['Stage_03a_original.py', 'Stage_03b_original.py', 'Stage_03b.py'])
def test_plotting_scripts_use_the_shared_calculations_without_starting_pipelines(filename):
    path = Path(__file__).parents[1] / 'scripts' / filename
    print(f'Reading from {path}...')
    tree = ast.parse(path.read_text())
    assert not any(isinstance(node, ast.FunctionDef) and node.name in {'mylog', 'percentage'} for node in tree.body)
    imports = [
        node for node in tree.body
        if isinstance(node, ast.ImportFrom) and node.module in {'src.Utils', 'src.Utils.calculations'}
    ]
    namespace = {}

    exec(compile(ast.Module(body=imports, type_ignores=[]), str(path), 'exec'), namespace)

    assert namespace['mylog'] is mylog
    assert namespace['percentage'] is namespace['calc'].percentage
    first, second = namespace['mylog']([1.0, 10.0, 100.0], [1.0, -1.0, 100.0])
    np.testing.assert_array_equal(first, [0.0, 2.0])
    np.testing.assert_array_equal(second, [0.0, 2.0])
    np.testing.assert_array_equal(namespace['percentage'](4, np.array([1, 2])), [25.0, 50.0])
