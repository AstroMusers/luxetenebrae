import ast
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.Utils.calculations import mylog


@pytest.mark.parametrize('verbosity', ['0', '1'])
def test_rename_runs_preserves_contents_and_uses_numeric_order(tmp_path, monkeypatch, capsys, verbosity):
    from src.Utils.utils import rename_runs

    for number in (10, 2):
        folder = tmp_path / f'COMPAS_Output_{number} (1)'
        folder.mkdir()
        (folder / 'result.txt').write_text(f'original run {number}', encoding='utf-8')
    unrelated = tmp_path / 'other_run'
    unrelated.mkdir()
    matching_file = tmp_path / 'COMPAS_Output_7 (1)'
    matching_file.write_text('not a run directory', encoding='utf-8')
    monkeypatch.setenv('TDPY_VERBOSITY', verbosity)

    rename_runs(100, tmp_path)

    assert (tmp_path / 'COMPAS_Output_100' / 'result.txt').read_text() == 'original run 2'
    assert (tmp_path / 'COMPAS_Output_101' / 'result.txt').read_text() == 'original run 10'
    assert unrelated.is_dir()
    assert matching_file.read_text() == 'not a run directory'
    assert not (tmp_path / 'COMPAS_Output_2 (1)').exists()
    assert not (tmp_path / 'COMPAS_Output_10 (1)').exists()
    expected = (f'Reading from {tmp_path}...\n'
                f'Writing to {tmp_path / "COMPAS_Output_100"}...\n'
                f'Writing to {tmp_path / "COMPAS_Output_101"}...\n')
    assert capsys.readouterr().out == (expected if verbosity == '1' else '')


@pytest.mark.parametrize('destination_kind', ['directory', 'file', 'broken_symlink'])
def test_rename_runs_rejects_all_destination_collisions_before_moving(tmp_path, destination_kind):
    from src.Utils.utils import rename_runs

    for number in (2, 10):
        folder = tmp_path / f'COMPAS_Output_{number} (1)'
        folder.mkdir()
        (folder / 'result.txt').write_text(f'original run {number}', encoding='utf-8')
    collision = tmp_path / 'COMPAS_Output_101'
    if destination_kind == 'directory':
        collision.mkdir()
    elif destination_kind == 'file':
        collision.write_text('existing destination', encoding='utf-8')
    else:
        collision.symlink_to(tmp_path / 'absent_target')
    original_entries = {entry.name: entry.lstat().st_ino for entry in tmp_path.iterdir()}

    with pytest.raises(FileExistsError, match='COMPAS_Output_101'):
        rename_runs(100, tmp_path)

    assert {entry.name: entry.lstat().st_ino for entry in tmp_path.iterdir()} == original_entries
    for number in (2, 10):
        assert (tmp_path / f'COMPAS_Output_{number} (1)' / 'result.txt').read_text() == f'original run {number}'
    if destination_kind == 'file':
        assert collision.read_text() == 'existing destination'


@pytest.mark.parametrize("value", [0.0, -0.012, 1.234, 9876.0])
def test_one_significant_figure_formatter_matches_python_format(value):
    from src.Utils import plots

    assert plots.format_to_1sf(value, 0) == f"{value:.1g}"


def test_kroupa_analysis_import_does_not_run_plots_or_print(capsys, monkeypatch):
    import matplotlib.pyplot as plt

    source = Path(__file__).parents[1] / 'scripts' / 'kroupa_dist_test.py'
    spec = importlib.util.spec_from_file_location('kroupa_import_check', source)
    module = importlib.util.module_from_spec(spec)
    with plt.rc_context():
        original_size = plt.rcParams['savefig.dpi']
        monkeypatch.setattr(plt, 'figure', lambda *args, **kwargs: pytest.fail('Analysis plotted during import'))
        spec.loader.exec_module(module)
        assert plt.rcParams['savefig.dpi'] == original_size
    assert capsys.readouterr().out == ''
    assert callable(module.main)
    assert module.kroupa_pdf(10.0) > 0


@pytest.mark.parametrize('x_labels', [None, ['First time [Myr]', 'Second time [Myr]']])
@pytest.mark.parametrize('panel_count', [1, 2])
def test_time_evolution_labels_axes_without_changing_series(x_labels, panel_count, tmp_path, monkeypatch):
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from src.Utils import plots

    times = np.array([0.0, 1.0, 2.0])  # [Myr]
    series = [np.array([1.0, 2.0, 3.0]), np.array([4.0, 3.0, 2.0])][:panel_count]
    original_times = times.copy()
    original_series = [values.copy() for values in series]
    figures = []
    monkeypatch.setattr(plots.plt, 'savefig', lambda *args, **kwargs: figures.append(plots.plt.gcf()))

    plots.plot_time_evolution(
        [times], [[values] for values in series], ['First quantity', 'Second quantity'][:panel_count],
        ['Model'], str(tmp_path), x_labels=None if x_labels is None else x_labels[:panel_count],
    )

    figure = figures[0]
    canvas = FigureCanvasAgg(figure)
    canvas.draw()
    if x_labels is None:
        assert figure._supxlabel.get_text() == 'Time [Myr]'
    else:
        assert [axis.get_xlabel() for axis in figure.axes] == x_labels[:panel_count]
    for axis, expected in zip(figure.axes, original_series):
        np.testing.assert_array_equal(axis.lines[0].get_xdata(), original_times)
        np.testing.assert_array_equal(axis.lines[0].get_ydata(), expected)
    np.testing.assert_array_equal(times, original_times)
    for actual, expected in zip(series, original_series):
        np.testing.assert_array_equal(actual, expected)
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize('rows,cols', [(1, 1), (1, 2), (2, 1), (2, 2)])
@pytest.mark.parametrize('label_axes', [False, True])
def test_property_histograms_support_all_subplot_shapes(rows, cols, label_axes, tmp_path, monkeypatch):
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from src.Utils import plots

    samples = [np.linspace(0.2, 1.2, 41) + 0.1 * index for index in range(rows * cols)]
    originals = [values.copy() for values in samples]
    labels = [f'Population {index + 1}' for index in range(rows * cols)]
    x_labels = [f'Mass parameter {index + 1} [M_sun]' for index in range(rows * cols)] if label_axes else None
    y_labels = [f'Population {index + 1} count' for index in range(rows * cols)] if label_axes else None
    figures = []
    monkeypatch.setattr(plots.plt, 'savefig', lambda *args, **kwargs: figures.append(plots.plt.gcf()))

    plots.plot_bhms_properties_histogram(
        rows, cols, samples, samples, y_labels, labels, str(tmp_path) + '/',
        x_labels=x_labels,
    )

    figure = figures[0]
    canvas = FigureCanvasAgg(figure)
    canvas.draw()
    assert len(figure.axes) == rows * cols
    for index, (axis, original, label) in enumerate(zip(figure.axes, originals, labels)):
        heights, edges = np.histogram(original, bins=20)
        np.testing.assert_array_equal([rectangle.get_height() for rectangle in axis.patches], heights)
        np.testing.assert_allclose([rectangle.get_x() for rectangle in axis.patches], edges[:-1])
        np.testing.assert_allclose([rectangle.get_width() for rectangle in axis.patches], np.diff(edges))
        assert axis.get_title() == label
        assert axis.get_xlabel() == (x_labels[index] if label_axes else '')
        assert axis.get_ylabel() == (y_labels[index] if label_axes else '')
    for actual, original in zip(samples, originals):
        np.testing.assert_array_equal(actual, original)
    assert list(tmp_path.iterdir()) == []


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
