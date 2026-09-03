import numpy as np
import pytest

from ruptura import (
    MIXTURE_FILE_CHEMISORPTION_COLUMNS,
    MIXTURE_FILE_COLUMNS,
    read_mixture_prediction,
)


def write(tmp_path, rows):
    path = tmp_path / "component_1_CO2.data"
    path.write_text(
        "# column 1: total pressure [Pa]\n"
        + "\n".join(" ".join(f"{value:.14g}" for value in row) for row in rows)
        + "\n",
        encoding="utf-8",
    )
    return path


PHYSICAL_ROWS = [
    [1.0e2, 0.10, 0.09, 1.0, 1.0, 1.0e2, 0.5],
    [1.0e5, 2.10, 1.90, 1.0, 1.0, 1.0e5, 3.5],
]
CHEMICAL_ROWS = [row + [2.5, 2.4, row[1] + 2.5, row[2] + 2.4] for row in PHYSICAL_ROWS]


def test_seven_column_output_keeps_the_physisorption_labels(tmp_path):
    result = read_mixture_prediction(write(tmp_path, PHYSICAL_ROWS))

    assert result.columns == MIXTURE_FILE_COLUMNS
    assert result.kind == "mixture_prediction_file"
    np.testing.assert_allclose(result.column("pressure"), [1.0e2, 1.0e5])
    np.testing.assert_allclose(result.column("mixture_physisorption_loading"), [0.09, 1.90])


def test_eleven_column_output_exposes_the_chemisorption_labels(tmp_path):
    result = read_mixture_prediction(write(tmp_path, CHEMICAL_ROWS))

    assert result.columns == MIXTURE_FILE_COLUMNS + MIXTURE_FILE_CHEMISORPTION_COLUMNS
    np.testing.assert_allclose(result.column("mixture_chemisorption_loading"), [2.4, 2.4])
    np.testing.assert_allclose(
        result.column("mixture_total_loading"),
        result.column("mixture_physisorption_loading") + result.column("mixture_chemisorption_loading"),
        rtol=1.0e-12,
        atol=1.0e-12,
    )


def test_unexpected_width_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="expected 7 or 11"):
        read_mixture_prediction(write(tmp_path, [[1.0, 2.0, 3.0]]))
