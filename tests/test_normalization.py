"""Tests for normalization module."""

import numpy as np
import pandas as pd
import pytest

from idp_interaction_map.normalization import (
    fitting_function,
    normalize_interaction_map,
)


def test_fitting_function():
    """Test power law fitting function."""
    df = pd.DataFrame({"distance": [1, 2, 5, 10]})
    result = fitting_function(df, a=13.12, b=-2.32)

    # Check it returns proper series
    assert isinstance(result, pd.Series)
    assert len(result) == 4

    # Check power law behavior (larger distance = smaller value)
    assert result.iloc[0] > result.iloc[1] > result.iloc[2] > result.iloc[3]


def test_normalize_interaction_map(sample_contact_data):
    """Test interaction map normalization."""
    result = normalize_interaction_map(sample_contact_data)

    # Check new columns are added
    assert "gs_standard" in result.columns
    assert "relative_strength" in result.columns
    assert "plot_value" in result.columns

    # Check dimensions
    assert len(result) == len(sample_contact_data)

    # Check plot values are in expected range
    assert result["plot_value"].isin([-2, -1, 0, 1, 2]).all()


def test_normalize_zero_contact_probability():
    """Test normalization handles zero contact probability."""
    df = pd.DataFrame(
        {
            "r_1": [1, 2],
            "r_2": [5, 6],
            "cont_prob": [0.0, 0.5],
            "distance": [4, 4],
        }
    )

    result = normalize_interaction_map(df)

    # Zero contact probability should give zero relative strength
    assert result.loc[0, "relative_strength"] == 0
    assert result.loc[1, "relative_strength"] != 0


def test_normalize_custom_cutoffs():
    """Test normalization with custom cutoff values."""
    df = pd.DataFrame(
        {
            "r_1": [1, 2, 3],
            "r_2": [5, 6, 7],
            "cont_prob": [0.8, 0.5, 0.2],
            "distance": [4, 4, 4],
        }
    )

    custom_cutoffs = (3, 1.5, -1.5, -3)
    custom_values = (3, 2, 0, -2, -3)

    result = normalize_interaction_map(
        df, inter_cutoff=custom_cutoffs, value_list=custom_values
    )

    # Check custom values are used
    assert set(result["plot_value"].unique()).issubset({-3, -2, 0, 2, 3})


def test_normalize_preserves_original_columns(sample_contact_data):
    """Test that normalization preserves original data columns."""
    original_columns = set(sample_contact_data.columns)
    result = normalize_interaction_map(sample_contact_data)

    # All original columns should still exist
    assert original_columns.issubset(set(result.columns))

    # Original data should be unchanged
    for col in original_columns:
        pd.testing.assert_series_equal(
            sample_contact_data[col], result[col], check_names=True
        )
