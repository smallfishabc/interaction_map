"""Normalization and interaction strength calculation module."""

import logging
from typing import List, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def fitting_function(df: pd.DataFrame, a: float, b: float) -> pd.Series:
    """
    Calculate standard contact probability curve.

    Uses power law relationship: P(d) = a * d^b

    Args:
        df: DataFrame containing 'distance' column
        a: Coefficient parameter
        b: Exponent parameter (typically negative)

    Returns:
        Series of fitted values
    """
    return a * df["distance"] ** b


def normalize_interaction_map(
    target_map: pd.DataFrame,
    a1: float = 13.12,
    b1: float = -2.32,
    inter_cutoff: Tuple[float, ...] = (2, 1, -1, -2),
    value_list: Tuple[int, ...] = (2, 1, 0, -1, -2),
) -> pd.DataFrame:
    """
    Normalize contact map and categorize interaction strengths.

    Compares observed contact probabilities against ideal polymer model
    to identify favorable/unfavorable interactions.

    Default parameters are for coarse-grained (CG) simulations.
    For all-atom simulations with CA atoms, use: a1=1.64, b1=-1.32

    Args:
        target_map: DataFrame with contact probability data
        a1: Coefficient for standard curve fitting (CG: 13.12, all-atom: 1.64)
        b1: Exponent for standard curve fitting (CG: -2.32, all-atom: -1.32)
        inter_cutoff: Thresholds for categorizing relative strengths
                      CG: (2, 1, -1, -2), all-atom: (1.5, 0.5, -1, -2)
        value_list: Integer labels for interaction categories

    Returns:
        DataFrame with added columns:
            - gs_standard: Expected contact probability from ideal polymer
            - relative_strength: Log ratio of observed/expected
            - plot_value: Categorical interaction strength
    """
    logger.info("Normalizing interaction map")

    # Calculate expected contact probability from ideal polymer model
    target_map["gs_standard"] = target_map.apply(fitting_function, axis=1, args=(a1, b1))

    # Calculate relative strength (log ratio)
    target_map["relative_strength"] = np.where(
        target_map["cont_prob"] == 0,
        0,
        np.log(target_map["cont_prob"] / target_map["gs_standard"]),
    )

    # Categorize interactions based on relative strength
    col = "relative_strength"
    conditions = [
        target_map[col] >= inter_cutoff[0],
        (target_map[col] >= inter_cutoff[1]) & (target_map[col] < inter_cutoff[0]),
        (target_map[col] < inter_cutoff[1]) & (target_map[col] > inter_cutoff[2]),
        (target_map[col] <= inter_cutoff[2]) & (target_map[col] > inter_cutoff[3]),
        target_map[col] <= inter_cutoff[3],
    ]

    target_map["plot_value"] = np.select(conditions, value_list, default=0)

    logger.info(
        f"Categorized {len(target_map)} interactions: "
        f"{(target_map['plot_value'] > 0).sum()} favorable, "
        f"{(target_map['plot_value'] < 0).sum()} unfavorable"
    )

    return target_map
