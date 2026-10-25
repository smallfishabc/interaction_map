# -*- coding: utf-8 -*-
"""
Created on Fri Oct 24 2026

Shape checks for the tables that move between stages.

These catch the failures that are silent rather than loud: a truncated CSV, a file written
by an older version with different columns, a contact probability that is not a
probability, a separation column that no longer agrees with the residue pair it was
derived from.

`validate` returns a list of problems rather than raising, so a caller can decide whether
a given problem is fatal. `require` raises.

@author: ShaharGroup-fyu
"""


class Schema(object):
    """
    Attributes:
        name      what the table is, used in messages
        columns   columns that must be present
        checks    (description, predicate) pairs run against the dataframe
    """

    def __init__(self, name, columns, checks=()):
        self.name = name
        self.columns = list(columns)
        self.checks = list(checks)


def _probabilities_in_range(df):
    return bool(((df['cont_prob'] >= 0.0) & (df['cont_prob'] <= 1.0)).all())


def _separation_matches_pair(df):
    return bool((df['distance'] == df['r_2'] - df['r_1']).all())


def _pairs_are_ordered(df):
    return bool((df['r_2'] > df['r_1']).all())


def _levels_are_known(df):
    return bool(df['plot_value'].isin([-2, -1, 0, 1, 2]).all())


def _no_missing_values(df):
    return bool(df.notna().all().all())


# The table ContactProbData produces and saves.
CONTACT_TABLE = Schema(
    'contact table',
    ['r_1', 'r_2', 'cont_prob', 'distance'],
    [('contact probabilities lie between 0 and 1', _probabilities_in_range),
     ('distance equals r_2 - r_1', _separation_matches_pair),
     ('r_2 is greater than r_1', _pairs_are_ordered),
     ('no missing values', _no_missing_values)])

# The same table after normalization.
INTERACTION_TABLE = Schema(
    'interaction table',
    ['r_1', 'r_2', 'cont_prob', 'distance', 'gs_standard', 'relative_strength', 'plot_value'],
    CONTACT_TABLE.checks + [('interaction levels are -2..2', _levels_are_known)])


def validate(df, schema):
    """
    :param df: dataframe to check
    :param schema: a Schema
    :return: list of problem descriptions, empty when the table is sound
    """
    problems = []
    missing = [c for c in schema.columns if c not in df.columns]
    if missing:
        problems.append('{0}: missing column(s) {1}'.format(schema.name, ', '.join(missing)))
        # The checks below index those columns, so stop here.
        return problems
    if len(df) == 0:
        problems.append('{0}: no rows'.format(schema.name))
        return problems
    for description, predicate in schema.checks:
        try:
            if not predicate(df):
                problems.append('{0}: {1} -- failed'.format(schema.name, description))
        except Exception as exc:
            problems.append('{0}: {1} -- could not be checked ({2})'.format(
                schema.name, description, exc))
    return problems


def require(df, schema):
    """Raise ValueError if the table does not satisfy the schema."""
    problems = validate(df, schema)
    if problems:
        raise ValueError('; '.join(problems))
    return df


def warn_if_invalid(df, schema):
    """Print each problem and return whether the table was sound."""
    problems = validate(df, schema)
    for problem in problems:
        print('WARNING: {0}'.format(problem))
    return not problems
