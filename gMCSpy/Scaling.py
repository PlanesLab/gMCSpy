import numpy as np
from cobra.util.array import create_stoichiometric_matrix


def scaleModel(cobraModel, passes=4):
    """Return a copy of the model whose stoichiometric matrix is equilibrated by powers of two.

    Rows (metabolites) and columns (reactions) are scaled so that the largest and smallest
    coefficients of each are centred around 1. This is a change of units: the same
    reactions can carry flux, so the gMCS are unchanged, and powers of two keep every
    coefficient exact in floating point. Flux bounds are rescaled with their reactions.
    Genome-scale models span coefficients from 1e-6 to ~1e2; unscaled, the gMCS dual
    could drift to |u| ~ 1e15 and accept certificates of gene sets that are not cuts.
    """
    S = create_stoichiometric_matrix(cobraModel, array_type="dok").tocoo()
    rows, cols, logs = S.row, S.col, np.log2(np.abs(S.data))
    rowExp = np.zeros(S.shape[0])
    colExp = np.zeros(S.shape[1])
    for _ in range(passes):
        rowExp -= np.round(_centre(logs + rowExp[rows] + colExp[cols], rows, S.shape[0]))
        colExp -= np.round(_centre(logs + rowExp[rows] + colExp[cols], cols, S.shape[1]))
    rowFactor, colFactor = np.exp2(rowExp), np.exp2(colExp)

    scaled = cobraModel.copy()
    metaboliteIndex = {m.id: i for i, m in enumerate(scaled.metabolites)}
    for j, reaction in enumerate(scaled.reactions):
        reaction.add_metabolites(
            {m: c * rowFactor[metaboliteIndex[m.id]] * colFactor[j] for m, c in reaction.metabolites.items()},
            combine=False,
        )
        # flux = colFactor * scaled flux, so the bounds divide by it
        reaction.bounds = (reaction.lower_bound / colFactor[j], reaction.upper_bound / colFactor[j])
    return scaled


def _centre(logs, index, size):
    high = np.full(size, -np.inf)
    low = np.full(size, np.inf)
    np.maximum.at(high, index, logs)
    np.minimum.at(low, index, logs)
    centre = (high + low) / 2
    centre[~np.isfinite(centre)] = 0.0
    return centre
