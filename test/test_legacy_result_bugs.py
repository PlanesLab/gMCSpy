# Regression tests for result-affecting bugs found while planning central-framework
# (KD-1, KD-2, KD-3 in central-framework/KNOWN_DIFFERENCES.md).
import cobra
import pytest

from gMCSpy.calculateGMCS import buildGMatrix, calculate_GPR_MCS, calculateGeneMCS, mergeIsforms

solvers = ["cplex", "gurobi"]


def isozyme_model(n_isozymes):
    """EX -> A -(R1: g1 or ... or gn)-> B -(R2: g0)-> C -> BIO. gMCS: {g0}, {g1..gn}."""
    model = cobra.Model(f"toy_isozymes{n_isozymes}")
    a, b, c = (cobra.Metabolite(m, compartment="c") for m in "ABC")
    ex = cobra.Reaction("EX_A", lower_bound=0, upper_bound=10)
    ex.add_metabolites({a: 1})
    r1 = cobra.Reaction("R1", lower_bound=0, upper_bound=1000)
    r1.add_metabolites({a: -1, b: 1})
    r2 = cobra.Reaction("R2", lower_bound=0, upper_bound=1000)
    r2.add_metabolites({b: -1, c: 1})
    bio = cobra.Reaction("BIO", lower_bound=0, upper_bound=1000)
    bio.add_metabolites({c: -1})
    model.add_reactions([ex, r1, r2, bio])
    r1.gene_reaction_rule = " or ".join(f"g{i}" for i in range(1, n_isozymes + 1))
    r2.gene_reaction_rule = "g0"
    model.objective = "BIO"
    return model


class TestMaxKOLengthHonoured:
    # KD-1: buildGMatrix overwrote maxKOLength with 3, dropping rows longer than 3 genes.
    def test_gmatrix_keeps_four_gene_row(self):
        gObject = buildGMatrix(isozyme_model(4), 4, solver="gurobi")
        assert frozenset({"g1", "g2", "g3", "g4"}) in gObject["gDict"]

    def test_gmatrix_still_filters_longer_rows(self):
        gObject = buildGMatrix(isozyme_model(4), 3, solver="gurobi")
        assert frozenset({"g1", "g2", "g3", "g4"}) not in gObject["gDict"]

    @pytest.mark.parametrize("solver", solvers)
    def test_four_gene_gmcs_found(self, solver):
        solutions = calculateGeneMCS(isozyme_model(4), maxKOLength=4, solver=solver)
        found = {frozenset(v["solution"]) for v in solutions.values()}
        assert found == {frozenset({"g0"}), frozenset({"g1", "g2", "g3", "g4"})}


class TestEmptyGMatrix:
    # KD-2: an empty G matrix crashed with IndexError while building the MILP.
    @pytest.mark.parametrize("solver", solvers)
    def test_returns_no_solutions(self, solver):
        model = isozyme_model(4)
        model.reactions.R2.gene_reaction_rule = ""
        assert calculateGeneMCS(model, maxKOLength=3, solver=solver) == {}


class TestNestedGPRLength:
    # KD-3: mixed AND/OR GPRs were solved with MAX_LENGTH_MCS=5, losing longer cuts.
    @staticmethod
    def cuts(rule, maxKOLength=1e6):
        model = cobra.Model("m")
        reaction = cobra.Reaction("R")
        model.add_reactions([reaction])
        reaction.gene_reaction_rule = rule
        gDict = {}
        calculate_GPR_MCS(gDict, reaction.gpr, [reaction.id], mergeIsforms(None), maxKOLength)
        return set(gDict)

    def test_six_gene_cuts(self):
        assert self.cuts("(a and b) or c or d or e or f or g") == {
            frozenset("acdefg"),
            frozenset("bcdefg"),
        }

    def test_short_cuts_unchanged(self):
        assert self.cuts("(a and b) or c or d", maxKOLength=3) == {frozenset("acd"), frozenset("bcd")}
