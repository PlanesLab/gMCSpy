# Regression tests for result-affecting bugs found while planning central-framework
# (KD-* in central-framework/KNOWN_DIFFERENCES.md). Expected gMCS come from a brute-force
# oracle: knock out every gene subset in cobra and keep the minimal sets that stop growth.
import itertools
import random

import cobra
import pytest

from gMCSpy.calculateGMCS import (
    buildGMatrix,
    calculate_GPR_MCS,
    calculateGeneMCS,
    mergeIsforms,
)
from gMCSpy.ProblemDefinitions import buildDictionaryFBAWithDeletions

solvers = ["cplex", "gurobi"]
# (label, kwargs) for every enumeration path that must match the oracle
enumerations = [
    ("gurobi staged-lazy", dict(solver="gurobi")),
    ("gurobi classic", dict(solver="gurobi", enumeration="classic")),
    ("cplex staged-lazy", dict(solver="cplex")),
    ("cplex classic", dict(solver="cplex", enumeration="classic")),
    ("scip staged-lazy", dict(solver="scip")),
    ("scip classic", dict(solver="scip", enumeration="classic")),
]


def build_model(reactions, model_id="toy"):
    """reactions: list of (id, {met: coef}, lb, ub, gpr or None); objective is BIO."""
    model = cobra.Model(model_id)
    metabolites = {}
    cobra_reactions = []
    for rid, stoichiometry, lb, ub, _ in reactions:
        reaction = cobra.Reaction(rid, lower_bound=lb, upper_bound=ub)
        for met in stoichiometry:
            metabolites.setdefault(met, cobra.Metabolite(met, compartment="c"))
        reaction.add_metabolites({metabolites[m]: c for m, c in stoichiometry.items()})
        cobra_reactions.append(reaction)
    model.add_reactions(cobra_reactions)
    for rid, _, _, _, gpr in reactions:
        if gpr:
            model.reactions.get_by_id(rid).gene_reaction_rule = gpr
    model.objective = "BIO"
    return model


def brute_force_gmcs(model, max_length, tol=1e-6):
    genes = sorted(g.id for g in model.genes)
    cuts = []
    for k in range(1, max_length + 1):
        for subset in itertools.combinations(genes, k):
            if any(cut <= set(subset) for cut in cuts):
                continue
            with model:
                for gene in subset:
                    model.genes.get_by_id(gene).knock_out()
                growth = model.slim_optimize(error_value=0.0)
            if growth <= tol:
                cuts.append(frozenset(subset))
    return set(cuts)


def gmcs_list(model, **kwargs):
    kwargs.setdefault("saveSolutions", False)
    return [frozenset(v["solution"]) for v in calculateGeneMCS(model, **kwargs).values()]


def isozyme_model(n_isozymes):
    """EX -> A -(R1: g1 or ... or gn)-> B -(R2: g0)-> C -> BIO. gMCS: {g0}, {g1..gn}."""
    return build_model([
        ("EX_A", {"A": 1}, 0, 10, None),
        ("R1", {"A": -1, "B": 1}, 0, 1000, " or ".join(f"g{i}" for i in range(1, n_isozymes + 1))),
        ("R2", {"B": -1, "C": 1}, 0, 1000, "g0"),
        ("BIO", {"C": -1}, 0, 1000, None),
    ], f"toy_isozymes{n_isozymes}")


def chains_model(k, n):
    """k parallel routes A -> B, each a chain of n reactions with distinct genes."""
    reactions = [("EX", {"A": 1}, 0, 1000, None), ("BIO", {"B": -1}, 0, 1000, None)]
    for i in range(k):
        previous = "A"
        for j in range(n):
            following = "B" if j == n - 1 else f"P{i}_{j}"
            reactions.append((f"R{i}_{j}", {previous: -1, following: 1}, 0, 1000, f"p{i}g{j}"))
            previous = following
    return build_model(reactions, f"chains{k}x{n}")


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
        found = set(gmcs_list(isozyme_model(4), maxKOLength=4, solver=solver))
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


class TestAgainstBruteForce:
    @pytest.mark.parametrize("label,kwargs", enumerations)
    @pytest.mark.parametrize("k,n", [(2, 4), (3, 3)])
    def test_deepest_order_complete(self, label, kwargs, k, n):
        # KD-8: gurobi's default earlyStop returned 1 of 16 / 1 of 27 gMCS here.
        model = chains_model(k, n)
        assert set(gmcs_list(model, maxKOLength=k, **kwargs)) == brute_force_gmcs(model, k)

    @pytest.mark.parametrize("label,kwargs", enumerations)
    @pytest.mark.parametrize("r2_bounds,r2_gpr", [((0, 0), "g2"), ((-1000, 0), "g2"), ((0, 0), None)])
    def test_reaction_bounds(self, label, kwargs, r2_bounds, r2_gpr):
        # KD-9: closed and reverse-only reactions were treated as forward-capable.
        model = build_model([
            ("EX", {"X": 1}, 0, 1000, None),
            ("R1", {"X": -1, "Y": 1}, 0, 1000, "g1"),
            ("R2", {"X": -1, "Y": 1}, *r2_bounds, r2_gpr),
            ("BIO", {"Y": -1}, 0, 1000, None),
        ])
        assert set(gmcs_list(model, maxKOLength=3, **kwargs)) == {frozenset({"g1"})}

    @pytest.mark.parametrize("label,kwargs", enumerations)
    def test_no_gene_supersets(self, label, kwargs):
        # KD-10: no-good cuts on rows let {g2,g4,g6,g7} through besides {g2,g6,g7}.
        model = build_model([
            ("EX", {"A": 1}, 0, 1000, None),
            ("R1", {"A": -1, "B": 1}, 0, 1000, "g7 or g2"),
            ("R6", {"A": -1, "B": 1}, 0, 1000, "(g4 and g7) or g6"),
            ("BIO", {"B": -1}, 0, 1000, None),
        ])
        assert set(gmcs_list(model, maxKOLength=4, **kwargs)) == {frozenset({"g2", "g6", "g7"})}

    @pytest.mark.parametrize("label,kwargs", enumerations)
    def test_gene_count_weights(self, label, kwargs):
        # KD-11: pairwise-only intersections double-counted A, hiding all 3 gMCS at length 3.
        model = build_model([
            ("EX", {"X": 1}, 0, 1000, None),
            ("R1", {"X": -1, "Y": 1}, 0, 1000, "A or B or C"),
            ("R2", {"Y": -1, "Z": 1}, 0, 1000, "A or B or D"),
            ("R3", {"Z": -1, "W": 1}, 0, 1000, "A or C or D"),
            ("BIO", {"W": -1}, 0, 1000, None),
        ])
        assert set(gmcs_list(model, maxKOLength=3, **kwargs)) == brute_force_gmcs(model, 3)

    @pytest.mark.parametrize("label,kwargs", enumerations)
    def test_no_duplicates(self, label, kwargs):
        # KD-13: CPLEX reported {g1,g2} twice (optional zero-weight row in the pool).
        model = build_model([
            ("EX", {"A": 1}, 0, 1000, None),
            ("Ra", {"A": -1, "B": 1}, 0, 1000, "g1"),
            ("Rb", {"A": -1, "B": 1}, 0, 1000, "g2"),
            ("Rc", {"A": -1, "C": 1}, 0, 1000, "g1 or g2"),
            ("SC", {"C": -1}, 0, 1000, None),
            ("BIO", {"B": -1}, 0, 1000, None),
        ])
        found = gmcs_list(model, maxKOLength=3, **kwargs)
        assert len(found) == len(set(found))
        assert set(found) == brute_force_gmcs(model, 3)

    @pytest.mark.parametrize("label,kwargs", enumerations)
    def test_objective_coefficient_not_one(self, label, kwargs):
        # KD-16: objectiveVector.index(1) crashed when the coefficient was not exactly 1.
        model = isozyme_model(2)
        model.objective = {model.reactions.BIO: 2.5}
        assert set(gmcs_list(model, maxKOLength=2, **kwargs)) == brute_force_gmcs(model, 2)


def random_model(seed, n_metabolites=6, n_reactions=10, n_genes=8):
    rng = random.Random(seed)
    genes = [f"g{i}" for i in range(n_genes)]
    mets = [f"M{i}" for i in range(n_metabolites)]

    def random_gpr():
        t = rng.random()
        if t < 0.1:
            return None
        if t < 0.25:
            return rng.choice(genes)
        if t < 0.75:
            return " or ".join(rng.sample(genes, rng.randint(2, 3)))
        a, b, c = rng.sample(genes, 3)
        return rng.choice([f"({a} and {b}) or {c}", f"({a} or {b}) and {c}"])

    reactions = [("EX", {"M0": 1}, 0, 1000, None)]
    for i in range(n_reactions):
        a, b = rng.sample(mets, 2)
        u = rng.random()
        lb, ub = (-1000, 1000) if u < 0.3 else (0, 0) if u < 0.4 else (-1000, 0) if u < 0.45 else (0, 1000)
        reactions.append((f"R{i}", {a: -1, b: 1}, lb, ub, random_gpr()))
    reactions.append(("BIO", {mets[-1]: -1}, 0, 1000, None))
    return build_model(reactions, f"rand{seed}")


@pytest.mark.parametrize("label,kwargs", enumerations)
@pytest.mark.parametrize("seed", range(40))
def test_random_models_match_brute_force(label, kwargs, seed):
    model = random_model(seed)
    if model.slim_optimize(error_value=0.0) < 1e-6:
        pytest.skip("model cannot grow")
    expected = brute_force_gmcs(model, 3)
    found = gmcs_list(random_model(seed), maxKOLength=3, **kwargs)
    assert len(found) == len(set(found))
    assert set(found) == expected


TARGET_MODELS = {
    # route1 g1 in parallel with route2 (g2 or g3), then essential g4
    "isozymes": [
        ("EX", {"X": 1}, 0, 1000, None),
        ("R1", {"X": -1, "Y": 1}, 0, 1000, "g1"),
        ("R2", {"X": -1, "Y": 1}, 0, 1000, "g2 or g3"),
        ("R3", {"Y": -1, "Z": 1}, 0, 1000, "g4"),
        ("BIO", {"Z": -1}, 0, 1000, None),
    ],
    # same with a complex on route2
    "complex": [
        ("EX", {"X": 1}, 0, 1000, None),
        ("R1", {"X": -1, "Y": 1}, 0, 1000, "g1"),
        ("R2", {"X": -1, "Y": 1}, 0, 1000, "g2 and g3"),
        ("R3", {"Y": -1, "Z": 1}, 0, 1000, "g4"),
        ("BIO", {"Z": -1}, 0, 1000, None),
    ],
}


class TestTargetKOs:
    # KD-12: the targeted formulation crashed, returned [] or non-minimal/duplicated sets.
    # Expected: exactly the gMCS that contain at least one target gene.
    @pytest.mark.parametrize("label,kwargs", enumerations)
    @pytest.mark.parametrize("model_name", sorted(TARGET_MODELS))
    @pytest.mark.parametrize("targets", [["g1"], ["g4"], ["g2"], ["g1", "g3"]])
    def test_matches_brute_force(self, label, kwargs, model_name, targets):
        model = build_model(TARGET_MODELS[model_name])
        expected = {cut for cut in brute_force_gmcs(model, 3) if cut & set(targets)}
        found = gmcs_list(build_model(TARGET_MODELS[model_name]), maxKOLength=3, targetKOs=targets, **kwargs)
        assert len(found) == len(set(found))
        assert set(found) == expected

    @pytest.mark.parametrize("label,kwargs", enumerations)
    @pytest.mark.parametrize("seed", range(15))
    def test_random_models(self, label, kwargs, seed):
        model = random_model(seed)
        if model.slim_optimize(error_value=0.0) < 1e-6:
            pytest.skip("model cannot grow")
        target = sorted(g.id for g in model.genes)[seed % len(model.genes)]
        expected = {cut for cut in brute_force_gmcs(model, 3) if target in cut}
        found = gmcs_list(random_model(seed), maxKOLength=3, targetKOs=[target], **kwargs)
        assert len(found) == len(set(found))
        assert set(found) == expected

    def test_unknown_gene_raises(self):
        with pytest.raises(ValueError):
            gmcs_list(build_model(TARGET_MODELS["isozymes"]), maxKOLength=3, solver="gurobi", targetKOs=["nope"])


def nutrient_model():
    """Two nutrients, either suffices. EX_A reversible, EX_B uptake-only."""
    return build_model([
        ("EX_A", {"A": -1}, -1000, 1000, None),
        ("EX_B", {"B": -1}, -1000, 0, None),
        ("RA", {"A": -1, "Y": 1}, 0, 1000, "g1"),
        ("RB", {"B": -1, "Y": 1}, 0, 1000, "g2"),
        ("BIO", {"Y": -1}, 0, 1000, None),
    ])


class TestNutrient:
    # KD-14: uptake-only exchanges got no artificial gene, an explicit exchange list was
    # ignored, and the caller's model was modified.
    @pytest.mark.parametrize("label,kwargs", enumerations)
    def test_all_exchanges(self, label, kwargs):
        found = set(gmcs_list(nutrient_model(), maxKOLength=3, isNutrient=True, **kwargs))
        assert found == {
            frozenset({"g1", "g2"}),
            frozenset({"g1", "B_ag"}),
            frozenset({"A_ag", "g2"}),
            frozenset({"A_ag", "B_ag"}),
        }

    def test_explicit_exchange_list(self):
        found = set(gmcs_list(nutrient_model(), maxKOLength=3, isNutrient=True,
                              exchangeReactions=["EX_A"], solver="gurobi"))
        assert found == {frozenset({"g1", "g2"}), frozenset({"A_ag", "g2"})}

    def test_caller_model_unchanged(self):
        model = nutrient_model()
        before = {r.id: (r.bounds, r.gene_reaction_rule) for r in model.reactions}
        gmcs_list(model, maxKOLength=3, isNutrient=True, solver="gurobi")
        assert {r.id: (r.bounds, r.gene_reaction_rule) for r in model.reactions} == before


class TestNestedGPRSolver:
    # KD-19: the nested MCS for mixed AND/OR GPRs ignored `solver` and always used gurobi.
    @pytest.mark.parametrize("solver", ["scip", "cplex"])
    def test_nested_mcs_uses_chosen_solver(self, solver, monkeypatch):
        import gMCSpy.calculateGMCS as module

        used = []
        original = module.calculateMCS

        def recording(*args, **kwargs):
            used.append(kwargs.get("solver"))
            return original(*args, **kwargs)

        monkeypatch.setattr(module, "calculateMCS", recording)
        model = build_model(TARGET_MODELS["complex"])
        model.reactions.R2.gene_reaction_rule = "(g2 and g3) or g5"
        found = set(gmcs_list(model, maxKOLength=3, solver=solver))
        assert used and set(used) == {solver}
        assert found == brute_force_gmcs(model, 3)


class TestNumericalScaling:
    # Unscaled genome-scale models (coefficients 1e-6..1e2) let the gMCS dual drift to
    # |u| ~ 1e15 and accept certificates of non-cuts; w is now fixed and S is equilibrated.
    @staticmethod
    def cofactor_model():
        """Biomass needs a cofactor at 1e-6 (yeast-GEM pattern), made by g3 or imported."""
        return build_model([
            ("EX_A", {"A": 1}, 0, 1000, None),
            ("R1", {"A": -1, "B": 1}, 0, 1000, "g1 or g2"),
            ("RC", {"A": -1, "C": 1}, 0, 1000, "g3"),
            ("BIO", {"B": -1, "C": -1e-6}, 0, 1000, None),
        ], "cofactor")

    def test_scaling_is_exact_and_narrows_coefficients(self):
        import numpy as np
        from cobra.util.array import create_stoichiometric_matrix
        from gMCSpy.Scaling import scaleModel

        model = self.cofactor_model()
        scaled = scaleModel(model)
        before = create_stoichiometric_matrix(model)
        after = create_stoichiometric_matrix(scaled)
        nonzero = before != 0
        ratio = after[nonzero] / before[nonzero]
        assert np.all(ratio > 0) and np.all(np.log2(ratio) == np.round(np.log2(ratio)))
        assert np.ptp(np.log2(np.abs(after[nonzero]))) < np.ptp(np.log2(np.abs(before[nonzero])))
        assert brute_force_gmcs(scaled, 3) == brute_force_gmcs(model, 3) == {
            frozenset({"g1", "g2"}), frozenset({"g3"})}

    @pytest.mark.parametrize("label,kwargs", enumerations)
    def test_cofactor_gmcs(self, label, kwargs):
        assert set(gmcs_list(self.cofactor_model(), maxKOLength=3, **kwargs)) == {
            frozenset({"g1", "g2"}), frozenset({"g3"})}

    def test_w_is_fixed(self):
        from gMCSpy.ProblemDefinitions import buildDictionaryMCSGeneProblem

        model = self.cofactor_model()
        gObject = buildGMatrix(model, 3, solver="gurobi")
        problem = buildDictionaryMCSGeneProblem(
            model, gObject["gDict"], gObject["gMatrix"], relationships=gObject["relationships"],
            maxKOLength=3, numberNewGenesByKO=gObject["numberNewGenesByKO"])
        w = [v for v in problem.getVariables().values() if v["name"].startswith("w_")]
        assert len(w) == 1 and w[0]["lb"] == w[0]["ub"] == 1


class TestSolverParameters:
    # KD-21: parameters that could not be set were dropped silently (Warning(...) was
    # constructed but never emitted); the CPLEX path also built code strings for exec.
    @pytest.mark.parametrize("solver,good,bad", [
        ("cplex", {"mip.pool.intensity": 4}, {"mip.pool.notAParameter": 1}),
        ("gurobi", {"MIPFocus": 1}, {"NotAParameter": 1}),
    ])
    def test_parameters_applied_or_reported(self, solver, good, bad):
        from gMCSpy.OptimizationProblem import OptimizationProblem

        problem = OptimizationProblem()
        problem.addVariable(0, 0, 1, vtype="binary")
        problem.addObjective([0], [1], "minimize")
        problem.setSolver(solver)
        with pytest.warns(UserWarning, match="could not be set"):
            optModel, _ = problem.interpretProblem(parameters={**good, **bad})
        if solver == "cplex":
            assert optModel.parameters.mip.pool.intensity.get() == 4
        else:
            assert optModel.Params.MIPFocus == 1


class TestIsoformSeparator:
    # KD-15: split(sep)[0] cut at the first separator, merging A_B_1 and A_C_1 into A.
    def test_rsplit(self):
        model = build_model([
            ("EX", {"A": 1}, 0, 1000, None),
            ("R1", {"A": -1, "B": 1}, 0, 1000, "A_B_1 or A_C_1"),
            ("BIO", {"B": -1}, 0, 1000, None),
        ])
        gObject = buildGMatrix(model, 3, isoformSeparator="_", solver="gurobi")
        assert frozenset({"A_B", "A_C"}) in gObject["gDict"]


class TestValidatorGPREvaluation:
    # KD-18: the FBA validator evaluated GPRs by substring replace + eval, so knocking out
    # g1 also rewrote g10 ("False0").
    def test_gene_name_prefix(self):
        model = build_model([
            ("EX", {"A": 1}, 0, 1000, None),
            ("R1", {"A": -1, "B": 1}, 0, 1000, "g1 or g10"),
            ("BIO", {"B": -1}, 0, 1000, None),
        ])
        problem = buildDictionaryFBAWithDeletions(model, ["g1"], None)
        problem.setSolver("gurobi")
        optModel, _ = problem.interpretProblem(parameters={"OutputFlag": 0})
        optModel.optimize()
        assert optModel.ObjVal > 1e-6  # g10 still catalyses R1
