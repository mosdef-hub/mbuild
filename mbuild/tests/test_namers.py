import itertools

import pytest

from mbuild.path.namers import (
    BeadNamer,
    ConstantNamer,
    CyclicNamer,
    GradientNamer,
    MarkovNamer,
    RandomNamer,
)


class TestBeadNamerCoerce:
    def test_coerce_string_returns_constant_namer(self):
        namer = BeadNamer.coerce("_A")
        assert isinstance(namer, ConstantNamer)

    def test_coerce_bead_namer_returns_same_object(self):
        original = ConstantNamer("_A")
        assert BeadNamer.coerce(original) is original

    def test_coerce_bad_type_raises(self):
        with pytest.raises(TypeError, match="bead_name must be a str or BeadNamer"):
            BeadNamer.coerce(42)

    def test_coerce_none_raises(self):
        with pytest.raises(TypeError):
            BeadNamer.coerce(None)

    def test_namers_are_iterators(self):
        for namer in [
            ConstantNamer("_A"),
            RandomNamer(["_A", "_B"], seed=0),
            CyclicNamer(["_A", "_B"]),
        ]:
            assert iter(namer) is namer


class TestConstantNamer:
    def test_always_returns_same_name(self):
        namer = ConstantNamer("_A")
        assert [next(namer) for _ in range(5)] == ["_A"] * 5

    def test_different_name(self):
        namer = ConstantNamer("_B")
        assert next(namer) == "_B"

    def test_repr(self):
        assert repr(ConstantNamer("_A")) == "ConstantNamer('_A')"

    def test_coerce_string_behaves_like_constant(self):
        namer = BeadNamer.coerce("_X")
        assert [next(namer) for _ in range(3)] == ["_X"] * 3


class TestRandomNamer:
    def test_draws_from_pool(self):
        namer = RandomNamer(["_A", "_B"], seed=0)
        results = {next(namer) for _ in range(50)}
        assert results <= {"_A", "_B"}

    def test_uniform_weights_use_both_names(self):
        namer = RandomNamer(["_A", "_B"], seed=42)
        results = [next(namer) for _ in range(100)]
        assert "_A" in results and "_B" in results

    def test_weight_one_always_picks_that_name(self):
        namer = RandomNamer(["_A", "_B"], weights=[1.0, 0.0], seed=0)
        assert all(next(namer) == "_A" for _ in range(20))

    def test_weights_are_normalized(self):
        namer = RandomNamer(["_A", "_B"], weights=[2, 8], seed=0)
        results = [next(namer) for _ in range(1000)]
        ratio = results.count("_B") / len(results)
        assert 0.75 < ratio < 0.85

    def test_seeded_reproducibility(self):
        seq_a = [next(RandomNamer(["_A", "_B"], seed=7)) for _ in range(20)]
        seq_b = [next(RandomNamer(["_A", "_B"], seed=7)) for _ in range(20)]
        assert seq_a == seq_b

    def test_different_seeds_differ(self):
        seq_a = [next(RandomNamer(["_A", "_B", "_C"], seed=1)) for _ in range(30)]
        seq_b = [next(RandomNamer(["_A", "_B", "_C"], seed=2)) for _ in range(30)]
        assert seq_a != seq_b

    def test_single_name_pool(self):
        namer = RandomNamer(["_A"], seed=0)
        assert [next(namer) for _ in range(5)] == ["_A"] * 5


class TestCyclicNamer:
    def test_flat_sequence_alternating(self):
        namer = CyclicNamer(["_A", "_B"])
        assert [next(namer) for _ in range(6)] == ["_A", "_B", "_A", "_B", "_A", "_B"]

    def test_flat_sequence_three(self):
        namer = CyclicNamer(["_A", "_B", "_C"])
        assert [next(namer) for _ in range(6)] == ["_A", "_B", "_C", "_A", "_B", "_C"]

    def test_block_sequence(self):
        namer = CyclicNamer([("_A", 3), ("_B", 2)])
        expected = ["_A", "_A", "_A", "_B", "_B", "_A", "_A", "_A", "_B", "_B"]
        assert [next(namer) for _ in range(10)] == expected

    def test_block_sequence_equal_blocks(self):
        namer = CyclicNamer([("_A", 5), ("_B", 5)])
        result = [next(namer) for _ in range(20)]
        assert result[:5] == ["_A"] * 5
        assert result[5:10] == ["_B"] * 5
        assert result[10:15] == ["_A"] * 5
        assert result[15:20] == ["_B"] * 5

    def test_mixed_flat_and_block(self):
        namer = CyclicNamer(["_A", ("_B", 2)])
        expected = ["_A", "_B", "_B", "_A", "_B", "_B"]
        assert [next(namer) for _ in range(6)] == expected

    def test_single_name_cycles_forever(self):
        namer = CyclicNamer(["_A"])
        assert [next(namer) for _ in range(5)] == ["_A"] * 5

    def test_empty_sequence_raises(self):
        with pytest.raises(ValueError, match="must not be empty"):
            CyclicNamer([])

    def test_never_exhausted(self):
        namer = CyclicNamer(["_A", "_B"])
        result = list(itertools.islice(namer, 10_000))
        assert len(result) == 10_000

    def test_repr(self):
        assert "CyclicNamer" in repr(CyclicNamer(["_A", "_B"]))


class TestRandomNamerStrict:
    def test_exact_composition_one_period(self):
        namer = RandomNamer(["_A", "_B"], weights=[1, 1], strict=True, seed=0)
        result = [next(namer) for _ in range(2)]
        assert sorted(result) == ["_A", "_B"]

    def test_exact_composition_over_many_periods(self):
        namer = RandomNamer(["_A", "_B"], weights=[1, 1], strict=True, seed=0)
        result = [next(namer) for _ in range(100)]
        assert result.count("_A") == result.count("_B") == 50

    def test_weighted_counts_exact(self):
        namer = RandomNamer(["_A", "_B"], weights=[1, 3], strict=True, seed=0)
        result = [next(namer) for _ in range(40)]
        assert result.count("_A") == 10
        assert result.count("_B") == 30

    def test_each_period_is_shuffled_independently(self):
        namer = RandomNamer(["_A", "_B"], weights=[1, 1], strict=True, seed=0)
        period_a = [next(namer) for _ in range(2)]
        period_b = [next(namer) for _ in range(2)]
        # Both periods have exact composition but order may differ
        assert sorted(period_a) == ["_A", "_B"]
        assert sorted(period_b) == ["_A", "_B"]

    def test_zero_count_raises(self):
        with pytest.raises(ValueError, match="counts must be >= 1"):
            RandomNamer(["_A", "_B"], weights=[0, 1], strict=True)

    def test_strict_seeded_reproducibility(self):
        a = [
            next(RandomNamer(["_A", "_B"], weights=[1, 3], strict=True, seed=7))
            for _ in range(20)
        ]
        b = [
            next(RandomNamer(["_A", "_B"], weights=[1, 3], strict=True, seed=7))
            for _ in range(20)
        ]
        assert a == b

    def test_strict_repr(self):
        assert "strict=True" in repr(RandomNamer(["_A", "_B"], strict=True))


class TestMarkovNamer:
    def test_draws_from_names(self):
        matrix = [[0.5, 0.5], [0.5, 0.5]]
        namer = MarkovNamer(["_A", "_B"], matrix, seed=0)
        result = {next(namer) for _ in range(50)}
        assert result <= {"_A", "_B"}

    def test_perfect_alternation(self):
        namer = MarkovNamer(["_A", "_B"], [[0, 1], [1, 0]], start="_A", seed=0)
        assert [next(namer) for _ in range(6)] == ["_A", "_B", "_A", "_B", "_A", "_B"]

    def test_absorbing_state(self):
        # Once in _B, stays in _B
        namer = MarkovNamer(["_A", "_B"], [[0, 1], [0, 1]], start="_A", seed=0)
        result = [next(namer) for _ in range(6)]
        assert result[0] == "_A"
        assert all(r == "_B" for r in result[1:])

    def test_row_normalization(self):
        # Unnormalized rows (raw counts) should be normalized automatically
        namer = MarkovNamer(["_A", "_B"], [[2, 2], [1, 3]], start="_A", seed=42)
        result = {next(namer) for _ in range(50)}
        assert result <= {"_A", "_B"}

    def test_start_by_name(self):
        namer = MarkovNamer(["_A", "_B"], [[0, 1], [1, 0]], start="_B", seed=0)
        assert next(namer) == "_B"

    def test_seeded_reproducibility(self):
        matrix = [[0.7, 0.3], [0.4, 0.6]]
        a = [next(MarkovNamer(["_A", "_B"], matrix, seed=7)) for _ in range(20)]
        b = [next(MarkovNamer(["_A", "_B"], matrix, seed=7)) for _ in range(20)]
        assert a == b

    def test_wrong_matrix_shape_raises(self):
        with pytest.raises(ValueError, match="transition_matrix must be"):
            MarkovNamer(["_A", "_B"], [[0.5, 0.3, 0.2], [0.5, 0.5, 0.0]])

    def test_three_state(self):
        matrix = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]  # A→B→C→A→...
        namer = MarkovNamer(["_A", "_B", "_C"], matrix, start="_A", seed=0)
        assert [next(namer) for _ in range(6)] == ["_A", "_B", "_C", "_A", "_B", "_C"]

    def test_repr(self):
        assert "MarkovNamer" in repr(
            MarkovNamer(["_A", "_B"], [[0.5, 0.5], [0.5, 0.5]])
        )


class TestGradientNamer:
    def test_draws_from_names(self):
        namer = GradientNamer(["_A", "_B"], [1, 0], [0, 1], n_steps=10, seed=0)
        result = {next(namer) for _ in range(10)}
        assert result <= {"_A", "_B"}

    def test_pure_start(self):
        # weights=[1, 0] at step 0 → must be _A
        namer = GradientNamer(["_A", "_B"], [1, 0], [0, 1], n_steps=10, seed=0)
        assert next(namer) == "_A"

    def test_pure_end(self):
        # weights=[0, 1] at step n_steps-1 → must be _B
        namer = GradientNamer(["_A", "_B"], [1, 0], [0, 1], n_steps=5, seed=0)
        result = [next(namer) for _ in range(5)]
        assert result[-1] == "_B"

    def test_composition_shifts(self):
        # Over many independent runs: early beads mostly _A, late mostly _B
        n_steps, n_trials = 100, 300
        first, last = [], []
        for i in range(n_trials):
            namer = GradientNamer(["_A", "_B"], [1, 0], [0, 1], n_steps=n_steps, seed=i)
            result = [next(namer) for _ in range(n_steps)]
            first.append(result[0])
            last.append(result[-1])
        assert first.count("_A") > n_trials * 0.95
        assert last.count("_B") > n_trials * 0.95

    def test_continues_past_n_steps(self):
        namer = GradientNamer(["_A", "_B"], [1, 0], [0, 1], n_steps=5, seed=0)
        for _ in range(5):
            next(namer)
        for _ in range(10):
            assert next(namer) in {"_A", "_B"}

    def test_seeded_reproducibility(self):
        a = [
            next(GradientNamer(["_A", "_B"], [1, 0], [0, 1], n_steps=10, seed=5))
            for _ in range(10)
        ]
        b = [
            next(GradientNamer(["_A", "_B"], [1, 0], [0, 1], n_steps=10, seed=5))
            for _ in range(10)
        ]
        assert a == b

    def test_repr(self):
        assert "GradientNamer" in repr(
            GradientNamer(["_A", "_B"], [1, 0], [0, 1], n_steps=10)
        )
