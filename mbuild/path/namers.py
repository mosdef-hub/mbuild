"""Bead name generators for use with hard_sphere_random_walk."""

import itertools
from abc import ABC, abstractmethod

import numpy as np


class BeadNamer(ABC):
    """Abstract base class for bead name generators.

    Subclasses must implement ``__next__``.

    Instances are iterators:
    ``next(namer)`` returns the name for the next accepted bead.

    Pass a plain string anywhere a ``BeadNamer`` is expected and
    ``BeadNamer.coerce()`` will wrap it in a ``ConstantNamer`` automatically.

    To build a custom namer, subclass ``BeadNamer`` and implement ``__next__``.
    """

    @abstractmethod
    def __next__(self) -> str: ...

    def __iter__(self):
        return self

    @classmethod
    def coerce(cls, value):
        """Return value unchanged if it is a BeadNamer; wrap strings in ConstantNamer."""
        if isinstance(value, BeadNamer):
            return value
        if isinstance(value, str):
            return ConstantNamer(value)
        raise TypeError(
            f"bead_name must be a str or BeadNamer, got {type(value).__name__!r}"
        )


class ConstantNamer(BeadNamer):
    """Always yields the same bead name."""

    def __init__(self, name: str):
        self.name = name

    def __next__(self) -> str:
        return self.name

    def __repr__(self):
        return f"ConstantNamer({self.name!r})"


class RandomNamer(BeadNamer):
    """Randomly draws bead names from a list, with optional weights.

    Parameters
    ----------
    names : list of str
        Pool of bead names to draw from.
    weights : list of float or int, optional, default None
        When ``strict=False``: sampling probabilities, normalized automatically.
        When ``strict=True``: integer counts per period (e.g. ``[1, 3]`` gives
        one ``_A`` and three ``_B`` per four-bead shuffle). Rounded to the
        nearest integer. Defaults to 1 of each name if not provided.
    strict : bool, default False
        If False, samples with replacement — composition is approximate.
        If True, samples without replacement from a shuffled pool that is
        rebuilt each period — composition is exact over every period.
    seed : int, optional, default None
        Random seed for reproducibility.
    """

    def __init__(self, names, weights=None, strict=False, seed=None):
        self.names = list(names)
        self.strict = strict
        self._rng = np.random.default_rng(seed)

        if strict:
            if weights is not None:
                counts = np.round(np.asarray(weights, dtype=float)).astype(int)
                if np.any(counts <= 0):
                    raise ValueError(
                        "All counts must be >= 1 when strict=True. "
                        "Pass weights as integer counts, e.g. weights=[1, 3]."
                    )
            else:
                counts = np.ones(len(self.names), dtype=int)
            self._pool_template = []
            for name, count in zip(self.names, counts):
                self._pool_template.extend([name] * int(count))
            self._pool = []
        else:
            if weights is not None:
                w = np.asarray(weights, dtype=float)
                self._weights = (w / w.sum()).tolist()
            else:
                self._weights = None

    def __next__(self) -> str:
        if self.strict:
            if not self._pool:
                self._pool = list(self._pool_template)
                self._rng.shuffle(self._pool)
            return self._pool.pop()
        idx = int(self._rng.choice(len(self.names), p=self._weights))
        return self.names[idx]

    def __repr__(self):
        if self.strict:
            return f"RandomNamer({self.names!r}, strict=True)"
        return f"RandomNamer({self.names!r}, weights={self._weights})"


def _flatten_sequence(sequence):
    """Expand a mixed list of names / (name, count) tuples into a flat list."""
    result = []
    for item in sequence:
        if isinstance(item, tuple):
            name, count = item
            result.extend([name] * int(count))
        else:
            result.append(item)
    return result


class MarkovNamer(BeadNamer):
    """Generates bead names from a first-order Markov chain.

    Each bead name is drawn based on the transition probabilities from the
    current state.

    Parameters
    ----------
    names : list of str
        Bead types (states in the chain).
    transition_matrix : array-like, shape (N, N)
        Entry ``[i][j]`` is the probability of transitioning from name ``i``
        to name ``j``. Rows are normalized automatically, so raw reactivity
        ratios can be passed directly.
    start : str or int, optional, default None
        Starting state, given as a name string or index. If None, a state is
        chosen uniformly at random.
    seed : int, optional, default None
        Random seed for reproducibility.

    Examples
    --------
    Perfect alternation:

    >>> namer = MarkovNamer(["_A", "_B"], [[0, 1], [1, 0]], start="_A", seed=0)
    >>> [next(namer) for _ in range(6)]
    ['_A', '_B', '_A', '_B', '_A', '_B']

    Blocky character (high self-transition probability):

    >>> namer = MarkovNamer(["_A", "_B"], [[0.9, 0.1], [0.1, 0.9]], start="_A", seed=0)
    """

    def __init__(self, names, transition_matrix, start=None, seed=None):
        self.names = list(names)
        matrix = np.asarray(transition_matrix, dtype=float)
        if matrix.shape != (len(self.names), len(self.names)):
            raise ValueError(
                f"transition_matrix must be ({len(self.names)}, {len(self.names)}), "
                f"got {matrix.shape}"
            )
        row_sums = matrix.sum(axis=1, keepdims=True)
        self._matrix = matrix / row_sums
        self._rng = np.random.default_rng(seed)

        if start is None:
            self._current = int(self._rng.integers(len(self.names)))
        elif isinstance(start, str):
            self._current = self.names.index(start)
        else:
            self._current = int(start)

    def __next__(self) -> str:
        name = self.names[self._current]
        self._current = int(
            self._rng.choice(len(self.names), p=self._matrix[self._current])
        )
        return name

    def __repr__(self):
        return f"MarkovNamer({self.names!r})"


class GradientNamer(BeadNamer):
    """Generates bead names with composition that shifts along the chain.

    Linearly interpolates sampling weights from ``start_weights`` to
    ``end_weights`` over ``n_steps``. After ``n_steps``, holds at
    ``end_weights`` so the walk can run longer without raising.

    Parameters
    ----------
    names : list of str
        Bead types to draw from.
    start_weights : list of float
        Sampling weights at position 0. Normalized automatically.
    end_weights : list of float
        Sampling weights at position ``n_steps - 1``. Normalized automatically.
    n_steps : int
        Number of steps over which the gradient is applied.
    seed : int, optional, default None
        Random seed for reproducibility.

    Examples
    --------
    Pure A at the start, pure B at the end:

    >>> namer = GradientNamer(["_A", "_B"], [1, 0], [0, 1], n_steps=100, seed=0)
    """

    def __init__(self, names, start_weights, end_weights, n_steps, seed=None):
        self.names = list(names)
        start = np.asarray(start_weights, dtype=float)
        end = np.asarray(end_weights, dtype=float)
        self._start = start / start.sum()
        self._end = end / end.sum()
        self._n_steps = int(n_steps)
        self._rng = np.random.default_rng(seed)
        self._step = 0

    def __next__(self) -> str:
        t = min(self._step / max(self._n_steps - 1, 1), 1.0)
        weights = (1 - t) * self._start + t * self._end
        idx = int(self._rng.choice(len(self.names), p=weights))
        self._step += 1
        return self.names[idx]

    def __repr__(self):
        return f"GradientNamer({self.names!r}, n_steps={self._n_steps})"


class CyclicNamer(BeadNamer):
    """Cycles through a sequence of bead names indefinitely.

    Accepts a flat list of names or a list of ``(name, count)`` tuples for
    block patterns.

    Parameters
    ----------
    sequence : list of str or list of (str, int) tuples
        Names to cycle. Tuples expand to blocks, e.g.
        ``[("_A", 5), ("_B", 5)]`` produces five ``_A`` then five ``_B``,
        repeating forever.

    Examples
    --------
    >>> list(itertools.islice(CyclicNamer(["_A", "_B"]), 6))
    ['_A', '_B', '_A', '_B', '_A', '_B']
    >>> list(itertools.islice(CyclicNamer([("_A", 2), ("_B", 3)]), 10))
    ['_A', '_A', '_B', '_B', '_B', '_A', '_A', '_B', '_B', '_B']
    """

    def __init__(self, sequence):
        self._flat = _flatten_sequence(sequence)
        if not self._flat:
            raise ValueError("CyclicNamer sequence must not be empty.")
        self._cycle = itertools.cycle(self._flat)

    def __next__(self) -> str:
        return next(self._cycle)

    def __repr__(self):
        return f"CyclicNamer({self._flat!r})"
