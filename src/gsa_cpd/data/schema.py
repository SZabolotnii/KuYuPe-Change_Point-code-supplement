"""Unified data schema for benchmark datasets."""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Optional


@dataclass
class AnnotatedChangepoint:
    """A known change-point with metadata."""
    index: int
    label: str = ""
    confidence: float = 1.0


@dataclass
class DatasetInfo:
    """Metadata about a benchmark dataset."""
    name: str
    domain: str = ""
    source: str = ""
    n_observations: int = 0
    kurtosis: float = 0.0
    skewness: float = 0.0


@dataclass
class TimeSeriesData:
    """Container for a time series with annotated change-points."""
    values: np.ndarray = field(default_factory=lambda: np.array([]))
    changepoints: List[AnnotatedChangepoint] = field(default_factory=list)
    info: DatasetInfo = field(default_factory=lambda: DatasetInfo(name=""))

    @property
    def n(self) -> int:
        return len(self.values)
