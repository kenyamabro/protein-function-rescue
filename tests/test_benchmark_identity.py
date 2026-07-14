from Bio.Align import PairwiseAligner, substitution_matrices

from benchmark.run_benchmark import _identity


def _aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = "local"
    return aligner


def test_pairwise_identity_is_symmetric() -> None:
    # Different lengths/repeats exercise the tied-local-alignment case that
    # previously made the cached matrix depend on dataset ordering.
    a = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"
    b = "QISFVKSHFSRQLQQQLEERLGLI"
    aligner = _aligner()
    assert _identity(aligner, a, b) == _identity(aligner, b, a)
