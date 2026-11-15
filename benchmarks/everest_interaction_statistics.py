"""Per-primary Everest observables for paired reference regression tests.

CPU collection only. Rows in the output arrays are indexed by interaction
code; columns by original primary particle_id (0..n-1).
"""
import numpy as np

# Stable InteractionRecord codes, part of the snapshot schema.
CONTINUOUS = frozenset((13, 100, 105))
RECORD_COLUMNS = ['particle_id_before', 'particle_id_after', 's_before', 's_after']


def interaction_snapshot(record, n, length, names):
    """Aggregate recorded events without retaining the entire record.

    names maps interaction codes to the existing human-readable labels.
    length is the element length, used only for a numerical roundoff tolerance.
    """
    nr = int(record._index.num_recorded)
    code = np.asarray(record._inter[:nr], dtype=np.int64)
    pid = np.asarray(record.particle_id_before[:nr], dtype=np.int64)
    if n < 2 or not np.isfinite(length) or length <= 0:
        raise ValueError('Expected n >= 2 and a positive finite element length')
    if np.any((pid < 0) | (pid >= n)):
        raise ValueError('Unexpected particle IDs in standalone Everest record')
    if np.any(code == 0):
        raise ValueError('Uninitialised interaction in recorded rows')
    codes = np.unique(code)
    counts = np.zeros((len(codes), n), dtype=np.uint32)
    length_codes = np.array(sorted(CONTINUOUS.intersection(codes)), dtype=np.int64)
    lengths = np.zeros((len(length_codes), n), dtype=np.float64)
    for i, k in enumerate(codes):
        mask = code == k
        cnt = np.bincount(pid[mask], minlength=n)
        if cnt.max(initial=0) > np.iinfo(np.uint32).max:
            raise OverflowError('Per-primary interaction count exceeds uint32')
        counts[i] = cnt
    for i, k in enumerate(length_codes):
        mask = code == k
        before = np.asarray(record.s_before[:nr])[mask]
        after = np.asarray(record.s_after[:nr])[mask]
        child = np.asarray(record.particle_id_after[:nr])[mask]
        if np.any(child != pid[mask]):
            raise ValueError(f'Incomplete or unexpected child for continuous event {k}')
        ds = after - before
        if (not np.all(np.isfinite(ds)) or np.any(ds < -1e-12 * length)):
            raise ValueError(f'Invalid longitudinal extent for continuous event {k}')
        lengths[i] = np.bincount(pid[mask], weights=np.maximum(ds, 0), minlength=n)
    return dict(event_code=codes,
                event_name=np.array([names[int(k)] for k in codes], dtype=str),
                event_n=counts.sum(axis=1, dtype=np.int64),
                event_count_by_particle=counts,
                event_length_code=length_codes,
                event_length_by_particle=lengths)


def read_interactions(snapshot, n):
    """Read and validate arrays; absence of a code means an all-zero row."""
    codes = snapshot['event_code']
    counts = snapshot['event_count_by_particle']
    length_codes = snapshot['event_length_code']
    lengths = snapshot['event_length_by_particle']
    names = snapshot['event_name']
    if (codes.ndim != 1 or length_codes.ndim != 1
            or names.shape != codes.shape or len(set(codes)) != len(codes)
            or len(set(length_codes)) != len(length_codes)
            or counts.shape != (len(codes), n)
            or lengths.shape != (len(length_codes), n)
            or not np.issubdtype(counts.dtype, np.integer)
            or np.any(counts < 0) or np.any(codes == 0)
            or not np.array_equal(counts.sum(axis=1, dtype=np.int64), snapshot['event_n'])
            or set(length_codes) != CONTINUOUS.intersection(codes)
            or not np.all(np.isfinite(lengths)) or np.any(lengths < 0)):
        raise ValueError('Invalid per-primary interaction snapshot')
    return (dict(zip(codes, counts)), dict(zip(length_codes, lengths)),
            dict(zip(codes, names)))


def paired_statistics(ref, new, *, relative_tol=0.05, sigma=5.0):
    """Mean difference and paired standard error over incident primaries.

    The tolerance is an engineering screening rule, not a calibrated
    hypothesis test. The arrays include zero-valued non-interacting primaries.
    """
    ref, new = np.asarray(ref, dtype=float), np.asarray(new, dtype=float)
    if ref.shape != new.shape or ref.ndim != 1 or len(ref) < 2:
        raise ValueError('Expected aligned per-primary vectors of length >= 2')
    if not np.all(np.isfinite(ref)) or not np.all(np.isfinite(new)):
        raise ValueError('Non-finite per-primary observable')
    d = new - ref
    stderr = d.std(ddof=1) / np.sqrt(len(d))
    mean_ref, mean_new = ref.mean(), new.mean()
    difference = abs(d.mean())
    return dict(ref=mean_ref, new=mean_new, difference=difference,
                stderr=stderr, limit=max(relative_tol * abs(mean_ref), sigma * stderr),
                affected_ref=int(np.count_nonzero(ref)),
                affected_new=int(np.count_nonzero(new)))


def compare_interactions(a, b, n, check, *, relative_tol=0.05, sigma=5.0):
    """Call the caller's existing check(label, value, limit) for observables."""
    ca, la, na = read_interactions(a, n)
    cb, lb, nb = read_interactions(b, n)
    zero = np.zeros(n)

    def compare_one(label, x, y):
        result = paired_statistics(x, y, relative_tol=relative_tol, sigma=sigma)
        print(f'    {label}: mean {result["ref"]:.6g} / {result["new"]:.6g}; '
              f'paired SE {result["stderr"]:.6g}; '
              f'nonzero primaries {result["affected_ref"]} / {result["affected_new"]}')
        check(label, result['difference'], result['limit'])
        if min(result['affected_ref'], result['affected_new']) < 100:
            print('    NOTE: low incidence; limited statistical sensitivity')

    for code in sorted(ca.keys() | cb.keys()):
        if code < 0:  # entry/exit bookkeeping, not scattering physics
            continue
        name = na.get(code, nb.get(code))
        x, y = ca.get(code, zero), cb.get(code, zero)
        if code in CONTINUOUS:
            # Segment count is diagnostic only; compare accumulated distance.
            u, v = la.get(code, zero), lb.get(code, zero)
            compare_one(f'{name} longitudinal extent / primary [m]', u, v)
            compare_one(f'{name} nonzero-length incidence', u > 0, v > 0)
        else:
            compare_one(f'{name} multiplicity / primary', x, y)
            compare_one(f'{name} incidence', x > 0, y > 0)
