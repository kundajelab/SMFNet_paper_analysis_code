#Run from SMFNet_paper_analysis_code github repo: python ./Processing_code/state_calling/call_states.py

import os
from pathlib import Path
import tempfile
import time

import h5py
import numpy as np
from numba import njit, prange, set_num_threads, config

#Set the input path and output prefix before running.
H5_IN = None
OUT_DIR = Path("output/state_calling")
OUT_STR = "sm"
OC_BORDER_VAL = 50
AC_BORDER_VAL = 20
OC_THRESH = 70.0
AC_THRESH = 80.0
WORKERS = min(8, int(os.getenv("SLURM_CPUS_PER_TASK", "1")))
CHUNK_READS = 2048
BLOCK_MB = 256
MATRIX_ORDER = "auto"
SUBGROUP = "predicted"
MAX_STRETCHES = 1000000
SEED = 0
MAKE_PLOTS = True
OVERWRITE = False


@njit(parallel=True, cache=True)
def stretches(matrix, positions, border, accessible=False):
    #Each row is a molecule; columns correspond to the measured positions.
    n_reads, n_positions = matrix.shape
    stretch_matrix = matrix.astype(np.float64)
    stretch_lengths = np.zeros((n_reads, n_positions), dtype=np.float64)
    for read_index in prange(n_reads):
        valid_columns = np.empty(n_positions, dtype=np.int64)
        n_valid = 0
        for column_index in range(n_positions):
            if matrix[read_index, column_index] != -1:
                valid_columns[n_valid] = column_index
                n_valid += 1
        #At least three measured sites are needed to estimate a stretch.
        if n_valid < 3:
            continue
        run_start = 0
        while run_start < n_valid:
            value = matrix[read_index, valid_columns[run_start]]
            #Intermediate codes distinguish short footprints from nucleosomes.
            if accessible:
                if value == -2:
                    value = 0
                elif value == -3 or value == -5:
                    value = 1
            run_end = run_start
            while run_end + 1 < n_valid:
                next_value = matrix[read_index, valid_columns[run_end + 1]]
                if accessible:
                    if next_value == -2:
                        next_value = 0
                    elif next_value == -3 or next_value == -5:
                        next_value = 1
                if next_value != value:
                    break
                run_end += 1
            #Extend a run to the midpoint of its flanking measured sites.
            if value == (0 if accessible else 1):
                first_column, last_column = (
                    valid_columns[run_start],
                    valid_columns[run_end],
                )
                left = (
                    positions[first_column]
                    if first_column == 0
                    else positions[first_column]
                    - (positions[first_column] - positions[first_column - 1]) / 2
                )
                right = (
                    positions[last_column]
                    if last_column == n_positions - 1
                    else positions[last_column]
                    + (positions[last_column + 1] - positions[last_column]) / 2
                )
                length = right - left + 1
                touches_missing_or_edge = (
                    first_column == 0
                    or matrix[read_index, first_column - 1] == -1
                    or last_column == n_positions - 1
                    or matrix[read_index, last_column + 1] == -1
                )
                trimmed = length < border and touches_missing_or_edge
                if not trimmed:
                    stretch_lengths[read_index, first_column] = length
                replacement = (-1 if accessible else -5) if trimmed else length
                for column_index in range(first_column, last_column + 1):
                    if not accessible or matrix[read_index, column_index] == 0:
                        stretch_matrix[read_index, column_index] = replacement
            run_start = run_end + 1
    flat_lengths = stretch_lengths.ravel()
    return stretch_matrix, flat_lengths[flat_lengths > 0]


def occupied_classes(matrix, threshold):
    positive = matrix > 0
    matrix[positive & (matrix < threshold)] = -2
    matrix[positive & (matrix >= threshold)] = -3
    return matrix


def final_classes(matrix, threshold):
    positive = matrix > 0
    matrix[positive & (matrix < threshold)] = 1
    matrix[positive & (matrix >= threshold)] = 0
    matrix[matrix == -2] = 2
    matrix[matrix == -3] = 3
    matrix[matrix == -5] = -1
    return matrix.astype(np.int8)


def regions(h5, subgroup):
    if "smf_mat" in h5 and "smf_pos" in h5:
        for name in h5["smf_mat"]:
            read_ids = h5["read_id"].get(name) if "read_id" in h5 else None
            yield name, h5["smf_mat"][name], h5["smf_pos"][name], read_ids
    elif "region_ids" in h5 and "regions" in h5:
        seen = set()
        for region_index, name in enumerate(h5["region_ids"].asstr()):
            if name in seen:
                raise ValueError(f"Duplicate region ID: {name}")
            seen.add(name)
            group = h5["regions"][f"{region_index:08d}"][subgroup]
            yield name, group["smf"], group["positions"], group.get("read_id")
    else:
        raise ValueError(
            "Expected smf_mat/smf_pos or region_ids/regions/<index>/predicted"
        )


def blocks(h5):
    for name, dataset, positions, read_ids in regions(h5, SUBGROUP):
        if "/" in name:
            raise ValueError(f"Region ID contains an HDF5 path separator: {name}")
        positions = np.asarray(positions, dtype=np.float64)
        if (
            positions.ndim != 1
            or not np.isfinite(positions).all()
            or np.any(np.diff(positions) <= 0)
        ):
            raise ValueError(f"{name}: positions must be finite, strictly increasing")
        if dataset.ndim != 2:
            raise ValueError(f"{name}: matrix must be two dimensional")
        transpose = MATRIX_ORDER == "positions-by-reads"
        if MATRIX_ORDER == "auto":
            transpose = (
                dataset.shape[1] != positions.size
                and dataset.shape[0] == positions.size
            )
        shape = dataset.shape[::-1] if transpose else dataset.shape
        if shape[1] != positions.size:
            raise ValueError(
                f"{name}: matrix {shape} does not match {positions.size} positions"
            )
        #Cap intermediate arrays to roughly the selected memory budget.
        block_size = min(
            CHUNK_READS,
            max(1, BLOCK_MB * 1024**2 // max(40 * positions.size, 1)),
        )
        if shape[0] == 0:
            yield name, shape, positions, read_ids, 0, np.empty(shape, dtype=np.int8)
        for start in range(0, shape[0], block_size):
            matrix = (
                dataset[:, start : start + block_size].T
                if transpose
                else dataset[start : start + block_size]
            )
            if not np.isin(matrix, [-1, 0, 1]).all():
                raise ValueError(
                    f"{name}: expected hard calls -1/0/1, found other values"
                )
            yield name, shape, positions, read_ids, start, np.ascontiguousarray(
                matrix, dtype=np.int8
            )


class Sample:

    def __init__(self, maximum, seed):
        self.maximum = maximum
        self.rng = np.random.default_rng(seed)
        self.values = np.empty(0)
        self.keys = np.empty(0)

    def add(self, values):
        if not values.size:
            return
        keys = self.rng.random(values.size)
        if self.keys.size == self.maximum:
            keep = keys < self.keys.max()
            keys, values = keys[keep], values[keep]
        self.keys = np.concatenate((self.keys, keys))
        self.values = np.concatenate((self.values, values))
        if self.values.size > self.maximum:
            take = np.argpartition(self.keys, self.maximum - 1)[: self.maximum]
            self.keys, self.values = self.keys[take], self.values[take]


def find_threshold(values):
    from scipy.stats import gaussian_kde

    if not values.size:
        raise ValueError(
            "No stretches for automatic threshold; supply a fixed threshold"
        )
    standard_deviation = np.std(values, ddof=1) if values.size > 1 else 0
    interquartile_range = np.subtract(*np.percentile(values, [75, 25]))
    bandwidth = (
        0.9
        * min(standard_deviation, interquartile_range / 1.34)
        * values.size ** (-1 / 5)
    )
    if not np.isfinite(bandwidth) or bandwidth <= 0:
        return float(np.median(values))
    kde = gaussian_kde(values)
    kde.set_bandwidth(bw_method=bandwidth / standard_deviation)
    grid = np.linspace(values.min(), values.max(), 512)
    density = kde(grid)
    turning_points = np.diff(np.sign(np.diff(density)))
    peaks, valleys = (
        np.where(turning_points == -2)[0] + 1,
        np.where(turning_points == 2)[0] + 1,
    )
    if len(peaks) < 2 or not len(valleys):
        return float(np.median(values))
    merged = [peaks[0]]
    for peak_index in peaks[1:]:
        if grid[peak_index] - grid[merged[-1]] > 10:
            merged.append(peak_index)
    peaks = np.array(merged)
    if len(peaks) < 2:
        return float(np.median(values))
    first_peak, second_peak = np.sort(peaks[np.argsort(density[peaks])[-2:]])
    valleys_between_peaks = valleys[(valleys > first_peak) & (valleys < second_peak)]
    return (
        float(grid[valleys_between_peaks[np.argmin(density[valleys_between_peaks])]])
        if valleys_between_peaks.size
        else float((grid[first_peak] + grid[second_peak]) / 2)
    )


def plot_density(values, threshold, path, label):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set(style="whitegrid")
    fig, ax = plt.subplots(figsize=(8, 4))
    if values.size > 1 and np.ptp(values) > 0:
        sns.kdeplot(values, ax=ax, color="green")
    ax.axvline(threshold, color="red", linestyle="--")
    ax.set(
        title=f"{label} - threshold: {threshold:.2f}",
        xlabel="Stretch length (bp)",
        ylabel="Density",
    )
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def threshold_value(value):
    if value is None or (isinstance(value, str) and value.lower() == "auto"):
        return None
    result = float(value)
    if not np.isfinite(result) or result < 0:
        raise ValueError("Threshold must be finite and nonnegative, or auto")
    return result


def main():
    if H5_IN is None:
        raise ValueError("Set H5_IN to the input reads H5 path")
    if MATRIX_ORDER not in {"auto", "reads-by-positions", "positions-by-reads"}:
        raise ValueError("Invalid MATRIX_ORDER")
    if min(WORKERS, CHUNK_READS, BLOCK_MB, MAX_STRETCHES) < 1:
        raise ValueError("workers, chunk sizes and max-stretches must be positive")
    if not all(
        np.isfinite(x) and x >= 0 for x in [OC_BORDER_VAL, AC_BORDER_VAL]
    ):
        raise ValueError("border thresholds must be finite and nonnegative")
    if Path(OUT_STR).name != OUT_STR or OUT_STR in ("", ".", ".."):
        raise ValueError("OUT_STR must be a filename prefix")
    set_num_threads(min(WORKERS, config.NUMBA_NUM_THREADS))
    output = Path(OUT_DIR)
    output.mkdir(parents=True, exist_ok=True)
    target = output / f"{OUT_STR}_smf_states.h5"
    if target.resolve() == Path(H5_IN).resolve():
        raise ValueError("output must differ from input")
    if target.exists() and not OVERWRITE:
        raise ValueError(f"{target} exists; set OVERWRITE = True")
    occupied_threshold = threshold_value(OC_THRESH)
    accessible_threshold = threshold_value(AC_THRESH)
    started = time.perf_counter()
    occupied_sample, accessible_sample = Sample(MAX_STRETCHES, SEED), Sample(
        MAX_STRETCHES, SEED + 1
    )
    #Automatic thresholds need a preliminary pass over sampled run lengths.
    with h5py.File(H5_IN, "r") as source:
        if occupied_threshold is None:
            for _, _, site_positions, _, _, matrix in blocks(source):
                _, lengths = stretches(matrix, site_positions, OC_BORDER_VAL)
                occupied_sample.add(lengths)
            occupied_threshold = find_threshold(occupied_sample.values)
        if accessible_threshold is None:
            for _, _, site_positions, _, _, matrix in blocks(source):
                occupied, _ = stretches(matrix, site_positions, OC_BORDER_VAL)
                _, lengths = stretches(
                    occupied_classes(occupied, occupied_threshold),
                    site_positions,
                    AC_BORDER_VAL,
                    True,
                )
                accessible_sample.add(lengths)
            accessible_threshold = find_threshold(accessible_sample.values)
        #Write to a temporary file so failed processing cannot publish partial output.
        fd, temporary = tempfile.mkstemp(
            prefix=target.name + ".", suffix=".partial", dir=output
        )
        os.close(fd)
        try:
            with h5py.File(temporary, "w") as destination:
                state_group, position_group = destination.create_group(
                    "smf_states"
                ), destination.create_group("smf_pos")
                destination.attrs.update(
                    oc_threshold=occupied_threshold,
                    ac_threshold=accessible_threshold,
                    oc_border=OC_BORDER_VAL,
                    ac_border=AC_BORDER_VAL,
                    matrix_order="reads_by_positions",
                    seed=SEED,
                )
                for name, shape, site_positions, read_ids, start, matrix in blocks(
                    source
                ):
                    occupied, lengths = stretches(
                        matrix, site_positions, OC_BORDER_VAL
                    )
                    if threshold_value(OC_THRESH) is not None and MAKE_PLOTS:
                        occupied_sample.add(lengths)
                    accessible, lengths = stretches(
                        occupied_classes(occupied, occupied_threshold),
                        site_positions,
                        AC_BORDER_VAL,
                        True,
                    )
                    if threshold_value(AC_THRESH) is not None and MAKE_PLOTS:
                        accessible_sample.add(lengths)
                    called = final_classes(accessible, accessible_threshold)
                    if name not in state_group:
                        kwargs = (
                            dict(compression="gzip", compression_opts=1, shuffle=True)
                            if all(shape)
                            else {}
                        )
                        state_group.create_dataset(
                            name, shape=shape, dtype="i1", **kwargs
                        )
                        position_group.create_dataset(name, data=site_positions)
                    state_group[name][start : start + len(matrix)] = called
            os.replace(temporary, target)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
    if MAKE_PLOTS:
        plot_density(
            occupied_sample.values,
            occupied_threshold,
            output / f"{OUT_STR}_states_oc_density.pdf",
            "Occupied",
        )
        plot_density(
            accessible_sample.values,
            accessible_threshold,
            output / f"{OUT_STR}_states_ac_density.pdf",
            "Accessible",
        )
    print(
        f"Wrote {target}; OC={occupied_threshold:g}, AC={accessible_threshold:g}; {time.perf_counter()-started:.2f}s",
        flush=True,
    )


if __name__ == "__main__":
    main()
