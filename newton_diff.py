import warnings

import numpy as np


def _float_array(values, name, allow_scalar=False):
    """Convert values to a finite float array.

    Parameters: values is array-like, name labels errors, and allow_scalar
    permits a scalar in addition to a one-dimensional array.
    Returns: A NumPy float array.
    Raises: TypeError for nonnumeric data and ValueError for invalid
    shape or nonfinite values.
    Example: _float_array([1, 2], "nodes")
    """
    try:
        array = np.asarray(values, dtype=float)
    except (TypeError, ValueError, OverflowError) as exc:
        raise TypeError(f"{name} must contain real numbers.") from exc

    if array.ndim != 1 and not (allow_scalar and array.ndim == 0):
        shape = "a scalar or a 1-D array" if allow_scalar else "a 1-D array"
        raise ValueError(f"{name} must be {shape}.")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values.")
    return array


def _sorted_data(x_nodes, y_nodes):
    """Sort paired node data by x.

    Parameters: x_nodes and y_nodes are paired 1-D numeric arrays.
    Returns: Two float arrays sorted by x, preserving each (x, y) pair.
    Raises: TypeError for nonnumeric data and ValueError for empty,
    mismatched, nonfinite, or wrongly shaped data.
    Example: _sorted_data([2, 1], [20, 10])
    """
    x = _float_array(x_nodes, "x_nodes")
    y = _float_array(y_nodes, "y_nodes")
    if x.size != y.size or x.size == 0:
        raise ValueError("x_nodes and y_nodes must have the same nonzero length.")
    order = np.argsort(x)
    return x[order], y[order]


def difference_table(y_nodes):
    """Build the full forward-difference table.

    Parameters: y_nodes is a nonempty 1-D array in node order.
    Returns: A list with table[0] = y and table[j][k] = Δ^j y_k.
    Raises: TypeError for nonnumeric data and ValueError for empty,
    nonfinite, or wrongly shaped data.
    Example: difference_table([1, 2, 4])[-1] gives array([1.]).
    """
    y = _float_array(y_nodes, "y_nodes")
    if y.size == 0:
        raise ValueError("y_nodes must contain at least one value.")

    table = [y.copy()]
    while table[-1].size > 1:
        table.append(np.diff(table[-1]))
    return table


def print_difference_table(x_nodes, y_nodes, decimals=6):
    """Print the forward differences in a triangular table.

    Parameters: x_nodes and y_nodes are paired 1-D data. decimals is
    the nonnegative integer number of decimal places to display.
    Returns: None.
    Raises: TypeError for nonnumeric data and ValueError for invalid
    precision, empty, mismatched, nonfinite, or wrongly shaped data.
    Example: print_difference_table([0, 1, 2], [1, 2, 4], decimals=2)
    """
    if not isinstance(decimals, (int, np.integer)) or decimals < 0:
        raise ValueError("decimals must be a nonnegative integer.")

    x, y = _sorted_data(x_nodes, y_nodes)
    table = difference_table(y)
    headings = ["x", "y"] + [f"Delta^{j} y" for j in range(1, len(table))]
    columns = [[f"{value:.{decimals}f}" for value in col] for col in [x] + table]
    widths = [max(len(h), max(map(len, col))) for h, col in zip(headings, columns)]

    print("  ".join(h.rjust(w) for h, w in zip(headings, widths)))
    for i in range(x.size):
        print("  ".join(col[i].rjust(w) for col, w in zip(columns, widths)
                        if i < len(col)))


def _uniform_spacing(x, rtol):
    """Check the gaps between sorted nodes.

    Parameters: x is a sorted, finite 1-D float array. rtol is the
    maximum permitted relative deviation from the mean gap.
    Returns: The mean gap as a float.
    Raises: ValueError for too few nodes, repeated nodes, unrepresentable
    spacing, or gaps outside the tolerance.
    Example: _uniform_spacing(np.array([0., 1., 2.]), 1e-9) gives 1.0.
    """
    if x.size < 2:
        raise ValueError("Newton interpolation requires at least two nodes.")

    gaps = np.diff(x)
    h = float(np.mean(gaps))
    if not np.isfinite(h):
        raise ValueError("Node spacing is too large to represent as a float.")
    if h <= 0:
        raise ValueError("Node spacing is zero. Remove repeated nodes before using Lagrange.")

    deviation = float(np.max(np.abs(gaps - h)))
    relative = deviation / h
    if np.any(gaps <= 0) or relative > rtol:
        raise ValueError(
            f"Largest spacing deviation is {deviation:.6g} "
            f"(relative {relative:.6g}, tolerance {rtol:.6g}). "
            "Use Lagrange interpolation with distinct nodes instead."
        )
    return h


def check_equal_spacing(x_nodes, rtol=1e-9):
    """Return the mean spacing if the sorted nodes are equally spaced.

    Parameters: x_nodes is 1-D numeric data. rtol is a finite,
    nonnegative relative tolerance.
    Returns: The mean spacing h as a float.
    Raises: TypeError for nonnumeric data and ValueError for invalid
    tolerance, invalid nodes, repeated nodes, or unequal spacing.
    Example: check_equal_spacing([2, 0, 1]) gives 1.0.
    """
    tolerance = _float_array(rtol, "rtol", allow_scalar=True)
    if tolerance.ndim != 0 or tolerance < 0:
        raise ValueError("rtol must be a finite, nonnegative scalar.")

    x = np.sort(_float_array(x_nodes, "x_nodes"))
    return _uniform_spacing(x, float(tolerance))


def _newton(x_nodes, y_nodes, x_eval, degree, backward):
    """Evaluate either Newton difference formula.

    Parameters: x_nodes and y_nodes are paired data, x_eval is a scalar
    or 1-D array, and degree is a nonnegative integer or None.
    backward selects the backward formula when True.
    Returns: A result dictionary in the format of the public functions.
    Raises: TypeError for nonnumeric data and ValueError for invalid
    data, spacing, or degree. Warns if any point extrapolates.
    Example: _newton([0, 1], [1, 3], 0.5, None, False)["value"] gives 2.0.
    """
    x, y = _sorted_data(x_nodes, y_nodes)
    h = _uniform_spacing(x, 1e-9)
    evaluation = _float_array(x_eval, "x_eval", allow_scalar=True)
    points = np.atleast_1d(evaluation)
    n = x.size - 1

    if degree is None:
        degree = n
    if (isinstance(degree, (bool, np.bool_))
            or not isinstance(degree, (int, np.integer)) or not 0 <= degree <= n):
        raise ValueError(f"degree must be an integer between 0 and {n}.")
    degree = int(degree)

    if np.any((points < x[0]) | (points > x[-1])):
        warnings.warn("Evaluation extends outside the node interval.", stacklevel=3)

    if backward:
        base = np.clip(np.searchsorted(x, points, side="left"), degree, n)
        sign = 1
        offset_key = "t"
    else:
        base = np.clip(np.searchsorted(x, points, side="right") - 1, 0, n - degree)
        sign = -1
        offset_key = "s"

    table = difference_table(y)
    offset = (points - x[base]) / h
    factor = np.ones_like(points)
    terms = [table[0][base]]
    value = terms[0].copy()

    for j in range(1, degree + 1):
        factor *= (offset + sign * (j - 1)) / j
        index = base - j if backward else base
        term = factor * table[j][index]
        terms.append(term)
        value += term

    estimate = None
    next_order = degree + 1
    if next_order <= n:
        index = base - next_order if backward else base
        available = (index >= 0) & (index < table[next_order].size)

        if np.any(available):
            factor *= (offset + sign * degree) / next_order
            estimate = np.full(points.shape, np.nan)
            estimate[available] = np.abs(
                factor[available] * table[next_order][index[available]]
            )

    if evaluation.ndim == 0:
        return {
            "value": float(value[0]),
            "degree": degree,
            "base_node": float(x[base[0]]),
            offset_key: float(offset[0]),
            "terms": [float(term[0]) for term in terms],
            "error_estimate": None if estimate is None else float(estimate[0]),
        }

    return {
        "value": value,
        "degree": degree,
        "base_node": x[base],
        offset_key: offset,
        "terms": terms,
        "error_estimate": estimate,
    }


def newton_forward(x_nodes, y_nodes, x_eval, degree=None):
    """Interpolate with Newton's forward-difference formula.

    Parameters: x_nodes and y_nodes are paired 1-D numeric data.
    x_eval is a float or 1-D array. degree is an integer from 0 to
    len(x_nodes) - 1, or None to use all nodes.
    Returns: A dict with value, degree, base_node, s, terms, and
    error_estimate. Pointwise fields are floats for scalar input.
    For arrays, terms is a list of arrays and missing estimates are
    NaN when some estimates exist. If all estimates are missing,
    error_estimate is None.
    Raises: TypeError for nonnumeric input and ValueError for invalid
    data, spacing, or degree. Warns when extrapolating.
    Example: newton_forward([0, 1], [1, 3], 0.5)["value"] gives 2.0.
    """
    return _newton(x_nodes, y_nodes, x_eval, degree, backward=False)


def newton_backward(x_nodes, y_nodes, x_eval, degree=None):
    """Interpolate with Newton's backward-difference formula.

    Parameters: x_nodes and y_nodes are paired 1-D numeric data.
    x_eval is a float or 1-D array. degree is a nonnegative integer
    up to len(x_nodes) - 1, or None to use all nodes.
    Returns: A dict as in newton_forward, with t in place of s.
    Raises: TypeError for nonnumeric input and ValueError for invalid
    data, spacing, or degree. Warns when extrapolating.
    Example: newton_backward([0, 1], [1, 3], 0.5)["value"] gives 2.0.
    """
    return _newton(x_nodes, y_nodes, x_eval, degree, backward=True)


if __name__ == "__main__":
    x = np.linspace(0, 1, 6)
    y = np.exp(x)

    print("Forward differences for exp(x)")
    print_difference_table(x, y)

    print("\nMethod     x       Degree      Value          True error    Next term")
    for label, method, point in [
        ("Forward", newton_forward, 0.05),
        ("Backward", newton_backward, 0.95),
    ]:
        for degree in range(1, 5):
            result = method(x, y, point, degree)
            error = abs(np.exp(point) - result["value"])
            print(f"{label:8} {point:6.2f} {degree:7d} "
                  f"{result['value']:15.10f} {error:13.6e} "
                  f"{result['error_estimate']:13.6e}")

    grid = np.linspace(0, 1, 200)
    for degree in (None, 2):
        forward = newton_forward(x, y, grid, degree)
        backward = newton_backward(x, y, grid, degree)
        gap = np.max(np.abs(forward["value"] - backward["value"]))
        print(f"Degree {forward['degree']}: max |forward - backward| = {gap:.6e}")
