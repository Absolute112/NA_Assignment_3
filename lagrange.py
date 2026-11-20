import math
import warnings
import numpy as np


def _numeric_array(values, name):
    try:
        arr = np.asarray(values)
    except Exception as exc:
        raise TypeError(f"{name} must contain numeric values") from exc
    if not np.issubdtype(arr.dtype, np.number):
        raise TypeError(f"{name} must contain numeric values")
    try:
        arr = arr.astype(float)
    except (TypeError, ValueError, OverflowError) as exc:
        raise TypeError(f"{name} must contain numeric values") from exc
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain finite values")
    return arr


def validate_nodes(x_nodes):
    """Validate interpolation nodes and reject repeated or near-repeated nodes.


    Parameters
    ----------
    x_nodes : array-like
        Numeric, finite, one-dimensional interpolation nodes.


    Returns
    -------
    numpy.ndarray
        Validated nodes as a float array.


    Raises
    ------
    TypeError
        If nodes are not numeric.
    ValueError
        If nodes are empty, not one-dimensional, non-finite, or coincident.

    Example
    -------
    >>> validate_nodes([0, 1])
    array([0., 1.])
    """
    nodes = _numeric_array(x_nodes, "x_nodes")
    if nodes.ndim != 1:
        raise ValueError("x_nodes must be one-dimensional")
    if nodes.size == 0:
        raise ValueError("At least one interpolation node is required")
    ordered = np.sort(nodes)
    for j in range(ordered.size - 1):
        a, b = ordered[j], ordered[j + 1]
        if abs(a - b) < 1e-12 * max(1.0, abs(a)):
            raise ValueError(f"node {b} appears twice")
    return nodes


def _restore_scalar(value, scalar_input):
    if scalar_input:
        return float(np.asarray(value))
    return value


def lagrange_basis(x_nodes, k, x_eval):
    """Evaluate the k-th Lagrange basis polynomial at evaluation points.


    Parameters
    ----------
    x_nodes : array-like
        Interpolation nodes.
    k : int
        Zero-based basis index.
    x_eval : number or array-like
        Evaluation point(s).

    Returns
    -------
    float or numpy.ndarray
        Basis values, returned as a float for scalar input.

    Raises
    ------
    TypeError
        If inputs are non-numeric or k is not an integer.
    ValueError
        If nodes are invalid or k is outside the valid range.

    Example
    -------
    >>> lagrange_basis([0, 1], 0, 0.25)
    0.75
    """
    nodes = validate_nodes(x_nodes)
    if isinstance(k, bool) or not isinstance(k, (int, np.integer)):
        raise TypeError("k must be an integer")
    if k < 0 or k >= nodes.size:
        raise ValueError("k must be between 0 and len(x_nodes)-1")
    x = _numeric_array(x_eval, "x_eval")
    scalar_input = x.ndim == 0
    result = np.ones_like(x, dtype=float)
    for i in range(nodes.size):
        if i != k:
            result *= (x - nodes[i]) / (nodes[k] - nodes[i])
    return _restore_scalar(result, scalar_input)


def lagrange_interpolate(x_nodes, y_nodes, x_eval, return_basis=False):
    """Evaluate the Lagrange interpolant and return metadata.


    Parameters
    ----------
    x_nodes, y_nodes : array-like
        Node locations and corresponding data values.
    x_eval : number or array-like
        Evaluation point(s).
    return_basis : bool, optional
        Include evaluated basis values in the returned dictionary if True.

    Returns
    -------
    dict
        Keys are value, degree, basis, and extrapolating. Scalar input produces
        a float value; basis is None unless requested.

    Raises
    ------
    TypeError
        If numeric inputs are invalid or return_basis is not boolean.
    ValueError
        If nodes are invalid or x_nodes and y_nodes have different lengths.

    Example
    -------
    >>> lagrange_interpolate([0, 1], [1, 3], 0.5)['value']
    2.0
    """
    nodes = validate_nodes(x_nodes)
    y = _numeric_array(y_nodes, "y_nodes")
    x = _numeric_array(x_eval, "x_eval")
    if y.ndim != 1:
        raise ValueError("y_nodes must be one-dimensional")
    if y.size != nodes.size:
        raise ValueError("x_nodes and y_nodes must have the same length")
    if not isinstance(return_basis, (bool, np.bool_)):
        raise TypeError("return_basis must be boolean")
    scalar_input = x.ndim == 0
    basis = []
    value = np.zeros_like(x, dtype=float)
    for k in range(nodes.size):
        Lk = np.ones_like(x, dtype=float)
        for i in range(nodes.size):
            if i != k:
                Lk *= (x - nodes[i]) / (nodes[k] - nodes[i])
        basis.append(Lk)
        value += y[k] * Lk
    extrapolating = bool(np.any((x < nodes.min()) | (x > nodes.max())))
    if extrapolating:
        warnings.warn("Evaluation outside the node interval (extrapolation)", RuntimeWarning)
    basis_out = np.asarray(basis) if return_basis else None
    if return_basis and scalar_input:
        basis_out = np.asarray(basis_out, dtype=float).reshape(nodes.size)
    return {
        "value": _restore_scalar(value, scalar_input),
        "degree": int(nodes.size - 1),
        "basis": basis_out,
        "extrapolating": extrapolating,
    }


def interpolation_error_bound(x_nodes, x_eval, deriv_bound):
    """Compute the classical Lagrange interpolation error bound.


    Parameters
    ----------
    x_nodes : array-like
        Interpolation nodes.
    x_eval : number or array-like
        Evaluation point(s).
    deriv_bound : number
        Bound on the absolute value of the (n+1)-th derivative.

    Returns
    -------
    float or numpy.ndarray
        M/(n+1)! times the absolute nodal product, scalar for scalar input.

    Raises
    ------
    TypeError
        If nodes, evaluation points, or derivative bound are not numeric.
    ValueError
        If nodes are invalid or deriv_bound is negative or non-scalar.

    Example
    -------
    >>> interpolation_error_bound([0, 1], 0.5, 2)
    0.125
    """
    nodes = validate_nodes(x_nodes)
    x = _numeric_array(x_eval, "x_eval")
    M = _numeric_array(deriv_bound, "deriv_bound")
    if M.size != 1 or M.ndim != 0:
        raise ValueError("deriv_bound must be a scalar")
    M = float(M)
    if M < 0:
        raise ValueError("deriv_bound must be non-negative")
    scalar_input = x.ndim == 0
    product = np.ones_like(x, dtype=float)
    for node in nodes:
        product *= x - node
    bound = M * np.abs(product) / math.factorial(nodes.size)
    return _restore_scalar(bound, scalar_input)


def _run_problem1_tests():
    """Run the Problem 1 verification tests and generate the Runge figure."""
    import matplotlib.pyplot as plt
    import os


    nodes4 = np.array([0.0, 1.5, 2.0, 4.2])
    matrix = np.column_stack([lagrange_basis(nodes4, k, nodes4) for k in range(4)])
    deviation = np.max(np.abs(matrix - np.eye(4)))
    print("Kronecker matrix:")
    print(matrix)
    print("Maximum deviation from identity:", deviation)
    assert np.allclose(matrix, np.eye(4), atol=1e-12)


    grid = np.linspace(0.0, 4.2, 501)
    basis_sum = sum(lagrange_basis(nodes4, k, grid) for k in range(4))
    basis_sum_deviation = np.max(np.abs(basis_sum - 1.0))
    print("Maximum basis-sum deviation from 1:", basis_sum_deviation)
    assert np.allclose(basis_sum, 1.0, atol=1e-12)


    f_cubic = lambda x: 2*x**3 - x + 5
    cubic_nodes = np.array([-1.0, -0.2, 0.4, 1.0])
    cubic_grid = np.linspace(-1.0, 1.0, 2001)
    cubic_result = lagrange_interpolate(cubic_nodes, f_cubic(cubic_nodes), cubic_grid)
    cubic_error = np.max(np.abs(cubic_result["value"] - f_cubic(cubic_grid)))
    print("Cubic exactness maximum error:", cubic_error)
    assert cubic_error < 1e-12


    f_runge = lambda x: 1.0 / (1.0 + 25.0*x**2)
    eval_grid = np.linspace(-1.0, 1.0, 2001)
    degrees = [5, 10, 15, 20]
    print("\nRunge maximum errors (2001-point grid)")
    print(" n | equally spaced | Chebyshev")
    print("---|---------------|----------")
    for n in degrees:
        equal_nodes = np.linspace(-1.0, 1.0, n + 1)
        i = np.arange(n + 1)
        cheb_nodes = np.cos((2*i + 1)*np.pi/(2*(n + 1)))
        p_equal = lagrange_interpolate(equal_nodes, f_runge(equal_nodes), eval_grid)["value"]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            p_cheb = lagrange_interpolate(cheb_nodes, f_runge(cheb_nodes), eval_grid)["value"]
        err_equal = np.max(np.abs(p_equal - f_runge(eval_grid)))
        err_cheb = np.max(np.abs(p_cheb - f_runge(eval_grid)))
        print(f"{n:2d} | {err_equal:.8e} | {err_cheb:.8e}")


    n = 10
    equal_nodes = np.linspace(-1.0, 1.0, n + 1)
    i = np.arange(n + 1)
    cheb_nodes = np.cos((2*i + 1)*np.pi/(2*(n + 1)))
    p_equal = lagrange_interpolate(equal_nodes, f_runge(equal_nodes), eval_grid)["value"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        p_cheb = lagrange_interpolate(cheb_nodes, f_runge(cheb_nodes), eval_grid)["value"]
    plt.figure(figsize=(9, 5.5))
    plt.plot(eval_grid, f_runge(eval_grid), label="Runge function")
    plt.plot(eval_grid, p_equal, label="Equally spaced, n=10")
    plt.plot(eval_grid, p_cheb, label="Chebyshev, n=10")
    plt.scatter(equal_nodes, f_runge(equal_nodes), marker="o", label="Equally spaced nodes")
    plt.scatter(cheb_nodes, f_runge(cheb_nodes), marker="x", label="Chebyshev nodes")
    plt.xlabel("x")
    plt.ylabel("f(x) / interpolant")
    plt.title("Runge Interpolation (n=10)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    os.makedirs("figures", exist_ok=True)
    plt.savefig(os.path.join("figures", "runge_n10.png"), dpi=200)
    plt.close()
    print("Saved plot to figures/runge_n10.png")


if __name__ == "__main__":
    _run_problem1_tests()
