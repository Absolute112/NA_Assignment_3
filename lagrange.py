import math
import warnings
import numpy as np


def _numeric_array(values, name):
    """Convert input to a numeric NumPy array, raising TypeError for non-numeric data."""
    try:
        arr = np.asarray(values)
    except Exception as exc:
        raise TypeError(f"{name} must contain numeric values") from exc

    if not np.issubdtype(arr.dtype, np.number):
        raise TypeError(f"{name} must contain numeric values")

    arr = arr.astype(float)
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain finite values")
    return arr



def _validate_nodes(x_nodes):
    """Validate numeric, finite, nonempty interpolation nodes and reject duplicates."""
    nodes = _numeric_array(x_nodes, "x_nodes")
    if nodes.ndim != 1:
        raise ValueError("x_nodes must be a one-dimensional array")
    if nodes.size == 0:
        raise ValueError("At least one interpolation node is required")

    sorted_nodes = np.sort(nodes)
    for j in range(len(sorted_nodes) - 1):
        a, b = sorted_nodes[j], sorted_nodes[j + 1]
        tol = 1e-12 * max(1.0, abs(a))
        if abs(a - b) < tol:
            raise ValueError(f"node {b} appears twice")
    return nodes



def lagrange_basis(x_nodes, k, x_eval):
    """Evaluate the k-th Lagrange basis polynomial at x_eval.



    Parameters
    ----------
    x_nodes : array-like of numeric values
        Interpolation nodes.
    k : int
        Index of the basis polynomial (zero-based).
    x_eval : numeric or array-like
        Locations at which to evaluate the basis polynomial.



    Returns
    -------
    dict
        Contains 'value', 'degree', 'basis', and 'extrapolating'.
        'value' is L_k(x_eval); 'basis' contains all basis values evaluated
        at x_eval, with shape (number_of_nodes, ...).



    Raises
    ------
    TypeError
        If nodes or evaluation points are not numeric, or k is not an integer.
    ValueError
        If nodes are empty, non-finite, or coincident.
    IndexError
        If k is outside the valid node-index range.



    Example
    -------
    >>> lagrange_basis([0, 1], 0, 0.25)['value']
    0.75
    """
    nodes = _validate_nodes(x_nodes)
    if isinstance(k, bool) or not isinstance(k, (int, np.integer)):
        raise TypeError("k must be an integer")
    if k < 0 or k >= len(nodes):
        raise IndexError("k is outside the node-index range")

    x = _numeric_array(x_eval, "x_eval")
    all_basis = []
    for j in range(len(nodes)):
        Lj = np.ones_like(x, dtype=float)
        for i in range(len(nodes)):
            if i != j:
                Lj *= (x - nodes[i]) / (nodes[j] - nodes[i])
        all_basis.append(Lj)

    extrap = bool(np.any((x < np.min(nodes)) | (x > np.max(nodes))))
    if extrap:
        warnings.warn("Evaluation outside the node interval (extrapolation)", RuntimeWarning)

    basis_values = np.asarray(all_basis)
    return {
        "value": basis_values[k],
        "degree": len(nodes) - 1,
        "basis": basis_values,
        "extrapolating": extrap,
    }



def lagrange_interpolate(x_nodes, y_nodes, x_eval, return_basis=False):
    """Evaluate the Lagrange interpolating polynomial at x_eval.



    Parameters
    ----------
    x_nodes, y_nodes : array-like of numeric values
        Node locations and corresponding function values.
    x_eval : numeric or array-like
        Locations at which to evaluate the interpolant.
    return_basis : bool, optional
        If True, include evaluated Lagrange basis functions in the result.



    Returns
    -------
    dict
        'value' contains the interpolant, 'degree' its degree (n_nodes - 1),
        'basis' contains basis evaluations when requested (otherwise None), and
        'extrapolating' indicates whether any evaluation is outside the node range.



    Raises
    ------
    TypeError
        If numerical inputs are not numeric or return_basis is not boolean.
    ValueError
        If nodes are empty, invalid, coincident, or y_nodes has wrong length.



    Example
    -------
    >>> lagrange_interpolate([0, 1], [1, 3], 0.5)['value']
    2.0
    """
    nodes = _validate_nodes(x_nodes)
    y = _numeric_array(y_nodes, "y_nodes")
    x = _numeric_array(x_eval, "x_eval")



    if y.ndim != 1:
        raise ValueError("y_nodes must be a one-dimensional array")
    if len(y) != len(nodes):
        raise ValueError("x_nodes and y_nodes must have the same length")
    if not isinstance(return_basis, (bool, np.bool_)):
        raise TypeError("return_basis must be boolean")



    basis = []
    for k in range(len(nodes)):
        Lk = np.ones_like(x, dtype=float)
        for i in range(len(nodes)):
            if i != k:
                Lk *= (x - nodes[i]) / (nodes[k] - nodes[i])
        basis.append(Lk)



    basis_array = np.asarray(basis)
    values = np.zeros_like(x, dtype=float)
    for k in range(len(nodes)):
        values += y[k] * basis_array[k]



    extrap = bool(np.any((x < np.min(nodes)) | (x > np.max(nodes))))
    if extrap:
        warnings.warn("Evaluation outside the node interval (extrapolation)", RuntimeWarning)



    return {
        "value": values,
        "degree": len(nodes) - 1,
        "basis": basis_array if return_basis else None,
        "extrapolating": extrap,
    }



def interpolation_error_bound(x_eval, x_nodes, M):
    """Compute M/(n+1)! * |product(x-x_i)| for the Lagrange error bound.



    Parameters
    ----------
    x_eval : numeric or array-like
        Evaluation points.
    x_nodes : array-like of numeric values
        Interpolation nodes.
    M : nonnegative numeric
        Bound on the absolute value of the (n+1)-st derivative.



    Returns
    -------
    ndarray or float
        The error bound at each evaluation point.
    """
    nodes = _validate_nodes(x_nodes)
    x = _numeric_array(x_eval, "x_eval")
    M_arr = _numeric_array(M, "M")
    if M_arr.size != 1:
        raise ValueError("M must be a scalar")
    M_value = float(M_arr)
    if M_value < 0:
        raise ValueError("M must be non-negative")

    degree = len(nodes) - 1
    product = np.ones_like(x, dtype=float)
    for node in nodes:
        product *= (x - node)
    return (M_value / math.factorial(degree + 1)) * np.abs(product)



def _run_problem1_tests():
    """Run Problem 1 checks, including the Runge comparison and plot."""
    import matplotlib.pyplot as plt



    # Kronecker property: L_k(x_j) = 1 if j=k, otherwise 0.
    nodes4 = np.array([-1.0, -0.25, 0.5, 1.0])
    matrix = np.column_stack([lagrange_basis(nodes4, k, nodes4)["value"] for k in range(4)])
    assert np.allclose(matrix, np.eye(4), atol=1e-12)



    # Partition of unity: sum_k L_k(x) = 1.
    grid = np.linspace(-1, 1, 501)
    basis_sum = sum(lagrange_basis(nodes4, k, grid)["value"] for k in range(4))
    assert np.allclose(basis_sum, 1.0, atol=1e-12)



    # Exactness for a polynomial of degree 3.
    f_cubic = lambda x: 2*x**3 - x + 5
    cubic_nodes = np.array([-1.0, -0.2, 0.4, 1.0])
    cubic_grid = np.linspace(-1, 1, 2001)
    cubic_result = lagrange_interpolate(cubic_nodes, f_cubic(cubic_nodes), cubic_grid)
    cubic_error = np.max(np.abs(cubic_result["value"] - f_cubic(cubic_grid)))
    print("Cubic exactness maximum error:", cubic_error)
    assert cubic_error < 1e-12



    # Runge function comparison for degrees n=5, 10, 15, 20.
    f_runge = lambda x: 1.0 / (1.0 + 25.0*x**2)
    eval_grid = np.linspace(-1, 1, 2001)
    degrees = [5, 10, 15, 20]
    errors_equal = []
    errors_cheb = []
    print("\nRunge max errors (2001-point grid)")
    print(" n | equally spaced | Chebyshev")
    print("---|---------------|----------")
    for n in degrees:
        equal_nodes = np.linspace(-1, 1, n + 1)
        i = np.arange(n + 1)
        cheb_nodes = np.cos((2*i + 1)*np.pi/(2*(n + 1)))
        p_equal = lagrange_interpolate(equal_nodes, f_runge(equal_nodes), eval_grid)["value"]
        p_cheb = lagrange_interpolate(cheb_nodes, f_runge(cheb_nodes), eval_grid)["value"]
        err_equal = np.max(np.abs(p_equal - f_runge(eval_grid)))
        err_cheb = np.max(np.abs(p_cheb - f_runge(eval_grid)))
        errors_equal.append(err_equal)
        errors_cheb.append(err_cheb)
        print(f"{n:2d} | {err_equal:.8e} | {err_cheb:.8e}")



    # Plot degree n=10, with both node sets marked.
    n = 10
    equal_nodes = np.linspace(-1, 1, n + 1)
    i = np.arange(n + 1)
    cheb_nodes = np.cos((2*i + 1)*np.pi/(2*(n + 1)))
    p_equal = lagrange_interpolate(equal_nodes, f_runge(equal_nodes), eval_grid)["value"]
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
    plt.savefig("runge_n10.png", dpi=200)
    plt.show()



if __name__ == "__main__":
    _run_problem1_tests()
