"""Hermite interpolation using repeated nodes."""
import warnings
import numpy as np
def hermite_coefficients(x_nodes, y_nodes, dydx_nodes):
    """Build the Hermite divided-difference table.
    Parameters: x_nodes, y_nodes, dydx_nodes are 1-D lists of x, y, and slope values.
    Returns: dict with doubled nodes 'z', diagonal 'coeffs', and full 'table'.
    Raises: TypeError for non-numbers; ValueError for bad lengths, non-finite
        values, or repeated x nodes.
    Example: hermite_coefficients([0, 1], [1, 2], [1, 1])['coeffs']
        gives array([1., 1., 0., 0.]).
    """
    try:
        x = np.asarray(x_nodes, dtype=float)
        y = np.asarray(y_nodes, dtype=float)
        dy = np.asarray(dydx_nodes, dtype=float)
    except (TypeError, ValueError, OverflowError):
        raise TypeError("x, y and dydx must be numbers")
    if x.ndim != 1 or y.ndim != 1 or dy.ndim != 1: raise ValueError("Inputs must be 1-D")
    if len(x) == 0: raise ValueError("At least one node is needed")
    if len(x) != len(y) or len(x) != len(dy): raise ValueError("Inputs must have the same length")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)) or not np.all(np.isfinite(dy)):
        raise ValueError("All values must be finite")
    sorted_x = np.sort(x)
    for i in range(1, len(x)):
        if abs(sorted_x[i] - sorted_x[i - 1]) <= 1e-12 * max(1, abs(sorted_x[i - 1])):
            raise ValueError("node " + str(sorted_x[i]) + " appears twice (or is too close)")
    z = np.repeat(x, 2)
    m = len(z)
    q = np.zeros((m, m))
    q[:, 0] = np.repeat(y, 2)
    for i in range(len(x)):
        q[2 * i + 1, 1] = dy[i]
        if i > 0: q[2 * i, 1] = (y[i] - y[i - 1]) / (x[i] - x[i - 1])
    for j in range(2, m):
        for i in range(j, m):
            q[i, j] = (q[i, j - 1] - q[i - 1, j - 1]) / (z[i] - z[i - j])
    return {"z": z, "coeffs": np.diag(q).copy(), "table": q}
def hermite_interpolate(x_nodes, y_nodes, dydx_nodes, x_eval, return_derivative=False):
    """Evaluate the Hermite polynomial and optionally its derivative.
    Parameters: x_nodes, y_nodes, dydx_nodes are 1-D lists of x, y, and slopes.
        x_eval is a number or 1-D array. return_derivative asks for H'(x).
    Returns: dict with 'value', 'derivative' (or None), 'degree', and
        'extrapolating'. A single x_eval gives float results.
    Raises: TypeError for non-numbers; ValueError for invalid data or x_eval.
    Example: hermite_interpolate([0, 1], [1, 2], [1, 1], 0.5)['value']
        gives 1.5.
    """
    data = hermite_coefficients(x_nodes, y_nodes, dydx_nodes)
    try:
        x = np.asarray(x_eval, dtype=float)
    except (TypeError, ValueError, OverflowError):
        raise TypeError("x_eval must be numeric")
    if x.ndim > 1 or not np.all(np.isfinite(x)): raise ValueError("x_eval must be finite and 1-D")
    scalar = x.ndim == 0
    x = np.atleast_1d(x)
    z = data["z"]
    c = data["coeffs"]
    p = np.full_like(x, c[-1])
    dp = np.zeros_like(x)
    for k in range(len(c) - 2, -1, -1):
        dp = p + (x - z[k]) * dp
        p = c[k] + (x - z[k]) * p
    outside = bool(np.any((x < np.min(z)) | (x > np.max(z))))
    if outside: warnings.warn("Extrapolating outside the nodes", UserWarning, stacklevel=2)
    value = float(p[0]) if scalar else p
    derivative = None
    if return_derivative: derivative = float(dp[0]) if scalar else dp
    return {"value": value, "derivative": derivative, "degree": len(c) - 1, "extrapolating": outside}
if __name__ == "__main__":
    x3 = np.array([0.0, 0.5, 1.0])
    x6 = np.linspace(0, 1, 6)
    grid = np.linspace(0, 1, 1001)
    try:
        from lagrange import lagrange_interpolate
    except ImportError:
        def lagrange_interpolate(x_nodes, y_nodes, x_eval):
            result = np.zeros_like(np.asarray(x_eval, dtype=float))
            for i in range(len(x_nodes)):
                term = np.ones_like(result)
                for j in range(len(x_nodes)):
                    if i != j: term = term * (x_eval - x_nodes[j]) / (x_nodes[i] - x_nodes[j])
                result = result + y_nodes[i] * term
            return {"value": result}
    print("Problem 3(i): exp(x) comparison")
    print("Method                 error at 0.25     max grid error")
    for name, nodes in [("Lagrange 3 nodes", x3), ("Hermite 3 nodes", x3), ("Lagrange 6 nodes", x6)]:
        if name == "Hermite 3 nodes":
            at_point = hermite_interpolate(nodes, np.exp(nodes), np.exp(nodes), 0.25)["value"]
            on_grid = hermite_interpolate(nodes, np.exp(nodes), np.exp(nodes), grid)["value"]
        else:
            at_point = lagrange_interpolate(nodes, np.exp(nodes), 0.25)["value"]
            on_grid = lagrange_interpolate(nodes, np.exp(nodes), grid)["value"]
        at_point = abs(float(at_point) - np.exp(0.25))
        max_error = np.max(abs(on_grid - np.exp(grid)))
        print(f"{name:22s} {at_point:14.6e} {max_error:18.6e}")
    check = hermite_interpolate(x3, np.exp(x3), np.exp(x3), x3, True)
    print("\nProblem 3(ii): errors at nodes")
    print("value error:", np.max(abs(check["value"] - np.exp(x3))))
    print("derivative error:", np.max(abs(check["derivative"] - np.exp(x3))))
    y5 = x3**5 - 3*x3**2 + 1
    dy5 = 5*x3**4 - 6*x3
    result5 = hermite_interpolate(x3, y5, dy5, grid)["value"]
    result6 = hermite_interpolate(x3, x3**6, 6*x3**5, grid)["value"]
    print("\nProblem 3(iii): maximum errors")
    print("x^5 - 3x^2 + 1:", np.max(abs(result5 - (grid**5 - 3*grid**2 + 1))))
    print("x^6:", np.max(abs(result6 - grid**6)))
