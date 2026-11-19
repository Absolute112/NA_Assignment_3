
import numpy as np
import warnings
import math


def validate_nodes(nodes):
    nodes = np.asarray(nodes, dtype=float)

    if len(nodes) == 0:
        raise ValueError("At least one node is required")

    if not np.all(np.isfinite(nodes)):
        raise ValueError("Nodes must be finite")

    sorted_nodes = np.sort(nodes)

    for i in range(len(sorted_nodes) - 1):
        a = sorted_nodes[i]
        b = sorted_nodes[i + 1]

        tolerance = 1e-12 * max(1.0, abs(a))

        if abs(a - b) < tolerance:
            raise ValueError(f"node {b} appears twice")

    return nodes


def lagrange_basis(x_eval, nodes, k):
    nodes = validate_nodes(nodes)

    if k < 0 or k >= len(nodes):
        raise IndexError("Invalid basis index k")

    x_eval = np.asarray(x_eval, dtype=float)

    if not np.all(np.isfinite(x_eval)):
        raise ValueError("x_eval must contain finite values")

    result = np.ones_like(x_eval, dtype=float)

    for i in range(len(nodes)):
        if i != k:
            result *= (x_eval - nodes[i]) / (nodes[k] - nodes[i])

    return result


def lagrange_interpolate(x_eval, nodes, y_values):
    nodes = validate_nodes(nodes)

    y_values = np.asarray(y_values, dtype=float)
    x_eval = np.asarray(x_eval, dtype=float)

    if len(nodes) != len(y_values):
        raise ValueError("Nodes and y_values must have the same length")

    if not np.all(np.isfinite(y_values)):
        raise ValueError("y_values must be finite")

    if not np.all(np.isfinite(x_eval)):
        raise ValueError("x_eval must be finite")

    result = np.zeros_like(x_eval, dtype=float)

    for k in range(len(nodes)):
        result += y_values[k] * lagrange_basis(x_eval, nodes, k)

    extrapolating = np.any(
        (x_eval < np.min(nodes)) | (x_eval > np.max(nodes))
    )

    if extrapolating:
        warnings.warn("Extrapolation is being performed", UserWarning)

    return result, extrapolating


def interpolation_error_bound(x_eval, nodes, M):
    nodes = validate_nodes(nodes)
    x_eval = np.asarray(x_eval, dtype=float)

    if not np.all(np.isfinite(x_eval)):
        raise ValueError("x_eval must contain finite values")

    if not np.isfinite(M) or M < 0:
        raise ValueError("M must be finite and non-negative")

    n = len(nodes) - 1

    product = np.ones_like(x_eval, dtype=float)

    for i in range(len(nodes)):
        product *= (x_eval - nodes[i])

    return (M / math.factorial(n + 1)) * np.abs(product)
