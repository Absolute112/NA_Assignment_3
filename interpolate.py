import warnings

import numpy as np

from lagrange import validate_nodes, lagrange_interpolate
from newton_differences import check_equal_spacing, newton_forward, newton_backward
from hermite import hermite_interpolate


def choose_method(x_nodes, y_nodes, x_eval, dydx_nodes=None, spacing_rtol=1e-9):
    """Decide which interpolation method to use for this data.

    Hermite if derivatives are given. Otherwise, if the nodes are equally
    spaced, Newton forward when x_eval is in the first half of the table and
    Newton backward when it's in the second half. Otherwise Lagrange.

    Parameters: x_nodes, y_nodes are the data, x_eval is a float or array,
        dydx_nodes is the derivative data or None, spacing_rtol is the
        tolerance for counting the nodes as equally spaced.
    Returns: dict with method, reason, equally_spaced, h, extrapolating.
    Raises: TypeError for non-numbers, ValueError for bad or repeated nodes
        or if x and y have different lengths.
    Example: choose_method([0, 1, 2, 3], [0, 1, 4, 9], 0.5)["method"]
        gives "newton_forward".
    """
    x = np.sort(validate_nodes(x_nodes))
    if len(x) != len(y_nodes):
        raise ValueError("x_nodes and y_nodes must be the same length")

    x_eval = np.asarray(x_eval, dtype=float)
    point = float(np.mean(x_eval))  
    n = len(x) - 1
    extrapolating = bool(np.any((x_eval < x[0]) | (x_eval > x[-1])))

    # check_equal_spacing raises an error if the gaps aren't equal
    try:
        h = check_equal_spacing(x, spacing_rtol)
        equally_spaced = True
    except ValueError:
        h = None
        equally_spaced = False

    if dydx_nodes is not None:
        method = "hermite"
        reason = "derivative values are given, so Hermite can use them"
    elif equally_spaced:
        s = (point - x[0]) / h
        if s <= n / 2:
            method = "newton_forward"
            reason = f"nodes equally spaced with h={h:g}; x_eval={point:g} is near the start of the table (s={s:g} of n={n})"
        else:
            method = "newton_backward"
            reason = f"nodes equally spaced with h={h:g}; x_eval={point:g} is near the end of the table (s={s:g} of n={n})"
    else:
        method = "lagrange"
        reason = "nodes are not equally spaced, so Lagrange is used"

    return {
        "method": method,
        "reason": reason,
        "equally_spaced": equally_spaced,
        "h": h,
        "extrapolating": extrapolating,
    }


def interpolate(x_nodes, y_nodes, x_eval, dydx_nodes=None, degree=None, verbose=False):
    """Pick a method with choose_method, then run it.

    Parameters: same as choose_method, plus degree (only used by the Newton
        formulas) and verbose (prints which method was picked and why).
    Returns: the chosen method's result dict, plus method_used and decision.
    Raises: whatever choose_method or the chosen method raises.
    Example: interpolate([0, 1, 2, 3], [0, 1, 4, 9], 1.5)["value"] gives 2.25.
    """
    decision = choose_method(x_nodes, y_nodes, x_eval, dydx_nodes)
    method = decision["method"]

    if verbose:
        print(f"Using {method}: {decision['reason']}")

    if method == "hermite":
        result = hermite_interpolate(x_nodes, y_nodes, dydx_nodes, x_eval)
    elif method == "newton_forward":
        result = newton_forward(x_nodes, y_nodes, x_eval, degree)
    elif method == "newton_backward":
        result = newton_backward(x_nodes, y_nodes, x_eval, degree)
    else:
        result = lagrange_interpolate(x_nodes, y_nodes, x_eval)

    result["method_used"] = method
    result["decision"] = decision
    return result


if __name__ == "__main__":
    # Test cases T1-T7 from the assignment, with y = sin(x)
    cases = [
        ("T1", [0, 1, 2, 3, 4, 5], False, 0.4),
        ("T2", [0, 1, 2, 3, 4, 5], False, 4.7),
        ("T3", [0, 0.3, 1.1, 2.0, 3.5], False, 1.5),
        ("T4", [0, 0.3, 1.1, 2.0, 3.5], True, 1.5),
        ("T5", [5, 4, 3, 2, 1, 0], False, 0.4),
        ("T6", [0, 0.1, 0.2000000001, 0.3], False, 0.15),
        ("T7", [0, 1, 2, 3, 4, 5], False, 6.5),
    ]

    warnings.simplefilter("ignore")

    print(f"{'Case':<5} {'Method':<16} {'Value':>11} {'sin(x)':>11} {'Error':>10}")
    for name, x, has_derivs, x_eval in cases:
        x = np.array(x, dtype=float)
        dydx = np.cos(x) if has_derivs else None
        r = interpolate(x, np.sin(x), x_eval, dydx_nodes=dydx)
        error = abs(r["value"] - np.sin(x_eval))
        print(f"{name:<5} {r['method_used']:<16} {r['value']:>11.6f} {np.sin(x_eval):>11.6f} {error:>10.2e}")
        print(f"      why: {r['decision']['reason']}")