import numpy as np
from lagrange import validate_nodes, lagrange_interpolate
from newton_diff import _float_array, check_equal_spacing, newton_forward, newton_backward
from hermite import hermite_interpolate


def choose_method(x_nodes, y_nodes, x_eval, dydx_nodes=None, spacing_rtol=1e-9):
    x = validate_nodes(x_nodes)
    y = _float_array(y_nodes, "y_nodes")
    if x.ndim != 1 or x.size != y.size:
        raise ValueError("x_nodes and y_nodes must be 1-D arrays of the same length")

    eval_arr = _float_array(x_eval, "x_eval", allow_scalar=True)
    xe = float(np.mean(eval_arr))
    xs = np.sort(x)
    n = xs.size - 1
    extrapolating = bool(np.any((eval_arr < xs[0]) | (eval_arr > xs[-1])))

    try:
        h = check_equal_spacing(xs, spacing_rtol) if xs.size >= 2 else None
    except ValueError:
        h = None

    equally = h is not None

    if dydx_nodes is not None:
        return {
            "method": "hermite",
            "reason": "derivative values f'(x_i) are supplied, so Hermite uses them",
            "equally_spaced": equally,
            "h": h,
            "extrapolating": extrapolating,
        }

    if equally:
        s = (xe - xs[0]) / h
        if s <= n / 2 + 1e-12:
            method, where = "newton_forward", "start"
        else:
            method, where = "newton_backward", "end"
        reason = (
            f"nodes equally spaced with h={h:g}; x_eval={xe:g} lies near the "
            f"{where} of the table (s={s:g} of n={n})"
        )
    else:
        method = "lagrange"
        reason = "nodes are not equally spaced, so only a general-purpose method applies"

    return {
        "method": method,
        "reason": reason,
        "equally_spaced": equally,
        "h": h,
        "extrapolating": extrapolating,
    }


def interpolate(x_nodes, y_nodes, x_eval, dydx_nodes=None, degree=None, verbose=False):
    decision = choose_method(x_nodes, y_nodes, x_eval, dydx_nodes)
    method = decision["method"]

    if verbose:
        print(f"{method}: {decision['reason']}")

    if method == "hermite":
        result = dict(hermite_interpolate(x_nodes, y_nodes, dydx_nodes, x_eval))
    elif method == "newton_forward":
        result = dict(newton_forward(x_nodes, y_nodes, x_eval, degree))
    elif method == "newton_backward":
        result = dict(newton_backward(x_nodes, y_nodes, x_eval, degree))
    else:
        val, extrap = lagrange_interpolate(x_eval, x_nodes, y_nodes)
        is_scalar = np.asarray(x_eval).ndim == 0
        result = {
            "value": float(val) if is_scalar else val,
            "degree": len(np.asarray(x_nodes)) - 1,
            "extrapolating": bool(extrap),
        }

    result["method_used"] = method
    result["decision"] = decision
    return result


if __name__ == "__main__":
    import warnings

    cases = [
        ("T1", [0, 1, 2, 3, 4, 5], False, 0.4),
        ("T2", [0, 1, 2, 3, 4, 5], False, 4.7),
        ("T3", [0, 0.3, 1.1, 2.0, 3.5], False, 1.5),
        ("T4", [0, 0.3, 1.1, 2.0, 3.5], True, 1.5),
        ("T5", [5, 4, 3, 2, 1, 0], False, 0.4),
        ("T6", [0, 0.1, 0.2000000001, 0.3], False, 0.15),
        ("T7", [0, 1, 2, 3, 4, 5], False, 6.5),
    ]

    print(f"{'case':<5} {'method chosen':<16} {'value':>12} {'sin(x_eval)':>12} {'abs error':>11}")
    for name, xn, use_d, xv in cases:
        xn = np.array(xn, dtype=float)
        d = np.cos(xn) if use_d else None
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = interpolate(xn, np.sin(xn), xv, dydx_nodes=d)
        print(f"{name:<5} {r['method_used']:<16} {r['value']:>12.8f} {np.sin(xv):>12.8f} {abs(r['value'] - np.sin(xv)):>11.2e}")
        print("      reason:", r["decision"]["reason"], "| extrapolating:", r["decision"]["extrapolating"])