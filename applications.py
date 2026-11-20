import os
import warnings

import numpy as np
import matplotlib.pyplot as plt

from lagrange import lagrange_interpolate
from newton_differences import print_difference_table
from hermite import hermite_interpolate
from interpolate import choose_method, interpolate

# hide the extrapolation warnings so the output is easier to read
warnings.simplefilter("ignore")
os.makedirs("figures", exist_ok=True)


def application_a():
    """Census data: difference table, then estimates for 1955 and 2015.

    Parameters: none, the census data is inside the function.
    Returns: nothing, it prints the table and the two estimates.
    Raises: nothing expected, errors would come from interpolate.
    Example: application_a() prints 1955 as about 166.0 million.
    """
    print("\nApplication A")
    years = np.array([1950, 1960, 1970, 1980, 1990, 2000, 2010, 2020])
    pop = np.array([151.3, 179.3, 203.3, 226.5, 248.7, 281.4, 308.7, 331.4])

    print_difference_table(years, pop, decimals=1)

    for year in [1955, 2015]:
        r = interpolate(years, pop, year, degree=3)
        print(year, round(r["value"], 2), "million, method:", r["method_used"],
              "error estimate:", round(r["error_estimate"], 4))


def application_b():
    """Steam table: estimates at 25 and 55 C with degrees 1 to 4, plus a plot.

    Parameters: none, the steam table data is inside the function.
    Returns: nothing, it prints the estimates and saves figures/steam_table.png.
    Raises: nothing expected, errors would come from interpolate.
    Example: application_b() prints P(25) with degree 1 as 3.2931 kPa.
    """
    print("\nApplication B")
    T = np.array([20, 30, 40, 50, 60])
    P = np.array([2.3392, 4.2469, 7.3851, 12.352, 19.947])
    true_values = {25: 3.1698, 55: 15.763}

    for temp in [25, 55]:
        print("T =", temp, "method:", choose_method(T, P, temp)["method"])
        for deg in [1, 2, 3, 4]:
            r = interpolate(T, P, temp, degree=deg)
            rel_error = abs(r["value"] - true_values[temp]) / true_values[temp]
            estimate = r["error_estimate"]
            if estimate is not None:
                estimate = round(estimate, 6)
            print("  degree", deg, "P =", round(r["value"], 4),
                  "rel error =", round(rel_error, 6), "error estimate =", estimate)

    T_fine = np.linspace(20, 60, 400)
    P_fine = lagrange_interpolate(T, P, T_fine)["value"]
    plt.figure()
    plt.plot(T_fine, P_fine, label="Degree 4 interpolant")
    plt.plot(T, P, "--", label="Linear interpolation")
    plt.plot(T, P, "ko", label="Table values")
    plt.plot([25, 55], [3.1698, 15.763], "rs", label="Reference values")
    plt.xlabel("Temperature (°C)")
    plt.ylabel("Saturation pressure (kPa)")
    plt.legend()
    plt.savefig("figures/steam_table.png")
    plt.close()


def application_c():
    """Drug concentration: degree 7 Lagrange, check for negatives, find the peak.

    Parameters: none, the blood sample data is inside the function.
    Returns: nothing, it prints the results and saves figures/pk_curve.png.
    Raises: nothing expected, errors would come from lagrange_interpolate.
    Example: application_c() prints C(3) as about 5.6694 mg/L.
    """
    print("\nApplication C")
    t = np.array([0, 0.5, 1, 2, 4, 6, 8, 12])
    C = np.array([0.000, 4.215, 6.033, 6.516, 4.543, 2.809, 1.709, 0.629])

    print("method:", choose_method(t, C, 3)["method"])
    print("C(3) =", round(lagrange_interpolate(t, C, 3)["value"], 4), "mg/L")

    t_fine = np.linspace(0, 12, 2000)
    C_fine = lagrange_interpolate(t, C, t_fine)["value"]

    negative = C_fine < 0
    if negative.any():
        print("goes negative between t =", round(t_fine[negative].min(), 2),
              "and", round(t_fine[negative].max(), 2), "h")
        print("lowest value:", round(C_fine.min(), 4), "mg/L at t =",
              round(t_fine[np.argmin(C_fine)], 2), "h")

    peak = np.argmax(C_fine)
    print("degree 7 peak:", round(C_fine[peak], 4), "mg/L at t =", round(t_fine[peak], 3), "h")

    # simpler estimate: quadratic through the 3 points around the peak
    t_local = np.linspace(1, 4, 2000)
    C_local = lagrange_interpolate([1, 2, 4], [6.033, 6.516, 4.543], t_local)["value"]
    peak = np.argmax(C_local)
    print("quadratic peak:", round(C_local[peak], 4), "mg/L at t =", round(t_local[peak], 3), "h")

    plt.figure()
    plt.plot(t_fine, C_fine, label="Degree 7 interpolant")
    plt.plot(t, C, "ko", label="Measured samples")
    plt.axhline(0, color="gray")
    plt.xlabel("Time (h)")
    plt.ylabel("Concentration (mg/L)")
    plt.legend()
    plt.savefig("figures/pk_curve.png")
    plt.close()


def application_d():
    """Speeding car: average speeds, then Hermite for position and speed.

    Parameters: none, the checkpoint data is inside the function.
    Returns: nothing, it prints the results and saves figures/speed_profile.png.
    Raises: nothing expected, errors would come from hermite_interpolate.
    Example: application_d() prints the speed at t = 33 s as about 28.97 m/s.
    """
    print("\nApplication D")
    t = np.array([0, 12, 25, 40, 52])
    s = np.array([0, 312, 655, 1062, 1378])
    v = np.array([24.5, 25.0, 24.0, 25.5, 24.0])
    limit = 100 / 3.6

    print("limit:", round(limit, 2), "m/s, highest recorded speed:", v.max())
    for i in range(4):
        average = (s[i + 1] - s[i]) / (t[i + 1] - t[i])
        print("average speed from", t[i], "to", t[i + 1], "s:", round(average, 2))

    print("method:", choose_method(t, s, 18, dydx_nodes=v)["method"])

    for time in [18, 33, 46]:
        r = hermite_interpolate(t, s, v, time, return_derivative=True)
        print("t =", time, "position =", round(r["value"], 2), "m, speed =",
              round(r["derivative"], 2), "m/s")

    t_fine = np.linspace(0, 52, 2000)
    speed = hermite_interpolate(t, s, v, t_fine, return_derivative=True)["derivative"]
    print("highest speed:", round(speed.max(), 2), "m/s at t =", round(t_fine[np.argmax(speed)], 2), "s")

    # go through the curve and print where it goes over the limit
    start = None
    for i in range(len(t_fine)):
        if speed[i] > limit and start is None:
            start = t_fine[i]
        elif speed[i] <= limit and start is not None:
            print("over the limit from", round(start, 2), "to", round(t_fine[i], 2), "s")
            start = None
    if start is not None:
        print("over the limit from", round(start, 2), "to", round(t_fine[-1], 2), "s")

    plt.figure()
    plt.plot(t_fine, speed, label="Hermite speed")
    plt.plot(t, v, "ko", label="Recorded speeds")
    plt.axhline(limit, color="red", linestyle="--", label="Speed limit")
    plt.xlabel("Time (s)")
    plt.ylabel("Speed (m/s)")
    plt.legend()
    plt.savefig("figures/speed_profile.png")
    plt.close()


if __name__ == "__main__":
    application_a()
    application_b()
    application_c()
    application_d()