[![](https://jitpack.io/v/micycle1/BetterBeziers.svg)](https://jitpack.io/#micycle1/BetterBeziers)

# BetterBeziers
High-precision utils for 2D Cubic Bezier Curves

## Overview

This library offers high-precision methods for computing the length of a cubic bezier and finding the time `t` at which a specific arc length occurs.

* Length computation uses _Legendre-Gauss quadrature_.
* Finding time `t` for a specified arc length uses a hybrid approach of Newton’s method and bisection.

Both of these algorithms are described in _Moving Along a Curve with Specified Speed_ [[1]].

Having these two functions, it becomes trivial to sample the curve (find `t`) for a distance (or fractional distance) along the curve, as well as subdividing the curve into equal-length segments – the library also offers methods for this.

[1]: https://www.geometrictools.com/Documentation/MovingAlongCurveSpecifiedSpeed.pdf
