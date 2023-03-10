[![](https://jitpack.io/v/micycle1/BetterBeziers.svg)](https://jitpack.io/#micycle1/BetterBeziers)

# BetterBeziers
High-precision utils for 2D Cubic Bezier Curves.

## Overview

BetterBeziers is a Java library that provides high-precision utilities for working with 2D Cubic Bezier Curves. The library offers two main functions: computing the length of a cubic bezier curve and finding the time `t` at which a specific arc length occurs.

The length computation function uses _Legendre-Gauss quadrature_, while the function for finding time `t` for a specified arc length uses a hybrid approach of Newton’s method and bisection. Both of these algorithms are described in the paper _Moving Along a Curve with Specified Speed_ [[1]].

With these two functions, it is easy to sample the curve (find t) for a given distance (or fractional distance) along the curve, as well as subdividing the curve into equal-length segments. BetterBeziers offers methods for all of these operations.

## Usage
To use BetterBeziers in your project, you can add it as a dependency using [JitPack](https://jitpack.io/#micycle1/BetterBeziers).

To compute the length of a cubic bezier curve, you can use the `getLength()` function. For example:

```
final double[][] controlPoints = {{0, 0}, {0.2, 0.3}, {0.7, 0.8}, {1, 1}};
CubicBezier bezier = new CubicBezier(controlPoints);
double length = bezier.getLength();
```

To find the time t at which a specific arc length occurs, you can use the `getCurveParameter()` function. For example:
```
final double[][] controlPoints = {{0, 0}, {0.2, 0.3}, {0.7, 0.8}, {1, 1}};
CubicBezier bezier = new CubicBezier(controlPoints);
arcLength = bezier.getLength()/3;
double t = bezier.getCurveParameter(arcLength);
```

[1]: https://www.geometrictools.com/Documentation/MovingAlongCurveSpecifiedSpeed.pdf
