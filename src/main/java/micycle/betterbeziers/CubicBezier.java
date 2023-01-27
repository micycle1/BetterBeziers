package micycle.betterbeziers;


import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;

/**
 * Better Cubic Bezier
 * 
 * BetterBeziers
 * 
 * @author Michael Carleton
 *
 */
public class CubicBezier {

	double[] p1, cp1, cp2, p2;
	final UnivariateFunction speedFunction;
	double length = Double.NaN;

	public CubicBezier(double[] p1, double[] cp1, double[] cp2, double[] p2) {
		this.p1 = p1;
		this.cp1 = cp1;
		this.cp2 = cp2;
		this.p2 = p2;

		speedFunction = makeSpeedFunction();
		init();
	}

	/**
	 * Computes length of this bezier via Legendre-Gauss quadrature (smartly
	 * integrating the curve to find its length). High precision. Approximate
	 * approach is better than other approximate approaches that simply lerp on the
	 * cached segment containing `t` because it also uses slope.
	 */
	public double getLength() {
		if (Double.isNaN(length)) {
			length = arcLength(1);
		}
		return length;
	}

	public double arcLength(double t) {
		return arcLength(0, t);
	}

	public double arcLength(double tStart, double tEnd) {
		/*
		 * The IterativeLegendreGaussIntegrator class is used to approximate the
		 * definite integral of the speed function from 0 to 1. The method returns an
		 * approximation of the definite integral of the speed function, which can be
		 * taken as a very good approximation of the length of the Bezier curve.
		 */
		int n = 5;
		IterativeLegendreGaussIntegrator integrator = new IterativeLegendreGaussIntegrator(n, 1e-3, 1e-5);
		double length = integrator.integrate(Integer.MAX_VALUE, speedFunction, tStart, tEnd);
		return length;
	}

	/**
	 * curvature/slope of the line tangent to curve at t
	 * 
	 * @param t
	 * @return
	 */
	public double getGradientAtParameter(double t) {
		return speedFunction.value(t);
	}

	public double[] getPointAtParameter(final double t) {
		final double t1 = 1.0f - t;
		final double x = p1[0] * t1 * t1 * t1 + 3 * cp1[0] * t * t1 * t1 + 3 * cp2[0] * t * t * t1 + p2[0] * t * t * t;
		final double y = p1[1] * t1 * t1 * t1 + 3 * cp1[1] * t * t1 * t1 + 3 * cp2[1] * t * t * t1 + p2[1] * t * t * t;
		return new double[] { x, y };
	}

	/**
	 * distance ratio
	 * 
	 * @param fraction
	 * @return
	 */
	public double[] getPointAtFraction(double fraction) {
		return getPointAtParameter(getCurveParameter(fraction * getLength()));
	}

	public double[] getPointAtLength(double length) {
		return getPointAtParameter(getCurveParameter(length));
	}

	private static final double CURVE_PARAMETER_EPSILON = 1E-2;
	private static final int CURVE_PARAMETER_MAX_ITER = 10;

	/**
	 * Finds time t for a specified arc length s. The algorithm is a hybrid of
	 * Newton’s method and bisection.
	 * 
	 * Given an arc length s, we want to know the time t at which this arc length
	 * occurs; that is, you first decide the distance the object should move along
	 * the curve and then figure out the time t at which that occurs. *
	 * <p>
	 * 
	 * @param s The arc length parameter. The value s is a measure of distance along
	 *          the curve. The arc-length parameterization of the curve is written
	 *          as X(s), where s ∈ [0, L] and L is the total length of the curve
	 * @return
	 */
	public double getCurveParameter(final double s) {
		if (s == 0) {
			return 0;
		}

		final double tmin = 0;
		final double tmax = 1;
		final double L = getLength();

		// Choose the initial guess for Newton’s method.
		double t = tmin + (tmax - tmin) * (s / L);
		double lower = tmin, upper = tmax;
		for (int iterate = 0; iterate < CURVE_PARAMETER_MAX_ITER; iterate++) {
			double F = arcLength(tmin, t) - s;
			if (Math.abs(F) <= CURVE_PARAMETER_EPSILON) {
				// |F(t)| is close enough to zero. Report t as the time at which length s is
				// attained.
				return t;
			}
			// Generate a candidate for Newton’s method.
			double DFDT = speedFunction.value(t); // F'(t) > 0 is guaranteed.
			double tCandidate = t - F / DFDT;
			// Update the root-bounding interval and test for containment of the candidate.
			if (F > 0) {
				// F > 0 implies tCandidate < t ≤ upper.
				upper = t;
				if (tCandidate <= lower) {
					// The candidate is outside the root-bounding interval, so use bisection
					// instead.
					t = (upper + lower) / 2;
				} else {
					// The candidate is in [lower,upper].
					t = tCandidate;
				}
			} else {
				// F < 0 implies lower ≤ t < tCandidate.
				lower = t;
				if (tCandidate >= upper) {
					// The candidate is outside the root-bounding interval, so use bisection
					// instead.
					t = (upper + lower) / 2;
				} else {
					// The candidate is in [lower,upper].
					t = tCandidate;
				}
			}
		}

		/*
		 * A root was not found according to the specified number of iterations and
		 * tolerance. You might want to increase imax or epsilon or the integration
		 * accuracy. However, in this application it is likely that the time values are
		 * oscillating because of the limited numerical precision of floating-point
		 * numbers. It is safe to report the last computed time as the t-value
		 * corresponding to the s-value.
		 */
		return t;
	}

	public double getCurveParameterApproximate(double s) {
		if (s <= 0) {
			return 0; // tmin
		}
		if (s >= getLength()) {
			return 1; // tmax
		}
		int i;
		for (i = 1; i < samples; i++) {
			if (s < sSample[i]) {
				break;
			}
		}
		double t = tSample[i - 1] + tsSlope[i] * (s - sSample[i - 1]);
		return t;
	}

	private int samples = 10;
	private double[] sSample, tSample, tsSlope;

	void init() {
		int n = samples;
		sSample = new double[n + 1];
		tSample = new double[n + 1];
		tsSlope = new double[n + 1];
		tsSlope[0] = 0; // tmin
		for (int i = 1; i < n; i++) {
			sSample[i] = getLength() * i / n;
			tSample[i] = getCurveParameter(sSample[i]);
			tsSlope[i] = (tSample[i] - tSample[i - 1]) / (sSample[i] - sSample[i - 1]);
		}
		sSample[n] = getLength();
		tSample[n] = 1; // tmax
		tsSlope[n] = (tSample[n] - tSample[n - 1]) / (sSample[n] - sSample[n - 1]);
	}
	
	public double[] subdivideByTime(double t) {
		return null;
	}

	/**
	 * Creates a speed function for this cubic bezier curve.
	 * <p>
	 * The speed is the magnitude of the first derivative of the Bezier curve (with
	 * respect to t), which is calculated using the standard Bezier curve equations.
	 * 
	 * @return
	 */
	private UnivariateFunction makeSpeedFunction() {
		return (double t) -> {
			double dx = 3 * (1 - t) * (1 - t) * (cp1[0] - p1[0]) + 6 * (1 - t) * t * (cp2[0] - cp1[0]) + 3 * t * t * (p2[0] - cp2[0]);
			double dy = 3 * (1 - t) * (1 - t) * (cp1[1] - p1[1]) + 6 * (1 - t) * t * (cp2[1] - cp1[1]) + 3 * t * t * (p2[1] - cp2[1]);
			return Math.sqrt(dx * dx + dy * dy);
		};
	}

}
