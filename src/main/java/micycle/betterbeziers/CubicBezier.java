package micycle.betterbeziers;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;

/**
 * Represents a 2D cubic Bezier curve defined by four control points.
 * <p>
 * Provides methods to calculate the curve's arclength, and sample point(s)
 * along the curve according to different parameter types.
 * <p>
 * The curve is defined by four control points: a start point (p1), an end point
 * (p2), and two control points (cp1 and cp2) that define the shape of the
 * curve. These control points are specified as arrays of two doubles that
 * represent the x and y coordinates of each point.
 * <p>
 * The length of the curve is calculated using numerical integration (the
 * <i>Gauss–Legendre</i> form) to approximate the integral of the gradient
 * function of the curve to a high precision.
 * <p>
 * The CubicBezier class is immutable; once a CubicBezier is instantiated, its
 * properties cannot be changed.
 * 
 * @author Michael Carleton
 */
public class CubicBezier {

	/**
	 * The maximum difference in measured vs user-specified arc length that is
	 * considered acceptable when checking for convergence during curve parameter
	 * computation.
	 */
	protected static final double CURVE_PARAMETER_EPSILON = 1E-3;
	protected static final int CURVE_PARAMETER_MAX_ITER = 10;

	/**
	 * The IterativeLegendreGaussIntegrator class is used to approximate the
	 * definite integral of the speed function from 0 to 1. The method returns an
	 * approximation of the definite integral of the speed function, which can be
	 * taken as a very good approximation of the length of the Bezier curve. For a
	 * given n (=6 here) the integrator can evaluate 2n-1 degree polynomials exactly
	 */
	protected static final IterativeLegendreGaussIntegrator integrator = new IterativeLegendreGaussIntegrator(6, 1e-3, 1e-5);

	/**
	 * The number of samples to use for approximating the curve parameter along the
	 * curve length. Increasing the number of samples will improve the accuracy of
	 * the approximation, but will also increase the memory and processing time
	 * requirements.
	 */
	protected static final int APPROXIMATION_SAMPLES = 50;
	private boolean approximationIsSetup = false;
	private double[] sSample, tSample, tsSlope;

	protected final double[] p1, cp1, cp2, p2;
	protected final UnivariateFunction speedFunction;
	protected double length = Double.NaN; // NaN marks yet to be computed and cached

	/**
	 * Constructs an instance of a {@code CubicBezier} using the given control
	 * points.
	 * 
	 * @param p1  first anchor point
	 * @param cp1 first control point
	 * @param cp2 second control point
	 * @param p2  second anchor point
	 */
	public CubicBezier(double[] p1, double[] cp1, double[] cp2, double[] p2) {
		this.p1 = p1;
		this.cp1 = cp1;
		this.cp2 = cp2;
		this.p2 = p2;

		speedFunction = makeSpeedFunction();
	}

	/**
	 * Constructs an instance of a {@code CubicBezier} using the given control
	 * points.
	 * 
	 * @param x1  the x-coordinate of the starting point
	 * @param y1  the y-coordinate of the starting point
	 * @param cx1 the x-coordinate of the first control point
	 * @param cy1 the y-coordinate of the first control point
	 * @param cx2 the x-coordinate of the second control point
	 * @param cy2 the y-coordinate of the second control point
	 * @param x2  the x-coordinate of the end point
	 * @param y2  the y-coordinate of the end point
	 */
	public CubicBezier(double x1, double y1, double cx1, double cy1, double cx2, double cy2, double x2, double y2) {
		this(new double[] { x1, y1 }, new double[] { cx1, cy1 }, new double[] { cx2, cy2 }, new double[] { x2, y2 });
	}

	/**
	 * Calculates the length of this bezier curve.
	 * 
	 * @return the length of the curve
	 * @see #getArcLength(double)
	 * @see #getArcLength(double, double)
	 */
	public double getLength() {
		if (Double.isNaN(length)) {
			length = getArcLength(1);
		}
		return length;
	}

	/**
	 * Calculates the arc length of the Bezier curve from the start point of the
	 * curve to the point on the curve specified by parameter value t.
	 * <p>
	 * The arc length is precisely approximated by numerically integrating the
	 * gradient function of the curve.
	 * 
	 * @param t The parameter value that specifies the point on the curve for which
	 *          the arc length is calculated. t ∈ [0...1].
	 * @return The arc length of the Bezier curve from the start point to the point
	 *         specified by t.
	 */
	public double getArcLength(double t) {
		return getArcLength(0, t);
	}

	/**
	 * Calculates the arc length of the Bezier curve between two bezier parameters:
	 * <code>tStart</code> and <code>tEnd</code>
	 * <p>
	 * The arc length is precisely approximated by numerically integrating the
	 * gradient function of the curve.
	 * 
	 * @param tStart The starting parameter value for calculating arc length.
	 * @param tEnd   The ending parameter value for calculating arc length.
	 * @return The arc length of the Bezier curve from tStart to tEnd.
	 */
	public double getArcLength(double tStart, double tEnd) {

		double length = integrator.integrate(Integer.MAX_VALUE, speedFunction, tStart, tEnd);
		return length;
	}

	/**
	 * 
	 * Returns the gradient of the Bezier curve at a specified bezier parameter.
	 * 
	 * @param t a value representing the parameter of the Bezier curve at which to
	 *          evaluate the gradient.
	 * @return a value representing the gradient of the Bezier curve at the
	 *         specified parameter.
	 */
	public double getGradientAtParameter(double t) {
		return speedFunction.value(t);
	}

	/**
	 * Returns the point on the Bezier curve corresponding to a specified parameter
	 * value.
	 * 
	 * @param t a value representing the parameter value of the Bezier curve at
	 *          which to evaluate the point. t ∈ [0...1].
	 * @return an array of size 2 representing the point on the Bezier curve
	 *         corresponding to the specified parameter value.
	 */

	public double[] getPointAtParameter(final double t) {
		final double t1 = 1.0f - t;
		final double x = p1[0] * t1 * t1 * t1 + 3 * cp1[0] * t * t1 * t1 + 3 * cp2[0] * t * t * t1 + p2[0] * t * t * t;
		final double y = p1[1] * t1 * t1 * t1 + 3 * cp1[1] * t * t1 * t1 + 3 * cp2[1] * t * t * t1 + p2[1] * t * t * t;
		return new double[] { x, y };
	}

	/**
	 * Returns the coordinates of the point on the cubic Bezier curve at the
	 * specified fraction of the curve's length.
	 * 
	 * @param fraction a value between 0 and 1, representing the fraction of the
	 *                 curve's length at which the point is to be evaluated.
	 * @return an array of doubles containing the x and y coordinates of the point
	 *         on the curve at the specified fraction of the curve's length.
	 */
	public double[] getPointAtFraction(double fraction) {
		return getPointAtParameter(getCurveParameter(fraction * getLength()));
	}

	/**
	 * Returns the coordinates of the point on the cubic Bezier curve at the
	 * specified distance along the curve.
	 * 
	 * @param length a double value representing the length/distance along the curve
	 *               at which the point is to be evaluated. length ∈ [0, L] where L
	 *               is the total length of the curve.
	 * @return an array of doubles containing the x and y coordinates of the point
	 *         on the curve at the specified length.
	 */
	public double[] getPointAtLength(double length) {
		return getPointAtParameter(getCurveParameter(length));
	}

	/**
	 * Samples <code>N</code> equidistant points along the cubic Bezier curve.
	 * <p>
	 * Note: the method will always return the last point on the Bezier curve.
	 * 
	 * @param numSamples an integer value representing the number of equidistant
	 *                   points to be sampled along the curve.
	 * @return a two-dimensional array of doubles containing the x and y coordinates
	 *         of the sampled points along the curve.
	 */
	public double[][] sampleEquidistantPoints(int numSamples) {
		numSamples = Math.max(2, numSamples);
		double[][] t = new double[numSamples][2];
		t[0] = p1.clone();
		t[numSamples - 1] = p2.clone();
		for (int i = 1; i < numSamples - 1; i++) {
			double fraction = (double) i / (numSamples - 1);
			t[i] = getPointAtFraction(fraction);
		}
		return t;
	}

	/**
	 * Samples and returns equidistant points along the cubic Bezier curve at a
	 * specified distance between points.
	 * <p>
	 * Note: the method will always return the last point on the Bezier curve.
	 * 
	 * @param distance a double value representing the distance between sampled
	 *                 points on the curve.
	 * @return a two-dimensional array of doubles containing the x and y coordinates
	 *         of the sampled points along the curve.
	 */
	public double[][] sampleEquidistantPoints(final double distance) {
		int numSamples = (int) Math.floor(getLength() / distance) + 2;
		numSamples = Math.max(2, numSamples);
		double[][] t = new double[numSamples][2];
		t[0] = p1.clone();
		t[numSamples - 1] = p2.clone();
		for (int i = 1; i < numSamples - 1; i++) {
			double length = i * distance;
			t[i] = getPointAtLength(length);
		}
		return t;
	}

	/**
	 * Finds the bezier parameter <code>t</code> for a specified arc length
	 * <code>s</code>. The algorithm is a hybrid of Newton’s method and bisection.
	 * 
	 * @param s The arc length parameter. The value s is a measure of distance along
	 *          the curve. s ∈ [0, L] where L is the total length of the curve.
	 * @return
	 */
	public double getCurveParameter(final double s) {
		/*
		 * The Newton-Raphson method is used to iteratively refine the estimate for the
		 * curve parameter t such that the arc length of the curve from the start point
		 * to t is as close to the target length s as possible. At each iteration, the
		 * value of F is computed as the difference between the actual length and the
		 * target length.
		 */
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
			double F = getArcLength(tmin, t) - s;
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

	/**
	 * Returns an approximate value of the curve parameter (<code>t</code>) for a
	 * given arc length (<code>s</code>) along the curve.
	 * 
	 * @param s the arc length along the curve to find the corresponding curve
	 *          parameter (t) for.
	 * @return an approximate value of the curve parameter (t) for the given arc
	 *         length (s).
	 */
	public double getCurveParameterApproximate(double s) {
		if (s <= 0) {
			return 0; // tmin
		}
		if (s >= getLength()) {
			return 1; // tmax
		}
		initApproximation();
		int i;
		for (i = 1; i < APPROXIMATION_SAMPLES; i++) {
			if (s < sSample[i]) {
				break;
			}
		}
		double t = tSample[i - 1] + tsSlope[i] * (s - sSample[i - 1]);
		return t;
	}

	/**
	 * Initialises the approximation data structures used to efficiently compute an
	 * approximate value of the curve parameter for a given arc length along the
	 * curve.
	 * <p>
	 * Note that this method only needs to be called once per instance of the class,
	 * as the approximation data structures are reused across calls to
	 * getCurveParameterApproximate.
	 */
	protected void initApproximation() {
		if (!approximationIsSetup) {
			int n = APPROXIMATION_SAMPLES;
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

			approximationIsSetup = true;
		}
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
