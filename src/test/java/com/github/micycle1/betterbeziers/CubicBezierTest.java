package com.github.micycle1.betterbeziers;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.BeforeAll;

public class CubicBezierTest {

	// A “straight‐line” cubic: p1=(0,0), cp1=(1/3,0), cp2=(2/3,0), p2=(1,0).
	// The curve is exactly the segment from (0,0) to (1,0), so its length is
	// exactly 1.0,
	// and for any t, pointAtParameter(t) == (t,0), arcLength(t) == t.
	private static CubicBezier straight;

	@BeforeAll
	public static void setup() {
		straight = new CubicBezier(new double[] { 0, 0 }, new double[] { 1.0 / 3.0, 0 }, new double[] { 2.0 / 3.0, 0 }, new double[] { 1, 0 });
	}

	@Test
	public void testLengthAndArcLength() {
		// Total length must be exactly 1.0
		double L = straight.getLength();
		assertEquals(1.0, L, 1e-9, "straight‐line total length");

		// arc length from 0 to t should equal t
		for (double t : new double[] { 0.0, 0.1, 0.5, 0.9, 1.0 }) {
			double s = straight.getArcLength(t);
			assertEquals(t, s, 1e-9, "arcLength(" + t + ")");
		}
	}

	@Test
	public void testGetPointAtParameter() {
		// getPointAtParameter(0) == p1, (1) == p2, and in between (t,0)
		assertArrayEquals(new double[] { 0, 0 }, straight.getPointAtParameter(0), 1e-9);
		assertArrayEquals(new double[] { 1, 0 }, straight.getPointAtParameter(1), 1e-9);
		for (double t : new double[] { 0.25, 0.5, 0.75 }) {
			double[] pt = straight.getPointAtParameter(t);
			assertEquals(t, pt[0], 1e-9);
			assertEquals(0.0, pt[1], 1e-9);
		}
	}

	@Test
	public void testGetCurveParameterExact() {
		double L = straight.getLength();
		// for s in [0,L], getCurveParameter(s) ≈ s/L = s
		for (double s : new double[] { 0, 0.1, 0.5, 0.9, 1.0 }) {
			double t = straight.getCurveParameter(s);
			assertEquals(s, t, CubicBezier.CURVE_PARAMETER_EPSILON, "exact getCurveParameter(" + s + ")");
		}
		// edge cases
		assertEquals(0.0, straight.getCurveParameter(0.0), 1e-12);
		assertEquals(1.0, straight.getCurveParameter(L), 1e-12);
	}

	@Test
	public void testGetCurveParameterApproximate() {
		double L = straight.getLength();
		// approximate version should be within about one slope‐interval tolerance
		for (double s : new double[] { 0.0, 0.123, 0.5, 0.876, 1.0 }) {
			double tex = straight.getCurveParameter(s);
			double tap = straight.getCurveParameterApproximate(s);
			assertEquals(tex, tap, 1e-3, "approximate vs exact t for s=" + s);
		}
		// call twice to ensure initApproximation is only one‐time
		double t1 = straight.getCurveParameterApproximate(0.42);
		double t2 = straight.getCurveParameterApproximate(0.42);
		assertEquals(t1, t2, 1e-9);
	}

	@Test
	public void testSampleEquidistantPointsByNumber() {
		// numSamples=5 ⇒ x coordinates should be 0,0.25,0.5,0.75,1
		double[][] pts = straight.sampleEquidistantPoints(5);
		assertEquals(5, pts.length);
		for (int i = 0; i < 5; i++) {
			assertEquals(i * 0.25, pts[i][0], 1e-9, "sample #" + i);
			assertEquals(0.0, pts[i][1], 1e-9);
		}
		// asking for 1 sample is clamped to at least 2
		double[][] two = straight.sampleEquidistantPoints(1);
		assertEquals(2, two.length);
		assertArrayEquals(new double[] { 0, 0 }, two[0], 1e-9);
		assertArrayEquals(new double[] { 1, 0 }, two[1], 1e-9);
	}

	@Test
	public void testSampleEquidistantPointsByDistance() {
		// distance = 0.3 ⇒ floor(1/0.3)=3 +2 ⇒ 5 points: at s=0,0.3,0.6,0.9,1.0
		double[][] pts = straight.sampleEquidistantPoints(0.3);
		assertEquals(5, pts.length);
		double[] expectedX = { 0.0, 0.3, 0.6, 0.9, 1.0 };
		for (int i = 0; i < pts.length; i++) {
			assertEquals(expectedX[i], pts[i][0], 1e-3, "distance sample #" + i);
			assertEquals(0.0, pts[i][1], 1e-9);
		}
		// distance > length ⇒ floor(1/2)=0+2 ⇒ 2 points
		double[][] two = straight.sampleEquidistantPoints(2.0);
		assertEquals(2, two.length);
	}
}