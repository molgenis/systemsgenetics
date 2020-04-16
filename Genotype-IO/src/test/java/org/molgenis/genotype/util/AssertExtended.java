/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.molgenis.genotype.util;

import org.testng.Assert;
import static org.testng.Assert.fail;

/**
 *
 * @author Patrick Deelen
 */
public class AssertExtended extends Assert {

	/**
	 * Asserts that two float arrays contain the same elements in the same
	 * order. If they do not, an AssertionError, with the given message, is
	 * thrown.
	 *
	 * @param actual the actual value
	 * @param expected the expected value
	 * @param delta the absolute tolerate value value between the actual and
	 * expected value
	 * @param message the assertion error message
	 */
	static public void assertEquals(final float[] actual, final float[] expected, float delta, final String message) {
		if (expected == actual) {
			return;
		}
		if (null == expected) {
			fail("expected a null array, but not null found. " + message);
		}
		if (null == actual) {
			fail("expected not null array, but null found. " + message);
		}

		org.testng.Assert.assertEquals(actual.length, expected.length, "arrays don't have the same size. " + message);

		for (int i = 0; i < expected.length; i++) {

			org.testng.Assert.assertEquals(actual[i], expected[i], delta, "arrays differ firstly at element [" + i + "]; "
					+ "expected value is <" + expected[i] + "> but was <"
					+ actual[i] + ">. "
					+ message);

		}
	}

	/**
	 * Asserts that two double arrays contain the same elements in the same
	 * order. If they do not, an AssertionError, with the given message, is
	 * thrown.
	 *
	 * @param actual the actual value
	 * @param expected the expected value
	 * @param delta the absolute tolerate value value between the actual and
	 * expected value
	 * @param message the assertion error message
	 */
	static public void assertEquals(final double[] actual, final double[] expected, double delta, final String message) {
		if (expected == actual) {
			return;
		}
		if (null == expected) {
			fail("expected a null array, but not null found. " + message);
		}
		if (null == actual) {
			fail("expected not null array, but null found. " + message);
		}

		org.testng.Assert.assertEquals(actual.length, expected.length, "arrays don't have the same size. " + message);

		for (int i = 0; i < expected.length; i++) {

			org.testng.Assert.assertEquals(actual[i], expected[i], delta, "arrays differ firstly at element [" + i + "]; "
					+ "expected value is <" + expected[i] + "> but was <"
					+ actual[i] + ">. "
					+ message);

		}
	}

	/**
	 * Asserts that two 2D float arrays contain the same elements in the same
	 * order. If they do not, an AssertionError, with the given message, is
	 * thrown.
	 *
	 * @param actual the actual value
	 * @param expected the expected value
	 * @param delta the absolute tolerate value value between the actual and
	 * expected value
	 * @param message the assertion error message
	 */
	static public void assertEquals(final float[][] actual, final float[][] expected, float delta, final String message) {
		if (expected == actual) {
			return;
		}
		if (null == expected) {
			fail("expected a null array, but not null found. " + message);
		}
		if (null == actual) {
			fail("expected not null array, but null found. " + message);
		}

		org.testng.Assert.assertEquals(actual.length, expected.length, "arrays don't have the same size. " + message);

		for (int i = 0; i < expected.length; i++) {

			assertEquals(actual[i], expected[i], delta, "arrays differ firstly at dimension 1 pos [" + i + "]; "
					+ message);

		}
	}

	/**
	 * Asserts that two 2D double arrays contain the same elements in the same
	 * order. If they do not, an AssertionError, with the given message, is
	 * thrown.
	 *
	 * @param actual the actual value
	 * @param expected the expected value
	 * @param delta the absolute tolerate value value between the actual and
	 * expected value
	 * @param message the assertion error message
	 */
	static public void assertEquals(final double[][] actual, final double[][] expected, double delta, final String message) {
		if (expected == actual) {
			return;
		}
		if (null == expected) {
			fail("expected a null array, but not null found. " + message);
		}
		if (null == actual) {
			fail("expected not null array, but null found. " + message);
		}

		org.testng.Assert.assertEquals(actual.length, expected.length, "arrays don't have the same size. " + message);

		for (int i = 0; i < expected.length; i++) {

			assertEquals(actual[i], expected[i], delta, "arrays differ firstly at dimension 1 pos [" + i + "]; "
					+ message);

		}
	}

	/**
	 * Asserts that two 3D double arrays contain the same elements in the same
	 * order. If they do not, an AssertionError, with the given message, is
	 * thrown.
	 *
	 * @param actual the actual value
	 * @param expected the expected value
	 * @param delta the absolute tolerate value value between the actual and
	 * expected value
	 * @param message the assertion error message
	 */
	static public void assertEquals(final double[][][] actual, final double[][][] expected, double delta, final String message) {
		if (expected == actual) {
			return;
		}
		if (null == expected) {
			fail("expected a null array, but not null found. " + message);
		}
		if (null == actual) {
			fail("expected not null array, but null found. " + message);
		}

		org.testng.Assert.assertEquals(actual.length, expected.length, "arrays don't have the same size. " + message);

		for (int i = 0; i < expected.length; i++) {

			assertEquals(actual[i], expected[i], delta, "arrays differ firstly at dimension 1 pos [" + i + "]; "
					+ message);

		}
	}
	
}
