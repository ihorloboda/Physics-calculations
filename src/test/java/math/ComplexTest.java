package math;

import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class ComplexTest {
    Complex c1;
    Complex c2;
    Complex c3;
    Complex c4;
    double d;

    double precision = 1E-10;

    @Before
    public void before() {
        c1 = new Complex(-2, 3);
        c2 = new Complex(3, -5);
        c3 = new Complex(1, Math.PI);
        c4 = new Complex(-1, 2);
        d = 2;
    }

    @Test
    public void getModulus() {
        double expected = 13;
        assertEquals(expected, c1.getSquaredModulus(), 1E-30);
    }

    @Test
    public void add() {
        Complex expected = new Complex(1, -2);
        assertEquals(expected, Complex.add(c1, c2));
    }

    @Test
    public void addMany() {
        Complex expected = new Complex(0, 0);
        assertEquals(expected, Complex.add(c1, c2, c4));
    }

    @Test
    public void subtractComplex() {
        Complex expected = new Complex(-5, 8);
        assertEquals(expected, Complex.subtract(c1, c2));
    }

    @Test
    public void substractDouble() {
        Complex expected = new Complex(4, -3);
        assertEquals(expected, Complex.subtract(d, c1));
    }

    @Test
    public void getNegative() {
        Complex expected = new Complex(2, -3);
        assertEquals(expected, Complex.negative(c1));
    }

    @Test
    public void multiplyComplex() {
        Complex expected = new Complex(9, 19);
        assertEquals(expected, Complex.multiply(c1, c2));
    }

    @Test
    public void multiplyMany() {
        Complex expected = new Complex(-47, -1);
        assertEquals(expected, Complex.multiply(c1, c2, c4));
    }

    @Test
    public void multiplyDouble() {
        Complex expected = new Complex(-4, 6);
        assertEquals(expected, Complex.multiply(d, c1));
    }

    @Test
    public void sqr() {
        Complex expected = new Complex(-3, -4);
        assertEquals(expected, Complex.sqr(c4));
    }

    @Test
    public void divideComplex() {
        Complex expected = new Complex(-21. / 34, -1. / 34);
        assertEquals(expected, Complex.divide(c1, c2));
    }

    @Test
    public void divideDouble() {
        Complex expected = new Complex(-1, 1.5);
        assertEquals(expected, Complex.divide(c1, d));
    }

    @Test
    public void exp() {
        Complex expected = new Complex(-Math.exp(1), 0);
        Complex actual = Complex.exp(c3);
        assertEquals(expected.getReal(), actual.getReal(), precision);
        assertEquals(expected.getImage(), actual.getImage(), precision);
    }
}
