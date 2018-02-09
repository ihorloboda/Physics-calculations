package problems.doublewell.structure;

import math.Complex;

import static java.lang.Math.*;
import static math.Complex.*;

import problems.QuantumProblem;

public class GrossPitaevskiiDoubleWellStructureNewtonMethodComplex extends QuantumProblem {
    private final String fileName = "output/dws_newton_complex.txt";

    //for sym A = B = 0.5;
    //for asym A = 0.6, B = 0.3
    //for antisym A = 2;
    private final double coeffA = 2;
    private final double coeffB = 0.3;
    private final boolean isAntisymmetric = true;
    private final double boundaryCoordinate = 10;
    private final int countOfPoints = 1001;
    private final int countOfWriting = 10;
    private final double widthOfPotential = 0.2;

    private Complex[] waveFunction;
    private final Complex[] initialWaveFunction;
    private Complex[] deltaWaveFunction;
    private final double[] potential;

    private final Complex a;
    private Complex[] b;
    private final Complex c;
    private Complex[] d;
    private Complex[] p;
    private Complex[] q;

    private final double dx;
    private final double dt = 0.1;

    public GrossPitaevskiiDoubleWellStructureNewtonMethodComplex() {
        waveFunction = new Complex[countOfPoints];
        initialWaveFunction = new Complex[countOfPoints];
        deltaWaveFunction = new Complex[countOfPoints];
        potential = new double[countOfPoints];
        a = new Complex(0, 0.5);
        c = a;
        b = new Complex[countOfPoints];
        d = new Complex[countOfPoints];
        p = new Complex[countOfPoints];
        q = new Complex[countOfPoints];
        dx = 2 * boundaryCoordinate / (countOfPoints - 1);
    }

    @Override
    protected void computeNormCondition() {
        double sum = 0;
        for (int i = 0; i < countOfPoints; i++) {
            sum += waveFunction[i].getSquaredModulus();
        }
        System.out.println("Norm: " + sum);
    }

    @Override
    public void setInitialConditions() {
        for (int i = 0; i < countOfPoints; i++) {
            double x = i * dx - boundaryCoordinate;
            potential[i] = -1 / (widthOfPotential * sqrt(PI)) * (exp(-pow((x + 1) / widthOfPotential, 2))
                    + exp(-pow((x - 1) / widthOfPotential, 2)));
            waveFunction[i] = initialFunction(x);
            initialWaveFunction[i] = waveFunction[i];
        }
    }

    private Complex initialFunction(double x) {
        return isAntisymmetric ? new Complex(coeffA * sin(x) / cosh(x), 0) :
                new Complex(coeffA / cosh(2 * (x - 1)) + coeffB / cosh(2 * (x + 1)), 0);
    }

    @Override
    public void computeIteration() {
        for (int i = 0; i < countOfWriting; i++) {
            computeCoefficients();
            computeWaveFunction();
        }
    }

    private void computeWaveFunction() {
        deltaWaveFunction[countOfPoints - 1] = ZERO;
        for (int i = countOfPoints - 2; i >= 0; i--) {
            deltaWaveFunction[i] = add(
                    multiply(
                            deltaWaveFunction[i + 1],
                            p[i]
                    ),
                    q[i]
            );
            waveFunction[i] = add(waveFunction[i], deltaWaveFunction[i]);
        }
    }

    private void computeCoefficients() {
        for (int j = 0; j < countOfWriting; j++) {
            p[0] = ZERO;
            q[0] = ZERO;
            for (int i = 1; i < countOfPoints - 1; i++) {
//                b[i] = computeB(i);
//                d[i] = computeD(i);
                Complex denominator = subtract(
                        b[i],
                        multiply(
                                c,
                                p[i - 1]
                        )
                );
                p[i] = computeP(denominator);
                q[i] = computeQ(i, denominator);
            }
        }
    }

    private Complex computeB(int i) {
        return add(
                multiply(
                        pow(dx, 2),
                        subtract(1 / dt, freeTermDerivative(i))
                ),
                EIN
        );
    }

    private Complex computeD(int i) {
        return add(
                multiply(
                        pow(dx, 2) / dt,
                        subtract(initialWaveFunction[i], waveFunction[i])
                ),
                multiply(
                        a,
                        add(waveFunction[i + 1], multiply(-2, waveFunction[i]), waveFunction[i - 1])
                ),
                multiply(
                        pow(dx, 2),
                        freeTermFunction(i)
                )
        );
    }

    private Complex freeTermFunction(int i) {
        return multiply(
                -potential[i] * waveFunction[i].getSquaredModulus(),
                multiply(EIN, waveFunction[i])
        );
    }

    private Complex freeTermDerivative(int i) {
        return multiply(
                multiply(-potential[i], EIN),
                add(waveFunction[i].getSquaredModulus(),
                        multiply(2, sqr(waveFunction[i])))
        );
    }

    private Complex computeP(Complex denominator) {
        return divide(
                a,
                denominator
        );
    }

    private Complex computeQ(int i, Complex denominator) {
        return divide(
                add(
                        d[i],
                        multiply(
                                c,
                                q[i - 1]
                        )
                ),
                denominator
        );
    }

    @Override
    public void writeToFile() {
        double[] modulus = new double[countOfPoints];
        double[] real = new double[countOfPoints];
        double[] image = new double[countOfPoints];
        for (int i = 0; i < countOfPoints; i++) {
            modulus[i] = waveFunction[i].getSquaredModulus();
            real[i] = waveFunction[i].getReal();
            image[i] = waveFunction[i].getImage();
        }
        writeToFile(fileName, -boundaryCoordinate, dx, countOfPoints, potential, real, image, modulus);
    }
}
