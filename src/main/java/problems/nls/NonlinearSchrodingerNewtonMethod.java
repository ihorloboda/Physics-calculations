package problems.nls;

import math.Complex;
import problems.QuantumProblem;

import static java.lang.Math.*;
import static math.Complex.*;

public class NonlinearSchrodingerNewtonMethod extends QuantumProblem {
    private final String fileName = "output/nls_newton_method.txt";

    private final int countOfPoints = 5000;
    private final int countOfWriting = 3;

    //TODO set these params
    //Equation coefficients
    private final double dispersionCoeff = 0.1; //b
    private final double nonlinearCoeff = 1; //g
    private final double nonlinearDefectCoeff = -0.26; // d
    private final double linearDefectCoeff = 0; //e

    //Inputs: soliton
    private final double boundaryCoordinate = 150;
    private final double relativeCoordinateOfWaveCenter = 0.4;
    private final double waveWidth = 5.75;
    private final double waveVelocity = 0.05;

    //Inputs: defect
    private final double relativeCoordinateOfDefect = 0.5;
    private final double defectWidth = 0.5;

    //For calculations
    private final double dx;
    private final double dt = 10;

    private final double coordinateOfWaveCenter;
    private final double coordinateOfDefect;

    private Complex[] waveFunction;
    private double waveInitialAmplitude;
    private double waveNumber;
    private double frequency;
    private boolean isBright;
    private Complex[] initialWaveFunction;
    private Complex[] deltaWaveFunction;

    private final double[] deltaFunction;
    private final double[] linearTerm;
    private final double[] nonlinearTerm;

    private Complex a;
    private Complex[] b;
    private Complex c;
    private Complex[] d;
    private Complex[] p;
    private Complex[] q;

    public NonlinearSchrodingerNewtonMethod() {
        dx = 2 * boundaryCoordinate / (countOfPoints - 1);
        coordinateOfWaveCenter = boundaryCoordinate * (2 * relativeCoordinateOfWaveCenter - 1);
        coordinateOfDefect = boundaryCoordinate * (2 * relativeCoordinateOfDefect - 1);
        waveNumber = waveVelocity / (2 * dispersionCoeff);
        isBright = dispersionCoeff / nonlinearCoeff >= 0 ? true : false;
        if (isBright) {
            waveInitialAmplitude = dispersionCoeff / (nonlinearCoeff * pow(waveWidth, 2));
            frequency = pow(waveVelocity, 2) / (4 * dispersionCoeff) - dispersionCoeff / pow(waveWidth, 2);
        } else {
            waveInitialAmplitude = -dispersionCoeff / (nonlinearCoeff * pow(waveWidth, 2));
            frequency = pow(waveVelocity, 2) / (4 * dispersionCoeff) + dispersionCoeff / pow(waveWidth, 2);
        }
//        waveInitialAmplitude = 1000;
        waveFunction = new Complex[countOfPoints];
        deltaWaveFunction = new Complex[countOfPoints];
        initialWaveFunction = new Complex[countOfPoints];
        deltaFunction = new double[countOfPoints];
        linearTerm = new double[countOfPoints];
        nonlinearTerm = new double[countOfPoints];
        a = multiply(dispersionCoeff, EIN);
        c = a;
        b = new Complex[countOfPoints];
        d = new Complex[countOfPoints];
        p = new Complex[countOfPoints];
        q = new Complex[countOfPoints];
    }

    @Override
    public void setInitialConditions() {
//        for (int i = 0; i < countOfPoints; i++) {
//            double x = i * dx - boundaryCoordinate;
//            setInitialWaveFunction(i, x);
//            deltaFunction[i] = 1 / (defectWidth * sqrt(PI)) * exp(-pow((x - coordinateOfDefect) / defectWidth, 2));
//            linearTerm[i] = linearDefectCoeff * deltaFunction[i];
//            nonlinearTerm[i] = 2 * nonlinearCoeff * (1 + nonlinearDefectCoeff * deltaFunction[i]);
//        }
    }

    private void setInitialWaveFunction(int i, double x) {
        double real = isBright ? waveInitialAmplitude / cosh((x - coordinateOfWaveCenter) / waveWidth) :
                waveInitialAmplitude * tanh((x - coordinateOfWaveCenter) / waveWidth);
        waveFunction[i] = multiply(real, exp(new Complex(0, waveNumber * x)));
        initialWaveFunction[i] = waveFunction[i];
    }

    @Override
    protected void computeIteration() {
        for (int j = 0; j < countOfWriting; j++) {
            computeCoefficients();
            computeWaveFunction();
        }
        initialWaveFunction = waveFunction;
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
        p[0] = ZERO;
        q[0] = ZERO;
        for (int i = 1; i < countOfPoints - 1; i++) {
//                b[i] = computeB(i);
//                d[i] = computeD(i);
//                Complex denominator = subtract(
//                        b[i],
//                        multiply(
//                                c,
//                                p[i - 1]
//                        )
//                );
//                p[i] = computeP(denominator);
//                q[i] = computeQ(i, denominator);
        }
    }

    private Complex computeB(int i) {
        return add(
                multiply(
                        pow(dx, 2),
                        subtract(1 / dt, freeTermDerivative(i))
                ),
                multiply(
                        2 * dispersionCoeff,
                        EIN
                )
        );
    }

    private Complex freeTermDerivative(int i) {
        return multiply(
                EIN,
                add(
                        linearTerm[i],
                        multiply(
                                nonlinearTerm[i],
                                add(
                                        waveFunction[i].getSquaredModulus(),
                                        multiply(2, sqr(waveFunction[i]))
                                ))
                )
        );
    }

    private Complex computeD(int i) {
        return add(
                multiply(
                        pow(dx, 2) / dt,
                        subtract(initialWaveFunction[i], waveFunction[i])
                ),
                multiply(
                        multiply(dispersionCoeff, EIN),
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
                nonlinearTerm[i] * waveFunction[i].getSquaredModulus(),
                multiply(
                        EIN,
                        waveFunction[i]
                )
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
    protected void computeNormCondition() {
        double sum = 0;
        for (int i = 0; i < countOfPoints; i++) {
            sum += waveFunction[i].getSquaredModulus();
        }
        System.out.println("Norm: " + sum);
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
        writeToFile(fileName, -boundaryCoordinate, dx, countOfPoints, deltaFunction, real, image, modulus);
    }
}
