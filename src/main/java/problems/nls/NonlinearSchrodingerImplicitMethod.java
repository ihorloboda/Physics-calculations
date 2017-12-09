package problems.nls;

import math.Complex;
import problems.QuantumProblem;

import static java.lang.Math.*;
import static math.Complex.*;

public class NonlinearSchrodingerImplicitMethod extends QuantumProblem {
    private final String fileName = "nlsComplex.txt";

    private final int countOfPoints = 5000;
    private final int countOfWriting = 10;

    //TODO set these params
    private final double dispersionCoeff = 2;
    private final double nonlinearCoeff = 1;
    private final double nonlinearDefectCoeff = 1;
    private final double linearDefectCoeff = 1;

    private final double borderCoordinate = 150;
    private final double relativeCoordinateOfWaveCenter = 0.5;
    private final double waveWidth = 30;
    private final double waveVelocity = 1;
    private final Complex leftBorderCondition = new Complex(); //zero
    private final Complex rightBorderCondition = new Complex(); //zero

    private final double relativeCoordinateOfDefect = 0.5;
    private final double defectWidth = 0.5;

    private final double coordinateOfWaveCenter;
    private final double coordinateOfDefect;

    private Complex[] waveFunction;
    private double[] deltaFunction;

    private Complex[] nonlinearTerm;
    private Complex[] linearTerm;
    private final Complex bTerm;
    private final double waveVector;

    private final double dx;
    private final double dt;

    private final Complex a;
    private Complex[] b;
    private final Complex c;
    private Complex[] d;
    private Complex[] p;
    private Complex[] q;

    public NonlinearSchrodingerImplicitMethod() {
        dx = 2 * borderCoordinate / countOfPoints;
        dt = 0; //TODO set dt
        coordinateOfWaveCenter = borderCoordinate * (2 * relativeCoordinateOfWaveCenter - 1);
        coordinateOfDefect = borderCoordinate * (2 * relativeCoordinateOfDefect - 1);
        a = new Complex(dispersionCoeff / (2 * pow(dx, 2)), 0);
        c = a;
        bTerm = new Complex(dispersionCoeff / pow(dx, 2), -dt);
        waveVector = waveVelocity / (2 * dispersionCoeff);
        waveFunction = new Complex[countOfPoints];
        deltaFunction = new double[countOfPoints];
        nonlinearTerm = new Complex[countOfPoints];
        linearTerm = new Complex[countOfPoints];
        b = new Complex[countOfPoints];
        d = new Complex[countOfPoints];
        p = new Complex[countOfPoints];
        q = new Complex[countOfPoints];
    }

    @Override
    public void setInitialConditions() {
        for (int i = 0; i < countOfPoints; i++) {
            double x = i * dx - borderCoordinate;
//            waveFunction[i] = exp(new Complex(-pow((x - coordinateOfWaveCenter) / waveWidth, 2),
//                    waveVector * x));
            waveFunction[i] = multiply(dispersionCoeff / (nonlinearCoeff * pow(waveWidth, 2)
                    * cosh(x / waveWidth)), exp(new Complex(0, waveVector * x)));
            deltaFunction[i] = 1 / (defectWidth * sqrt(PI)) * exp(-pow((x - coordinateOfDefect) / defectWidth, 2));
            linearTerm[i] = new Complex(linearDefectCoeff * deltaFunction[i] / 2, 0);
        }
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
        double[] waveFunctionSquaredModulus = new double[countOfPoints];
        for (int i = 0; i < countOfPoints; i++) {
            waveFunctionSquaredModulus[i] = waveFunction[i].getSquaredModulus();
        }
        writeToFile(fileName, -borderCoordinate, dx, countOfPoints, waveFunctionSquaredModulus, deltaFunction);
    }

    @Override
    protected void computeIteration() {
        for (int j = 0; j < countOfWriting; j++) {
            p[0] = new Complex(-1, 0);
            q[0] = multiply(2, leftBorderCondition);
            for (int i = 1; i < countOfPoints - 1; i++) {
                nonlinearTerm[i] = computeNonlinearTerm(i);
                b[i] = computeB(i);
                d[i] = computeD(i);
                p[i] = computeP(i);
                q[i] = computeQ(i);
            }
            waveFunction[countOfPoints - 1] = divide(
                    subtract(
                            multiply(2, rightBorderCondition),
                            q[countOfPoints - 2]),
                    add(1, p[countOfPoints - 2])
            );
            for (int i = countOfPoints - 2; i >= 0; i--) {
                waveFunction[i] = add(
                        multiply(
                                waveFunction[i + 1],
                                p[i]
                        ),
                        q[i]
                );
            }
        }
    }

    private Complex computeNonlinearTerm(int i) {
        return new Complex(nonlinearCoeff * (1 + nonlinearDefectCoeff * deltaFunction[i])
                * waveFunction[i].getSquaredModulus(), 0);
    }

    private Complex computeB(int i) {
        return add(bTerm, nonlinearTerm[i], linearTerm[i]);
    }

    private Complex computeD(int i) {
        return add(
                multiply(
                        waveFunction[i],
                        add(nonlinearTerm[i], linearTerm[i], new Complex(0, 1 / dt))),
                multiply(
                        dispersionCoeff / (2 * pow(dx, 2)),
                        add(waveFunction[i + 1], multiply(-2, waveFunction[i]), waveFunction[i - 1])
                ));
    }

    private Complex computeP(int i) {
        return divide(
                a,
                subtract(
                        b[i],
                        multiply(
                                c,
                                p[i - 1]
                        )
                ));
    }

    private Complex computeQ(int i) {
        return divide(
                add(
                        d[i],
                        multiply(
                                c,
                                q[i - 1]
                        )
                ),
                subtract(
                        b[i],
                        multiply(
                                c,
                                p[i - 1]
                        )
                )
        );
    }
}
