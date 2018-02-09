package problems.nls;

import math.Complex;
import problems.QuantumProblem;

import static math.Complex.*;
import static java.lang.Math.*;

public abstract class NonlinearSchrodinger extends QuantumProblem {
    protected final String fileName;

    protected final int countOfPoints = (int) pow(2, 14);
    protected final int countOfIterations = 4;
    protected final int countOfWriting = 100;

    //Equation coefficients
    protected final double dispersionCoeff = 0.1; //b
    protected final double nonlinearCoeff = 1; //g
    protected final double nonlinearDefectCoeff = -0.26; // d
    protected final double linearDefectCoeff = 0; //e

    //Inputs: soliton
    protected final double boundaryCoordinate = 150;
    protected final double relativeCoordinateOfWaveCenter = 0.3;
    protected final double waveWidth = 20;
    protected final double waveVelocity = 5;
    protected double waveInitialAmplitude = 100; //use it if g = 0
    protected final Complex leftBoundaryCondition = ZERO;
    protected final Complex rightBoundaryCondition = ZERO;

    //Inputs: defect
    protected final double relativeCoordinateOfDefect = 0.5;
    protected final double defectWidth = 0.01;

    //For calculations
    protected final double dx;
    protected final double dt = 0.001;
    protected double t;

    protected final double coordinateOfWaveCenter;
    protected final double coordinateOfDefect;

    protected Complex[] waveFunction;
    protected double waveVector;
    protected double frequency;
    protected boolean isBright;

    protected double[] deltaFunction;
    protected double[] linearTerm;
    protected double[] nonlinearTerm;

    protected Complex a;
    protected Complex[] b;
    protected Complex c;
    protected Complex[] d;
    protected Complex[] p;
    protected Complex[] q;

    protected NonlinearSchrodinger(String fileName) {
        this.fileName = fileName;
        waveVector = waveVelocity / (2 * dispersionCoeff);
        dx = 2 * boundaryCoordinate / (countOfPoints - 1);
//        dt = 0.01;
        coordinateOfWaveCenter = boundaryCoordinate * (2 * relativeCoordinateOfWaveCenter - 1);
        coordinateOfDefect = boundaryCoordinate * (2 * relativeCoordinateOfDefect - 1);
        if (nonlinearCoeff != 0) {
            isBright = dispersionCoeff / nonlinearCoeff >= 0 ? true : false;
            if (isBright) {
                waveInitialAmplitude = dispersionCoeff / (nonlinearCoeff * pow(waveWidth, 2));
                frequency = pow(waveVelocity, 2) / (4 * dispersionCoeff) - dispersionCoeff / pow(waveWidth, 2);
            } else {
                waveInitialAmplitude = -dispersionCoeff / (nonlinearCoeff * pow(waveWidth, 2));
                frequency = pow(waveVelocity, 2) / (4 * dispersionCoeff) + dispersionCoeff / pow(waveWidth, 2);
            }
        }
        waveFunction = new Complex[countOfPoints];
        deltaFunction = new double[countOfPoints];
        linearTerm = new double[countOfPoints];
        nonlinearTerm = new double[countOfPoints];
        b = new Complex[countOfPoints];
        d = new Complex[countOfPoints];
        p = new Complex[countOfPoints];
        q = new Complex[countOfPoints];
        t = 0;
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
            setInitialWaveFunction(i, x);
            deltaFunction[i] = 1 / (defectWidth * sqrt(PI)) * exp(-pow((x - coordinateOfDefect) / defectWidth, 2));
//            deltaFunction[i] = -100 * exp(-pow((x - coordinateOfDefect) / defectWidth, 2));
            linearTerm[i] = linearDefectCoeff * deltaFunction[i];
            nonlinearTerm[i] = 2 * nonlinearCoeff * (1 + nonlinearDefectCoeff * deltaFunction[i]);
        }
    }

    private void setInitialWaveFunction(int i, double x) {
        double real = isBright ? waveInitialAmplitude / cosh((x - coordinateOfWaveCenter) / waveWidth) :
                waveInitialAmplitude * tanh((x - coordinateOfWaveCenter) / waveWidth);
        waveFunction[i] = multiply(real, exp(new Complex(0, waveVector * (x - coordinateOfWaveCenter))));
    }

    @Override
    public void writeToFile() {
        double[] modulus = new double[countOfPoints];
        double[] real = new double[countOfPoints];
        double[] image = new double[countOfPoints];
        double[] analyticalSolution = new double[countOfPoints];
        for (int i = 0; i < countOfPoints; i++) {
            double x = i * dx - boundaryCoordinate;
            real[i] = waveFunction[i].getReal();
            image[i] = waveFunction[i].getImage();
            modulus[i] = waveFunction[i].getSquaredModulus();
            analyticalSolution[i] = pow(waveInitialAmplitude /
                    cosh((x - coordinateOfWaveCenter - waveVelocity * t) / waveWidth), 2);
        }
        writeToFile(fileName, -boundaryCoordinate, dx, countOfPoints, nonlinearTerm, real, image, modulus, analyticalSolution);
    }

    protected void computeCoefficients() {
        setInitialPQ();
        for (int i = 1; i < countOfPoints - 1; i++) {
            b[i] = computeB(i);
            d[i] = computeD(i);
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

    protected abstract void setInitialPQ();

    protected abstract Complex computeB(int i);

    protected abstract Complex computeD(int i);

    protected Complex computeP(Complex denominator) {
        return divide(
                a,
                denominator
        );
    }

    protected Complex computeQ(int i, Complex denominator) {
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

    protected abstract void computeWaveFunction();
}
