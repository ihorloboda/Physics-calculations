package problems.doublewell.structure;

import math.Complex;
import problems.QuantumProblem;

import static java.lang.Math.*;
import static math.Complex.*;

public abstract class GrossPitaevskiiDoubleWellStructure extends QuantumProblem {
    protected final String fileName;

    //for sym A = B = 0.5;
    //for asym A = 0.6, B = 0.3
    //for antisym A = 2;
    protected final double coeffA = 2;
    protected final double coeffB = 0.3;
    protected final boolean isAntisymmetric = true;
    protected final double borderCoordinate = 10;
    protected final int countOfPoints = 1001;
    protected final int countOfWriting = 1;
    protected final double chemicalPotential = -0.075;
    protected final double widthOfPotential = 0.2;
    protected final double leftBorderCondition = 0;
    protected final double rightBorderCondition = 0;

    protected double[] function;
    protected Complex[] waveFunction;
    protected final double[] potential;

    protected double a;
    protected final double b;
    protected double c;
    protected double[] d;
    protected double[] p;
    protected double[] q;

    protected final double dx;
    protected final double dt = 10;
    protected double t;

    protected GrossPitaevskiiDoubleWellStructure(String fileName) {
        this.fileName = fileName;
        potential = new double[countOfPoints];
        function = new double[countOfPoints];
        waveFunction = new Complex[countOfPoints];
        dx = 2 * borderCoordinate / countOfPoints;
        t = 0;
        a = dt / pow(dx, 2);
        b = 2 * dt / pow(dx, 2) + 1;
        c = dt / pow(dx, 2);
        d = new double[countOfPoints];
        p = new double[countOfPoints];
        q = new double[countOfPoints];
    }

    @Override
    public void writeToFile() {
        double[] real = new double[countOfPoints];
        double[] image = new double[countOfPoints];
        double[] modulus = new double[countOfPoints];
        for (int i = 0; i < countOfPoints; i++) {
            real[i] = waveFunction[i].getReal();
            image[i] = waveFunction[i].getImage();
            modulus[i] = waveFunction[i].getSquaredModulus();
        }
        writeToFile(fileName, -borderCoordinate, dx, countOfPoints, potential, function, real, image, modulus);
    }

    @Override
    public void setInitialConditions() {
        for (int i = 0; i < countOfPoints; i++) {
            double x = i * dx - borderCoordinate;
            potential[i] = -1 / (widthOfPotential * sqrt(PI)) * (exp(-pow((x + 1) / widthOfPotential, 2))
                    + exp(-pow((x - 1) / widthOfPotential, 2)));
            function[i] = initialFunction(x);
        }
        computeWaveFunction();
    }

    @Override
    public void computeIteration() {
        for (int i = 0; i < countOfWriting; i++) {
            t += dt;
            computeCoefficients();
            computeFunction();
        }
    }

    protected abstract void computeCoefficients();

    protected abstract double freeTermFunction(int i);

    protected void computeFunction() {
        function[countOfPoints - 1] = (2 * rightBorderCondition - q[countOfPoints - 2]) / (1 + p[countOfPoints - 2]);
        for (int i = countOfPoints - 2; i >= 0; i--) {
            function[i] = p[i] * function[i + 1] + q[i];
        }
    }

    protected void computeWaveFunction() {
        for (int i = 0; i < countOfPoints; i++) {
            waveFunction[i] = multiply(function[i], exp(new Complex(0, -chemicalPotential * t)));
        }
    }

    @Override
    protected void computeNormCondition() {
        double sum = 0;
        for (int i = 0; i < countOfPoints; i++) {
            sum += pow(function[i], 2);
        }
        System.out.println("Norm: " + sum);
        System.out.println("t = " + t);
    }

    protected double initialFunction(double x) {
        return isAntisymmetric ? coeffA * sin(x) / cosh(x) : coeffA / cosh(2 * (x - 1)) + coeffB / cosh(2 * (x + 1));
    }
}
