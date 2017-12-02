package problems.gross.pitaevskii.doublewell.structure;

import problems.QuantumProblem;

import static java.lang.Math.*;

public abstract class GrossPitaevskiiDoubleWellStructure extends QuantumProblem {
    protected final String fileName;

    //for sym A = B = 0.5;
    //for asym A = 0.6, B = 0.3
    //for antisym A = 2;
    protected final double coeffA = 0.6;
    protected final double coeffB = 0.3;
    protected final double borderCoordinate = 10;
    protected final int countOfPoints = 1001;
    protected final int countOfWriting = 1;
    protected final double chemicalPotential = -0.075;
    protected final double widthOfPotential = 0.2;
    protected final double leftBorderCondition = 0;
    protected final double rightBorderCondition = 0;

    protected double[] waveFunction;
    protected final double[] potential;

    protected double a;
    protected final double b;
    protected double c;
    protected double[] d;
    protected double[] p;
    protected double[] q;

    protected final double dx;
    protected final double dt = 0.1;

    protected GrossPitaevskiiDoubleWellStructure(String fileName) {
        this.fileName = fileName;
        potential = new double[countOfPoints];
        waveFunction = new double[countOfPoints];
        dx = 2 * borderCoordinate / countOfPoints;
        a = dt / pow(dx, 2);
        b = 2 * dt / pow(dx, 2) + 1;
        c = dt / pow(dx, 2);
        d = new double[countOfPoints];
        p = new double[countOfPoints];
        q = new double[countOfPoints];
    }

    @Override
    public void writeToFile() {
        writeToFile(fileName, -borderCoordinate, dx, countOfPoints, potential, waveFunction);
    }

    @Override
    public void setInitialConditions() {
        for (int i = 0; i < countOfPoints; i++) {
            double x = i * dx - borderCoordinate;
            potential[i] = -1 / (widthOfPotential * sqrt(PI)) * (exp(-pow((x + 1) / widthOfPotential, 2))
                    + exp(-pow((x - 1) / widthOfPotential, 2)));
            waveFunction[i] = initialFunction(x);
        }
    }

    @Override
    public void computeIteration() {
        for (int i = 0; i < countOfWriting; i++) {
            computeCoefficients();
            computeWaveFunction();
        }
    }

    protected abstract void computeCoefficients();

    protected abstract double freeTermFunction(int i);

    protected void computeWaveFunction() {
        waveFunction[countOfPoints - 1] = (2 * rightBorderCondition - q[countOfPoints - 2]) / (1 + p[countOfPoints - 2]);
        for (int i = countOfPoints - 2; i >= 0; i--) {
            waveFunction[i] = p[i] * waveFunction[i + 1] + q[i];
        }
    }

    @Override
    protected void computeNormCondition() {
        double sum = 0;
        for (int i = 0; i < countOfPoints; i++) {
            sum += pow(waveFunction[i], 2);
        }
        System.out.println("Norm: " + sum);
    }

    protected double initialFunction(double x) {
        return coeffA / cosh(2 * (x - 1)) + coeffB / cosh(2 * (x + 1));
//        return coeffA*sin(x)/cosh(x);
    }
}
