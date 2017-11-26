package problems.gross.pitaevskii.doublewell.structure;

import problems.QuantumProblem;

import static java.lang.Math.*;

public abstract class GrossPitaevskiiDoubleWellStructure extends QuantumProblem {
    protected final String fileName;

    protected final double borderCoordinate = 30;
    protected final int countOfPoints = 1001;
    protected final int countOfWriting = 1;
    protected final double chemicalPotential = -0.075;
    protected final double widthOfPotential = 0.02;
    protected final double coeffA = 0.5;
    protected final double coeffB = 1;
    protected final double leftBorderCondition;
    protected final double rightBorderCondition;

    protected double[] waveFunction;
    protected final double[] potential;

    protected final double dx;
    protected final double dt = 0.1;

    protected GrossPitaevskiiDoubleWellStructure(String fileName) {
        this.fileName = fileName;
        potential = new double[countOfPoints];
        leftBorderCondition = initialFunction(-borderCoordinate);
        rightBorderCondition = initialFunction(borderCoordinate);
        dx = 2 * borderCoordinate / countOfPoints;
    }

    @Override
    public void writeToFile() {
        writeToFile(fileName, -borderCoordinate, dx, countOfPoints, waveFunction);
    }

    @Override
    public void setInitialConditions() {
        for (int i = 0; i < countOfPoints; i++) {
            double x = i * dx - borderCoordinate;
            potential[i] = 2 / (widthOfPotential * sqrt(PI)) * (exp(-pow((x + 1) / widthOfPotential, 2))
                    + exp(-pow((x - 1) / widthOfPotential, 2)));
            waveFunction[i] = initialFunction(x);
        }
    }

    protected void computeNormCondition() {
        double sum = 0;
        for (int i = 0; i < countOfPoints; i++) {
            sum += pow(waveFunction[i], 2);
        }
        System.out.println("Norm: " + sum);
    }

    protected double initialFunction(double x) {
        return coeffA / cosh(2 * (x - 1)) + coeffB / cosh(2 * (x + 1));
    }

    protected abstract void computeIteration();
}
