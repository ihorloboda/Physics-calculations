package problems.diffusion;

import problems.Problem;

//this problem is not working
public class NonlinearDiffusion extends Problem {
    private final String fileName = "output/nld.txt";

    private final int countOfPoints = 3000;
    private final int countOfWriting = 10;
    private final double coordinateOfFunctionCenter = 50;
    private final double length = 150;
    private final double borderConditionAtTheBeginning = 100;
    private final double borderConditionAtTheEnd = 0;
    private final double functionMaxValue = 200;
    private final double widthOfFunction = 10;

    private final double dx;
    private final double dt = 0.1; //t relax ~  L^2/2D

    private double[] function;
    private double[] a;
    private double[] b;
    private double[] c;
    private double[] d;
    private double[] p;
    private double[] q;

    private final double diffFactor;

    public NonlinearDiffusion() {
        dx = length / countOfPoints;
        function = new double[countOfPoints];
        a = new double[countOfPoints];
        b = new double[countOfPoints];
        c = new double[countOfPoints];
        d = new double[countOfPoints];
        p = new double[countOfPoints];
        q = new double[countOfPoints];
        diffFactor = dt / (Math.pow(dx, 2));
    }

    @Override
    public void setInitialConditions() {
        for (int i = 0; i < countOfPoints; i++) {
            function[i] = functionMaxValue * Math.exp(-Math.pow((i * dx - coordinateOfFunctionCenter) / widthOfFunction, 2));
        }
    }

    @Override
    public void writeToFile() {
        writeToFile(fileName, 0, dx, countOfPoints, function);
    }

    protected void computeIteration() {
        for (int j = 0; j < countOfWriting; j++) {
            p[0] = -1;
            q[0] = 2 * borderConditionAtTheBeginning;
            for (int i = 1; i < countOfPoints - 1; i++) {
                a[i] = calculateA(i);
                c[i] = calculateC(i);
                b[i] = calculateB(i);
                d[i] = calculateD(i);
                p[i] = calculateP(i);
                q[i] = calculateQ(i);
            }
            function[countOfPoints - 1] = (q[countOfPoints - 2] + dx * borderConditionAtTheEnd) / (1 - p[countOfPoints - 2]);
            for (int i = countOfPoints - 2; i >= 0; i--) {
                function[i] = p[i] * function[i + 1] + q[i];
            }
        }
    }

    private double calculateA(int i) {
        double functionAverageValue = (function[i] + function[i - 1]) / 2;
        double diffusionCoeff = (Math.pow(functionAverageValue, 2) + 1) / (functionAverageValue + 1);
        return diffFactor * diffusionCoeff;
    }

    private double calculateB(int i) {
        return 1 + a[i] + c[i];
    }

    private double calculateC(int i) {
        double functionAverageValue = (function[i + 1] + function[i]) / 2;
        double diffusionCoeff = (Math.pow(functionAverageValue, 2) + 1) / (functionAverageValue + 1);
        return diffFactor * diffusionCoeff;
    }

    private double calculateD(int i) {
        return function[i];
    }

    private double calculateP(int i) {
        return a[i] / (b[i] - c[i] * p[i - 1]);
    }

    private double calculateQ(int i) {
        return (d[i] + c[i] * q[i - 1]) / (b[i] - c[i] * p[i - 1]);
    }
}
