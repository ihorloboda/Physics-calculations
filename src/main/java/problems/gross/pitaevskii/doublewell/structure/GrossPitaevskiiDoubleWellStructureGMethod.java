package problems.gross.pitaevskii.doublewell.structure;

import static java.lang.Math.pow;

//this problem is not working
public class GrossPitaevskiiDoubleWellStructureGMethod extends GrossPitaevskiiDoubleWellStructure {
    private final double a;
    private final double c;
    private double[] b;
    private double[] d;
    private double[] p;
    private double[] q;

    private final double bTerm;

    public GrossPitaevskiiDoubleWellStructureGMethod() {
        super("GrossPitaevskiiDoubleWellStructureGMethod.txt");
        a = dt / (2 * dx);
        c = a;
        waveFunction = new double[countOfPoints];
        b = new double[countOfPoints];
        d = new double[countOfPoints];
        p = new double[countOfPoints];
        q = new double[countOfPoints];
        bTerm = 1 - chemicalPotential + 2 * a;
    }

    @Override
    protected void computeIteration() {
        for (int j = 0; j < countOfWriting; j++) {
            p[0] = -1;
            q[0] = 2 * leftBorderCondition;
            for (int i = 1; i < countOfPoints - 1; i++) {
                double nonlinearTerm = potential[i] * pow(waveFunction[i], 2) / 2;
                double denominator = (b[i] - c * p[i - 1]);
                b[i] = computeB(nonlinearTerm);
                d[i] = computeD(i, nonlinearTerm);
                p[i] = computeP(i, denominator);
                q[i] = computeQ(i, denominator);
            }
            waveFunction[countOfPoints - 1] = (2 * rightBorderCondition - q[countOfPoints - 2]) / (1 + p[countOfPoints - 2]);
            for (int i = countOfPoints - 2; i >= 0; i--) {
                waveFunction[i] = p[i] * waveFunction[i + 1] + q[i];
            }
        }
    }

    private double computeB(double nonlinearTerm) {
        return bTerm - nonlinearTerm;
    }

    private double computeD(int i, double nonlinearTerm) {
        return 1 + a * (waveFunction[i + 1] - 2 * waveFunction[i] + waveFunction[i - 1])
                + (chemicalPotential + nonlinearTerm) * waveFunction[i];
    }

    private double computeP(int i, double denominator) {
        return a / denominator;
    }

    private double computeQ(int i, double denominator) {
        return (d[i] + c * q[i - 1]) / denominator;
    }
}
