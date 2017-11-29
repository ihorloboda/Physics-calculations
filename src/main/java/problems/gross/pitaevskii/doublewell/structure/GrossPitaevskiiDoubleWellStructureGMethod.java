package problems.gross.pitaevskii.doublewell.structure;

import static java.lang.Math.abs;
import static java.lang.Math.pow;

//doesn't work (tends to zero)
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
        a = dt / (2 * pow(dx, 2));
        c = a;
        bTerm = 1 - chemicalPotential + 2 * a;
        waveFunction = new double[countOfPoints];
        b = new double[countOfPoints];
        d = new double[countOfPoints];
        p = new double[countOfPoints];
        q = new double[countOfPoints];
    }

    @Override
    protected void computeCoefficients() {
        p[0] = -1;
        q[0] = 2 * leftBorderCondition;
        for (int i = 1; i < countOfPoints - 1; i++) {
            double nonlinearTerm = freeTermFunction(i);
            b[i] = computeB(nonlinearTerm);
            d[i] = computeD(i, nonlinearTerm);
            double denominator = (b[i] - c * p[i - 1]);
            p[i] = computeP(denominator);
            if (abs(p[i]) > 1) {
                throw new IllegalArgumentException("p[" + i + "] > 1");
            }
            q[i] = computeQ(i, denominator);
        }
    }

    private double computeB(double nonlinearTerm) {
        return bTerm - nonlinearTerm;
    }

    private double computeD(int i, double nonlinearTerm) {
        return 1 + a * (waveFunction[i + 1] - 2 * waveFunction[i] + waveFunction[i - 1])
                + (chemicalPotential + nonlinearTerm) * waveFunction[i];
    }

    private double computeP(double denominator) {
        return a / denominator;
    }

    private double computeQ(int i, double denominator) {
        return (d[i] + c * q[i - 1]) / denominator;
    }

    @Override
    protected double freeTermFunction(int i) {
        return -2 * potential[i] * pow(waveFunction[i], 2);
    }
}
