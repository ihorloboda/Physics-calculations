package problems.gross.pitaevskii.doublewell.structure;

import static java.lang.Math.*;

public class GrossPitaevskiiDoubleWellStructureNewtonMethod extends GrossPitaevskiiDoubleWellStructure {
    private double[] b;
    private double[] deltaWaveFunction;

    public GrossPitaevskiiDoubleWellStructureNewtonMethod() {
        super("GrossPitaevskiiDoubleWellStructureNewtonMethod.txt");
        b = new double[countOfPoints];
        deltaWaveFunction = new double[countOfPoints];
        a = 1;
        c = 1;
    }

    @Override
    protected void computeCoefficients() {
        p[0] = -1;
        q[0] = 2 * leftBorderCondition;
        for (int i = 1; i < countOfPoints - 1; i++) {
            double freeTerm = freeTermFunction(i);
            b[i] = computeB(i, freeTerm);
            d[i] = computeD(i, freeTerm);
            double denominator = (b[i] - c * p[i - 1]);
            p[i] = computeP(denominator);
            if (abs(p[i]) > 1) {
                throw new IllegalArgumentException("p[" + i + "] > 1");
            }
            q[i] = computeQ(i, denominator);
        }
    }

    //with df/dt
//    private double computeB(int i) {
//        return 2 * (1 + pow(dx, 2) * (1 / dt + 3 * potential[i] * pow(waveFunction[i], 2) - chemicalPotential));
//    }

    //with df/dt
//    private double computeD(int i) {
//        return waveFunction[i + 1] - 2 * waveFunction[i] + waveFunction[i - 1] + pow(dx, 2) * freeTermFunction(i);
//    }

    private double computeB(int i, double freeTerm) {
        return 2 + pow(dx, 2) * freeTerm;
    }

    private double computeD(int i, double freeTerm) {
        return waveFunction[i + 1] - 2 * waveFunction[i] + waveFunction[i - 1] - pow(dx, 2) * freeTerm;
    }


    private double computeP(double denominator) {
        return a / denominator;
    }

    private double computeQ(int i, double denominator) {
        return (d[i] + c * q[i - 1]) / denominator;
    }

    @Override
    protected void computeWaveFunction() {
        deltaWaveFunction[countOfPoints - 1] = (2 * rightBorderCondition - q[countOfPoints - 2]) / (1 + p[countOfPoints - 2]);
        for (int i = countOfPoints - 2; i >= 0; i--) {
            deltaWaveFunction[i] = p[i] * deltaWaveFunction[i + 1] + q[i];
        }
        for (int i = 0; i < countOfPoints; i++) {
            waveFunction[i] += deltaWaveFunction[i];
        }
    }

    @Override
    protected double freeTermFunction(int i) {
        return -2 * (chemicalPotential - potential[i] * pow(waveFunction[i], 2)) * waveFunction[i];
    }
}
