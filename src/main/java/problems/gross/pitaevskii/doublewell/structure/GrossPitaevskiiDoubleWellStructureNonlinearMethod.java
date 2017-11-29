package problems.gross.pitaevskii.doublewell.structure;

import static java.lang.Math.*;

//may be it works, but use 2-3 iterations (see Kalitkin p.387)
public class GrossPitaevskiiDoubleWellStructureNonlinearMethod extends GrossPitaevskiiDoubleWellStructure {
    private final double b;

    public GrossPitaevskiiDoubleWellStructureNonlinearMethod() {
        super("GrossPitaevskiiDoubleWellStructureNonlinearMethod.txt");
        b = 2 * a + 1;
    }

    @Override
    protected void computeCoefficients() {
        p[0] = -1;
        q[0] = 2 * leftBorderCondition;
        for (int i = 1; i < countOfPoints - 1; i++) {
            d[i] = computeD(i);
            double denominator = (b - c * p[i - 1]);
            p[i] = computeP(denominator);
            if (abs(p[i]) > 1) {
                throw new IllegalArgumentException("p[" + i + "] > 1");
            }
            q[i] = computeQ(i, denominator);
        }
    }

    private double computeD(int i) {
        return dt * freeTermFunction(i) + waveFunction[i];
    }

    private double computeP(double denominator) {
        return a / denominator;
    }

    private double computeQ(int i, double denominator) {
        return (d[i] + c * q[i - 1]) / denominator;
    }

    @Override
    protected double freeTermFunction(int i) {
        //there is no difference between these two approximations
//        double approximation = (waveFunction[i + 1] + 2 * waveFunction[i] + waveFunction[i - 1]) / 4;
        double approximation = (waveFunction[i + 1] + waveFunction[i]) / 2;
        return approximation * (2 * chemicalPotential - potential[i] * pow(approximation, 2));
    }
}
