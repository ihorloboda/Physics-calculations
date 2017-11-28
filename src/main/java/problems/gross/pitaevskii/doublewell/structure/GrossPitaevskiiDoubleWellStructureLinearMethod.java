package problems.gross.pitaevskii.doublewell.structure;

import static java.lang.Math.pow;

//doesn't work
public class GrossPitaevskiiDoubleWellStructureLinearMethod extends GrossPitaevskiiDoubleWellStructure {
    public GrossPitaevskiiDoubleWellStructureLinearMethod() {
        super("GrossPitaevskiiDoubleWellStructureLinearMethod.txt");
    }

    @Override
    protected void computeCoefficients() {
        p[0] = -1;
        q[0] = 2 * leftBorderCondition;
        for (int i = 1; i < countOfPoints - 1; i++) {
            double denominator = b - c * p[i - 1];
            d[i] = computeD(i);
            p[i] = computeP(i, denominator);
            q[i] = computeQ(i, denominator);
        }
    }

    private double computeD(int i) {
        return waveFunction[i] + freeTermFunction(i);
    }

    private double computeP(int i, double denominator) {
        return a / denominator;
    }

    private double computeQ(int i, double denominator) {
        return (d[i] + c * q[i - 1]) / denominator;
    }

//    @Override
//    protected double freeTermFunction(int i) {
//        double function = 0.25 * (waveFunction[i + 1] + 2 * waveFunction[i] + waveFunction[i - 1]);
//        return 2 * chemicalPotential * function + potential[i] * pow(function, 3);
//    }

    @Override
    protected double freeTermFunction(int i) {
        double function = 0.5 * (waveFunction[i + 1] + waveFunction[i]);
        return 2 * chemicalPotential * function + potential[i] * pow(function, 3);
    }
}
