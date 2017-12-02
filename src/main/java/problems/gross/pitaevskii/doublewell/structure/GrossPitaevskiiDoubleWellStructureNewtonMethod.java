package problems.gross.pitaevskii.doublewell.structure;

import static java.lang.Math.*;

//it works
public class GrossPitaevskiiDoubleWellStructureNewtonMethod extends GrossPitaevskiiDoubleWellStructure {
    private double[] b;

    public GrossPitaevskiiDoubleWellStructureNewtonMethod() {
        super("GrossPitaevskiiDoubleWellStructureNewtonMethod.txt");
        b = new double[countOfPoints];
        a = 1;
        c = 1;
    }

    @Override
    protected void computeCoefficients() {
        p[0] = 0;
        q[0] = leftBorderCondition;
        for (int i = 1; i < countOfPoints - 1; i++) {
            double freeTerm = freeTermFunction(i);
            double freeTermDerivative = freeTermFunctionDerivative(i);
            b[i] = computeB(freeTermDerivative);
            d[i] = computeD(i, freeTerm, freeTermDerivative);
            double denominator = b[i] - c * p[i - 1];
            p[i] = computeP(denominator);
            q[i] = computeQ(i, denominator);
        }
    }

    private double computeB(double freeTermDerivative) {
        return 2 + pow(dx, 2) * freeTermDerivative;
    }

    private double computeD(int i, double freeTerm, double freeTermDerivative) {
        return pow(dx, 2) * (-freeTerm + freeTermDerivative * waveFunction[i]);
    }


    private double computeP(double denominator) {
        return a / denominator;
    }

    private double computeQ(int i, double denominator) {
        return (d[i] + c * q[i - 1]) / denominator;
    }

    @Override
    protected void computeWaveFunction() {
        waveFunction[countOfPoints - 1] = rightBorderCondition;
        for (int i = countOfPoints - 2; i >= 0; i--) {
            waveFunction[i] = p[i] * waveFunction[i + 1] + q[i];
        }
    }

    @Override
    protected double freeTermFunction(int i) {
        return -2 * (chemicalPotential - potential[i] * pow(waveFunction[i], 2)) * waveFunction[i];
    }

    private double freeTermFunctionDerivative(int i) {
        return -2 * chemicalPotential + 6 * potential[i] * pow(waveFunction[i], 2);
    }
}
