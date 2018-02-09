package problems.nls;

import math.Complex;

import static java.lang.Math.*;
import static math.Complex.*;

public class NonlinearSchrodingerNonlinearMethod extends NonlinearSchrodinger {
    private Complex[] initialWaveFunction;
    private Complex[] prevWaveFunction;
    private final Complex b;

    public NonlinearSchrodingerNonlinearMethod() {
        super("output/nls_nonlinear.txt");
        initialWaveFunction = new Complex[countOfPoints];
        prevWaveFunction = new Complex[countOfPoints];
        a = multiply(dispersionCoeff * dt / pow(dx, 2), EIN);
        b = add(1, multiply(2, a));
        c = a;
    }

    @Override
    protected void computeIteration() {
        for (int i = 0; i < countOfWriting; i++) {
            initialWaveFunction = waveFunction;
            prevWaveFunction = waveFunction;
            for (int j = 0; j < countOfIterations; j++) {
                computeCoefficients();
                computeWaveFunction();
                prevWaveFunction = waveFunction;
            }
            t += dt;
        }
        System.out.println("t = " + t);
    }

    @Override
    protected void computeCoefficients() {
        setInitialPQ();
        for (int i = 1; i < countOfPoints - 1; i++) {
            d[i] = computeD(i);
            Complex denominator = subtract(
                    b,
                    multiply(
                            c,
                            p[i - 1]
                    )
            );
            p[i] = computeP(denominator);
            q[i] = computeQ(i, denominator);
        }
    }

    @Override
    protected Complex computeB(int i) {
        return null;
    }

    @Override
    protected Complex computeD(int i) {
        return add(
                multiply(dt, freeTermFunction(i)),
                initialWaveFunction[i]
        );
    }

    private Complex freeTermFunction(int i) {
        return multiply(
                nonlinearTerm[i] * prevWaveFunction[i].getSquaredModulus() + linearTerm[i],
                multiply(EIN, initialWaveFunction[i])
        );
    }

    @Override
    protected void computeWaveFunction() {
        waveFunction[countOfPoints - 1] = rightBoundaryCondition;
        for (int i = countOfPoints - 2; i >= 0; i--) {
            waveFunction[i] = add(
                    multiply(
                            waveFunction[i + 1],
                            p[i]
                    ),
                    q[i]
            );
        }
    }

    @Override
    protected void setInitialPQ() {
        p[0] = ZERO;
        q[0] = leftBoundaryCondition;
    }
}
