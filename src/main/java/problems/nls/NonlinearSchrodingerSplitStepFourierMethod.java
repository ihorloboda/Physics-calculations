package problems.nls;

import math.Complex;

import static java.lang.Math.*;
import static math.Complex.*;
import static math.FastFourierTransform.*;

public class NonlinearSchrodingerSplitStepFourierMethod extends NonlinearSchrodinger {
    private Complex[] image;

    public NonlinearSchrodingerSplitStepFourierMethod() {
        super("output/nls_fft.txt");
        image = new Complex[countOfPoints];
    }

    @Override
    protected void computeIteration() {
        for (int j = 0; j < countOfWriting; j++) {
            for (int i = 0; i < countOfPoints; i++) {
                waveFunction[i] = multiply(
                        waveFunction[i],
                        exp(multiply(
                                dt,
                                multiply(EIN, nonlinearOperator(i))
                        ))
                );
            }
            image = fftShift(fft(waveFunction));
            for (int i = 0; i < countOfPoints; i++) {
                double waveNumber = 2 * PI * (i - countOfPoints / 2) / (2 * boundaryCoordinate);
                image[i] = multiply(
                        image[i],
                        exp(multiply(
                                -dt * dispersionCoeff * pow(2 * PI, 3 / 2) * pow(waveNumber, 2),
                                EIN
                                )
                        )
                );
            }
            waveFunction = ifft(fftShift(image));
            t += dt;
        }
        System.out.println("t = " + t);
    }

    private Complex nonlinearOperator(int i) {
        return new Complex(2 * nonlinearCoeff * (1 + nonlinearDefectCoeff * deltaFunction[i]) * waveFunction[i].getSquaredModulus()
                + linearDefectCoeff * deltaFunction[i], 0);
    }

    @Override
    protected void setInitialPQ() {
        throw new UnsupportedOperationException();
    }

    @Override
    protected Complex computeB(int i) {
        throw new UnsupportedOperationException();
    }

    @Override
    protected Complex computeD(int i) {
        throw new UnsupportedOperationException();
    }

    @Override
    protected void computeWaveFunction() {
        throw new UnsupportedOperationException();
    }
}
