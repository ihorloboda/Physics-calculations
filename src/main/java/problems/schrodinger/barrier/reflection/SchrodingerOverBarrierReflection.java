package problems.schrodinger.barrier.reflection;

import math.Complex;
import problems.Problem;

import java.io.*;

public class SchrodingerOverBarrierReflection implements Problem {
    private final String fileName = "SchrodingerOverBarrierReflection.txt";

    private final int countOfPoints = 20000;
    private final int countOfWriting = 1000;
    private final double widthOfWavePackage = 5;
    private final double widthOfPotential = 15;
    private final double heightOfPotential = 23;
    private final double relativeCoordinateOfPackageCenter = 0.1;
    private final double relativeCoordinateOfPotentialCenter = 0.5;
    private final double waveVector = 7;

    private final double dx;
    private final double dt;
    private final double coordinateOfPackageCenter;
    private final double coordinateOFPotentialCenter;

    private Complex[] waveFunction;
    private double[] potential;

    private Complex a;
    private Complex[] b;
    private Complex c;
    private Complex[] d;
    private Complex[] p;
    private Complex[] q;

    public SchrodingerOverBarrierReflection() {
        dx = Math.PI / (50 * waveVector);
        dt = dx / waveVector;
        coordinateOfPackageCenter = (countOfPoints - 1) * dx * relativeCoordinateOfPackageCenter;
        coordinateOFPotentialCenter = (countOfPoints - 1) * dx * relativeCoordinateOfPotentialCenter;
        waveFunction = new Complex[countOfPoints];
        potential = new double[countOfPoints];
        a = new Complex(0, dt / (4 * Math.pow(dx, 2)));
        c = new Complex(0, dt / (4 * Math.pow(dx, 2)));
        b = new Complex[countOfPoints];
        d = new Complex[countOfPoints];
        p = new Complex[countOfPoints];
        q = new Complex[countOfPoints];
    }

    @Override
    public void setInitialConditions() {
        for (int i = 0; i < countOfPoints; i++) {
            double x = i * dx;
            potential[i] = heightOfPotential * Math.exp(-Math.pow((x - coordinateOFPotentialCenter) / widthOfPotential, 2));
            waveFunction[i] = Complex.exp(new Complex(-Math.pow((x - coordinateOfPackageCenter) / widthOfWavePackage, 2),
                    waveVector * x));
        }
    }

    @Override
    public void calculate() {
        setInitialConditions();
        while (true) {
            calculateWaveFunction();
            calculateNormCondition();
            writeToFile();
            try {
                System.in.read();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    @Override
    public void writeToFile() {
        File file = new File(fileName);
        try (OutputStream out = new FileOutputStream(file);
             Writer writer = new OutputStreamWriter(out)) {
            for (int i = 0; i < countOfPoints; i++) {
                writer.append(String.valueOf(i * dx));
                writer.append(" ");
                writer.append(String.valueOf(waveFunction[i].getReal()));
                writer.append(" ");
                writer.append(String.valueOf(potential[i]));
                writer.append(System.lineSeparator());
            }
        } catch (IOException e) {
            System.out.println("You have not access to file");
        }
    }

    private void calculateWaveFunction() {
        for (int j = 0; j < countOfWriting; j++) {
            p[0] = new Complex(-1, 0);
            q[0] = new Complex(0, 0);
            for (int i = 1; i < countOfPoints - 1; i++) {
                b[i] = calculateB(i);
                d[i] = calculateD(i);
                p[i] = calculateP(i);
                q[i] = calculateQ(i);
            }
            waveFunction[countOfPoints - 1] = Complex.divide(
                    Complex.getNegative(q[countOfPoints - 2]),
                    Complex.add(1, p[countOfPoints - 2])
            );
            for (int i = countOfPoints - 2; i >= 0; i--) {
                waveFunction[i] = Complex.add(
                        Complex.multiply(
                                waveFunction[i + 1],
                                p[i]
                        ),
                        q[i]
                );
            }
        }
    }

    private void calculateNormCondition() {
        double totalIntensity = 0;
        double sumBeforeBarrier = 0;
        double sumAfterBarrier = 0;
        for (int i = 0; i < countOfPoints; i++) {
            double waveFunctionModulus = waveFunction[i].getSquaredModulus();
            totalIntensity += waveFunctionModulus;
            if (i * dx < coordinateOFPotentialCenter) {
                sumBeforeBarrier += waveFunctionModulus;
            } else {
                sumAfterBarrier += waveFunctionModulus;
            }
        }
        System.out.println("Total intensity = " + totalIntensity);
        System.out.println("Reflected = " + sumBeforeBarrier / totalIntensity * 100 + " %");
        System.out.println("Transmitted = " + sumAfterBarrier / totalIntensity * 100 + " %");
        System.out.println(System.lineSeparator());
    }

    private Complex calculateB(int i) {
        return new Complex(1, dt / 2 * (1 / Math.pow(dx, 2) + potential[i]));
    }

    private Complex calculateD(int i) {
        return Complex.add(
                Complex.multiply(
                        a,
                        Complex.add(waveFunction[i - 1], waveFunction[i + 1])
                ),
                Complex.multiply(
                        Complex.subtract(new Complex(2, 0), b[i]),
                        waveFunction[i]
                )
        );
    }

    private Complex calculateP(int i) {
        return Complex.divide(
                a,
                Complex.subtract(
                        b[i],
                        Complex.multiply(
                                c,
                                p[i - 1]
                        )
                ));
    }

    private Complex calculateQ(int i) {
        return Complex.divide(
                Complex.add(
                        d[i],
                        Complex.multiply(
                                c,
                                q[i - 1]
                        )
                ),
                Complex.subtract(
                        b[i],
                        Complex.multiply(
                                c,
                                p[i - 1]
                        )
                )
        );
    }
}
