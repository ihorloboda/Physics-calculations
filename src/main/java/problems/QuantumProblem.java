package problems;

import java.io.IOException;

public abstract class QuantumProblem extends Problem {
    @Override
    public void compute() {
        setInitialConditions();
        writeToFile();
        computeNormCondition();
        while (true) {
            try {
                System.out.println("Press enter to do next iteration");
                System.in.read();
            } catch (IOException e) {
                e.printStackTrace();
            }
            computeIteration();
            writeToFile();
            computeNormCondition();
        }
    }

    protected abstract void computeNormCondition();
}
