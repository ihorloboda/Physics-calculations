package problems;

import java.io.*;

public abstract class Problem {
    public abstract void setInitialConditions();

    public void compute() {
        setInitialConditions();
        writeToFile();
        while (true) {
            try {
                System.out.println("Press enter to do next iteration");
                System.in.read();
            } catch (IOException e) {
                e.printStackTrace();
            }
            computeIteration();
            writeToFile();
        }
    }

    public abstract void writeToFile();

    protected void writeToFile(String fileName, double from, double dx, int countOfPoints, double[]... data) {
        File file = new File(fileName);
        try (OutputStream out = new FileOutputStream(file);
             Writer writer = new OutputStreamWriter(out)) {
            for (int i = 0; i < countOfPoints; i++) {
                writer.append(String.valueOf(i * dx + from));
                writer.append(" ");
                for (int j = 0; j < data.length; j++) {
                    writer.append(String.valueOf(data[j][i]));
                    writer.append(" ");
                }
                writer.append(System.lineSeparator());
            }
        } catch (IOException e) {
            System.out.println("You have not access to file");
        }
    }

    protected abstract void computeIteration();
}
