package problems;

public interface Problem {
    void setInitialConditions();

    void calculate();

    void writeToFile();
}
