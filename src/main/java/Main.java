import problems.Problem;
import problems.nls.NonlinearSchrodingerImplicitMethod;

import java.io.IOException;

public class Main {
    public static void main(String[] args) throws IOException, InterruptedException {
        Problem problem = new NonlinearSchrodingerImplicitMethod();
        problem.compute();
    }
}
