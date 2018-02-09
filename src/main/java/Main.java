import problems.Problem;
import problems.nls.NonlinearSchrodingerSplitStepFourierMethod;

import java.io.IOException;

public class Main {
    public static void main(String[] args) throws IOException, InterruptedException {
//        Problem problem = new OldNonlinearSchrodingerCrankNicolsonMethod();
//        Problem problem = new OldNonlinearSchrodingerImplicitMethod();
        Problem problem = new NonlinearSchrodingerSplitStepFourierMethod();
//        Problem problem = new NonlinearSchrodingerNewtonMethod();
//        Problem problem = new SchrodingerOverBarrierReflection();
        problem.compute();
    }
}
