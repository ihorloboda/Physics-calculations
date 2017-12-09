import problems.Problem;
import problems.doublewell.structure.GrossPitaevskiiDoubleWellStructure;
import problems.doublewell.structure.GrossPitaevskiiDoubleWellStructureNewtonMethod;
import problems.nls.NonlinearSchrodingerImplicitMethod;

import java.io.IOException;

public class Main {
    public static void main(String[] args) throws IOException, InterruptedException {
        Problem problem = new GrossPitaevskiiDoubleWellStructureNewtonMethod();
        problem.compute();
    }
}
