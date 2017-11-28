import problems.Problem;
import problems.gross.pitaevskii.doublewell.structure.GrossPitaevskiiDoubleWellStructureGMethod;
import problems.gross.pitaevskii.doublewell.structure.GrossPitaevskiiDoubleWellStructureLinearMethod;
import problems.gross.pitaevskii.doublewell.structure.GrossPitaevskiiDoubleWellStructureNewtonMethod;
import problems.gross.pitaevskii.doublewell.structure.GrossPitaevskiiDoubleWellStructureNonlinearMethod;
import problems.schrodinger.barrier.reflection.SchrodingerOverBarrierReflection;

import java.io.IOException;

public class Main {
    public static void main(String[] args) throws IOException, InterruptedException {
        Problem problem = new GrossPitaevskiiDoubleWellStructureNewtonMethod();
        problem.compute();
    }
}
