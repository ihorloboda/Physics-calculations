import problems.Problem;
import problems.diffusion.NonLinearDiffusion;

import java.io.IOException;

public class Main {
    public static void main(String[] args) throws IOException, InterruptedException {
        Problem problem = new NonLinearDiffusion();
//        Problem problem= new SchrodingerOverBarrierReflection();
        problem.calculate();
    }
}