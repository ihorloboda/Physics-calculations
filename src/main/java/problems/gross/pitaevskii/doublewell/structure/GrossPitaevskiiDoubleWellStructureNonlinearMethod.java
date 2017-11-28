package problems.gross.pitaevskii.doublewell.structure;

//TODO save prevIteration
public class GrossPitaevskiiDoubleWellStructureNonlinearMethod extends GrossPitaevskiiDoubleWellStructure {
    private double[] prevWaveFunction;

    public GrossPitaevskiiDoubleWellStructureNonlinearMethod(String fileName) {
        super("GrossPitaevskiiDoubleWellStructureNonlinearMethod.txt");
        prevWaveFunction = new double[countOfPoints];
    }

    @Override
    protected void computeCoefficients() {

    }

    @Override
    protected double freeTermFunction(int i) {
        return 0;
    }
}
