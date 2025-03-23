package NumericalMethod;
public interface Function {
    double[] calculateValue(double x, double[] y, int index);
    // x:volume of the reactor.
    // y: the variables in the reactor that require for an ordinary differential equation to be solved
    // index: indicates the iteration number.
    //This method is called in the RK4 class allowing for all the ODEs that call on it in the PFR class to be calculated by the solver
}
