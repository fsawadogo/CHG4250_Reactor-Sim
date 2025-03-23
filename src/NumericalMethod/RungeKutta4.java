package NumericalMethod;

public class RungeKutta4 {
    public static double [] rungeKutta4Integrate (double x_0, double x_f, double [] y_0, double delx, int maxIt, Function f) {
        double x = x_0;
        double[] y = new double[y_0.length];
        System.arraycopy(y_0, 0, y, 0, y_0.length);
        int i = 0;
double[] rhs;
        while (x <= x_f && i <= maxIt) {
            double[] k1 = new double[y.length];
            double[] k2 = new double[y.length];
            double[] k3 = new double[y.length];
            double[] k4 = new double[y.length];

            for (int j = 0; j < y.length; j++) {
                rhs=f.calculateValue(x, y, j);
                k1[j] = delx*rhs[j];
            double[] y2 = new double[y.length];
            for (int h=0;h<y.length;h++){
                y2[h]=y[h]+k1[h]/2.0; }
            rhs=f.calculateValue(x + delx / 2.0, y2, j);
            k2[j] = delx *rhs[j];
            double [] y3 = new double[y.length];
            for (int h=0;h<y.length;h++){
                y3[h]=y[h]+k2[h]/2.0;
            }
            rhs=f.calculateValue(x+delx/2.0,y3,j);
            k3[j]=delx*rhs[j];
            double [] y4= new double[y.length];
            for (int h=0;h<y.length;h++){
                y4[h]=y[h]+k3[h];
            }
            rhs=f.calculateValue(x+delx,y4,j);
            k4[j]=delx*rhs[j];
        }
            for (int j=0;j<y.length;j++){
                y[j]=y[j]+(k1[j]+2*k2[j]+2*k3[j]+k4[j])/6.0;
    }
            x+=delx;
            i++;
            }
        if(i > maxIt) throw new ArithmeticException("endpoint not reached within maximum specified iterations in Euler integrate method");
        return y;
    }
}
