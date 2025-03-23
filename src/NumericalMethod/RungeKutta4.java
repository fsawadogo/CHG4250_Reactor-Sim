package NumericalMethod;

public class RungeKutta4 {
	public static double[] rungeKutta4Integrate(double x_0, double x_f, double[] y_0, double delx, int maxIt, Function f) {
	    // Initialize position and state
	    double x = x_0;
	    double[] y = new double[y_0.length];
	    System.arraycopy(y_0, 0, y, 0, y_0.length);
	    int i = 0;
	    
	    // Backup state in case integration step fails
	    double[] y_backup = new double[y.length];
	    
	    while (x < x_f && i < maxIt) {
	        // Save current state as backup
	        System.arraycopy(y, 0, y_backup, 0, y.length);
	        
	        // Calculate appropriate step size (don't go beyond x_f)
	        double stepSize = Math.min(delx, x_f - x);
	        
	        try {
	            // Initialize RK4 coefficient arrays
	            double[] k1 = new double[y.length];
	            double[] k2 = new double[y.length];
	            double[] k3 = new double[y.length];
	            double[] k4 = new double[y.length];
	            double[] rhs;
	            
	            // First RK4 stage
	            for (int j = 0; j < y.length; j++) {
	                rhs = f.calculateValue(x, y, j);
	                k1[j] = stepSize * rhs[j];
	            }
	            
	            // Second RK4 stage
	            double[] y2 = new double[y.length];
	            for (int h = 0; h < y.length; h++) {
	                y2[h] = y[h] + k1[h] / 2.0;
	            }
	            
	            for (int j = 0; j < y.length; j++) {
	                rhs = f.calculateValue(x + stepSize / 2.0, y2, j);
	                k2[j] = stepSize * rhs[j];
	            }
	            
	            // Third RK4 stage
	            double[] y3 = new double[y.length];
	            for (int h = 0; h < y.length; h++) {
	                y3[h] = y[h] + k2[h] / 2.0;
	            }
	            
	            for (int j = 0; j < y.length; j++) {
	                rhs = f.calculateValue(x + stepSize / 2.0, y3, j);
	                k3[j] = stepSize * rhs[j];
	            }
	            
	            // Fourth RK4 stage
	            double[] y4 = new double[y.length];
	            for (int h = 0; h < y.length; h++) {
	                y4[h] = y[h] + k3[h];
	            }
	            
	            for (int j = 0; j < y.length; j++) {
	                rhs = f.calculateValue(x + stepSize, y4, j);
	                k4[j] = stepSize * rhs[j];
	            }
	            
	            // Calculate new state using RK4 formula
	            for (int j = 0; j < y.length; j++) {
	                y[j] = y[j] + (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) / 6.0;
	            }
	            
	            // Ensure physical variables stay within reasonable bounds
	            // Index 0 is temperature, should be positive
	            if (y[0] <= 0 || y[0] > 1000) {
	                throw new IllegalArgumentException("Temperature out of physical range: " + y[0] + " K");
	            }
	            
	            // Index 1 is pressure, should be positive
	            if (y[1] < 0.1) {
	                y[1] = 0.1;  // Minimum pressure of 0.1 bar
	                System.out.println("Warning: Pressure reached minimum bound at x = " + x);
	            }
	            
	            // Indices 2-7 are flow rates, should be non-negative
	            for (int j = 2; j <= 7; j++) {
	                if (y[j] < 0) {
	                    y[j] = 0;  // Enforce non-negative flow rates
	                    System.out.println("Warning: Flow rate " + j + " reached zero at x = " + x);
	                }
	            }
	            
	            // Advance position
	            x += stepSize;
	            i++;
	            
	        } catch (Exception e) {
	            // If integration step fails, try a smaller step size
	            System.err.println("Integration error at x = " + x + ": " + e.getMessage());
	            
	            // Restore state from backup
	            System.arraycopy(y_backup, 0, y, 0, y.length);
	            
	            // Reduce step size by half
	            delx = delx / 2.0;
	            System.out.println("Reducing step size to " + delx);
	            
	            // If step size becomes too small, abort
	            if (delx < 1e-10) {
	                System.err.println("Step size too small, stopping integration at x = " + x);
	                break;
	            }
	            
	            // Try again with smaller step (don't increment counter)
	            continue;
	        }
	    }
	    
	    // Check if we completed integration or hit max iterations
	    if (i >= maxIt && x < x_f) {
	        System.err.println("Warning: Maximum iterations reached before endpoint. Reached x = " + 
	                         x + " of desired x_f = " + x_f);
	    }
	    
	    return y;
	}
}
