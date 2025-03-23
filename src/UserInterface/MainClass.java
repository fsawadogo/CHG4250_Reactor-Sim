package UserInterface;

import NumericalMethod.RungeKutta4;
import Reactions.PrimaryReaction;
import Reactions.SecondaryReaction;
import Reactor.PFR;
import fileio.CSVFile;

import java.awt.*;
import java.util.ArrayList;
import java.util.List;

public class MainClass {
    public static final double M_EDC = 83.5; // molecular weight of EDC(kg/kmol)
    public static final double M_NH3 = 17.031; // molecular weight of NH3(kg/kmol)
    public static final double M_EDA = 117.26; // molecular weight of EDA (kg/kmol)
    public static final double M_HCl = 36.46; // molecular weight of HCl(kg/kmol)
    public static final double M_VC = 62.498; // molecular weight of VC(kg/kmol)
    public static final double M_NH4Cl = 53.491; // molecular weight of NH4Cl(kg/kmol)
    public static double X = 0.987;// final conversion in the reactor
    public static double I = 10;// inlet ratio between the NH3 and EDC inlet 30:1
    public static final double F_EDC0 = (42000 * 1000 * (1. / 365) * (1. / 24) * (1. / 60) * 0.962);


    //missing heat capacities for the temperature range given can i extrapolate
public static void main(String[] args) {
    double x_f = 7.74;//volume of the reactor in L
    double delx = 0.0774; //increment for 100 iterations
    System.out.println("The initial flow rate of the limiting reactant is " + F_EDC0);

    FileReader.readContents();
    if (FileReader.g_opConditions == null) {
        throw new IllegalArgumentException("The operating conditions text file is null");
    }
    
    // Initial flow rates
    double F_NH30 = (F_EDC0 * I);
    double F_HCl0 = 0.0;
    double F_NH4Cl0 = 0.0;
    double F_VC0 = 0.0;
    double F_EDA0 = 0.0;
    double Ta = 298;

    // Initialize reactors
    PrimaryReaction primaryReaction = new PrimaryReaction(PrimaryReaction.k1_0);
    SecondaryReaction secondaryReaction = new SecondaryReaction(SecondaryReaction.k2_0);
    PFR pfr = new PFR(primaryReaction, secondaryReaction);

    // Initial state
    double[] initialState = {
        FileReader.g_opConditions[0], // Initial Temperature
        FileReader.g_opConditions[1], // Initial Pressure
        F_EDC0,                       // Initial EDC flow 
        F_NH30,                       // Initial NH3 flow
        F_EDA0,                       // Initial EDA flow
        F_HCl0,                       // Initial HCl flow
        F_VC0,                        // Initial VC flow
        F_NH4Cl0,                     // Initial NH4Cl flow
        Ta,                           // Initial coolant temperature
    };

    List<double[]> resultsList = new ArrayList<>();
    
    // Add initial state to results
    double[] initialRow = new double[10];
    System.arraycopy(initialState, 0, initialRow, 0, 9);
    initialRow[9] = 0.0; // Initial conversion is 0
    resultsList.add(initialRow);
    
    // Debug output for initial state
    double X_initial = 0.0;
    System.out.printf("V: %.3f, T: %.2f, P: %.5f, F_EDC: %.5f, F_NH3: %.5f, F_EDA: %.5f, FHCl: %.5f, FVC: %.5f, F_NH4Cl: %.5f, X: %.5f%n", 
                     0.0, initialState[0], initialState[1], initialState[2], initialState[3], 
                     initialState[4], initialState[5], initialState[6], initialState[7], X_initial);

    // Integration loop
    double x = 0.0;
    double[] y = initialState.clone();
    
    while (x < x_f) {
        // Calculate next position
        double nextX = Math.min(x + delx, x_f);
        double stepSize = nextX - x;
        
        try {
            // Integrate from current position to next position
            // IMPORTANT: Pass current state y as starting point, not initialState
            y = RungeKutta4.rungeKutta4Integrate(x, nextX, y, stepSize/10, 100, pfr);
            
            // Calculate conversion at this point
            double X = (F_EDC0 - y[2]) / F_EDC0;
            if (X > 1.0) X = 1.0; // Cap conversion at 100%
            
            // Print current state
            System.out.printf("V: %.3f, T: %.2f, P: %.5f, F_EDC: %.5f, F_NH3: %.5f, F_EDA: %.5f, FHCl: %.5f, FVC: %.5f, F_NH4Cl: %.5f, X: %.5f%n", 
                             nextX, y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7], X);
            
            // Store current state
            double[] resultRow = new double[10];
            System.arraycopy(y, 0, resultRow, 0, 9);
            resultRow[9] = X;
            resultsList.add(resultRow);
            
            // Move to next position
            x = nextX;
            
        } catch (Exception e) {
            System.err.println("Error during integration: " + e.getMessage());
            // Reduce step size and try again
            delx = delx / 2.0;
            if (delx < 1e-6) {
                System.err.println("Step size too small, stopping integration");
                break;
            }
            System.out.println("Reducing step size to " + delx);
            // Don't update x, try again with smaller step
        }
    }
    
    // Convert results to array and write to CSV
    double[][] results = resultsList.toArray(new double[0][]);
    String[] headers = {"V (m3)", "T (K)", "P (atm)", "F_EDC (kmol/h)", "F_NH3 (kmol/h)", 
                       "F_EDA (kmol/h)", "F_HCl (kmol/h)", "F_VC (kmol/h)", "F_NH4Cl (kmol/h)", "X"};
    CSVFile.writeData("EDA_PFR_Results.csv", headers, results);
    System.out.println("Results have been exported successfully to EDA_PFR_Results.csv");
	}
}
