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
        int maxIt = 100; // Maximum Iterations
        System.out.println("The initial flow rate of the limiting reactant is"+F_EDC0);


        FileReader.readContents();
        if (FileReader.g_opConditions == null) {
            throw new IllegalArgumentException("The operating conditions text file is null");
        }
        // The initial flow rates of all components are calculated to ensure the annual production is attained

        double F_NH30 = (F_EDC0 * I);
        double F_HCl0 = 0.0;
        double F_NH4Cl0 = 0.0;
        double F_VC0 = 0.0;
        double F_EDA0 = 0.0;
        double Ta = 298;
        double C_EDC0 = (F_EDC0/(F_NH30+F_EDC0));
        double C_NH30 = (F_NH30/(F_NH30+F_EDC0));

        System.out.println("the initial flow rate of EDC" + F_EDC0);


        // Initialize reactors
        PrimaryReaction primaryReaction = new PrimaryReaction(PrimaryReaction.k1_0);
        SecondaryReaction secondaryReaction = new SecondaryReaction(SecondaryReaction.k2_0);
        PFR PFR = new PFR(primaryReaction, secondaryReaction);

        // Define initial state (8 elements)
        double[] initialState = {
                FileReader.g_opConditions[0], // Initial Temperature
                FileReader.g_opConditions[1], // Initial Pressure
                F_EDC0, //
                F_NH30,                         // FB0
                F_EDA0,                         // Initial Flow Rate
                F_HCl0,                         // Initial Flow Rate
                F_VC0,                          // Initial Flow Rate of Vinyl Chloride
                F_NH4Cl0,                       // Initial Flow Rate of Ammonium Chloride
                Ta,                             //Initial temperature of coolant



        };

        List<double[]> ResultsList = new ArrayList<>();

        double x = 0.0;
        double[] y = initialState.clone();
        while (x <= x_f) {
            double F_T = y[2] + y[3] + y[4] + y[5] + y[6] + y[7];
            double X = (F_EDC0 - y[2]) / F_EDC0;
            double C_EDC = y[2] / F_T;
            double C_NH3 = y[3] / F_T;
            System.out.println("the initial concentration of EDC and Ammonia are"+" "+C_EDC+C_NH3);
            // Store results for this step
            ResultsList.add(new double[]{x, y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7], y[8]});


            // Print the step
            System.out.printf("V: %.3f, T: %.2f, P: %.5f, F_EDC: %.5f, F_NH3: %.5f, F_EDA: %.5f, FHCl: %.5f, FVC: %.5f, F_NH4Cl: %.5f, X: %.5f%n", x, y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7], X);
            //Perform one integration step
            y = RungeKutta4.rungeKutta4Integrate(x, 7.74, initialState, delx, 100, PFR);
            x += delx;// Increment x
        }
        double [][] results = ResultsList.toArray(new double [0][]);
        //Write results to CSV
        String[] headers = {"V (m3)", "T (K)", "P (atm)", "F_EDC (kmol/h)", "F_NH3 (kmol/h)", "F_EDA (kmol/h)", "F_HCl (kmol/h)", "F_VC (kmol/h)", "F_NH4Cl (kmol/h)","X"};
        CSVFile.writeData("EDA_PFR_Results", headers, results);
        System.out.println("Results have been exported successfully to the EDA_PFR_Results.csv");
    }
}
