package UserInterface;
import java.io.BufferedReader;
import java.io.IOException;

public class FileReader {

    public static double[] g_opConditions = new double[3];

    public static void readContents() {


        try {
            //BufferedReader Constructor is the more efficient F/IO
            BufferedReader reader = new BufferedReader(new java.io.FileReader("OperationalConditions.txt"));
            String line;
            int i = 0;

            while ((line = reader.readLine()) != null && i < g_opConditions.length)// The loop runs until the g_opConditions array is filled and contains 6 objects
            //the values on each line of the text file is assigned to one element in the array
            {
                String datapoint = line.split("=")[1].strip().trim();// this allows for the reader to only save the string placed after the = sign in the textfile and to remove anything before it
                //the trim.() removes any blank space that could be an issue later in the code
                g_opConditions[i] = Double.parseDouble(datapoint);// This parses the string input in the text file into a double value to store the values in the double array g_opConditions
                i++;
            }
            reader.close();//closing the reader  ensures that the reader does not run after it is not needed anymore, this release the BufferedReader constructor and reduces memory use and potential leaks1
            // Printing the initial operating conditions, allowing for the user to notice mistake in the text file
            System.out.println("Your current input values are as follow:");
            System.out.println("Initial Temperature = " + g_opConditions[0]+"K");//K
            System.out.println("Initial Pressure = " + g_opConditions[1]+"bar");// atm
            System.out.println("Annual Production Rate of EDA = " + g_opConditions[2]+"tons/year");// in tons/year
        } catch (IOException ex) {
            throw new RuntimeException(ex);
        }

    }
}
