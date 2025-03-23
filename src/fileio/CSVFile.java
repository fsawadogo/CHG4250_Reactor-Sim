package fileio;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

public class CSVFile {
    public static void writeData(String filename, String[] headers, double[][] data) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filename))) {
            // Write headers
            writer.write(String.join(",", headers));
            writer.newLine();

            // Write data rows
            for (double[] row : data) {
                for (int i = 0; i < row.length; i++) {
                    writer.write(row[i] + (i < row.length - 1 ? "," : ""));
                }
                writer.newLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}