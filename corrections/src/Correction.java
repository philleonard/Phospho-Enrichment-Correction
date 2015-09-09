/**
 * Created by ppl on 15/06/15.
 */
import java.io.*;
import java.util.ArrayList;
import java.util.List;

public class Correction implements Runnable{
    static double PE_THRESH = 75.00;
    static int isoforms = 0;
    static String destination = "";
    static String newDestinationOver = "";
    static String newDestinationUnder = "";
    static MainPanel mp;
    static boolean remUn;
    static boolean remDup;
    static boolean remDupSeq;

    public Correction(String dest, String newDest, String newDestUnd, boolean remDup, boolean remUn, boolean remDupSeq, double thresh, MainPanel mp){
        destination = dest;
        newDestinationOver = newDest;
        newDestinationUnder = newDestUnd;
        this.remDup = remDup;
        this.remUn = remUn;
        this.mp = mp;
        this.remDupSeq = remDupSeq;
        PE_THRESH = thresh;

    }

    public static void fix() {
        isoforms = 0;
        System.out.printf("==========================\nPhospho Enrichment Correction Program v0.1\n");
        System.out.printf("Opening %s\n", destination);
        System.out.printf("Reading in data...");
        String [] raw = readIn(destination);
        System.out.printf(" done\n");

        System.out.printf("Correcting data\n");
        String [] dataFinal = corrections(raw);

        if (dataFinal != null) {
            System.out.printf("Writing out corrections...");
            writeOut(dataFinal);
            System.out.printf(" done\n");
            System.out.printf("Wrote over %.2f corrections to %s\n", PE_THRESH, newDestinationOver);
            System.out.printf("Wrote under %.2f corrections to %s\n", PE_THRESH, newDestinationUnder);
        }
        else
            System.out.printf("Not writing corrections\n");


        System.out.printf("Finished!\n==========================\n");
        mp.getRunButton().setEnabled(true);
    }

    public static String[] corrections(String [] raw) {
        //Check columns here, and return null if incorrect or missing
        String[] metaData = getMeta(raw);
        if (metaData != null) {

            //Remove Unused...
            if (remUn)
                raw = removeUnused(getOffset("Quan Usage", metaData), raw);

            //Remove Duplicate Scans
            if (remDup)
                raw = removeDuplicateScans(getOffset("Scan", metaData), raw);

            //Capitalise sequence Data
            raw = capitalise(getOffset("Sequence", metaData), raw);

            //Repairing the probabilities in Sequence
            raw = repairProb(getOffset("Sequence", metaData), getOffset("phosphoRS Site Probabilities", metaData), raw);

            if (remDupSeq)
                raw = removeDuplicateSeq(raw, getOffset("Sequence", metaData), getOffset("Heavy/Light", metaData));
        }
        else
            return null;
        return raw;
    }

    public static String [] removeDuplicateSeq(String [] raw, int seqOff, int heavyOff) {
        System.out.printf("Removing duplicate sequence, choosing closest weight to 1...");
        int dups = 0;
        String [] line;
        int notNumber = 0;
        int i = 0;
        while (i < raw.length) {
            if (!raw[i].equals("")) {
                line = raw[i].split(",");
                int j = i + 1;
                boolean foundDup = false;
                try {
                    double minWeight = Math.abs(1 - Double.parseDouble(line[heavyOff]));
                    if (j < raw.length) {
                        String[] nextLine = raw[j].split(",");
                        int keep = i;
                        boolean cont = false;
                        try {
                            cont = line[seqOff].equals(nextLine[seqOff]);
                        } catch (Exception e) {}
                        while ((cont || raw[j].equals("")) && j < raw.length) {
                            try {
                                cont = line[seqOff].equals(nextLine[seqOff]);
                            } catch (Exception e) {}
                            if (!raw[j].equals("")) {
                                double newMinWeight = 0.0;
                                try {
                                    newMinWeight = Math.abs(Double.parseDouble(nextLine[heavyOff]));
                                } catch (ArrayIndexOutOfBoundsException aob) {
                                    //raw[j] = "";
                                }
                                if (newMinWeight < minWeight) {
                                    minWeight = newMinWeight;
                                    keep = j;
                                }
                                j++;
                                nextLine = raw[j].split(",");
                                foundDup = true;
                                if (j > raw.length)
                                    cont = false;
                            }
                            else
                                j++;

                            if (!cont)
                                break;
                        }
                        int newI = i;

                        while (newI < j && foundDup) {
                            if (newI != keep) {
                                raw[newI] = "";
                            }
                            else
                                dups++;
                            newI++;
                        }
                    }
                } catch (NumberFormatException e) { notNumber++; }
            }
            i++;
        }
        System.out.printf(" done\nRemoved %d duplicate sequences\n", dups);
        System.out.printf("Got %d number format exceptions\n", notNumber);
        return raw;
    }

    public static String [] repairProb(int seqOff, int probOff, String [] raw) {
        String [] line;
        String sequence;
        System.out.printf("Repairing probabilities over %.2f... ", PE_THRESH);
        for (int i = 1; i < raw.length; i++) {
            if (!raw[i].equals("")) {
                line = raw[i].split(",");
                if (!line[probOff].equals("")) {
                    sequence = line[seqOff];
                    double[] probs = getProbabilities(line[probOff], sequence.length());
                    if (probs != null) {
                        for (int j = 0; j < probs.length; j++) {
                            if (probs[j] >= PE_THRESH) {
                                StringBuilder strBuilder = new StringBuilder(sequence);
                                strBuilder.setCharAt(j, Character.toLowerCase(strBuilder.charAt(j)));
                                sequence = strBuilder.toString();
                            }
                        }
                        line[seqOff] = sequence;
                        String finalLine = "";
                        for (int x = 0; x < line.length; x++) {
                            if (x == line.length)
                                finalLine = finalLine + line[x];
                            else
                                finalLine = finalLine + line[x] + ",";
                        }
                        raw[i] = finalLine;
                    }
                }
            }
        }
        System.out.printf(" done\n");
        if (isoforms > 0)
            System.out.printf("Found too many isoforms (x%d)\n", isoforms);
        return raw;
    }

    public static double [] getProbabilities(String probString, int length) {
        String [] probStringSp = probString.split(";");
        double [] probs = new double[length];
        for (int i = 0; i < length; i++) {
            probs[i] = 0.0d;
        }
        for (int j = 0; j < probStringSp.length; j++) {
            String[] singleProb = probStringSp[j].split(":");
            int index = getIndex(singleProb[0]);
            double prob;
            try {
                prob = Double.parseDouble(singleProb[1]);
            } catch (NumberFormatException nfe) {
                System.out.printf("Warning: Failed to get correct probability...");
                prob = 0.0d;
            } catch (ArrayIndexOutOfBoundsException e) {
                if (singleProb[0].equals("Too many isoforms"))
                    isoforms++;
                return null;
            }
            if (index >= 0)
                probs[index] = prob;
        }
        return probs;
    }

    public static int getIndex(String indexStr) {
        indexStr = indexStr.replaceAll("\\D+","");
        int index = -1;
        try {
            index = Integer.parseInt(indexStr);
        } catch (NumberFormatException nfe) {

        }
        return index - 1; //Align it to array
    }

    public static String [] capitalise(int offset, String[] raw) {
        String [] line;
        System.out.printf("Normalising sequence data...");
        for (int i = 1; i < raw.length; i++) {
            if (!raw[i].equals("")) {
                line = raw[i].split(",");
                String sequence = line[offset];
                StringBuilder strBuilder = new StringBuilder(sequence);
                for (int j = 0; j < sequence.length(); j++) {
                    strBuilder.setCharAt(j, Character.toUpperCase(strBuilder.charAt(j)));
                }
                line[offset] = strBuilder.toString();
                String finalLine = "";
                for (int x = 0; x < line.length; x++) {
                    if (x == line.length)
                        finalLine = finalLine + line[x];
                    else
                        finalLine = finalLine + line[x] + ",";
                }
                raw[i] = finalLine;
            }
        }
        System.out.printf(" done\n");
        return raw;
    }
    public static String [] removeUnused(int offset, String [] raw) {
        int unusedCnt = 0;
        String [] line;
        for (int i = 1; i < raw.length; i++) {
            if (!raw[i].equals("")) {
                line = raw[i].split(",");
                if (line[offset].equals("Not Used")) {
                    unusedCnt++;
                    raw[i] = "";
                }
            }
        }
        System.out.printf("Found %d unused\n", unusedCnt);
        return raw;
    }

    public static String [] removeDuplicateScans(int offset, String [] raw) {
        String [] line;
        List<Integer> scans = new ArrayList<Integer>();
        int scanNo, dupCnt = 0;
        for (int i = 1; i < raw.length; i++) {
            if (!raw[i].equals("")) {
                line = raw[i].split(",");
                try {
                    scanNo = Integer.parseInt(line[offset]);
                } catch (NumberFormatException e) {
                    System.out.printf("Detected commas in data, please remove.\n");
                    return null;
                }
                if (scans.contains(scanNo)) {
                    dupCnt++;
                    raw[i] = "";
                } else
                    scans.add(scanNo);
            }
        }
        System.out.printf("Found %d duplicates thereafter\n", dupCnt);
        return raw;
    }

    public static int getOffset(String column, String[] metaData) {
        int x = 0;
        for (x = 0; x < metaData.length; x++) {
            if (metaData[x].contains(column) || metaData[x].equals(column))
                break;
        }
        return x;
    }
    public static String[] getMeta(String [] raw) {
        String [] metaData = raw[0].split(",");
        boolean correct = true;
        String missing = "";
        if (!checkColumn("Heavy/Light", metaData)) {
            correct = false;
            missing = missing + "Heavy/Light ";
        }
        if (!checkColumn("Sequence", metaData)) {
            correct = false;
            missing = missing + "Sequence ";
        }
        if (!checkColumn("phosphoRS Site Probabilities", metaData)) {
            correct = false;
            missing = missing + "phosphoRS Site Probabilities ";
        }
        if (!checkColumn("Quan Usage", metaData)) {
            correct = false;
            missing = missing + "Quan Usage ";
        }
        if (!checkColumn("First Scan", metaData)) {
            missing = missing + "First Scan ";
            if (!checkColumn("Last Scan", metaData)) {
                correct = false;
                missing = missing + "Last Scan ";
            }
        }
        if (correct)
            return metaData;
        else {
            System.out.printf("Can't find the following columns: %s\n", missing);
            return null;
        }
    }

    public static boolean checkColumn(String column, String[] metaData) {
        for (int i = 0; i < metaData.length; i++) {
            if (metaData[i].equals(column))
                return true;
        }
        return false;
    }

    public static void writeOut(String [] dataFinal) {
        File underFile = new File(newDestinationUnder);
        underFile.getParentFile().mkdirs();

        File overFile = new File(newDestinationOver);
        overFile.getParentFile().mkdirs();

        PrintWriter writeUnder = null;
        try {
            writeUnder = new PrintWriter(underFile);
        } catch (IOException e) {
            e.printStackTrace();
        }

        PrintWriter writeOver = null;
        try {
            writeOver = new PrintWriter(overFile);
        } catch (IOException e) {
            e.printStackTrace();
        }

        int offset = getOffset("Sequence", getMeta(dataFinal));

        for (int i = 0; i < dataFinal.length; i++) {
            if (!dataFinal[i].equals("")) {
                String[] line = dataFinal[i].split(",");
                String seq = line[offset];
                boolean over = false;
                for (int j = 0; j < seq.length(); j++) {
                    if (Character.isLowerCase(seq.charAt(j)))
                        over = true;
                }

                //Here choose between under and over
                if (over)
                    writeOver.println(dataFinal[i]);
                else
                    writeUnder.println(dataFinal[i]);
            }
        }
        writeUnder.close();
        writeOver.close();
    }

    public static String[] readIn(String destination) {
        BufferedReader br = null;
        int line = 0;
        try {
            String sCurrentLine;
            br = new BufferedReader(new FileReader(destination));
            while ((sCurrentLine = br.readLine()) != null) {
                line++;
            }
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            try {
                if (br != null)br.close();
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }
        if (line == 0)
            return null;

        String [] raw = new String[line];
        br = null;
        int i = 0;
        try {
            String sCurrentLine;
            br = new BufferedReader(new FileReader(destination));
            while ((sCurrentLine = br.readLine()) != null) {
                raw[i] = sCurrentLine;
                i++;
            }
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            try {
                if (br != null)br.close();
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }
        try {
            br.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return raw;
    }

    @Override
    public void run() {
        fix();
    }
}
