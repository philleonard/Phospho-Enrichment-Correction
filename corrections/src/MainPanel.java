//import build.tools.javazic.Main;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

/**
 * Created by ppl on 15/06/15.
 */
public class MainPanel extends JFrame {
    private JButton runButton;
    private JTextField textField1;
    private JButton openButton;
    private JTextPane textPane1;
    private JPanel mainPanel;
    private JFormattedTextField a7500FormattedTextField;
    private JCheckBox removeUnusedCheckBox;
    private JCheckBox removeDupsCheckBox;
    private JCheckBox deleteDupSeqCheckBox;
    static MainPanel frame;
    final JFileChooser fc = new JFileChooser();
    double thresh = -1;

    static String destination = "";
    static String newDestination = "";

    public String getNewFile(boolean over) {
        String str = textField1.getText().toString();
        str = str.substring(0, str.lastIndexOf('.'));
        if (over)
            str = str + "_fixed_over" + thresh + ".csv";
        else
            str = str + "_fixed_under" + thresh + ".csv";
        return str;
    }

    public MainPanel() {
        super("Phospho Enrichment Sequence Corrector v0.1");
        setContentPane(mainPanel);
        pack();
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setBounds(100, 100, 594, 320);

        runButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent actionEvent) {

                boolean remDups = removeDupsCheckBox.isSelected();
                boolean remUn = removeUnusedCheckBox.isSelected();
                boolean remDupSeq = deleteDupSeqCheckBox.isSelected();
                thresh = -1;
                try {
                    thresh = Double.parseDouble(a7500FormattedTextField.getText().toString());
                } catch (NumberFormatException e) {
                    System.err.printf("Threshold not number\n");
                }
                if (thresh >= 0.0 && thresh <= 100.00) {
                    if (!textField1.getText().toString().equals("")) {
                        File file = new File(textField1.getText().toString());
                        if (file.exists()) {
                            Thread cl = new Thread(new Correction(textField1.getText().toString(), getNewFile(true), getNewFile(false), remDups, remUn, remDupSeq, thresh, frame));
                            cl.start();
                            runButton.setEnabled(false);
                        }
                        else
                            System.err.printf("File does not exist\n");
                    }
                    else
                        System.err.printf("No file specified\n");
                }
                else
                    System.err.printf("Threshold not in range 0.00 to 100.00\n");
            }
        });

        openButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent actionEvent) {
                int returnVal = fc.showOpenDialog(MainPanel.this);

                if (returnVal == JFileChooser.APPROVE_OPTION) {
                    File file = fc.getSelectedFile();
                    String extension = Utils.getExtension(file);
                    if (extension != null) {
                        if (extension.equals(Utils.csv)) {
                            textField1.setText(file.getAbsoluteFile().toString());
                            System.err.printf("Is it sorted alphabetically on sequence?\n");
                        }
                        else
                            System.err.printf("Not CSV file\n");
                    }
                    //System.out.printf("\n%s", file.getAbsoluteFile());

                } else {
                    System.out.printf("No file chosen\n");
                }

            }
        });

        MessageConsole mc = new MessageConsole(textPane1);
        mc.redirectOut();
        mc.redirectErr(Color.RED, null);
        mc.setMessageLines(100);
    }

    public static void main(String[] args) {
        EventQueue.invokeLater(new Runnable() {
            public void run() {
                try {
                    frame = new MainPanel();
                    frame.setVisible(true);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        });
        /*Correction co = new Correction();
        co.fix(destination, newDestination);*/
    }


    public JButton getRunButton() {
        return runButton;
    }
}
