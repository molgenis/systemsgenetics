/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.console;

import umcg.genetica.util.RunTimer;

/**
 * @author harmjan
 */
public class ProgressBar {

    private RunTimer timer;
    private long iterations = 0;
    private long maxIterations = 0;
    private int printEvery = 0;
    private int width = 25;
    private int maxwidth = 0;
    private boolean isLong = false;

    public ProgressBar(long length) {
        iterations = 0;

        isLong = true;
        maxIterations = length;
        printEvery = (int) Math.ceil((double) length / width);
        timer = new RunTimer();

        System.out.println("\nProgress:");

        String out = "|";
        for (int i = 0; i < width - 1; i++) {
            out += "-";
        }
        out += "| Waiting for update...\r";
        maxwidth = out.length();
        System.out.print(out);
    }

    public ProgressBar(long length, String title) {
        iterations = 0;
        maxIterations = length;
        printEvery = (int) Math.ceil((double) length / width);
        timer = new RunTimer();

        System.out.println("\n" + title);

        String out = "|";
        for (int i = 0; i < width - 1; i++) {
            out += "-";
        }
        out += "| Waiting for update...\r";

        maxwidth = out.length();
        System.out.print(out);
    }

    public void iterate() {
        iterations++;
        if (iterations % printEvery == 0) {
            print();
        }
    }

    public void set(long num) {
        iterations = num;
        print();
    }

    public void print() {

        if (printEvery > 0) {
            StringBuilder out = new StringBuilder(maxwidth);
            int numToPrint = (int) Math.ceil(iterations / (double) printEvery);
            if (numToPrint > width) {
                numToPrint = width;
            }
            out.append("|");
            for (int i = 0; i < numToPrint; i++) {
                out.append("#");
            }
            for (int i = 0; i < width - numToPrint - 1; i++) {
                out.append("-");
            }
            out.append("| ");
            int perc = (int) Math.ceil((double) iterations / maxIterations * 100);

            out.append(perc).append("% - T: ").append(timer.getTimeDesc());

            long diff = timer.getTimeDiff() / 1000000000;
            double timePerIter = (double) diff / iterations;
            double timeLeft = timePerIter * (maxIterations - iterations);
            String strTimeLeft = timer.getTimeDesc(((long) timeLeft) * 1000000000);
            out.append(" T-: ").append(strTimeLeft);
            out.append(" #: ").append(iterations).append("/").append(maxIterations);


            int length = out.toString().length();
            if (length < maxwidth) {
                int difflen = maxwidth - out.length();
                if (difflen > 0) {
                    int i = 0;
                    while (i < difflen) {
                        out.append(" ");
                        i++;
                    }
                }
            } else {
                maxwidth = length;
            }


            out.append("\r");
            System.out.print(out.toString());
        }

    }

    public void close() {
        iterations = maxIterations;
        print();
        System.out.println("");

    }

    public synchronized void iterateSynched() {
        iterations++;
        if (iterations % printEvery == 0) {
            print();
        }
    }
}
