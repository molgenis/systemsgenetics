/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.util;

/**
 *
 * @author harm-jan
 */
public class RunTimer {

    private long start = System.nanoTime();
    private long prevstart;
    private Granularity g = null;

    public enum Granularity {

        MICROS, SEC, MILLIS, MINS, HOURS
    };

    public RunTimer(Granularity g) {
        prevstart = start;
        this.g = g;
    }

    public RunTimer() {
        prevstart = start;
    }

    public String getTimeDesc() {
        long stopTime = System.nanoTime();
        long diff = stopTime - start;
        return getTimeDesc(diff);
    }

    public void start() {
        prevstart = System.nanoTime();

    }

    public String stop() {
        long currtime = System.nanoTime();
        long diff = currtime - prevstart;
        return getTimeDesc(diff);
    }

    public String getTimeDesc(long timeDiff) {
        String leftTime = "";
        String pastTime = "";
        double surplus = 0.0;

        if (g == null) {

            if (timeDiff < 1000) {
                pastTime = timeDiff + " ns";
            } else if (timeDiff < 1000000) {
                surplus = timeDiff % 1000000;
                surplus /= 1000;
                pastTime = (int) Math.floor((double) timeDiff / 1000000) + " micros " + (int) Math.floor(surplus) + " ns";
            } else if (timeDiff < 1000000000) {
                surplus = timeDiff % 1000000000;
                surplus /= 1000000;
                pastTime = (int) Math.floor((double) timeDiff / 1000000000) + " ms " + (int) Math.floor(surplus) + " micros";
            } else if (timeDiff > 1000000000) {
                int time = (int) Math.ceil((double) (timeDiff) / 1000000000);
                if (time < 60) {
                    pastTime = "" + (int) Math.floor(time) + "s";
                } else if (time >= 60 && time < 3600) {
                    surplus = time % 60;
                    time /= 60;
                    pastTime = "" + (int) Math.floor(time) + "m " + (int) Math.floor(surplus) + " s";
                } else if (time >= 3600) {
                    surplus = time % 3600;
                    surplus /= 60;
                    time /= 3600;
                    pastTime = "" + (int) Math.floor(time) + "h " + (int) Math.floor(surplus) + " m";
                }
            }
        } else {
            // display time in smallest unit required.
            switch (g) {
                case HOURS:
                    break;
                case MINS:
                    break;
                case SEC:
                    break;
                case MILLIS:
                    break;
                case MICROS:
                    break;
            }
        }
        return pastTime;

    }

    public long getTimeDiff() {

        return System.nanoTime() - prevstart;
    }
}
