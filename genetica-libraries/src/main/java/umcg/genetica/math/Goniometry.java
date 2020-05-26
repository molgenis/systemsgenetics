/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.math;

import umcg.genetica.containers.Pair;

/**
 *
 * @author harmjan
 */
public class Goniometry {
    /**
     * @param r radius of the circle
     * @param originX the origin (center of the circle) X coordinate
     * @param originY the origin (center of the circle) Y coordinate
     * @param angle the angle (in degrees)
     * @return Pair<Integer, Integer> position in pixels (x,y)
     */
    public static Pair<Integer, Integer> calcPosOnCircle(double r, double originX, double originY, double angle) {

        // x = cx + r * cos(a)
        // y = cy + r * sin(a)
        double radians = angle * Math.PI / 180;
        Integer x = (int) Math.floor(originX + (r * Math.cos(radians)));
        Integer y = (int) Math.floor(originY + (r * Math.sin(radians)));

        return new Pair<Integer, Integer>(x, y);
    }
}
