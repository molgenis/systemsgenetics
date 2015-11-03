/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.chrContacts;

/**
 *
 * @author MaKKie_Admin
 */
public class ChrContact implements Comparable<ChrContact> {

    final int chrLocationSmaller;
    final int chrLocationLarger;
    final double contactValue;

    public ChrContact(int chrLocSmal, int chrLocLarge, double contactVal) {
        this.chrLocationLarger = chrLocLarge;
        this.chrLocationSmaller = chrLocSmal;
        this.contactValue = contactVal;
    }

    @Override
    public int compareTo(ChrContact other) {
        if (other.getChrLocationSmaller() > this.chrLocationSmaller) {
            return -1;
        } else if (other.getChrLocationSmaller() < this.chrLocationSmaller) {
            return 1;
        } else {
            if (other.getChrLocationLarger() > this.chrLocationLarger) {
                return -1;
            } else if (other.getChrLocationLarger() < this.chrLocationLarger) {
                return 1;
            } else {
                return 0;
            }
        }
    }

    public int getChrLocationSmaller() {
        return chrLocationSmaller;
    }

    public int getChrLocationLarger() {
        return chrLocationLarger;
    }

    public double getContactValue() {
        return contactValue;
    }
}
