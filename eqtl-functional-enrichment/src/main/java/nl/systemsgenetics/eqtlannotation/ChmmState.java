/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.eqtlannotation;

/**
 *
 * @author MarcJan
 */
class ChmmState {
    private final long start;
    private final long stop;

    ChmmState(String[] parts) {
        start = Long.parseLong(parts[1]);
        stop = Long.parseLong(parts[2]);
    }

    public long getStart() {
        return start;
    }

    public long getStop() {
        return stop;
    }
    
}
