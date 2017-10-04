package deconvolution;

import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.Callable;

// from: http://stackoverflow.com/a/4010275/651779


import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;


public class Parallel {
    public static ExecutorService forPool;
    public Parallel(int numCores){
    	final int NUM_CORES = numCores;
    	forPool = Executors.newFixedThreadPool(NUM_CORES * 2, new NamedThreadFactory("Parallel.For"));
    }
    public <T> void For(final Iterable<T> elements, final Operation<T> operation) {
        try {
            // invokeAll blocks for us until all submitted tasks in the call complete
            forPool.invokeAll(createCallables(elements, operation));
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

    public static <T> Collection<Callable<Void>> createCallables(final Iterable<T> elements, final Operation<T> operation) {
        List<Callable<Void>> callables = new LinkedList<Callable<Void>>();
        for (final T elem : elements) {
            callables.add(new Callable<Void>() {
                @Override
                public Void call() {
                    operation.perform(elem);
                    return null;
                }
            });
        }

        return callables;
    }

    public static interface Operation<T> {
        public void perform(T pParameter);
    }
}