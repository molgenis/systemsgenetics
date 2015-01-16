package eqtlmappingpipeline.ase;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

/**
 *
 * @author Patrick Deelen
 */
public class AseCalculator implements Runnable {

	private final AseVariant[] aseVariants;
	private final AtomicInteger counter;	

	public AseCalculator(AseVariant[] aseVariants, AtomicInteger counter) {
		this.aseVariants = aseVariants;
		this.counter = counter;
	}
	
	@Override
	public void run() {

		int i;
		while( (i = counter.getAndIncrement()) < aseVariants.length){
			aseVariants[i].calculateStatistics();
		}
	}

	public static void startAseCalculators(AseVariant[] aseVariants, int threadCount) {

		List<Thread> threads = new ArrayList<Thread>(threadCount);
		final Ase.ThreadErrorHandler threadErrorHandler = new Ase.ThreadErrorHandler();

		final AtomicInteger count = new AtomicInteger(0);
		
		for (int i = 0; i < threadCount; ++i) {

			AseCalculator aseCalculator = new AseCalculator(aseVariants, count);
			Thread worker = new Thread(aseCalculator);
			worker.setUncaughtExceptionHandler(threadErrorHandler);
			worker.start();
			threads.add(worker);

		}

		int nextReport = 1000;
		boolean running;
		do {
			running = false;
			for (Thread thread : threads) {
				if (thread.isAlive()) {
					running = true;
				}
			}

			int currentProgress = count.get();

			if (currentProgress > nextReport) {
				System.out.println("Calculated " + currentProgress + " / " + aseVariants.length + " ASE variants");
				nextReport += 1000;
			}


			try {
				Thread.sleep(500);
			} catch (InterruptedException ex) {
			}

		} while (running);

	}
}
