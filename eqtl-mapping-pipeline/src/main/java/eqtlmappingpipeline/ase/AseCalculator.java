package eqtlmappingpipeline.ase;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author Patrick Deelen
 */
public class AseCalculator implements Runnable {

	private final int first;
	private final int last;
	private final AseVariant[] aseVariants;
	private int i;

	public AseCalculator(int first, int last, AseVariant[] aseVariants) {
		this.first = first;
		this.last = last;
		this.aseVariants = aseVariants;
	}

	@Override
	public void run() {

		for (i = first; i < last; ++i) {
			aseVariants[i].calculateStatistics();
		}
	}

	public int getProgress() {
		return i - first;
	}

	public static void startAseCalculators(AseVariant[] aseVariants, int threadCount) {

		List<Thread> threads = new ArrayList<Thread>(threadCount);
		List<AseCalculator> aseCalculators = new ArrayList<AseCalculator>(threadCount);
		final Ase.ThreadErrorHandler threadErrorHandler = new Ase.ThreadErrorHandler();

		int tasksPerThread = aseVariants.length / threadCount;
		if (aseVariants.length % threadCount != 0) {
			tasksPerThread++;
		}

		for (int i = 0; i < threadCount; ++i) {

			int first = i * tasksPerThread;
			int last = first + tasksPerThread;
			if (last > aseVariants.length) {
				last = aseVariants.length;
			}

			AseCalculator aseCalculator = new AseCalculator(first, last, aseVariants);
			Thread worker = new Thread(aseCalculator);
			worker.setUncaughtExceptionHandler(threadErrorHandler);
			worker.start();
			threads.add(worker);
			aseCalculators.add(aseCalculator);

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

			int currentProgress = 0;
			for (AseCalculator aseCalculator : aseCalculators) {
				currentProgress += aseCalculator.getProgress();
			}

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
