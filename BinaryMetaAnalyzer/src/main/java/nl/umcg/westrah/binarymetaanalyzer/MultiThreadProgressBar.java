package nl.umcg.westrah.binarymetaanalyzer;

import com.google.common.base.Strings;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.util.RunTimer;

import java.text.NumberFormat;

public class MultiThreadProgressBar {
	
	private final RunTimer timer;
	int[] ctrs;
	int[] subtasks;
	boolean[] complete;
	NumberFormat defaultFormat = NumberFormat.getPercentInstance();
	
	public MultiThreadProgressBar(int tasks) {
		this.ctrs = new int[tasks];
		this.subtasks = new int[tasks];
		for (int q = 0; q < tasks; q++) {
			subtasks[q] = -1;
		}
		complete = new boolean[tasks];
		defaultFormat.setMinimumFractionDigits(1);
		timer = new RunTimer();
		timer.start();
	}
	
	public void setSubtasks(int task, int subtaskNr) {
		subtasks[task] = subtaskNr;
	}
	
	int maxlen = 0;
	
	public synchronized void display() {
		StringBuilder outln = new StringBuilder().append("\rProgress = ");
		
		for (int q = 0; q < ctrs.length; q++) {
			if (q == 0) {
				outln.append(" t" + q + ": ");
			} else {
				outln.append(", t" + q + ": ");
			}
			
			if (subtasks[q] == -1) {
				outln.append("NO START");
			} else if (complete[q]) {
				outln.append("COMPLETE");
			} else {
				double perc = (double) ctrs[q] / subtasks[q];
				outln.append(defaultFormat.format(perc));
			}
		}
		
		// calculate approximate ETA
		int totalApproximateTasks = 0;
		// calculate average nr of subtasks
		double nrSubtasksTotalRegistered = 0;
		double nrSubtasksComplete = 0;
		int nrtasksStarted = 0;
		for (int q = 0; q < ctrs.length; q++) {
			if (subtasks[q] > -1) {
				nrSubtasksTotalRegistered += subtasks[q];
				nrSubtasksComplete += ctrs[q];
				nrtasksStarted++;
			}
		}
		
		if (nrtasksStarted > 0) {
			double totalExpectedTasks = ctrs.length * (nrSubtasksTotalRegistered / nrtasksStarted);
			long timediff = timer.getTimeDiff() / 1000000000;
			double timePerIter = (double) timediff / nrSubtasksComplete;
			double timeLeft = timePerIter * (totalExpectedTasks - nrSubtasksComplete);
			double averagePerc = nrSubtasksComplete / totalExpectedTasks;
			outln.append(" | ").append(defaultFormat.format(averagePerc)).append(" | T: ").append(timer.getTimeDesc());
			if (timeLeft > 0) {
				String strTimeLeft = timer.getTimeDesc(((long) timeLeft) * 1000000000);
				outln.append(" T-: ").append(strTimeLeft);
			}
			
		}
		String outstr = outln.toString();
		if (outstr.length() > maxlen) {
			maxlen = outstr.length();
		} else {
			int diff = maxlen - outstr.length();
			String whitespace = Strings.repeat(" ", diff);
			outln.append(whitespace);
		}
		System.out.print(outln.toString());
		
	}
	
	public void iterate(int taskid) {
		ctrs[taskid]++;
	}
	
	public void set(int taskid, int snp) {
		ctrs[taskid] = snp;
	}
	
	public boolean allCompleted() {
		for (boolean b : complete) {
			if (!b) {
				return false;
			}
		}
		return true;
	}
	
	public void complete(int taskid) {
		complete[taskid] = true;
	}
}
