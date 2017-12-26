package nl.umcg.westrah.binarymetaanalyzer;

import java.text.NumberFormat;

public class MultiThreadProgressBar {
	
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
	}
	
	public void setSubtasks(int task, int subtaskNr) {
		subtasks[task] = subtaskNr;
	}
	
	
	public synchronized void display() {
		String outln = "\rProgress = ";
	
		for (int q = 0; q < ctrs.length; q++) {
			if (q == 0) {
				outln += " t" + q + ": ";
			} else {
				outln += ", t" + q + ": ";
			}
			
			if (subtasks[q] == -1) {
				outln += "NO START";
			} else if (complete[q]) {
				outln += "COMPLETE";
			} else {
				double perc = (double) ctrs[q] / subtasks[q];
				outln += defaultFormat.format(perc);
			}
		}
		
		// calculate approximate ETA
		int totalApproximateTasks = 0;
		
		System.out.print(outln);
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
