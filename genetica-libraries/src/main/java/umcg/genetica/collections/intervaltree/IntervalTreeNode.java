package umcg.genetica.collections.intervaltree;

import java.util.ArrayList;
import java.util.Arrays;

/**
 *
 * @author Patrick Deelen
 */
class IntervalTreeNode<E extends Range> {

	private final int centerPoint;
	private final IntervalTreeNode left;
	private final IntervalTreeNode right;
	private final E[] startSorted;//sorted small to large
	private final E[] endSorted;//sorted large to small
	private static final RangeStartComparator startComparator = new RangeStartComparator();
	private static final RangeEndReverseComparator endReverseComparator = new RangeEndReverseComparator();

	protected IntervalTreeNode(int centerPoint, IntervalTreeNode left, IntervalTreeNode right, E[] elements) {
		this.centerPoint = centerPoint;
		this.left = left;
		this.right = right;

		this.startSorted = elements;
		this.endSorted = Arrays.copyOf(elements, elements.length);

		Arrays.sort(startSorted, startComparator);
		Arrays.sort(endSorted, endReverseComparator);

	}

	protected void queryNode(final ArrayList<E> results, final int query) {

		if (query == centerPoint) {
			results.ensureCapacity(results.size() + startSorted.length);
			for (E e : startSorted) {
				results.add(e);
			}
		} else if (query < centerPoint) {
			for (E e : startSorted) {
				if (e.getStart() <= query) {
					results.add(e);
				} else {
					break;
				}
			}
			if (left != null) {
				left.queryNode(results, query);
			}
		} else {
			//query > centerPoint
			for (E e : endSorted) {
				if (e.getEnd() >= query) {
					results.add(e);
				} else {
					break;
				}
			}
			if (right != null) {
				right.queryNode(results, query);
			}
		}

	}
}
