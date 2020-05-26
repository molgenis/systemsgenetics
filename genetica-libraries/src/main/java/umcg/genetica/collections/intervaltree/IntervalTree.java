package umcg.genetica.collections.intervaltree;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

/**
 *
 * @author Patrick Deelen
 */
public class IntervalTree<E extends Range> {

	private final IntervalTreeNode<E> rootNode;
	private final int size;
	private final Class<E> classE;

	public IntervalTree(final E[] elements, Class<E> classE) {

		//This does not make a copy of the array. 
		this(Arrays.asList(elements), classE);

	}

	public IntervalTree(final List<E> elements, Class<E> classE) {

		for (E e : elements) {
			if (e.getStart() > e.getEnd()) {
				throw new IllegalArgumentException("Start is greater than end. This is not allowed in the interval tree.");
			}
		}
		this.classE = classE;

		rootNode = createNode(elements);
		size = elements.size();
	}

	public ArrayList<E> getElementsOverlappingQuery(final int query) {
		final ArrayList<E> results = new ArrayList<E>();

		rootNode.queryNode(results, query);

		return results;

	}

	;
	
	
	private IntervalTreeNode<E> createNode(final List<E> elements) {

		final ArrayList<E> left = new ArrayList<E>();
		final ArrayList<E> right = new ArrayList<E>();
		final ArrayList<E> center = new ArrayList<E>();

		final E centerElement = elements.get(elements.size()/2);
		final int centerPoint = centerElement.getStart();

//		System.out.println("--------");
//		System.out.println("Elements " + elements.size());
//		System.out.println("Center point " + centerPoint);
		
		for (E element : elements) {
			if (element.getEnd() < centerPoint) {
				left.add(element);
			} else if (element.getStart() > centerPoint) {
				right.add(element);
			} else {
				center.add(element);
			}
		}
		
//		System.out.println("Left size " + left.size());
//		System.out.println("Right size " + right.size());
//		System.out.println("Center size " + center.size());
//		System.out.println("--------");
		
		final IntervalTreeNode<E> leftNode = left.isEmpty() ? null : createNode(left);
		final IntervalTreeNode<E> rightNode = right.isEmpty() ? null : createNode(right);

		//We copy ArrayList of center into array since we will have many nodes and need this list twice
		//Due to different sorting 
		@SuppressWarnings("unchecked")
		E[] centerArray = (E[]) Array.newInstance(classE, center.size());
		return new IntervalTreeNode<E>(centerPoint, leftNode, rightNode, center.toArray(centerArray));

	}

	public int size() {
		return size;
	}

	public boolean isEmpty() {
		return size == 0;
	}

}
