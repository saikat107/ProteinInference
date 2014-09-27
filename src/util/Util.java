package util;


import java.util.HashSet;
import java.util.Set;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author define
 */
public class Util {
        public static final int num_gen = 100;
        public static final int num_off = 100;
        public static final int max_unimproved = num_gen/5;
        
        
        public static  int union(Set<Integer> setA, Set<Integer> setB) {
            Set<Integer> tmp = new HashSet<Integer>(setA);
            tmp.addAll(setB);
            return tmp.size();
        }

  public static int intersection(Set<Integer> setA, Set<Integer> setB) {
    Set<Integer> tmp = new HashSet<Integer>();
    for (Integer x : setA)
      if (setB.contains(x))
        tmp.add(x);
    return tmp.size();
  }

  public static  int difference(Set<Integer> setA, Set<Integer> setB) {
    Set<Integer> tmp = new HashSet<Integer>(setA);
    tmp.removeAll(setB);
    return tmp.size();
  }
        

}
