package protinInference;




import java.util.Random;


/**
 * Artificial Intelligence A Modern Approach (3rd Edition): Figure 4.5, page
 * 126.<br>
 * <br>
 *
 * <pre>
 * function SIMULATED-ANNEALING(problem, schedule) returns a solution state
 *                    
 *   current &lt;- MAKE-NODE(problem.INITIAL-STATE)
 *   for t = 1 to INFINITY do
 *     T &lt;- schedule(t)
 *     if T = 0 then return current
 *     next &lt;- a randomly selected successor of current
 *     /\E &lt;- next.VALUE - current.value
 *     if /\E &gt; 0 then current &lt;- next
 *     else current &lt;- next only with probability e&circ;(/\E/T)
 * </pre>
 *
 * Figure 4.5 The simulated annealing search algorithm, a version of stochastic
 * hill climbing where some downhill moves are allowed. Downhill moves are
 * accepted readily early in the annealing schedule and then less often as time
 * goes on. The schedule input determines the value of the temperature T as a
 * function of time.
 *
 * @author Ravi Mohan
 * @author Mike Stampone
 */


public class SimulatedAnnealing  
{
        ProteinInference pI;
        private Scheduler scheduler;
        private Random rand;
        
        public SimulatedAnnealing(ProteinInference p) 
        {         
                pI = p;
        	scheduler = new Scheduler();
        	rand = new Random();
        }

        
        


    	public void improve(int [] initialNode) 
    	{   		
    		int []current;
                int []next;
    		current = new int[initialNode.length];
                System.arraycopy(initialNode, 0, current, 0, initialNode.length);
    		next = new int[initialNode.length];
                
            int timeStep = 0;
            
            do{
                    //System.out.println(pI.evaluateFitness(current, true));
	            double temperature = scheduler.getTemp(timeStep);
	            timeStep++;
                    if (temperature <= 0.00000000001){
	                    break;
	            }
	            System.arraycopy(current, 0, next, 0, initialNode.length);
                    applyMutation(next);
                    double currentFitness = pI.evaluateFitness(current, true);
                    double nextFitness = pI.evaluateFitness(next, true);
	            
                    double deltaE = nextFitness - currentFitness;
	
	            if (shouldAccept(temperature, deltaE)){
	                    current = next;
	            }
	            
            }while(true);
            if(pI.evaluateFitness(initialNode, true)< pI.evaluateFitness(next, true))
            	System.arraycopy(next, 0, initialNode, 0, initialNode.length);

    	}
		
        // if /\E > 0 then current <- next
        // else current <- next only with probability e^(/\E/T)
        private boolean shouldAccept(double temperature, double deltaE) 
        {
                return (deltaE > 0.0)
                                || (rand.nextDouble() <= probabilityOfAcceptance(
                                                temperature, deltaE));
        }

        
        public double probabilityOfAcceptance(double temperature, double deltaE) {
                return Math.exp(deltaE / temperature);
        }

		void applyMutation(int [] offspring)
		{
			
			int index = rand.nextInt(offspring.length);
                        offspring[index] = 1-offspring[index];

		}
}

class Scheduler 
{

        private final int k, limit;

        private final double lam;

        public Scheduler(int k, double lam, int limit) {
                this.k = k;
                this.lam = lam;
                this.limit = limit;
        }

        public Scheduler() 
        {
                this.k = 20;
                this.lam = 0.045;
                this.limit = 750;
        }

        public double getTemp(int t) {
                if (t < limit)
                {
                        double res = k * Math.exp((-1) * lam * t);
                        return res;
                } 
                else 
                {
                        return 0.0;
                }
        }
}



