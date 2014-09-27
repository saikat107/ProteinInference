package protinInference;


import util.Configuration;
import util.Util;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;

import java.io.FileReader;
import java.io.FileWriter;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;
import java.util.Scanner;
import java.util.StringTokenizer;
import java.util.logging.Level;
import java.util.logging.Logger;

public final class ProteinInference {
    static Scanner parameterScan;
    static Scanner parameterScan2;
    static PrintWriter parameterWriter;
    
    
    ArrayList<Integer> candidateProteins;
    ArrayList<Integer> testPeptides;
    int [][]population;
    double []fitness;
    int populationSize=100;
    static double beta=0.23;
    static double epsilon = 0.55;
    ArrayList<Integer>peptides;
    int total_iteration=0;

    ArrayList<Integer> actualProteinsList = new ArrayList<Integer>();
    Random randomGenerator = new Random();
    int [][]parentOffspringPopulation = new int[populationSize+Util.num_off][];
    SimulatedAnnealing sim;
    HashMap<Integer, ArrayList> protenPeptideMap = new HashMap<Integer, ArrayList>();
    ProteinInference(){
        sim = new SimulatedAnnealing(this);
        try{
                 
            FileWriter outputWriter = new FileWriter(Configuration.outputData,true);
            BufferedWriter outputBufferWriter = new BufferedWriter(outputWriter);
            
            FileWriter exampleWriter = new FileWriter(Configuration.outputIntermediate,true);
            BufferedWriter exampleBufferWriter = new BufferedWriter(exampleWriter);
            outputBufferWriter.flush();
            for(int test_case=1;test_case<=1;test_case++){

                double avgPrecision=0;
                double avgRecall=0;
                double avgFmeans=0;
                String testString;
                double avgAccuracy=0.0;
                double avgRedundancy=0.0;
                double avgConsistancy=0.0;
                double avgTP=0;
                double avgFP=0;
                double avgFN=0;
                FileReader testReader = new FileReader(Configuration.inputData);
                BufferedReader testBufferReader=new BufferedReader(testReader);
                int testExampleNo=0;
                while(true){
                    actualProteinsList = new ArrayList<Integer>();
                    testString = testBufferReader.readLine();
                    if(testString==null) break;

                    StringTokenizer tokens = new StringTokenizer(testString,":");
                    String ActualProteinsString=tokens.nextToken();
                    StringTokenizer actualProteinTokens = new StringTokenizer(ActualProteinsString,",");
                    while(actualProteinTokens.hasMoreTokens()){
                        actualProteinsList.add(Integer.parseInt(actualProteinTokens.nextToken()));
                    }
                    String PeptideSequence=tokens.nextToken();
                    int []result=this.geneticAlgorithm(PeptideSequence);                    
                    
                    for(int aa=0;aa<this.population.length;aa++){
                        result=this.population[aa];
                    ArrayList<Integer> inferredProteinsList = new ArrayList<Integer>();
                    for(int i=0; i<result.length; i++){
                        if(result[i]==1){
                            inferredProteinsList.add(this.candidateProteins.get(i));
                        }
                    }

                    int captured=0, notCaptured=0, extraCaptured=0;

   
                    for(int i=0; i < inferredProteinsList.size(); i++){
                        if(actualProteinsList.contains(inferredProteinsList.get(i))){
                            captured++;
                        }
                        else{
                            extraCaptured++;
                        }
                    }
                    for(int i=0; i < actualProteinsList.size(); i++){
                        if(inferredProteinsList.contains(actualProteinsList.get(i))==false){
                            notCaptured++;
                        }
                    }
                    testExampleNo++;
                    double Accuracy=captured/(double) (captured+notCaptured);
                    double Redundancy=extraCaptured/(double) (captured+extraCaptured);
                    double Consistancy=captured/(double) (captured+notCaptured+extraCaptured);
                    double reCall=Accuracy;
                    double precision=captured/(double)(captured+extraCaptured);
                    double fMeans=2*precision*reCall/(precision+reCall);
                    avgPrecision+=precision;
                    avgRecall+=reCall;
                    avgFmeans+=fMeans;
                    avgTP+=captured;
                    avgFN+=notCaptured;
                    avgFP+=extraCaptured;
                    System.out.println("\nTest Example: "+aa+" TP->"+captured+" FP->"+extraCaptured+" FN->"+notCaptured+" Precision->"+String.format("%.2f", precision)+" Recall->"+String.format("%.2f",reCall)+" FMeans->"+String.format("%.2f", fMeans));
                    
                    exampleBufferWriter.write(testExampleNo+","+captured+","+extraCaptured+","+notCaptured+","+String.format("%.2f", precision)+","+String.format("%.2f",reCall)+","+String.format("%.2f", fMeans)+"\n");
                    exampleBufferWriter.flush();
                    avgAccuracy+=Accuracy;
                    avgRedundancy+=Redundancy;
                    avgConsistancy+=Consistancy;
                    }
                }                
                avgPrecision/=testExampleNo;
                avgRecall/=testExampleNo;
                avgFmeans/=testExampleNo;
                avgAccuracy/=testExampleNo;
                avgRedundancy/=testExampleNo;
                avgConsistancy/=testExampleNo;
                avgTP/=testExampleNo;
                avgFN/=testExampleNo;
                avgFP/=testExampleNo;
                System.out.println(beta+"  TP->"+String.format("%.5f", avgTP)+" FP->"+String.format("%.5f",avgFP)+" FN->"+String.format("%.5f",avgFN)+" Precision->"+String.format("%.5f", avgPrecision)+" Recall->"+String.format("%.5f",avgRecall)+" FMeans->"+String.format("%.5f", avgFmeans));
                outputBufferWriter.write(total_iteration+","+epsilon+","+String.format("%.2f", avgTP)+","+String.format("%.2f",avgFP)+","+String.format("%.2f",avgFN)+","+String.format("%.2f", avgPrecision)+","+String.format("%.2f",avgRecall)+","+String.format("%.2f", avgFmeans)+"\n");
                parameterWriter.println(epsilon+","+beta+","+avgTP+","+avgFP+","+avgFN+","+avgPrecision+","+avgRecall+","+avgFmeans);
                parameterWriter.flush();
                testBufferReader.close();
                outputBufferWriter.flush();
            }
            outputBufferWriter.close();
            exampleBufferWriter.close();

        }catch(Exception e){
            e.printStackTrace();
        }    
    }

    

    int[] geneticAlgorithm(String testString){
        
        int iteration=0;
        this.candidateProteins = new ArrayList<Integer>();
        this.testPeptides = new ArrayList<Integer>();
        try{
            FileReader proteinPeptideRelationReader 
                    = new FileReader(Configuration.proteinPeptideMap);
            
            BufferedReader proteinPeptideRelationBufferReader
                    =new BufferedReader(proteinPeptideRelationReader);
            
            StringTokenizer tokens = new StringTokenizer(testString,",");
            
            while(tokens.hasMoreTokens()){
                testPeptides.add(Integer.parseInt(tokens.nextToken()));
            }

            while(true){
                    String s=proteinPeptideRelationBufferReader.readLine();
                    if(s==null) break;
                    tokens = new StringTokenizer(s,", ");
                    int protein=Integer.parseInt(tokens.nextToken());
                    int peptide=Integer.parseInt(tokens.nextToken());
                    ArrayList<Integer> peptides = protenPeptideMap.get(protein);
                    if(peptides==null){
                        peptides = new ArrayList<Integer>();
                    }
                    peptides.add(peptide);
                    protenPeptideMap.put(protein, peptides);
                    if(this.testPeptides.contains(peptide)){
                        if(this.candidateProteins.contains(protein)==false){
                            this.candidateProteins.add(protein);
                        }
                    }
            }
            proteinPeptideRelationBufferReader.close();
            proteinPeptideRelationReader.close();

        }catch(Exception e){
            e.printStackTrace();
        }
        
        for(int i=0;i<candidateProteins.size();i++){
            System.out.print(candidateProteins.get(i)+",");
        }
        this.population=new int[this.populationSize][];
        this.fitness=new double[this.populationSize];

        
        peptides=new ArrayList<Integer>();
        StringTokenizer tokens=new StringTokenizer(testString, ", ");
        while(tokens.hasMoreTokens()){
            peptides.add(Integer.parseInt(tokens.nextToken()));
        }
        for(int i=0; i<this.populationSize; i++)
        {
            this.population[i]=new int[this.candidateProteins.size()];

            for(int j=0; j < this.candidateProteins.size(); j++)
            {
                if(randomGenerator.nextGaussian()>0.5)
                {
                    this.population[i][j]=1;
                }
            }
        }        
        for(int i=0; i<this.population.length; i++){
            this.fitness[i]=evaluateFitness(this.population[i],true);
        }
        for(int i=0;i<populationSize+Util.num_off;i++){
            parentOffspringPopulation[i] = new int[population[0].length];
        }
        
        for(int generation = 0; generation<Util.num_gen ; generation++ ){
            System.out.print("\n"+generation+" "+fitness[0]);            
            sortByFitness();
            for(int i=0;i<populationSize;i++){
                parentOffspringPopulation[i] = population[i];
            }
            int offSpring = 0;
            while(offSpring < Util.num_off){                
                int firstSeletion = this.fussSelection();
                int secondSelection;
                do{
                    secondSelection = this.rouletteWheelSelection();
                } while (firstSeletion == secondSelection);
                for(int i=0;i<population[firstSeletion].length;i++) {
                    parentOffspringPopulation[offSpring+populationSize][i] 
                            = population[firstSeletion][i];
                }
                for(int i=0;i<population[secondSelection].length;i++) {
                    parentOffspringPopulation[populationSize+offSpring+1][i] 
                            = population[secondSelection][i];
                }
                crossOver(parentOffspringPopulation[populationSize+offSpring], 
                        parentOffspringPopulation[populationSize+offSpring+1]);
                
                mutation(parentOffspringPopulation[
                        populationSize+offSpring], 0.7);
                mutation(parentOffspringPopulation[
                        populationSize+offSpring+1], 0.7);
                offSpring+=2;                
            }
            sortBy(parentOffspringPopulation);
            
            int []forImprovement = new int[parentOffspringPopulation.length/10];
            int a=0;
            boolean visited[]=new boolean[forImprovement.length];
            for(int i=0;i<forImprovement.length;i++)visited[i]=false;
            for(a=0;a<forImprovement.length/10;a++) {
                forImprovement[a]=a;
                visited[a] = true;
            }                       
            for(;a<forImprovement.length;){
                int index = randomGenerator.nextInt(forImprovement.length-a);
                forImprovement[a]=index;
                visited[index] = true;
                a++;                
            }
            for(int im=0;im<forImprovement.length;im++){
                sim.improve(parentOffspringPopulation[forImprovement[im]]);
            }
            int length = population[0].length;
            int elitSolution = (int)(populationSize * 0.1);
            int i=0;
            for(i=0;i<elitSolution;i++){
                System.arraycopy(parentOffspringPopulation[i], 0, population[i], 0, length);
                fitness[i] = evaluateFitness(population[i], true);
            }      
            for(;i<populationSize;i++){
                int index = randomGenerator.nextInt(populationSize-i)+i;
                System.arraycopy(parentOffspringPopulation[index], 0, population[i], 0, length);
                fitness[i] = evaluateFitness(population[i], true);
            }
        }
        total_iteration=iteration;
        return this.population[0];
    }

    double evaluateFitness(int a[],boolean flag){
        HashSet<Integer> allPeptideSet = new HashSet<Integer>();
        int proteinCount = 0;
        for(int i=0;i<a.length;i++){
            if(a[i]==1){
                proteinCount++;
                int actualProtein = candidateProteins.get(i);
                ArrayList<Integer> peptidesOfThisProtein = protenPeptideMap.get(actualProtein);
                if(peptidesOfThisProtein!=null){
                    allPeptideSet.addAll(peptidesOfThisProtein);
                }
            }
        }
        Iterator<Integer>iterator = testPeptides.iterator();
        int tp=0,fp=0,fn=0;
        HashSet<Integer> testPeptideSet = new HashSet<Integer>(testPeptides);
        tp = Util.intersection(testPeptideSet, allPeptideSet);
        fp = Util.difference(allPeptideSet, testPeptideSet);
        fn = Util.difference(testPeptideSet, allPeptideSet);
        double pr = (double)tp/(tp+fp);
        double rc = (double)tp/(tp+fn);
        return 1.0/(beta*(1/pr)+(1-beta)*(1/rc))*Math.pow((double)proteinCount, 2*epsilon-1 );
    }


    void sortBy(int [][]popul){
        double []fit = new double[popul.length];
        for(int i=0;i<fit.length;i++){
            fit[i] = this.evaluateFitness(popul[i], true);
        }
        for(int i=0; i< popul.length ; i++){
            for(int j=popul.length-1; j>i ; j--){
                if(fit[j] > fit[j-1]){
                    int []temp = popul[j];
                    popul[j]=popul[j-1];
                    popul[j-1]=temp;
                    double tempFitness=fit[j];
                    fit[j]=fit[j-1];
                    fit[j-1]=tempFitness;
                }
            }
        }
        
    }

    void sortByFitness(){
        for(int i=0; i< this.population.length ; i++)
        {
            for(int j=this.population.length-1; j>i ; j--)
            {
                if(this.fitness[j] > this.fitness[j-1])
                {
                    int []temp=this.population[j];
                    this.population[j]=this.population[j-1];
                    this.population[j-1]=temp;
                    double tempFitness=this.fitness[j];
                    this.fitness[j]=this.fitness[j-1];
                    this.fitness[j-1]=tempFitness;
                }
            }
        }
        
    }


    int rouletteWheelSelection(){
            double []relativeFitness=new double[this.populationSize];
            double sum=0;
            for (int i=0; i<this.populationSize;i++)
            {
                sum=sum+this.fitness[i];
            }

            for (int i=0; i<this.populationSize;i++)
            {
                relativeFitness[i]=this.fitness[i]/sum;
            }

            double p=randomGenerator.nextGaussian();
            sum=0;

            for (int i=0; i<populationSize;i++)
            {
                sum+=relativeFitness[i];
                if(sum>p)
                {
                        return i;
                }
            }
            return 0;
    }
    
    int fussSelection(){
        double maxFitness = this.fitness[0];
        double minFitness = this.fitness[this.populationSize-1];
        //System.out.println(minFitness+" "+maxFitness);
        double p = (randomGenerator.nextDouble()*(maxFitness-minFitness))+minFitness;
        for(int i=populationSize-1;i>=0;i--){
            if(fitness[i]>p){
                return i;
            }
        }
        return 0;
    }



    void crossOver(int []firstChromosome,int []secondChromosome){
        int firstCrossoverPoint=randomGenerator.nextInt(this.population[0].length);
        int secondCrossoverPoint=randomGenerator.nextInt(this.population[0].length);
        if (firstCrossoverPoint > secondCrossoverPoint){
                int temp=firstCrossoverPoint;
                firstCrossoverPoint=secondCrossoverPoint;
                secondCrossoverPoint=temp;
        }
                    
        for(int i=firstCrossoverPoint; i<=secondCrossoverPoint; i++){
                int temp=firstChromosome[i];
                firstChromosome[i]=secondChromosome[i];
                secondChromosome[i]=temp;
        }
    }


    void mutation(int []chromosome,double prob){
        double p=randomGenerator.nextGaussian();
        if(p<prob){
            int mutationPoint=randomGenerator.nextInt(this.population[0].length);
            if(chromosome[mutationPoint]==0){
                chromosome[mutationPoint]=1;
            }
            else{
                chromosome[mutationPoint]=0;
            }
        }
    }

    
    public void hillClimbing(int []node){
        int []current = new int[node.length];
        int []next = new int[node.length];
        for(int i=0;i<node.length;i++){
            current[i]=node[i];
            next[i]=node[i];
        }
        double currentFitness = evaluateFitness(current, true);
        double nextFitness =evaluateFitness(next, true);
        
        int unimproved = 0;
        while(unimproved<node.length/4){
            System.out.println(currentFitness);
            mutation(next,1.0);
            nextFitness = evaluateFitness(next, true);
            if(nextFitness>currentFitness){
                currentFitness = nextFitness;
                System.arraycopy(next, 0, current, 0, current.length);
                unimproved=0;
            }
            else if(currentFitness-nextFitness <=0.00000000000000001){
                unimproved++;
            }
            else{
                System.arraycopy(current, 0, next, 0, current.length);
                unimproved++;
            }            
        }
        System.arraycopy(current, 0, node, 0, node.length);
        
    }
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args)
    {
        FileWriter timer = null;
        try {
            parameterScan = new Scanner(new File("parameters/epsilon.txt"));
            parameterWriter = new PrintWriter(new File("parameters/epsilonBetaOutput.csv"));
            while(parameterScan.hasNextLine()){
                parameterScan2 = new Scanner(new File("parameters/beta.txt"));
                while(parameterScan2.hasNextLine()){
                    System.out.println(epsilon+" "+beta);
                    new ProteinInference();
                }
                parameterScan2.close();;
            }
            
            parameterWriter.close();
           
        } catch (Exception ex) {
            Logger.getLogger(ProteinInference.class.getName()).log(Level.SEVERE, null, ex);
        } 
    }
}