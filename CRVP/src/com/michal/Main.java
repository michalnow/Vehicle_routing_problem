package com.michal;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

public class Main {

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);
        System.out.println("Enter starting population number: ");
        int populationQuantity = scanner.nextInt();
        System.out.println("Enter vehicle capacity: ");
        int vehicleCapacity = scanner.nextInt();
        System.out.println("Enter vehicle quantity: ");
        int vehicleQuantity = scanner.nextInt();
        System.out.println("Enter mutation rate (1 - 100): ");
        int mutationRate = scanner.nextInt(); // 1 / 2
        System.out.println("Enter crossover rate (1 - 100): ");
        int crossoverRate = scanner.nextInt(); // 1 / 2
        System.out.println("Enter number of Generations: ");
        int maxGeneration = scanner.nextInt();
        System.out.println("Number of chromosomes to be selected: ");
        int truncationNumber = scanner.nextInt();
        int numberOfCities = 0;
        //Path path = Paths.get("./data.txt");
        Path path = Paths.get("../../../data.txt");
        Map<Integer, City> cities = new TreeMap<Integer, City>();

        try {
            Files.lines(path).forEach(line ->
            {
                List<Double> coords = new ArrayList<Double>();
                String arr[] = line.split((";"));
                coords.add(Double.parseDouble(arr[2]));
                coords.add(Double.parseDouble(arr[3]));
                City city = new City(coords, Integer.parseInt(arr[4]),arr[1]);
                cities.put(Integer.parseInt(arr[0]),city);
            });
            numberOfCities = cities.size() - 1;
            System.out.println(numberOfCities);
        }
        catch (IOException e)
        {
            e.printStackTrace();
        }

        List<List<Double>> distanceMatrix = calculateDistanceMatrix(cities);
        List<Chromosome> chromosomes = generatePopulation(populationQuantity, vehicleQuantity, numberOfCities);

        int currentGeneration = 0;
        while(maxGeneration != currentGeneration ) {
            currentGeneration++;
            for (int i = 0; i < chromosomes.size(); i++) {
                countFitnessLvlForChromosome(chromosomes.get(i), distanceMatrix, vehicleCapacity, cities);
            }
            sort(chromosomes);

            List<Chromosome> truncatedPopulation = truncationSelection(chromosomes, truncationNumber);

            List<Chromosome> newChildren = crossoverOX(truncatedPopulation, vehicleQuantity, numberOfCities, crossoverRate);
            mutateChildren(newChildren, numberOfCities, mutationRate, vehicleQuantity);

            removeLastNChromosomes(chromosomes,(int) (Math.random() * chromosomes.size() * (1 / 10.0)));
            addChildrenToPopulation(chromosomes, newChildren);
            for (int i = 0; i < chromosomes.size(); i++) {
                countFitnessLvlForChromosome(chromosomes.get(i), distanceMatrix, vehicleCapacity, cities);
            }
            System.out.println("Generation " + currentGeneration);
            System.out.println("Population size " + chromosomes.size());
            sort(chromosomes);
            System.out.println("Best Solution= " + getFittestChromosome(chromosomes).getFitnessLvl());

            for(Gene gene: chromosomes.get(0).getGenes()){
                System.out.print(gene.getRoute() + ", ");
            }
            System.out.println();

            System.out.println("Best route: ");
            printRoute(chromosomes.get(0),cities);

        }
    }

    static double calculateEuclideanDistance(double x1, double x2, double y1, double y2){
        return Math.sqrt(Math.pow(x2 - x1,2) + Math.pow(y2 - y1,2)) * 73;
    }

    static List<List<Double>> calculateDistanceMatrix(Map<Integer, City> cities){

        Set<Integer> citiesSet = cities.keySet();
        List<List<Double>> distanceMatrix = new ArrayList<>();

        for(Integer city: citiesSet){

            List<Double> tmpDistance = new ArrayList<Double>();
            for(Integer city2: citiesSet){

                if(!city.equals(city2)){
                    List<Double> firstCityCoords = cities.get((city)).getCoords();
                    List<Double> secondCityCoords = cities.get((city2)).getCoords();

                    tmpDistance.add(calculateEuclideanDistance(firstCityCoords.get(0), secondCityCoords.get(0),
                            firstCityCoords.get(1), secondCityCoords.get(1)));
                }else{
                    tmpDistance.add(0.0);
                }
            }
            distanceMatrix.add(tmpDistance);
        }

        return distanceMatrix;
    }

    static List<Chromosome> generatePopulation(int populationQuantity , int vehicleQuantity, int cityQuantity){

        List<Chromosome> chromosomes = new ArrayList<Chromosome>();
        List<Integer> cityNumbers = new ArrayList<>();

        for(int i = 1; i <= cityQuantity; i++)
            cityNumbers.add(i);

        for(int i = 0 ; i < populationQuantity ; i++) {
            Collections.shuffle(cityNumbers);

            chromosomes.add(createChromosome(cityNumbers, vehicleQuantity, cityQuantity));
        }

        return chromosomes;
    }

    static void printRoute(Chromosome chromosome, Map<Integer,City> cities){
        Set<Integer> citiesSet = cities.keySet();
        int vehicle = 0;

        for(Gene gene: chromosome.getGenes()){
            vehicle++;
            System.out.print("Vehicle" + vehicle + " route = [ " + cities.get(0).getName() + " -> ");
            for(int i = 0; i < gene.getRoute().size();i++) {
                if(i == gene.getRoute().size() - 1){
                    System.out.print(cities.get(gene.getRoute().get(i)).getName() + " -> Kraków]");
                }else {
                    System.out.print(cities.get(gene.getRoute().get(i)).getName() + " -> ");
                }
            }
            System.out.println();
        }
    }

    //   car1     car2 ...
    //[ [route],[route], ... ]
    static Chromosome createChromosome(List<Integer> numbersOfCity, int vehicleQuantity, int cityQuantity){

        List<Gene> genesTmp = new ArrayList<>();
        Chromosome chromosome = new Chromosome();
        int numberOfCitiesInRoute = cityQuantity / vehicleQuantity;

        for(int k = 0 ; k < vehicleQuantity ; k++) {
            List<Integer> genes = new ArrayList<>();
            for (int j = 0; j < numberOfCitiesInRoute; j++) {
                if(j == numberOfCitiesInRoute - 1 && k == vehicleQuantity - 1){
                    for(int l = 0 ; l < cityQuantity % vehicleQuantity; l++){
                        genes.add(numbersOfCity.get(0));
                        int num = numbersOfCity.remove(0);
                        numbersOfCity.add(num);
                    }
                }
                genes.add(numbersOfCity.get(0));
                int num = numbersOfCity.remove(0);
                numbersOfCity.add(num);
            }

            Gene gene = new Gene();
            gene.setRoute(genes);
            genesTmp.add(gene);
        }

        chromosome.setGenes(genesTmp);
        return chromosome;
    }

    //calculating fitness lvl
    static void countFitnessLvlForChromosome(Chromosome chromosome, List<List<Double>> distanceMatrix, int vehicleCapacity, Map<Integer, City> cities){

        double currentDistance = 0;
        int currentDemand = 0;
        for(Gene gene: chromosome.getGenes())
        {
            //przyjazd do pierwszego miasta
            //currentDistance += distanceMatrix.get(0).get(gene.getRoute().get(0));

            int tmpVehicleCapacity = vehicleCapacity;
            int listSize = gene.getRoute().size() - 1;

            for(int i = 0 ; i < listSize; i++)
            {
                //przyjazd do miasta, dla i = 0, przyjazd z bazy Kraków
                if (i == 0) {
                    currentDistance += distanceMatrix.get(0).get(gene.getRoute().get(0));
                } else {
                    currentDistance += distanceMatrix.get(gene.getRoute().get(i - 1)).get(gene.getRoute().get(i));
                }

                currentDemand = cities.get(i).getDemand();

                while(currentDemand > 0)
                {
                    if (tmpVehicleCapacity - currentDemand < 0)
                    {
                        currentDemand -= tmpVehicleCapacity;
                        //powrót do bazy z miasta itego
                        currentDistance += distanceMatrix.get(gene.getRoute().get(i)).get(0);
                        //wyjazd z krakowa do miasta itego i wzięcie pozostałego towaru
                        currentDistance += distanceMatrix.get(0).get(gene.getRoute().get(i));
                        tmpVehicleCapacity =  vehicleCapacity;

                    } else
                    {
                        tmpVehicleCapacity -= currentDemand;
                        currentDemand = 0;
                    }
                }
            }

            //powrót do bazy
            currentDistance += distanceMatrix.get(gene.getRoute().get(listSize)).get(0);
        }
        chromosome.setFitnessLvl(currentDistance);
    }

    static Chromosome getFittestChromosome(List<Chromosome> chromosomes){
        return chromosomes.get(0);
    }

    static void addChildrenToPopulation(List<Chromosome> population, List<Chromosome> newChildren){
        population.addAll(newChildren);
    }

    static void removeLastNChromosomes(List<Chromosome> population, int N){
        for(int i = 0 ; i < N ; i++){
            int indexToRemove = (int) ((population.size() - N) + Math.random() * N);
            population.remove(indexToRemove);
        }
    }

    static void sort(List<Chromosome> chromosomes){
        chromosomes.sort(Comparator.comparing(Chromosome::getFitnessLvl));
    }

    static List<Chromosome> truncationSelection(List<Chromosome> chromosomes, int limit){
        return chromosomes.stream().limit(limit).collect(Collectors.toList());
    }

    static List<Chromosome> crossoverOX(List<Chromosome> chromosomes, int vehicleQuantity, int cityQuantity, int crossoverRate){
        Collections.shuffle(chromosomes);
        List<Chromosome> childrenChromosome = new ArrayList<Chromosome>();


        for (int i = 0; i < chromosomes.size(); i += 2){
            if(Math.random() < crossoverRate/100.0) {
                List<Integer> parent1 = new ArrayList<>();
                List<Integer> parent2 = new ArrayList<>();

                for (Gene gene : chromosomes.get(i).getGenes())
                    for (Integer route : gene.getRoute())
                        parent1.add(route);

                for (Gene gene : chromosomes.get(i + 1).getGenes())
                    for (Integer route : gene.getRoute())
                        parent2.add(route);

                List<Integer> child1 = new ArrayList<>(Collections.nCopies(30, 0));
                List<Integer> child2 = new ArrayList<>(Collections.nCopies(30, 0));

                Random rand = new Random();
                int n1 = rand.nextInt(parent1.size());
                int n2 = rand.nextInt(parent1.size() - 1);

                int start = Math.min(n1, n2);
                int end = Math.max(n1, n2);

                List<Integer> parent1GenesToCopy = new ArrayList<>(parent1.subList(start, end));
                List<Integer> parent2GenesToCopy = new ArrayList<>(parent2.subList(start, end));

                child1.addAll(start, parent1GenesToCopy);
                child2.addAll(start, parent2GenesToCopy);

                for (int j = 0; j <= parent1GenesToCopy.size() - 1; j++) {
                    parent1.remove(parent1.indexOf(parent2GenesToCopy.get(j)));
                    parent2.remove(parent2.indexOf(parent1GenesToCopy.get(j)));
                }

                for (int z = 0; z < parent1.size(); z++) {
                    child1.set(child1.indexOf(0), parent2.get(z));
                    child2.set(child2.indexOf(0), parent1.get(z));
                }

                Chromosome child1C = createChromosome(child1, vehicleQuantity, cityQuantity);
                Chromosome child2C = createChromosome(child2, vehicleQuantity, cityQuantity);

                childrenChromosome.add(child1C);
                childrenChromosome.add(child2C);
            }
        }

        return childrenChromosome;
    }

    static List<Chromosome> mutateChildren(List<Chromosome> children, int numberOfCities, int mutationRate, int vehicleQuantity){

        List<Chromosome> chromosomes = new ArrayList<Chromosome>();

        for (int i = 0; i < children.size(); i++) {
            List<Integer> currentChromosome = new ArrayList<>();
            for (Gene gene : children.get(i).getGenes())
                for (Integer route : gene.getRoute())
                    currentChromosome.add(route);

            for(Integer city: currentChromosome) {
                if (Math.random() < (mutationRate / 100.0)) {
                    int tmpIndex = currentChromosome.indexOf(city);
                    int tmpValue = city;
                    int indexToSwap = (int) (Math.random() * numberOfCities);
                    int valueToSwap = currentChromosome.get(indexToSwap);
                    currentChromosome.set(tmpIndex, valueToSwap);
                    currentChromosome.set(indexToSwap, tmpValue);
                }
            }

            Chromosome afterMutation = createChromosome(currentChromosome, vehicleQuantity, numberOfCities);
            chromosomes.add(afterMutation);
        }

        return chromosomes;
    }

}

class City {
    private String name;
    private List<Double> coords;
    private int demand;

    City(List<Double> coords, int demand, String name){
        this.coords = coords;
        this.demand = demand;
        this.name = name;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public List<Double> getCoords() {
        return coords;
    }

    public void setCoords(List<Double> coords) {
        this.coords = coords;
    }

    public int getDemand() {
        return demand;
    }

    public void setDemand(int demand) {
        this.demand = demand;
    }
}
class Gene {
    List<Integer> route;

    public Gene() { }

    public Gene(List<Integer> route) {
        this.route = route;
    }

    public List<Integer> getRoute() {
        return route;
    }

    public void setRoute(List<Integer> route) {
        this.route = route;
    }
}
class Chromosome {
    List<Gene> genes;
    double fitnessLvl;

    public Chromosome() { }

    public Chromosome(List<Gene> genes) {
        this.genes = genes;
    }

    public List<Gene> getGenes() {
        return genes;
    }

    public void setGenes(List<Gene> genes) {
        this.genes = genes;
    }

    public double getFitnessLvl() {
        return fitnessLvl;
    }

    public void setFitnessLvl(double fitnessLvl) {
        this.fitnessLvl = fitnessLvl;
    }
}