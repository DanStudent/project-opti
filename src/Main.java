import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.Random;

public class Main {
    static int[] head;
    static int[] succ;
    static int edges;
    static int vertices;
    static int[][] adjMat; // la matrice d'adjacence utile pour calculer la modularité

    static int[] create_tab_from_string(String str) {
        LinkedList<Integer> fifo = new LinkedList<Integer>();
        int i = 0;
        while (i < str.length()) {
            String parsedNumber = "";
            while (i < str.length() && str.charAt(i) != ' ') {
                parsedNumber += str.charAt(i);
                i++;
            }
            if(!parsedNumber.equals("")){
                fifo.add(Integer.parseInt(parsedNumber));
            }
            i++;
        }
        int[] tab = new int[fifo.size()];
        for (int j = 0; j < tab.length; j++) {
            tab[j] = fifo.removeFirst();
        }
        return tab;
    } // permet de creer un tableau d'entier à partir d'un string

    static void parse(String filename) {
        try {
            File file = new File(filename);
            String path = file.getAbsolutePath().substring(0, file.getAbsolutePath().length() - filename.length()); // permet de prendre le chemin correct du fichier
            FileReader fileReader = new FileReader(path + "src" + File.separator + filename);
            BufferedReader bufferedReader = new BufferedReader(fileReader);
            edges = Integer.parseInt(bufferedReader.readLine().strip());
            vertices = Integer.parseInt(bufferedReader.readLine().strip());
            head = create_tab_from_string(bufferedReader.readLine());
            succ = create_tab_from_string(bufferedReader.readLine());

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    } // s'occupe d'extraire les informations du fichier text

    static int density(int x) {
        return head[sommet(x) + 1] - head[sommet(x)];
    } //retourne la densité d'un sommet

    static int[][] init_tab(int[] head) {
        int[][] adj = new int[head.length - 1][head.length - 1];
        for (int i = 0; i < adj.length; i++) {
            for (int j = 0; j < adj[i].length; j++) {
                adj[i][j] = 0;
            }
        }
        return adj;
    } // permet d'initialiser une matrice contenant que des zeros

    static int[][] create_adjMatrix(int[] head, int[] succ) {
        int[][] adj = init_tab(head);
        for (int i = 0; i < head.length - 1; i++) {
            for (int j = head[i] - 1; j < head[i + 1] - 1; j++) {
                adj[i][succ[j] - 1] = 1;
            }
        }
        return adj;
    }

    static double modularity(ArrayList<ArrayList<Integer>> sets) {
        double res = 0;
        for (int k = 0; k < sets.size(); k++) {
            ArrayList<Integer> set = sets.get(k);
            for (int i = 0; i < set.size(); i++) {
                for (int j = 0; j < set.size(); j++) {
                    double numerator = (density(set.get(i)) * density(set.get(j)));
                    double secondTerm = numerator / (2*edges);
                    res += (adjMat[sommet(set.get(i))][sommet(set.get(j))] - secondTerm);
                }
            }
        }
        return res / (2*edges);
    }

    static int sommet(int x) {
        return x - 1;
    }

    static Cost evaluate_costs(ArrayList<Integer> candidates, ArrayList<ArrayList<Integer>> sols) {
        double costs[] = new double[vertices];
        ArrayList<ArrayList<Integer>> potentialSols[] = new ArrayList[vertices];

        for (int vertice : candidates) {
            if (sols.isEmpty()) {
                ArrayList<Integer> newSet = new ArrayList<>();
                newSet.add(vertice);
                sols.add(newSet);
                costs[sommet(vertice)] = modularity(sols);
                potentialSols[sommet(vertice)] = new ArrayList<>(sols);
                sols.remove(0);
            } else {
                ArrayList<Integer> newSet = new ArrayList<>();
                newSet.add(vertice);
                sols.add(newSet);
                costs[sommet(vertice)] = modularity(sols);
                potentialSols[sommet(vertice)] = new ArrayList<>(sols);
                sols.remove(sols.size() - 1);
                for (int i = 0; i < sols.size(); i++) {
                    ArrayList<Integer> set = sols.get(i);
                    set.add(vertice);
                    double cost = modularity(sols);
                    if (cost > costs[sommet(vertice)]) {
                        costs[sommet(vertice)] = cost;
                        potentialSols[sommet(vertice)] = copy(sols);

                    }

                    set.remove(set.size() - 1);
                }
            }
        }
        return new Cost(costs, potentialSols);
    }
    static ArrayList<ArrayList<Integer>> copy(ArrayList<ArrayList<Integer>> a){
        ArrayList<ArrayList<Integer>> newArrayList = new ArrayList<>();
        for (int i = 0; i < a.size(); i++) {
            ArrayList<Integer> b = new ArrayList<>();
            b.addAll(a.get(i));
            newArrayList.add(b);
        }
        return newArrayList;

    }

    static ArrayList<ArrayList<Integer>> greedy_randomized_construction(double alpha) {
        ArrayList<ArrayList<Integer>> solution = new ArrayList<>();
        ArrayList<Integer> usedCandidates = new ArrayList<>();
        ArrayList<Integer> candidates = init_candidates(vertices, usedCandidates);
        Cost costs = evaluate_costs(candidates, solution);
        while (!candidates.isEmpty()) {
            Double[] newCosts = Arrays.stream(costs.cost).boxed().toArray(Double[]::new);
            double c_min = getC_min(newCosts);
            double c_max = getC_max(newCosts);
            ArrayList<Integer> RCL = create_RCL(c_min, c_max, alpha, candidates, newCosts);
            Random random = new Random();
            int candidate_chosen = RCL.get(random.nextInt(RCL.size()));
            solution = costs.sols[sommet(candidate_chosen)];
            usedCandidates.add(candidate_chosen);
            candidates = init_candidates(vertices, usedCandidates);
            costs = evaluate_costs(candidates, solution);
        }
        return solution;

    }

    static ArrayList<Integer> create_RCL(double c_min, double c_max, double alpha, ArrayList<Integer> candidates, Double[] costs) {

        ArrayList<Integer> RCL = new ArrayList<>();
        System.out.println("c_min : "+c_min);
        System.out.println("c_max : "+c_max);
        System.out.println("borne inferieur : " + (c_max - alpha * (c_max - c_min)));
        for (int vertices : candidates) {
            if (costs[sommet(vertices)] != 0) {
                if (costs[sommet(vertices)] >= c_max - alpha * (c_max - c_min)) {
                    RCL.add(vertices);
                }
            }

        }
        return RCL;
    }

    static double getC_min(Double[] costs) {
        double min = Double.POSITIVE_INFINITY;
        for (int i = 0; i < costs.length; i++) {
            if (costs[i] != 0) {
                if (costs[i] < min) {
                    min = costs[i];
                }
            }
        }
        return min;
    }

    static double getC_max(Double[] costs) {
        double max = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < costs.length; i++) {
            if (costs[i] != 0) {
                if (costs[i] > max) {
                    max = costs[i];
                }
            }
        }
        return max;
    }

    static ArrayList<Integer> init_candidates(int n, ArrayList<Integer> usedCandidates) {
        ArrayList<Integer> candidates = new ArrayList<>();
        for (int i = 1; i <= n; i++) {
            if (!(usedCandidates.contains(i))) {
                candidates.add(i);
            }
        }
        return candidates;
    }

    static boolean feasability(ArrayList<ArrayList<Integer>> sols){
        for (int i = 0; i < sols.size(); i++) {
            ArrayList<Integer> set = sols.get(i);
            for (int j = 0; j < set.size(); j++) {
                int vertice = set.get(j);
                if(!linked(vertice, set)){
                    return false;
                }
            }
        }
        return true;
    }
    static boolean linked(int vertice, ArrayList<Integer> set){
        boolean value = false;
        for (int sommet : set) {
            if(sommet!=vertice){
                if(adjMat[sommet(vertice)][sommet(sommet)] == 1){
                    value = true;
                }
            }
        }
        return value;
    }
    static void repair_solution(ArrayList<ArrayList<Integer>> sols){
        for (ArrayList<Integer> set : sols) {
            if(set.size()>1){
                for (int i = 0; i < set.size(); i++) {
                    int vertice = set.get(i);
                    if(!linked(vertice, set)){
                        set.remove(vertice);
                        ArrayList<Integer> newSet = new ArrayList<>();
                        newSet.add(vertice);
                        sols.add(set);
                    }
                }


            }
        }
    }
    static void display(int[] tab) {
        for (int x : tab
        ) {
            System.out.print(x + " ");
        }
        System.out.print("\n");
    } // afficher un tableau 1d

    static void display(double[] tab) {
        for (double x : tab
        ) {
            System.out.print(x + " ");
        }
        System.out.print("\n");
    }

    static void display2d(int[][] tab) {
        System.out.println(Arrays.deepToString(tab).replace("], ", "]\n"));
    } // affiche une matrice
    static void displaySolutions(ArrayList<ArrayList<Integer>>[] tab){
        for (int i = 0; i < tab.length; i++) {
            System.out.println(tab[i]);
        }
    }

    public static void main(String[] args) {
        parse("File5.txt");
        adjMat = create_adjMatrix(head, succ);
        display2d(adjMat);
        ArrayList<ArrayList<Integer>> solution = greedy_randomized_construction(1);
         if (!feasability(solution)){
             System.out.println("before : "+ solution);
             repair_solution(solution);
         }
         System.out.println(solution);
         System.out.println(modularity(solution));


    }
}
