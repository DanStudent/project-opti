import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.Random;

public class Main {
    static int[] create_tab_from_string(String str) {
        LinkedList<Integer> fifo = new LinkedList<Integer>();
        int i = 0;
        while (i < str.length()) {
            String parsedNumber = "";
            while (i < str.length() && str.charAt(i) != ' ') {
                parsedNumber += str.charAt(i);
                i++;
            }
            if (!parsedNumber.equals("")) {
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

    static Datas parse(String filename) {
        try {
            File file = new File(filename);
            String path = file.getAbsolutePath().substring(0, file.getAbsolutePath().length() - filename.length()); // permet de prendre le chemin correct du fichier
            FileReader fileReader = new FileReader(path + "src" + File.separator + filename);
            BufferedReader bufferedReader = new BufferedReader(fileReader);
            int edges = Integer.parseInt(bufferedReader.readLine().strip());
            int vertices = Integer.parseInt(bufferedReader.readLine().strip());
            int[] head = create_tab_from_string(bufferedReader.readLine());
            int[] succ = create_tab_from_string(bufferedReader.readLine());
            return new Datas(head, succ, edges, vertices, create_adjMatrix(head, succ));

        } catch (FileNotFoundException e) {
            e.printStackTrace();
            return null;
        } catch (IOException e) {
            e.printStackTrace();
            return null;
        }
    } // s'occupe d'extraire les informations du fichier text

    static int density(int x, int[] head) {
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

    static double modularity(ArrayList<ArrayList<Integer>> sets, Datas datas) {
        double res = 0;
        for (int k = 0; k < sets.size(); k++) {
            ArrayList<Integer> set = sets.get(k);
            for (int i = 0; i < set.size(); i++) {
                for (int j = 0; j < set.size(); j++) {
                    double numerator = (density(set.get(i), datas.head) * density(set.get(j), datas.head));
                    double secondTerm = numerator / (2 * datas.edges);
                    res += (datas.adjMat[sommet(set.get(i))][sommet(set.get(j))] - secondTerm);
                }
            }
        }
        return res / (2 * datas.edges);
    }

    static int sommet(int x) {
        return x - 1;
    }

    static Cost evaluate_costs(ArrayList<Integer> candidates, ArrayList<ArrayList<Integer>> sols, Datas datas) {
        double costs[] = new double[datas.vertices];
        ArrayList<ArrayList<Integer>> potentialSols[] = new ArrayList[datas.vertices];

        for (int vertice : candidates) {
            if (sols.isEmpty()) {
                ArrayList<Integer> newSet = new ArrayList<>();
                newSet.add(vertice);
                sols.add(newSet);
                costs[sommet(vertice)] = modularity(sols, datas);
                potentialSols[sommet(vertice)] = new ArrayList<>(sols);
                sols.remove(0);
            } else {
                ArrayList<Integer> newSet = new ArrayList<>();
                newSet.add(vertice);
                sols.add(newSet);
                costs[sommet(vertice)] = modularity(sols, datas);
                potentialSols[sommet(vertice)] = new ArrayList<>(sols);
                sols.remove(sols.size() - 1);
                for (int i = 0; i < sols.size(); i++) {
                    ArrayList<Integer> set = sols.get(i);
                    set.add(vertice);
                    double cost = modularity(sols, datas);
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

    static ArrayList<ArrayList<Integer>> copy(ArrayList<ArrayList<Integer>> a) {
        ArrayList<ArrayList<Integer>> newArrayList = new ArrayList<>();
        for (int i = 0; i < a.size(); i++) {
            ArrayList<Integer> b = new ArrayList<>();
            b.addAll(a.get(i));
            newArrayList.add(b);
        }
        return newArrayList;

    }

    static ArrayList<ArrayList<Integer>> greedy_randomized_construction(double alpha, Datas datas) {
        ArrayList<ArrayList<Integer>> solution = new ArrayList<>();
        ArrayList<Integer> usedCandidates = new ArrayList<>();
        ArrayList<Integer> candidates = init_candidates(datas.vertices, usedCandidates);
        Cost costs = evaluate_costs(candidates, solution, datas);
        while (!candidates.isEmpty()) {
            Double[] newCosts = Arrays.stream(costs.cost).boxed().toArray(Double[]::new);
            double c_min = getC_min(newCosts);
            double c_max = getC_max(newCosts);
            ArrayList<Integer> RCL = create_RCL(c_min, c_max, alpha, candidates, newCosts);
            Random random = new Random();
            int candidate_chosen = RCL.get(random.nextInt(RCL.size()));
            solution = costs.sols[sommet(candidate_chosen)];
            usedCandidates.add(candidate_chosen);
            candidates = init_candidates(datas.vertices, usedCandidates);
            costs = evaluate_costs(candidates, solution, datas);
        }
        return solution;

    }

    static ArrayList<Integer> create_RCL(double c_min, double c_max, double alpha, ArrayList<Integer> candidates, Double[] costs) {

        ArrayList<Integer> RCL = new ArrayList<>();
        /*System.out.println("c_min : " + c_min);
        System.out.println("c_max : " + c_max);
        System.out.println("borne inferieur : " + (c_max - alpha * (c_max - c_min)));*/
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

    static boolean feasability(ArrayList<ArrayList<Integer>> sols, Datas datas) {
        for (int i = 0; i < sols.size(); i++) {
            ArrayList<Integer> set = sols.get(i);
            for (int j = 0; j < set.size(); j++) {
                int vertice = set.get(j);
                if (!linked(vertice, set, datas)) {
                    return false;
                }
            }
        }
        return true;
    }

    static boolean linked(int vertice, ArrayList<Integer> set, Datas datas) {
        boolean value = false;
        for (int sommet : set) {
            if (sommet != vertice) {
                if (datas.adjMat[sommet(vertice)][sommet(sommet)] == 1) {
                    value = true;
                }
            }
        }
        return value;
    }

    static ArrayList<ArrayList<ArrayList<Integer>>> get_neighborhood(ArrayList<ArrayList<Integer>> solution, Datas datas) {
        ArrayList<ArrayList<ArrayList<Integer>>> neighborhood = new ArrayList<>();
        Random random = new Random();
        int index = random.nextInt(solution.size());
        ArrayList<Integer> set = solution.get(index);
        if (set.size() == 1) {
            ArrayList<ArrayList<Integer>> copyOfSolution = copy(solution);
            int vertice = set.remove(0);
            copyOfSolution.remove(set);
            for (int i = 0; i < copyOfSolution.size(); i++) {
                ArrayList<ArrayList<Integer>> neighbor = copy(copyOfSolution);
                neighbor.get(i).add(vertice);
                if (feasability(neighbor, datas)) {
                    neighborhood.add(neighbor);
                }
            }
        } else {
            for (int i = 0; i < set.size(); i++) {
                int vertice = set.remove(i);
                for (int j = 0; j < solution.size(); j++) {
                    if (j != index) {
                        ArrayList<ArrayList<Integer>> neighbor = copy(solution);
                        neighbor.get(j).add(vertice);
                        if (feasability(neighbor, datas)) {
                            neighborhood.add(neighbor);
                        }

                    }
                }
                set.add(vertice);
            }
        }
        return neighborhood;
    }

    static ArrayList<ArrayList<Integer>> local_search(ArrayList<ArrayList<Integer>> solution, Datas datas) {
        ArrayList<ArrayList<ArrayList<Integer>>> neighborhood = get_neighborhood(solution, datas);
        ArrayList<ArrayList<Integer>> bestSol = solution;
        double modularity = modularity(solution, datas);
        boolean solutionChanged = true;
        while (solutionChanged) {
            solutionChanged = false;
            for (int i = 0; i < neighborhood.size(); i++) {
                if (modularity(neighborhood.get(i), datas) > modularity) {
                    bestSol = neighborhood.get(i);
                    modularity = modularity(neighborhood.get(i), datas);
                    solutionChanged = true;
                }
            }
        }
        return bestSol;
    }

    static ArrayList<ArrayList<Integer>> grasp(double alpha, int maxIteration, Datas datas) {
        ArrayList<ArrayList<Integer>> bestSolution = null;
        for (int i = 0; i < maxIteration; i++) {
            ArrayList<ArrayList<Integer>> solution;
            solution = greedy_randomized_construction(alpha, datas);
            if (feasability(solution, datas)) {
                solution = local_search(solution, datas);
                if (bestSolution == null) {
                    bestSolution = solution;
                } else {
                    if (modularity(solution, datas) > modularity(bestSolution, datas)) {
                        bestSolution = solution;
                    }
                }
            }
        }
        return bestSolution;
    }

    static void display2d(int[][] tab) {
        System.out.println(Arrays.deepToString(tab).replace("], ", "]\n"));
    } // affiche une matrice

    public static void main(String[] args) {
        /// PLUS IL Y A DE SOMMETS, PLUS IL EST INTERESSANT DE PRENDRE UN ALPHA GRAND
        Datas file1 = parse("File1.txt");
        Datas file2 = parse("File2.txt");
        Datas file3 = parse("File3.txt");
        Datas file4 = parse("File4.txt");
        Datas file5 = parse("File5.txt");
        Thread thread1 = new Thread() {
            @Override
            public void run() {
                ArrayList<ArrayList<Integer>> solution = grasp(0.5, 100, file1);
                System.out.println("Thread 1 over");
                System.out.println("La solution est : " + solution);
                System.out.println("La modularite est de : " + modularity(solution, file1));
            }
        };
        Thread thread2 = new Thread() {
            @Override
            public void run() {
                ArrayList<ArrayList<Integer>> solution = grasp(0.5, 100, file2);
                System.out.println("Thread 2 over");
                System.out.println("La solution est : " + solution);
                System.out.println("La modularite est de : " + modularity(solution, file2));
            }
        };
        Thread thread3 = new Thread() {
            @Override
            public void run() {
                ArrayList<ArrayList<Integer>> solution = grasp(0.5, 100, file3);
                System.out.println("Thread 3 over");
                System.out.println("La solution est : " + solution);
                System.out.println("La modularite est de : " + modularity(solution, file3));
            }
        };
        Thread thread4 = new Thread() {
            @Override
            public void run() {
                ArrayList<ArrayList<Integer>> solution = grasp(0.5, 100, file4);
                System.out.println("Thread 4 over");
                System.out.println("La solution est : " + solution);
                System.out.println("La modularite est de : " + modularity(solution, file4));
            }
        };
        Thread thread5 = new Thread() {
            @Override
            public void run() {
                ArrayList<ArrayList<Integer>> solution = grasp(1, 100, file5);
                System.out.println("Thread 5 over");
                System.out.println("La solution est : " + solution);
                System.out.println("La modularite est de : " + modularity(solution, file5));
            }
        };
        thread1.start();
        thread2.start();
        thread3.start();
        thread4.start();
        thread5.start();
    }
}
