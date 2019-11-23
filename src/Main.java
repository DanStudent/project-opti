import java.io.*;
import java.util.Arrays;
import java.util.LinkedList;

public class Main {
    static int []head;
    static int []succ;
    static int edges;
    static int vertices;
    static int[][] adjMat; // la matrice d'adjacence utile pour calculer la modularité
    public static int[] create_tab_from_string(String str){
        LinkedList<Integer> fifo = new LinkedList<Integer>();
        int i = 0;
        while (i< str.length()){
            String parsedNumber = "";
            while (i<str.length() && str.charAt(i) != ' '){
                parsedNumber += str.charAt(i);
                i++;
            }
            fifo.add(Integer.parseInt(parsedNumber));
            i++;
        }
        int tab[] = new int[fifo.size()];
        for (int j = 0; j < tab.length; j++) {
            tab[j] = fifo.removeFirst();
        }
        return tab;
    } // permet de creer un tableau d'entier à partir d'un string
    public static void parse(String filename){
        try {
            File file = new File(filename);
            String path = file.getAbsolutePath().substring(0, file.getAbsolutePath().length()-filename.length()); // permet de prendre le chemin correct du fichier
            FileReader fileReader = new FileReader(path+"src\\"+filename);
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
    public static int density(int x){
        return head[sommet(x)+1]-head[sommet(x)];
    } //retourne la densité d'un sommet
    public static int[][] init_tab(int []head){
        int adj[][] = new int[head.length-1][head.length-1];
        for (int i = 0; i < adj.length; i++) {
            for (int j = 0; j < adj[i].length; j++) {
                adj[i][j] = 0;
            }
        }
        return adj;
    } // permet d'initialiser une matrice contenant que des zeros
    public static int[][] create_adjMatrix(int []head, int []succ){
        int adj[][] = init_tab(head);
        for (int i = 0; i < head.length-1; i++) {
            for (int j = head[i]-1; j <head[i+1]-1 ; j++) {
                adj[i][succ[j]-1]=1;
            }
        }
        return adj;
    }
    public static double modularity(int [][]sets){
        double res = 0;
        for (int k = 0; k < sets.length; k++) {
            int []set = sets[k];
            for (int i = 0; i < set.length; i++) {
                for (int j = 0; j < set.length; j++) {
                    double numerator = (density(set[i]) * density(set[j]));
                    double secondTerm = numerator/edges;
                    res += (adjMat[sommet(set[i])][sommet(set[j])] - secondTerm);
                }
            }
        }
        return res/edges;
    }
    public static int sommet(int x){
        return x-1;
    }
    public static void display(int tab[]){
        for (int x:tab
             ) {
            System.out.print(x + " ");
        }
        System.out.print("\n");
    } // afficher un tableau 1d
    public static void display2d(int tab[][]){
        System.out.println(Arrays.deepToString(tab).replace("], ", "]\n"));;
    } // affiche une matrice
    public static void main(String[] args) {
        parse("test.txt");
        adjMat = create_adjMatrix(head, succ);
        display2d(adjMat);
        int sets[][]= new int[6][1];
        sets[0] = new int[]{1};
        sets[1] = new int[]{2};
        sets[2] = new int[]{3};
        sets[3] = new int[]{4};
        sets[4] = new int[]{5};
        sets[5] = new int[]{6};
        System.out.println(modularity(sets));
    }
}
