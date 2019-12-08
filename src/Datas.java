public class Datas {
    int[] head;
    int[] succ;
    int edges;
    int vertices;
    int[][] adjMat;

    public Datas(int[] head, int[] succ, int edges, int vertices, int[][] adjMat) {
        this.head = head;
        this.succ = succ;
        this.edges = edges;
        this.vertices = vertices;
        this.adjMat = adjMat;
    }
}
