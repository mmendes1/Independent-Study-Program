import java.util.Random;

public class Matrix {

    private int width, height;
    private int[][] matrix;

    public Matrix () {}

    public Matrix (int width, int height) {
        int[][] matrix = new int[width][height];
            this.width = width;
            this.height = height;
            this.matrix = matrix;
    }

    public void fillMatrix() {
        Random rand = new Random();

        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                int random = rand.nextInt(1000+ 1);
                //System.out.println("Value: " + random);
                matrix[i][j] = random;
            }
        }
        printMatrix(matrix, "Matrix 1");
    }

    public void printMatrix(int[][] matrix, String matrixID) {
        System.out.println(matrixID);
        System.out.println("-----------------------------------------------------------------------");
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                System.out.print(matrix[i][j] + " ");
            }
            System.out.println("");
        }
    }
}
