/**
 * Main class that performs the matrix calculation
 * and outputs the relevant data
 * Author: Mike Mendes
 */
import java.io.IOException;
import java.nio.file.Paths;
import java.util.Random;
import java.util.Scanner;
import java.io.File;

public class MatrixCalc {
    private boolean matrix = false;
    private static int size = 16, userIn, ceiling = 20, runCeiling = 15;
    private static int[][] matrix1;
    private static int[][] matrix2;
    private static long durationBasic, durationStrassen, durationLaderman;

    public static void main(String[] args) {
        Scanner input = new Scanner(System.in);
        MatrixCalc calc = new MatrixCalc();

        do{
            System.out.println("\nPlease selected one of the following options...\n------------------------------------------------------\n    1. Set size of matrices\n    2. Run standard time trial\n    3. Run Basic\n    4. Run Strassen's\n    5. Run Laderman's\n    6. Set Random Value Range for Matrices\n    7. Set Amount of Runs for Time Trial\n    8. Toggle Matrix mode\n    0. Exit");
            userIn = input.nextInt();

            switch (userIn)
            {
                case 1:
                    System.out.println("Please enter the size of the matrices.");
                    userIn = input.nextInt();
                        size = userIn;
                        break;
                case 2:
                    calc.timeTrial();
                        break;
                case 3:
                    calc.createMatrices(true, true);
                    calc.printMatrix(calc.basicMult(matrix1, matrix2), "Result", size);
                        break;
                case 4:
                    calc.createMatrices(true, false);
                    calc.printMatrix(calc.strassenMult(matrix1, matrix2), "Result", size);
                        break;
                case 5:
                    calc.createMatrices(false, false);
                    calc.printMatrix(calc.ladermanMult(matrix1, matrix2), "Result", size);
                        break;
                case 6:
                   System.out.println("Please enter a ceiling value for random values.");
                   ceiling = input.nextInt();
                        break;
                case 7:
                    System.out.println("Please enter the amount of runs to be done per time trial.");
                    runCeiling = input.nextInt();
                        break;
                case 8:
                    System.out.println("Matrix mode has been toogled " + calc.toggleMatrixMode());
                        break;
                case 0:
                    System.out.println("Goodbye!");
                        return;
                default: break;
            }
        }while(userIn != 0);
    
    }

    private String toggleMatrixMode () 
    { 
        if(matrix) {
            matrix = false;
            return "off."; 
        } else if(!matrix) {
            matrix = true;
            return "on.";
        } else return "ERROR: Did not toggle.";  
    }

    public void createMatrices(boolean pow2, boolean basic) {
        
        Random rand = new Random();

        matrix1 = new int[size][size];
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                int random = rand.nextInt(1 + ceiling) - 1;
                    matrix1[i][j] = random;
            }
        }
        if(!basic) matrix1 = appendMatrix(matrix1, size, pow2);
        if(matrix) printMatrix(matrix1, "Matrix 1", size);

        matrix2 = new int[size][size];
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                int random = 1 + rand.nextInt(1 + ceiling) - 1;
                matrix2[i][j] = random;
            }
        }
        if(!basic) matrix2 = appendMatrix(matrix2, size, pow2);
        if(matrix) printMatrix(matrix2, "Matrix 2", size);
    }

    public int[][] basicMult (int[][] A, int[][] B) {
        int[][] product = new int[A.length][A.length];
        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j < A.length; j++) {
                for (int k = 0; k < A.length; k++) {
                    product[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        return product;
    }

    //Using matrix1 as A and matrix2 as B to be closer to proof
    public int[][] strassenMult(int[][] A, int[][] B) {
        int n = A.length;
        int nhalf = n/2;
        int[][] result = new int[n][n];

        if (n == 1)
            result[0][0] = A[0][0] * B[0][0];

        else if (n == 2) {
            int[][] C = new int[2][2];

            //Strassen's method for 2x2;

            int m1 = (A[0][0] + A[1][1]) * (B[0][0] + B[1][1]);
            int m2 = (A[1][0] + A[1][1]) * B[0][0];
            int m3 = A[0][0] * (B[0][1] - B[1][1]);
            int m4 = A[1][1] * (B[1][0] - B[0][0]);
            int m5 = (A[0][0] + A[0][1]) * B[1][1];
            int m6 = (A[1][0] - A[0][0]) * (B[0][0] + B[0][1]);
            int m7 = (A[0][1] - A[1][1]) * (B[1][0] + B[1][1]);

            C[0][0] = m1 + m4 - m5 + m7;
            C[0][1] = m3 + m5;
            C[1][0] = m2 + m4;
            C[1][1] = m1 + m3 - m2 + m6;

            combine(C, result, 0, 0);
        }

        else {
            int[][] A11 = new int[nhalf][nhalf];
            int[][] A12 = new int[nhalf][nhalf];
            int[][] A21 = new int[nhalf][nhalf];
            int[][] A22 = new int[nhalf][nhalf];
            int[][] B11 = new int[nhalf][nhalf];
            int[][] B12 = new int[nhalf][nhalf];
            int[][] B21 = new int[nhalf][nhalf];
            int[][] B22 = new int[nhalf][nhalf];

            split(A, A11, 0, 0);
            split(A, A12, 0, nhalf);
            split(A, A21, nhalf, 0);
            split(A, A22, nhalf, nhalf);

            split(B, B11, 0, 0);
            split(B, B12, 0, nhalf);
            split(B, B21, nhalf, 0);
            split(B, B22, nhalf, nhalf);

            int[][] M1 = strassenMult(plus(A11, A22), plus(B11, B22));
            int[][] M2 = strassenMult(plus(A21, A22), B11);
            int[][] M3 = strassenMult(A11, minus(B12, B22));
            int[][] M4 = strassenMult(A22, minus(B21, B11));
            int[][] M5 = strassenMult(plus(A11, A12), B22);
            int[][] M6 = strassenMult(minus(A21, A11), plus(B11, B12));
            int[][] M7 = strassenMult(minus(A12, A22), plus(B21, B22));

            int[][] C11 = plus(minus(plus(M1, M4), M5), M7);
            int[][] C12 = plus(M3, M5);
            int[][] C21 = plus(M2, M4);
            int[][] C22 = plus(minus(plus(M1, M3), M2), M6);

            combine(C11, result, 0, 0);
            combine(C12, result, 0, nhalf);
            combine(C21, result, nhalf, 0);
            combine(C22, result, nhalf, nhalf);
        }
        return result;
    }

    public int[][] ladermanMult(int[][] A, int[][]B) {
        int n = A.length;
        int[][] result = new int[n][n];
        int nthird = n/3, n2thirds = 2 * nthird;

        if (n == 1) result[0][0] = A[0][0] * B[0][0];
        else if (n == 3) {
            int[][] C = new int[3][3];

            int m1 = (A[0][0] + A[0][1] + A[0][2] - A[1][0] - A[1][1] - A[2][1] - A[2][2])*B[1][1];
            int m2 = (A[0][0] - A[1][0])*(-B[0][1] + B[1][1]);
            int m3 = A[1][1] * (-B[0][0] + B[0][1] +B[1][0] - B[1][1] - B[1][2] - B[2][0] + B[2][2]);
            int m4 = (-A[0][0] + A[1][0] + A[1][1]) * (B[0][0] - B[0][1] + B[1][1]);
            int m5 = (A[1][0] + A[1][1]) * (-B[0][0] + B[0][1]);
            int m6 = A[0][0] * B[0][0];
            int m7 = (-A[0][0] + A[2][0] + A[2][1]) * (B[0][0] - B[0][2] + B[1][2]);
            int m8 = (-A[0][0] + A[2][0]) * (B[0][2] - B[1][2]);
            int m9 = (A[2][0] + A[2][1]) * (-B[0][0] + B[0][2]);
            int m10 = (A[0][0] + A[0][1] + A[0][2] - A[1][1] - A[1][2] - A[2][0] - A[2][1]) * B[1][2];
            int m11 = A[2][1] * (-B[0][0] + B[0][2] + B[1][0] - B[1][1] -B[1][2] - B[2][0] + B[2][1]);
            int m12 = (-A[0][2] + A[2][1] + A[2][2]) * (B[1][1] + B[2][0] - B[2][1]);
            int m13 = (A[0][2] - A[2][2]) * (B[1][1] - B[2][1]);
            int m14 = A[0][2] * B[2][0];
            int m15 = (A[2][1] + A[2][2]) * (-B[2][0] + B[2][1]);
            int m16 = (-A[0][2] + A[1][1] + A[1][2]) * (B[1][2] + B[2][0] - B[2][2]);
            int m17 = (A[0][2] - A[1][2]) * (B[1][2] - B[2][2]);
            int m18 = (A[1][1] + A[1][2]) * (-B[2][0] + B[2][2]);
            int m19 = A[0][1] * B[1][0];
            int m20 = A[1][2] * B[2][1];
            int m21 = A[1][0] * B[0][2];
            int m22 = A[2][0] * B[0][1];
            int m23 = A[2][2] * B[2][2];

            C[0][0] = m6 + m14 + m19;
            C[0][1] = m1 + m4 + m5 + m6 + m12 + m14 + m15; 
            C[0][2] = m6 + m7 + m9 + m10 + m14 + m16 + m18;
            C[1][0] = m2 + m3 + m4 + m6 + m14 + m16 + m17;
            C[1][1] = m2 + m4 + m5 + m6 + m20; 
            C[1][2] = m14 + m16 + m17 + m18 + m21;
            C[2][0] = m6 + m7 + m8 + m11 + m12 + m13 + m14;
            C[2][1] = m12 + m13 + m14 + m15 + m22; 
            C[2][2] = m6 + m7 + m8 + m9 + m23;

            combine(C, result, 0, 0);
        } else {
            int[][] A11 = new int[nthird][nthird];
            int[][] A12 = new int[nthird][nthird];
            int[][] A13 = new int[nthird][nthird];
            int[][] A21 = new int[nthird][nthird];
            int[][] A22 = new int[nthird][nthird];
            int[][] A23 = new int[nthird][nthird];
            int[][] A31 = new int[nthird][nthird];
            int[][] A32 = new int[nthird][nthird];
            int[][] A33 = new int[nthird][nthird];
            int[][] B11 = new int[nthird][nthird];
            int[][] B12 = new int[nthird][nthird];
            int[][] B13 = new int[nthird][nthird];
            int[][] B21 = new int[nthird][nthird];
            int[][] B22 = new int[nthird][nthird];
            int[][] B23 = new int[nthird][nthird];
            int[][] B31 = new int[nthird][nthird];
            int[][] B32 = new int[nthird][nthird];
            int[][] B33 = new int[nthird][nthird];

            //nthrid and n2thirds for positions, try on 3x3 and 9x9
            split(A, A11, 0, 0);
            split(A, A12, 0, nthird);
            split(A, A13, 0, n2thirds);
            split(A, A21, nthird, 0);
            split(A, A22, nthird, nthird);
            split(A, A23, nthird, n2thirds);
            split(A, A31, n2thirds, 0);
            split(A, A32, n2thirds, nthird);
            split(A, A33, n2thirds, n2thirds);

            split(B, B11, 0, 0);
            split(B, B12, 0, nthird);
            split(B, B13, 0, n2thirds);
            split(B, B21, nthird, 0);
            split(B, B22, nthird, nthird);
            split(B, B23, nthird, n2thirds);
            split(B, B31, n2thirds, 0);
            split(B, B32, n2thirds, nthird);
            split(B, B33, n2thirds, n2thirds);

            int[][] M1 = ladermanMult(minus(minus(minus(minus(plus(plus(A11, A12), A13), A21), A22), A32), A33), B22);
            int[][] M2 = ladermanMult(minus(A11, A21), minus(B22, B12));
            int[][] M3 = ladermanMult(A22, plus(minus(minus(minus(plus(minus(B12, B11), B21), B22), B23), B31), B33));
            int[][] M4 = ladermanMult(plus(minus(A21, A11), A22), plus(minus(B11, B12), B22));
            int[][] M5 = ladermanMult(plus(A21, A22), minus(B12, B11));
            int[][] M6 = ladermanMult(A11, B11);
            int[][] M7 = ladermanMult(plus(minus(A31, A11), A32), plus(minus(B11, B13), B23));
            int[][] M8 = ladermanMult(minus(A31, A11), minus(B13, B23));
            int[][] M9 = ladermanMult(plus(A31, A32), minus(B13, B11));
            int[][] M10 = ladermanMult(minus(minus(minus(minus(plus(plus(A11, A12), A13), A22), A23), A31), A32), B23);
            int[][] M11 = ladermanMult(A32, plus(minus(minus(minus(plus(minus(B13, B11), B21), B22), B23), B31), B32));
            int[][] M12 = ladermanMult(plus(minus(A32, A13), A33), minus(plus(B22, B31), B32));
            int[][] M13 = ladermanMult(minus(A13, A33), minus(B22, B32));
            int[][] M14 = ladermanMult(A13, B31);
            int[][] M15 = ladermanMult(plus(A32, A33), minus(B32, B31));
            int[][] M16 = ladermanMult(plus(minus(A22, A13), A23), minus(plus(B23, B31), B33));
            int[][] M17 = ladermanMult(minus(A13, A23), minus(B23, B33));
            int[][] M18 = ladermanMult(plus(A22, A23), minus(B33, B31));
            int[][] M19 = ladermanMult(A12, B21);
            int[][] M20 = ladermanMult(A23, B32);
            int[][] M21 = ladermanMult(A21, B13);
            int[][] M22 = ladermanMult(A31, B12);
            int[][] M23 = ladermanMult(A33, B33);

            int[][] C11 = plus(plus(M6, M14), M19);
            int[][] C12 = plus(plus(plus(plus(plus(plus(M1, M4), M5), M6), M12), M14), M15);
            int[][] C13 = plus(plus(plus(plus(plus(plus(M6, M7), M9), M10), M14), M16), M18);
            int[][] C21 = plus(plus(plus(plus(plus(plus(M2, M3), M4), M6), M14), M16), M17);
            int[][] C22 = plus(plus(plus(plus(M2, M4), M5), M6), M20);
            int[][] C23 = plus(plus(plus(plus(M14, M16), M17), M18), M21);
            int[][] C31 = plus(plus(plus(plus(plus(plus(M6, M7), M8), M11), M12), M13), M14);
            int[][] C32 = plus(plus(plus(plus(M12, M13), M14), M15), M22);
            int[][] C33 = plus(plus(plus(plus(M6, M7), M8), M9), M23);

            combine(C11, result, 0, 0);

            combine(C12, result, 0, nthird);
            combine(C13, result, 0, n2thirds);
            combine(C21, result, nthird, 0);
            combine(C22, result, nthird, nthird);
            combine(C23, result, nthird, n2thirds);
            combine(C31, result, n2thirds, 0);
            combine(C32, result, n2thirds, nthird);
            combine(C33, result, n2thirds, n2thirds);
        }
        if(matrix) if(result.length > size ) printMatrix(result, "Result", size);
        return result;
    }

    public void split (int[][] matrix, int[][] splitMatrix, int i2Pos, int j2Pos) {
            for(int i1 = 0, i2 = i2Pos; i1 < splitMatrix.length; i1++, i2++)
                for(int j1 = 0, j2 = j2Pos; j1 < splitMatrix.length; j1++, j2++)
                    splitMatrix[i1][j1] = matrix[i2][j2];
    }

    public void combine(int[][] subMatrix, int[][] result, int i2Pos, int j2Pos)
    {
        for(int i1 = 0, i2 = i2Pos; i1 < subMatrix.length; i1++, i2++)
            for(int j1 = 0, j2 = j2Pos; j1 < subMatrix.length; j1++, j2++)
                result[i2][j2] = subMatrix[i1][j1];
    }

    public int[][] plus (int[][] matrix1, int[][] matrix2) {
        int[][] temp = new int[matrix1.length][matrix1.length];
        for(int i = 0; i < matrix1.length; i++) {
            for (int j = 0; j < matrix1.length; j++) {
                temp[i][j] = matrix1[i][j] + matrix2[i][j];
            }
        }
        return temp;
    }

    public int[][] minus (int[][] matrix1, int[][] matrix2) {
        int[][] temp = new int[matrix1.length][matrix1.length];
        for(int i = 0; i < matrix1.length; i++) {
            for (int j = 0; j < matrix1.length; j++) {
                temp[i][j] = matrix1[i][j] - matrix2[i][j];
            }
        }
        return temp;
    }

    // Function to check if x is power of 2
    boolean checkPow(int n, int pow)
    {
        if(n==0)
            return false;
        return (int)(Math.ceil((Math.log(n) / Math.log(pow)))) ==
                (int)(Math.floor(((Math.log(n) / Math.log(pow)))));
    }

    private void fixResult(int[][] result) {
        for (int i = 0; i < result.length; i++) {
            for (int j = 1; j < result.length; j++) {
                if (i % 2 == 0 && j % 2 != 0) {
                    int swap = result[i][j];
                    result[i][j] = result[i + 1][j - 1];
                    result[i + 1][j - 1] = swap;
                }
            }
        }
    }

    void timeTrial () {
        int x = 0;
        int[][] result;
        MatrixCalc calc = new MatrixCalc();

        while (x < runCeiling) {
            calc.createMatrices(true, true);
            long basicStartTime = System.nanoTime();
            result = calc.basicMult(matrix1, matrix2);
                if(matrix)calc.printMatrix(result, "\nBasic Matrix", size);
                durationBasic += basicStartTime - System.nanoTime();
            
            calc.createMatrices(true, false);
            long strassenStartTime = System.nanoTime();
                result = calc.strassenMult(matrix1, matrix2);
                    if(matrix) calc.printMatrix(result, "Strassen Matrix", size);
                    durationStrassen += strassenStartTime - System.nanoTime();

            calc.createMatrices(false, false);
            long ladermanStartTime = System.nanoTime();
                result = calc.ladermanMult(matrix1, matrix2);
                if(matrix) calc.printMatrix(result, "\nLaderman Matrix", size);
                durationLaderman += ladermanStartTime - System.nanoTime();
            x++;
        }
        durationBasic *= (-1);
        durationStrassen *= (-1);
        durationLaderman *= (-1);
            System.out.println("\nBasic Method\n-----------------\nTotal runs: " + x + "\nRuntime = " + (double)durationBasic/1000000000 + " s, " + durationBasic/1000000 + " ms");
            System.out.println("\nStrassen's Method\n-----------------\nTotal runs: " + x + "\nRuntime = " + (double)durationStrassen/1000000000 + " s, " + durationStrassen/1000000 + " ms");
            System.out.println("\nLaderman's Method\n-----------------\nTotal runs: " + x + "\nRuntime = " + (double)durationLaderman/1000000000 + " s, " + durationLaderman/1000000 + " ms\n\n");
    }

    int[][] appendMatrix(int[][] matrix, int size, boolean pow2) {
        int count = 0;
        int newSize = 0, comp, power;
        int[] twoPowArr = new int[8], threePowArr = new int[6];
        
        if(pow2) {
            newSize = 2;

            twoPowArr[0] = Math.abs(2 - size);
            twoPowArr[1] = Math.abs(4 - size);
            twoPowArr[2] = Math.abs(8 - size);
            twoPowArr[3] = Math.abs(16 - size);
            twoPowArr[4] = Math.abs(32 - size);
            twoPowArr[5] = Math.abs(64 - size);
            twoPowArr[6] = Math.abs(128 - size);
            twoPowArr[7] = Math.abs(256 - size);
            
            comp = twoPowArr[0];
            while (count < twoPowArr.length) {
                if (twoPowArr[count] < comp) {
                    comp = twoPowArr[count];
                    power = count + 1;
                        newSize = (int)Math.pow(2.0, (double)count + 1);
                }
                //System.out.println("Val: " + powArr[count] +" Count: " + count + "\nNew Size: " + newSize);
                count++;
            }
        if (size > newSize) newSize *= 2;
        } else if(!pow2) {
            newSize = 3;
            threePowArr[0] = Math.abs(3 - size);
            threePowArr[1] = Math.abs(9 - size);
            threePowArr[2] = Math.abs(27 - size);
            threePowArr[3] = Math.abs(81 - size);
            threePowArr[4] = Math.abs(243 - size);
            threePowArr[5] = Math.abs(729 - size);
            
            comp = threePowArr[0];
            while (count < threePowArr.length) {
                if (threePowArr[count] < comp) {
                    comp = threePowArr[count];
                    power = count + 1;
                        newSize = (int)Math.pow(3.0, (double)count + 1);
                }
                //System.out.println("Val: " + powArr[count] +" Count: " + count + "\nNew Size: " + newSize);
                count++;
            }
        if (size > newSize) newSize *= 3;
        }
            //System.out.println("Size: " + size + "\nNew Size: " + newSize);

            int[][] append = new int [newSize][newSize];
            for (int i = 0; i < newSize; i++) {
                for (int j = 0; j < newSize; j++) {
                    append[i][j] = 0;
                }
            }

            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    append[i][j] = matrix[i][j];
                }
            }
            //printMatrix(append, "Appended Matrix: ", newSize);

        return append;
    }

    public void printMatrix(int[][] matrix, String matrixID, int size) {
        System.out.println(matrixID);
        System.out.println("-----------------------------------------------------------------------");
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                System.out.print(matrix[i][j] + " ");
            }
            System.out.println("");
        }
    }
}

/**public void outputToFile() {
 try {
 File myObj = new File("output.txt");
 if (myObj.createNewFile()) {
 System.out.println("File created: " + myObj.getName());
 } else {
 System.out.println("File already exists.");
 }
 } catch (IOException e) {
 System.out.println("An error occurred.");
 e.printStackTrace();
 }
             printMatrix(A11, "A11", A11.length);
            printMatrix(A12, "A12", A12.length);
            printMatrix(A13, "A13", A13.length);
            printMatrix(A21, "A21", A21.length);
            printMatrix(A22, "A22", A22.length);
            printMatrix(A23, "A23", A23.length);
            printMatrix(A31, "A31", A31.length);
            printMatrix(A32, "A32", A32.length);
            printMatrix(A33, "A33", A33.length);

            
    /** FOR NEXT TIME
     *     Changing matrix testing loop
     *        let m = matrix size
     *        m    Basic         Strassen       Laderman
     *        -----------------------------------------------
     *        2    100(2x2)      100(2x2)       N/A
     *        ...
     *        10   100 (10x10)   100 (10x10)    100 (10x10)
     *        20   100 (20x20)   100 (20x20)    100 (20x20)
     *        30   etc...
     *        ...
     *        100
     *        -----------------------------------------------
     *        Keep track of amount of time to do 100 cases
     *        '
     *        ALSO FOR SPLITTING MATRICES
     *        test with 3x3 matricies, HAVE TO MAKE MATRIX CREATION WORK FOR 3X3 AND POWER CHECKING ETC..
     */

    //CHECK RANDOM 10X10 MATRIX TO SEE IF CALCULATIONS ARE CORRECT