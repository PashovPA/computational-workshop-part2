package com.pashovpa;

import java.util.Arrays;

public class FirstTask {

  public double[] solveGauss(double[][] matrixInitAB) {
    int height = matrixInitAB.length;
    int width = matrixInitAB[0].length;
    double[][] matrixAB = new double[height][width];
    for (int i = 0; i < height; i++) {
      for (int j = 0; j < width; j++) {
        matrixAB[i][j] = matrixInitAB[i][j];
      }
    }

    // Прямой ход
    for (int k = 0; k < height; k++) {
      int maxLine = k;
      double maxAbsNum = 0;

      // Поиск главного элемента по стобцу
      for (int j = k; j < height; j++) {
        if (maxAbsNum < Math.abs(matrixAB[j][k])) {
          maxAbsNum = matrixAB[j][k];
          maxLine = j;
        }
      }

      // Перестановка строк
      double[] temp = matrixAB[k];
      matrixAB[k] = matrixAB[maxLine];
      matrixAB[maxLine] = temp;

      // Рассматриваем случай, когда главный элемент не равен 0
      if (Math.abs(matrixAB[k][k]) > 0) {

        // Преобразование главного элемента к единице, разделив его строку на его коэффициент
        double weightMainElem = matrixAB[k][k];
        for (int i = 0; i < width; i++) {
          matrixAB[k][i] /= weightMainElem;
        }

        // Получение нуля под главным элементом
        for (int i = k + 1; i < height; i++) {
          double weightColumn = matrixAB[i][k];
          for (int j = k; j < width; j++) {
            matrixAB[i][j] -= matrixAB[k][j] * weightColumn;
          }
        }
      }
    }

    // Обратный ход
    double[] result = new double[width - 1];
    Arrays.fill(result, 0);
    for (int i = height - 1; i >= 0; i--) {
      double tempResult = 0;
      for (int j = i + 1; j < height; j++) {
        tempResult += matrixAB[i][j] * result[j];
      }
      result[i] = matrixAB[i][height] - tempResult;
    }
    return result;
  }

  public double[][] inverseJordan(double[][] matrixInitA) {
    int height = matrixInitA.length;
    int width = matrixInitA[0].length;
    double[][] matrixAE = new double[height][2 * width];
    int heightAE = matrixAE.length;
    int widthAE = matrixAE[0].length;
    for (int i = 0; i < heightAE; i++) {
      Arrays.fill(matrixAE[i], 0);
      matrixAE[i][i + width] = 1;
      for (int j = 0; j < width; j++) {
        matrixAE[i][j] = matrixInitA[i][j];
      }
    }

    // Прямой ход
    for (int k = 0; k < heightAE; k++) {
      int maxLine = k;
      double maxAbsNum = 0;

      // Поиск главного элемента по стобцу
      for (int j = k; j < heightAE; j++) {
        if (maxAbsNum < Math.abs(matrixAE[j][k])) {
          maxAbsNum = matrixAE[j][k];
          maxLine = j;
        }
      }

      // Перестановка строк
      double[] temp = matrixAE[k];
      matrixAE[k] = matrixAE[maxLine];
      matrixAE[maxLine] = temp;

      // Рассматриваем случай, когда главный элемент не равен 0
      if (Math.abs(matrixAE[k][k]) > 0) {

        // Преобразование главного элемента к единице, разделив его строку на его коэффициент
        double weightMainElem = matrixAE[k][k];
        for (int i = 0; i < widthAE; i++) {
          matrixAE[k][i] /= weightMainElem;
        }

        // Получение нуля под и над главным элементом
        for (int i = 0; i < heightAE; i++) {
          if (i == k) continue;
          double weightColumn = matrixAE[i][k];
          for (int j = k; j < widthAE; j++) {
            matrixAE[i][j] -= matrixAE[k][j] * weightColumn;
          }
        }
      }
    }
    double[][] result = new double[height][width];
    for (int i = 0; i < height; i++) {
      for (int j = 0; j < width; j++) {
        result[i][j] = matrixAE[i][j + width];
      }
    }

    return result;
  }

  public double detLU(double[][] matrixInitA) {
    int height = matrixInitA.length;
    int width = matrixInitA[0].length;
    double[][] matrixL = new double[height][width];
    double[][] matrixU = new double[height][width];

    for (int i = 0; i < height; i++) {
      for (int j = 0; j < height; j++) {
        if (i > j) {
          matrixU[i][j] = 0;
        } else if (i == j) {
          matrixU[i][j] = 1;
        } else {
          matrixL[i][j] = 0;
        }
      }
    }

    double tempL, tempU;
    for (int i = 0; i < height; i++) {
      for (int j = i; j < height; j++) {
        tempL = 0;
        for (int k = 0; k < i; k++) {
          tempL += matrixL[j][k] * matrixU[k][i];
        }
        matrixL[j][i] = matrixInitA[j][i] - tempL;
      }

      for (int j = i; j < height; j++) {
        tempU = 0;
        for (int k = 0; k < i; k++) {
          tempU += matrixL[i][k] * matrixU[k][j];
        }
        matrixU[i][j] = (matrixInitA[i][j] - tempU) / matrixL[i][i];
      }
    }

    double result = 1;
    for (int i = 0; i < height; i++) {
      result *= matrixL[i][i];
    }
    return result;
  }

  public double cond(double[][] matrixInitA) {
    int height = matrixInitA.length;
    int width = matrixInitA[0].length;
    double[][] matrixA = new double[height][width];
    for (int i = 0; i < height; i++) {
      for (int j = 0; j < width; j++) {
        matrixA[i][j] = matrixInitA[i][j];
      }
    }
    double[][] matrixAInverse = inverseJordan(matrixA);

    return normMatrix(matrixA) * normMatrix(matrixAInverse);
  }

  public double fault(double[] result1, double[] result2) {
    double[] diffResult = new double[result1.length];
    for (int i = 0; i < result1.length; i++) {
      diffResult[i] = result2[i] - result1[i];
    }
    return normVector(diffResult) / normVector(result1);
  }

  public double normMatrix(double[][] matrix) {
    double temp = 0;
    for (int i = 0; i < matrix.length; i++) {
      for (int j = 0; j < matrix[0].length; j++) {
        temp += Math.pow(matrix[i][j], 2);
      }
    }
    return Math.sqrt(temp);
  }

  public double normVector(double[] vector) {
    double temp = 0;
    for (int i = 0; i < vector.length; i++) {
      temp += Math.pow(vector[i], 2);
    }
    return Math.sqrt(temp);
  }

  public static void main(String[] args) {
    System.out.println("Задание 1.1. Вариант 5.");
    System.out.println("===================================");
    double[][] matrixFirstA = {{8.673134, 1.041039, -2.677712}, {1.041039, 6.586211, 0.623016}, {-2.677712, 0.623016, 5.225935}};
    double[][] matrixFirstAB = {{8.673134, 1.041039, -2.677712, -1.289879}, {1.041039, 6.586211, 0.623016, 4.020225}, {-2.677712, 0.623016, 5.225935, 5.269671}};
    double[][] matrixSecondA = {{-401.46, 200.18}, {1201.08, -601.64}};
    double[][] matrixSecondAB1 = {{-401.46, 200.18, 200}, {1201.08, -601.64, -600}};
    double[][] matrixSecondAB2 = {{-401.46, 200.18, 199}, {1201.08, -601.64, -601}};
    FirstTask test = new FirstTask();
    System.out.println("Решение СЛАУ методом Гаусса с выбором главного жлемента по столбцу: ");
    var resultGauss = test.solveGauss(matrixFirstAB);
    for (int i = 0; i < resultGauss.length; i++) {
      System.out.println("x" + (i + 1) + " = " + resultGauss[i]);
    }
    System.out.println("===================================");

    System.out.println("Обратная матрица методом Жордана: ");
    var resultJordan = test.inverseJordan(matrixFirstA);
    for (int i = 0; i < resultJordan.length; i++) {
      for (int j = 0; j < resultJordan[0].length; j++) {
        System.out.print(resultJordan[i][j] + " ");
      }
      System.out.println();
    }
    System.out.println("===================================");

    System.out.println("Определитель матрицы методом LU - разложения: ");
    System.out.println("det(A) = " + test.detLU(matrixFirstA));
    System.out.println("===================================");

    System.out.println("Задание 1.2. Вариант 5.");
    System.out.println("Решение СЛАУ Ax = b: ");
    var result1 = test.solveGauss(matrixSecondAB1);
    for (int i = 0; i < result1.length; i++) {
      System.out.println("x" + (i + 1) + " = " + result1[i]);
    }
    System.out.println("===================================");

    System.out.println("Решение СЛАУ Ax' = b': ");
    var result2 = test.solveGauss(matrixSecondAB2);
    for (int i = 0; i < result2.length; i++) {
      System.out.println("x" + (i + 1) + " = " + result2[i]);
    }
    System.out.println("===================================");

    System.out.println("Число обусловленности матрицы A: ");
    double condA = test.cond(matrixSecondA);
    System.out.println("cond(A) = " + condA);
    System.out.println("Фактическая относительная погрешность: ");
    System.out.println("бx = " + test.fault(result1, result2));

    double[] deltaB = {-1, -1};
    double[] B1 = {200, -600};
    double[] B2 = {199, - 601};
    System.out.println("Оценка фактической относительной погрешности: ");
    System.out.println("бx <= " + condA * (test.fault(B1, B2)));
    System.out.println("===================================");
  }
}
