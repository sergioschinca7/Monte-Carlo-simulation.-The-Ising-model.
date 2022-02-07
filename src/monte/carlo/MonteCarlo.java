package monte.carlo;

import java.io.*;

public class MonteCarlo {
    private static final int J = 1;
    private static int L = 20; 
    private static int E;
    private static int M;
    private static double aceptacion;
    private static int red[][] = new int[L][L];

    public static void main(String[] args) throws IOException {

        float E, M;
        //rango de temperatura
        double T_max = 5.0;
        double T_min = 0.01;
        double T_step = 0.01;
        //sampleo
        int L = 20;
        int pasos_termalizacion = 4 * L * L;
        int n_sample = 100;
        int dist_sample = 2 * L * L;
        int simu = 1;
        
        File file = new File("termalizacion.dat");
        File file2 = new File("Output.dat");
        FileWriter fw = new FileWriter(file);
        FileWriter fw1 = new FileWriter(file2);
        PrintWriter pw = new PrintWriter(fw);
        PrintWriter pw1 = new PrintWriter(fw1);
        pw.println("# Termalizacion: " + "J = "+ J+"," +"\t"+ "pasos_termalizacion = "+pasos_termalizacion+"," +"\t"+ "n_sample = "+
                n_sample+"," +"\t"+ "dist_sample = "+dist_sample);
        //pw.println("L" + "T" + "simu" + "mean_M" + "mean_E" + "mean_M2" + "mean_E2" + "Aceptacion");
        pw1.println("L\t" + "T\t" + "simu\t" + "mean_M\t"+"mean_E\t"+ "mean_M2\t"+"mean_E2\t" + "Aceptacion" );

        for (int i = 0; i < simu; i++) {
            double p = 0.5;
            llenar(p);
            E = calcular_E_inicial();
            M = calcular_M_inicial();
            //consideramos diferentes temperaturas
            int cant_T = (int) ((T_max - T_min) / T_step + 1);
            for (int i_T = 0; i_T < cant_T; i_T++) {
                double T = T_max - T_min * T_step;
                //termalizamos
                for (int paso = 0; paso < pasos_termalizacion; paso++) {
                    metropolis(T);
                }
                //Sampleo
                double mean_M = 0; // magnetizacion media
                double mean_E = 0; // energia media
                double mean_M2 = 0; // magnetizacion cuadrada media
                double mean_E2 = 0; // energia cuadrada media
                aceptacion = 0; // prop de pasos de MM aceptados

                for (int sample = 0; sample < n_sample; sample++) {
                    //descorrelacionamos las muestras
                    for (int paso = 0; paso < dist_sample; paso++) {
                        metropolis(T);
                    }
                    mean_M = mean_M + M;
                    mean_E = mean_E + E;
                    mean_M2 = mean_M2 + (M * M);
                    mean_E2 = mean_E2 + (E * E);
                    pw.println("simu\t" + simu+"\t" + "T\t" + T+"\t" + "sample\t" + sample);
                }
                mean_M = mean_M / (double) n_sample;
                mean_E = mean_E / (double) n_sample;
                mean_M2 = mean_M2 / (double) n_sample;
                mean_E2 = mean_E2 / (double) n_sample;
                aceptacion = aceptacion / (double) (n_sample * dist_sample);

                pw1.println(L +"\t" + T +"\t" + simu +"\t" + mean_M +"\t" + mean_M2 +"\t" + mean_E +"\t" + mean_E2 +"\t" + aceptacion);
            }
        }
        pw.close();
        pw1.close();
    } // fin main

    public static float calcular_E_inicial() {
        //E = -J sumatoria (si*sj)
        float E = 0;
        int columna_derecha, fila_abajo;
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                fila_abajo = (i + 1) % L;
                columna_derecha = (j + 1) % L;

                E = E + red[i][j] * (red[i][columna_derecha] + red[fila_abajo][j]);
            }
        }
        E = -J * E;
        return E;
    }

    public static float calcular_M_inicial() {
        //M = sum (s_i)
        float resultado = 0;
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                resultado = resultado + red[i][j];
            }
        }
        return resultado;
    }

    public static void metropolis(double T) {
        int fila = (int) ((1 + Math.random() * 21) % L);
        int columna = (int) ((1 + Math.random() * 21) % L);
        flip(fila, columna, T);
    }

    public static void flip(int fila, int columna, double T) {
        int fila_abajo = (fila + 1) % L;
        int fila_arriba = (fila - 1 + L) % L;
        int columna_derecha = (columna + 1) % L;
        int columna_izquierda = (columna - 1 + L) % L;

        int suma_vecinos = red[fila_abajo][columna] + red[fila_arriba][columna]
                + red[fila][columna_derecha] + red[fila][columna_izquierda];
        int dE = 2 * J * red[fila][columna] * suma_vecinos;

        if (dE < 0) {
            red[fila][columna] = -red[fila][columna];
            //actualizo valores de magnetizacion y energia
            M = M + 2 * red[fila][columna];
            E = E + dE;
            aceptacion++;
        } else {
            double exponente = (double) (dE / T);
            double probabilidad = Math.exp(exponente);
            double ratio = Math.random();
            if (ratio <= probabilidad) {
                red[fila][columna] = -red[fila][columna];
                //como se acepta el paso se actualzan valores de mag y energia
                M = M + 2 * red[fila][columna];
                E = E + dE;
                aceptacion++;

            }
        }
    }

    public static void llenar(double p) {
        double ratio;
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                ratio = Math.random();
                if (ratio < p) {
                    red[i][j] = 1;
                } else {
                    red[i][j] = -1;
                }
            }
        }
    }

    public static void mostrar() {
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                System.out.print(red[i][j]);
            }
            System.out.println("");
        }

    }
} //fin class
