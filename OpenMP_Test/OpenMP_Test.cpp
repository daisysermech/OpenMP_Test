#include <iostream>
#include <omp.h>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <chrono>
#include <fstream>
//#pragma omp parallel

long long GCD(long long a, long long b)
{
    while (b != 0)
    {
        long long r = a % b;
        a = b;
        b = r;
    }
    return a;
}

bool is_sqrt(long long a)
{
    long long i = sqrt(a); 
    if (i * i == a) return true; else
    return false;
}

long long Mod(long long a, long long b)
{
    if (b < 0) return -Mod(-a, -b);
    long long res = a % b;
    return res < 0 ? res + b : res;
}

long long GCD_Ext(long long a, long long b, long long& x, long long& y)
{
    if (a == 0) { x = 0; y = 1; return b; }
    long long x1 = 0, y1 = 0;
    long long d = GCD_Ext(b % a, a, x1, y1);
    x = y1 - (b / a) * x1;
    y = x1;
    return d;
}

long long find_x(long long* equations, long long N, long long& modulus)
{
    long long modulo = 1;

    #pragma omp parallel for
    for (long long i = 0; i < N; i++)
    {
        if (equations[i * 3] != 1)
        {
            if (equations[i * 3 + 1] >= equations[i * 3 + 2])
                equations[i * 3 + 1] = Mod(equations[i * 3 + 1], equations[i * 3 + 2]);
            long long inv1, inv2;
            GCD_Ext(equations[i * 3], equations[i * 3 + 2], inv1, inv2);
            equations[i * 3 + 1] = Mod(equations[i * 3 + 1] * Mod(inv1, equations[i * 3 + 2]), equations[i * 3 + 2]);
        }
        modulo *= equations[i * 3 + 2];
    }
    modulus = modulo;

    long long* inv = new long long[N];

    #pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        long long x, y;
        GCD_Ext(modulo / equations[i * 3 + 2], equations[i * 3 + 2], x, y);
        x = Mod(x, equations[i * 3 + 2]);
        inv[i] = x;
    }

    long long res = 0;

    #pragma omp parallel for
    for (int i = 0; i < N; i++)
        res += Mod(equations[i * 3 + 1] * inv[i] * (modulo / equations[i * 3 + 2]), modulo);
    return Mod(res, modulo);
}

bool check(long long* equations, long long N)
{
    //#pragma omp parallel for collapse(2)
        /*for (long long i = 0; i < N - 1; i++)
            for (long long j = i + 1; j < N; j++)
                if (GCD(equations[3 * j + 2], equations[3 * i + 2]) != 1)
                    return false;*/
    bool res = true;
#pragma omp parallel for
    for (long long ij = 0; ij < N * (N - 1); ij++)
    {
        if (!res) continue;
        if ((ij / (N - 1) != (ij % (N - 1))) && (GCD(equations[3 * (ij / (N - 1)) + 2], equations[3 * (ij % (N - 1)) + 2]) != 1))
            res = false;

    }
    return res;
}

int main(int argc, char* argv[])
{
    printf("Enter the number of threads:\n");
    int Thread_Num;
    std::cin >> Thread_Num;
    omp_set_num_threads(Thread_Num);
    long long N, modulus, x;
    double time;
    long long* equations;

    std::cout << "Enter n:\n";
    std::cin >> N;

    equations = (long long*)malloc((N * 3) * sizeof(long long));
    char y = NULL;
    std::cout << "File input or manual? f/m\n";
    std::cin >> y;
    if (y == 'm')
    {
        for (int i = 0; i < N; i++)
        {
            std::cout << "Enter " << i + 1 << "th equation coefficients:\n";
            std::cin >> equations[i * 3] >> equations[i * 3 + 1] >> equations[i * 3 + 2];
        }
    }
    else
        if (y == 'f')
        {
            std::ifstream in;
            std::string file;
            while (!in.is_open())
            {
                std::cout << "Enter file name.\n";
                std::cin >> file;
                in.open("C:\\temp\\" + file);
            }

            for (int i = 0; i < N; i++)
            {
                long long temp1, temp2, temp3;
                in >> temp1 >> temp2 >> temp3;
                equations[3 * i] = temp1;
                equations[3 * i + 1] = temp2;
                equations[3 * i + 2] = temp3;
                //in >> equations[3 * i] >> equations[3 * i + 1] >> equations[3 * i + 2];
            }

            in.close();
        }
        else
        {
            std::cout << "Wrong input, program shall exit.\n";
            return -2;
        }

    auto begin = std::chrono::high_resolution_clock::now();

    if (!check(equations, N))
        return -1;
    x = find_x(equations, N, modulus);

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);

    std::cout << "Res is " << x << " mod " << modulus << " - calculated for " << elapsed.count() * 1e-9 << "\n";

    delete[] equations;
}