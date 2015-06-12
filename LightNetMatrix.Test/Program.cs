using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace LightNetMatrix.Test
{
    class Program
    {
        static void Main(string[] args)
        {
            TestInverse(300);

            Console.ReadKey();
        }

        static void TestInverse(int n)
        {
            var m = RandomMatrix(n, n);
            for (int i = 0; i < n; i++)
            {
                m[i, i] = 0.0;
            }

            var sp = System.Diagnostics.Stopwatch.StartNew();
            m.Determinant();
            sp.Stop();
            var invm = m.Inverse();
            Console.WriteLine("Inverse of {0}x{0} matrix tooks {1} milisecs", n, sp.ElapsedMilliseconds);
            
            

            var t2 = (m * m.Inverse() - Matrix.Eye(n)).coreArray.Select(i => Math.Abs(i)).Max();

            Console.WriteLine("Maximum residual: {0}", t2);
        }


        static Matrix RandomMatrix(int m, int n)
        {
            var buf = new Matrix(m, n);

            var rnd = new Random(0);


            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    buf[i, j] = rnd.NextDouble() * 100;
                }
            }

            return buf;
        }

    }
}
