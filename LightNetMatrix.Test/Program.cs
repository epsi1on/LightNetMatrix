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
            TestInverse(1000);

            TestLargeMultiply(1200, 1200, 1200, true);

            TestInitiation(20000, 20000);

            Console.ReadKey();
        }

        static void TestInverse(int n)
        {
            Matrix m = RandomMatrix(n, n);

            //for (int i = 0; i < n; i++)
            //    m[i, i] = 0.0;


            var sp = System.Diagnostics.Stopwatch.StartNew();
            var invm = m.Inverse();
            var trM = m.Transpose();

            Console.WriteLine("Inverse of {0}x{0} matrix tooks {1} milisecs", n, sp.ElapsedMilliseconds);

            sp.Restart();
            var expectedEye = m * invm;
            Console.WriteLine("multiply of two {0}x{0} matrixes tooks {1} milisecs", n, sp.ElapsedMilliseconds);

            sp.Restart();
            var expectedEye2 = invm * m;
            Console.WriteLine("multiply of two {0}x{0} matrixes tooks {1} milisecs", n, sp.ElapsedMilliseconds);

            sp.Stop();

            var z = expectedEye - Matrix.Eye(n);//expected zero
            var z2 = Matrix.Eye(n) - expectedEye;//expected zero

            var err = double.MinValue;

            int maxI = -1, maxJ = -1;

            for (var i = 0; i < n; i++)
                for (var j = 0; j < n; j++)
                {
                    if (Math.Abs(z[i, j]) > err)
                        err = Math.Abs(z[maxI = i, maxJ = j]);

                    if (Math.Abs(z2[i, j]) > err)
                        err = Math.Abs(z2[maxI = i, maxJ = j]);
                }
                    

            //var t2 = (m * m.Inverse() - Matrix.Eye(n)).coreArray.Select(i => i.Select(j =>) Math.Abs(i)).Max();

            Console.WriteLine("Maximum residual of M.inv * M : {0} at row:{1}, col:{2}", err, maxI, maxJ);

            sp.Restart();
        }

        static void TestLargeMultiply(int n,int m,int l,bool tryNormal)
        {
            Matrix mtx1 = RandomMatrix(n, m);
            Matrix mtx2 = RandomMatrix(m, l);

            //for (int i = 0; i < n; i++)
            //    m[i, i] = 0.0;

            var sp1 = 0.0;
            var sp2 = 0.0;

            var sp = System.Diagnostics.Stopwatch.StartNew();

            if(tryNormal)
            {
                var mlp1 = Matrix.Multiply(mtx1, mtx2);
                Console.WriteLine($"multiply of {n}x{m} matrix by a {m}x{l} matrix tooks {sp1 = sp.ElapsedMilliseconds} milisecs");
                sp.Restart();
            }
            

            var mlp2 = Matrix.FastMultiply(mtx1, mtx2);
            Console.WriteLine($"fast multiply of {n}x{m} matrix by a {m}x{l} matrix tooks {sp2 = sp.ElapsedMilliseconds} milisecs");

            if (tryNormal)
            {
                var ratio = sp1 / sp2;

                Console.WriteLine("fast multiply is {0:00} times faster than normal multiply!", ratio);
            }
            
            
        }

        static void TestInitiation(int n,int m)
        {
            long memory = GC.GetTotalMemory(true);

            Console.WriteLine($"memory in use {memory.ToByteSize()}");

            

            var mtx = new Matrix(n, m);
            Console.WriteLine($"matrix with size {n}x{m} initiated successfully");
            memory = GC.GetTotalMemory(true);
            Console.WriteLine($"memory in use {memory.ToByteSize()}");
        }

        static void TestCode()
        {
            Matrix m = RandomMatrix(5, 6);

            m[3, 4] = 1.0;
            

            Matrix inv = m.Inverse();
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
