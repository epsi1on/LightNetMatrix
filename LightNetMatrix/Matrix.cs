using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Runtime.Serialization;
using System.Text;
using System.Xml.Serialization;
using System.Linq;

namespace LightNetMatrix
{
    /// <summary>
    /// Represents a two dimensional dense matrix of real (double precision) elements.
    /// </summary>
    [DebuggerDisplay("Matrix {RowCount} x {ColumnCount}")]
    //[Serializable]
    [CLSCompliant(true)]
    public class Matrix
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="Matrix"/> class with specified dimensions.
        /// </summary>
        /// <param name="rowCount">The row count.</param>
        /// <param name="columnCount">The column count.</param>
        public Matrix(int rowCount, int columnCount)
        {
            this.rowCount = rowCount;
            this.columnCount = columnCount;

            this.coreArray = new double[rowCount][];

            for (var i = 0; i < rowCount; i++)
                this.coreArray[i] = new double[columnCount];
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="Matrix"/> class of specified dimension.
        /// </summary>
        /// <param name="dimension">The dimension.</param>
        public Matrix(int dimension) : this(dimension, dimension)
        {
        }

        private Matrix()
        {

        }

        public double[][] coreArray;

        /// <summary>
        /// Gets the number of rows of the matrix.
        /// </summary>

        public int RowCount
        {
            get { return rowCount; }
        }

        /// <summary>
        /// Gets the number of columns of the matrix.
        /// </summary>
        public int ColumnCount
        {
            get { return columnCount; }
        }

        /// <summary>
        /// Number of rows of the matrix.
        /// </summary>
        [XmlElement(ElementName = "RowCount")]
        private int rowCount;

        /// <summary>
        /// Number of columns of the matrix.
        /// </summary>
        [XmlElement(ElementName = "ColumnCount")]
        private int columnCount;

        /// <summary>
        /// Gets or sets the specified member.
        /// </summary>
        /// <value>
        /// The <see cref="System.Double"/>.
        /// </value>
        /// <param name="row">The row (zero based).</param>
        /// <param name="column">The column (zero based).</param>
        /// <returns></returns>
        [System.Runtime.CompilerServices.IndexerName("TheMember")]
        public double this[int row, int column]
        {
            get
            {
                if (row >= this.RowCount || column >= this.ColumnCount)
                    throw new Exception("Invalid column or row specified");

                //return this.coreArray[column * this.rowCount + row];
                return this.coreArray[row][column];
            }

            set
            {
                if (row >= this.RowCount || column >= this.ColumnCount)
                    throw new Exception("Invalid column or row specified");

                this.coreArray[row][column] = value;
            }
        }

        /// <summary>
        /// Multiplies the specified matrices.
        /// </summary>
        /// <param name="m1">The m1.</param>
        /// <param name="m2">The m2.</param>
        /// <returns>m1 * m2</returns>
        public static Matrix Multiply(Matrix m1, Matrix m2)
        {
            if (m1.ColumnCount != m2.RowCount)
                throw new InvalidOperationException("No consistent dimensions");

            var res = new Matrix(m1.RowCount, m2.ColumnCount);

            var m1a = m1.coreArray;
            var m2a = m2.coreArray;
            var resa = res.coreArray;

            var v1 = m1.rowCount;
            var v2 = m2.columnCount;
            var v3 = m1.columnCount;

            for (int i = 0; i < v1; i++)
            {
                var resi = resa[i];
                var m1i = m1a[i];

                for (int j = 0; j < v2; j++)
                    for (int k = 0; k < v3; k++)
                    {
                        //res.coreArray[j * res.rowCount + i] +=
                        resi[j] +=
                            //m1.coreArray[k * m1.rowCount + i] *
                            m1i[k] *
                            //m2.coreArray[j * m2.rowCount + k];
                            m2a[k][j];
                    }

            }



            return res;
        }

        /// <summary>
        /// Throws an exception with specified <see cref="message"/> if <see cref="condition"/> is true.
        /// </summary>
        /// <param name="condition">the condition.</param>
        /// <param name="message">The message.</param>
        /// <exception cref="System.Exception">If conditions met!</exception>
        private static void ThrowIf(bool condition, string message)
        {
            if (condition)
                throw new Exception(message);
        }

        /// <summary>
        /// To2s the d double array.
        /// </summary>
        /// <param name="mtx">The MTX.</param>
        /// <returns></returns>
        public static double[,] To2DDoubleArray(Matrix mtx)
        {
            var buf = new double[mtx.RowCount, mtx.ColumnCount];

            for (int i = 0; i < mtx.RowCount; i++)
                for (int j = 0; j < mtx.ColumnCount; j++)
                    //buf[i, j] = mtx.coreArray[j * mtx.RowCount + i];
                    buf[i, j] = mtx.coreArray[i][j];

            return buf;
        }

        /// <summary>
        /// Creates a new Identity matrix.
        /// </summary>
        /// <param name="n">The matrix dimension.</param>
        /// <returns>Eye matrix of size <see cref="n"/></returns>
        public static Matrix Eye(int n)
        {
            var buf = new Matrix(n, n);

            for (int i = 0; i < n; i++)
                //buf.coreArray[i * n + i] = 1.0;
                buf.coreArray[i][i] = 1.0;

            return buf;
        }

        /// <summary>
        /// Creates <see cref="m"/> by <see cref="n"/> matrix filled with zeros.
        /// </summary>
        /// <param name="m">Number of rows.</param>
        /// <param name="n">Number of columns.</param>
        /// <returns><see cref="m"/>x<see cref="n"/> matrix filled with zeros.</returns>
        public static Matrix Zeros(int m, int n)
        {
            return new Matrix(m, n);
        }

        /// <summary>
        /// Creates n by n matrix filled with zeros.
        /// </summary>       
        /// <param name="n">Number of rows and columns, resp.</param>
        /// <returns>n x n matrix filled with zeros.</returns>
        public static Matrix Zeros(int n)
        {
            return new Matrix(n);
        }

        /// <summary>
        /// Creates row by n matrix filled with ones.
        /// </summary>
        /// <param name="m">Number of rows.</param>
        /// <param name="n">Number of columns.</param>
        /// <returns>m x n matrix filled with zeros.</returns>
        public static Matrix Ones(int m, int n)
        {
            var buf = new Matrix(m, n);

            //for (int i = 0; i < m*n; i++)
            //buf.coreArray[i] = 1.0;
            var arr = buf.coreArray;

            for (int i = 0; i < m; i++)
                for (int j = 0; j < arr[i].Length; j++)
                    arr[i][j] = 1.0;

            return buf;
        }

        /// <summary>
        /// Creates n by n matrix filled with ones.
        /// </summary>        
        /// <param name="n">Number of columns.</param>
        /// <returns>n x n matrix filled with ones.</returns>        
        public static Matrix Ones(int n)
        {
            var buf = new Matrix(n);

            //for (int i = 0; i < n * n; i++)
            //    buf.coreArray[i] = 1.0;

            var arr = buf.coreArray;

            for (int i = 0; i < n; i++)
                for (int j = 0; j < arr[i].Length; j++)
                    arr[i][j] = 1.0;

            return buf;
        }

        /// <summary>
        /// Computes product of main diagonal entries.
        /// </summary>
        /// <returns>Product of diagonal elements</returns>
        public double DiagProd()
        {
            var buf = 1.0;
            int dim = System.Math.Min(this.rowCount, this.columnCount);

            for (int i = 0; i < dim; i++)
            {
                buf *= this.coreArray[i][i];
            }

            return buf;
        }

        /// <summary>
        /// Implements the operator *.
        /// </summary>
        /// <param name="m1">The m1.</param>
        /// <param name="m2">The m2.</param>
        /// <returns>
        /// The result of the operator.
        /// </returns>
        public static Matrix operator *(Matrix m1, Matrix m2)
        {
            return Matrix.Multiply(m1, m2);
        }

        /// <summary>
        /// Implements the operator *.
        /// </summary>
        /// <param name="coeff">The coeff.</param>
        /// <param name="mat">The mat.</param>
        /// <returns>
        /// The result of the operator.
        /// </returns>
        public static Matrix operator *(double coeff, Matrix mat)
        {
            //var newMat = new double[mat.RowCount*mat.ColumnCount];


            //for (int i = 0; i < newMat.Length; i++)
            //    newMat[i] = coeff*mat.coreArray[i];

            var buf = mat.Clone();

            for (int i = 0; i < buf.rowCount; i++)
                for (int j = 0; j < buf.coreArray[i].Length; j++)
                    buf.coreArray[i][j] *= coeff;

            return buf;

        }

        /// <summary>
        /// Implements the operator *.
        /// </summary>
        /// <param name="mat">The mat.</param>
        /// <param name="coeff">The coeff.</param>
        /// <returns>
        /// The result of the operator.
        /// </returns>
        public static Matrix operator *(Matrix mat, double coeff)
        {
            //var newMat = new double[mat.RowCount*mat.ColumnCount];


            //for (int i = 0; i < newMat.Length; i++)
            //{
            //    newMat[i] = coeff*mat.coreArray[i];
            //}

            //var buf = new Matrix(mat.RowCount, mat.ColumnCount);
            //buf.coreArray = newMat;

            //return buf;

            var buf = mat.Clone();

            for (int i = 0; i < buf.rowCount; i++)
                for (int j = 0; j < buf.coreArray[i].Length; j++)
                    buf.coreArray[i][j] *= coeff;

            return buf;
        }

        /// <summary>
        /// Implements the operator -.
        /// </summary>
        /// <param name="mat">The mat.</param>
        /// <returns>
        /// The result of the operator.
        /// </returns>
        public static Matrix operator -(Matrix mat)
        {
            var buf = mat.Clone();// new Matrix(mat.RowCount, mat.ColumnCount);


            for (int i = 0; i < buf.rowCount; i++)
                for (int j = 0; j < buf.coreArray[i].Length; j++)
                    buf.coreArray[i][j] *= -1;

            return buf;
        }

        /// <summary>
        /// Implements the operator +.
        /// </summary>
        /// <param name="mat1">The mat1.</param>
        /// <param name="mat2">The mat2.</param>
        /// <returns>
        /// The result of the operator.
        /// </returns>
        public static Matrix operator +(Matrix mat1, Matrix mat2)
        {
            ThrowIf(mat1.RowCount != mat2.RowCount || mat1.ColumnCount != mat2.ColumnCount,
                "Inconsistent matrix sizes");

            var buf = new Matrix(mat1.RowCount, mat1.ColumnCount);

            var n = mat1.rowCount;
            var m = mat1.columnCount;

            //for (int i = 0; i < buf.coreArray.Length; i++)
            for (int i = 0; i < n; i++)
                for (int j = 0; j < m; j++)
                {
                    //buf.coreArray[i] = mat1.coreArray[i] + mat2.coreArray[i];
                    buf.coreArray[i][j] = mat1.coreArray[i][j] + mat2.coreArray[i][j];
                }

            return buf;
        }

        /// <summary>
        /// Implements the operator -.
        /// </summary>
        /// <param name="mat1">The mat1.</param>
        /// <param name="mat2">The mat2.</param>
        /// <returns>
        /// The result of the operator.
        /// </returns>
        public static Matrix operator -(Matrix mat1, Matrix mat2)
        {
            ThrowIf(mat1.RowCount != mat2.RowCount || mat1.ColumnCount != mat2.ColumnCount,
                "Inconsistent matrix sizes");

            var buf = new Matrix(mat1.RowCount, mat1.ColumnCount);

            var n = mat1.rowCount;
            var m = mat1.columnCount;

            //for (int i = 0; i < buf.coreArray.Length; i++)
            for (int i = 0; i < n; i++)
                for (int j = 0; j < m; j++)
                {
                    //buf.coreArray[i] = mat1.coreArray[i] + mat2.coreArray[i];
                    buf.coreArray[i][j] = mat1.coreArray[i][j] - mat2.coreArray[i][j];
                }

            return buf;
        }

        /// <summary>
        /// Implements the operator ==.
        /// </summary>
        /// <param name="left">The left.</param>
        /// <param name="right">The right.</param>
        /// <returns>
        /// The result of the operator.
        /// </returns>
        public static bool operator ==(Matrix left, Matrix right)
        {
            return Equals(left, right);
        }

        /// <summary>
        /// Implements the operator !=.
        /// </summary>
        /// <param name="left">The left.</param>
        /// <param name="right">The right.</param>
        /// <returns>
        /// The result of the operator.
        /// </returns>
        public static bool operator !=(Matrix left, Matrix right)
        {
            return !Equals(left, right);
        }

        /// <summary>
        /// Sets the specified row of matrix with defined values.
        /// </summary>
        /// <param name="i">The row number.</param>
        /// <param name="values">The values.</param>
        /// <exception cref="System.ArgumentOutOfRangeException">values</exception>
        public void SetRow(int i, params double[] values)
        {
            if (values.Length != this.ColumnCount)
                throw new ArgumentOutOfRangeException("values");

            Array.Copy(values, this.coreArray[i], values.Length);

            //for (int j = 0; j < this.ColumnCount; j++)
            {
                //this.coreArray[j * this.RowCount + i] = values[j];
                //this.coreArray[i][j] = values[j];
            }
        }

        /// <summary>
        /// Sets the specified column of matrix with defined values.
        /// </summary>
        /// <param name="j">The column number.</param>
        /// <param name="values">The values.</param>
        /// <exception cref="System.ArgumentOutOfRangeException">values</exception>
        public void SetColumn(int j, params double[] values)
        {
            if (values.Length != this.RowCount)
                throw new ArgumentOutOfRangeException("values");


            for (int i = 0; i < this.RowCount; i++)
            {
                //this.coreArray[j * this.RowCount + i] = values[i];
                this.coreArray[i][j] = values[i];
            }
        }

        /// <summary>
        /// Provides a shallow copy of this matrix.
        /// </summary>
        /// <returns>Clone of this matrix</returns>
        public Matrix Clone()
        {
            var buf = new Matrix(this.RowCount, this.ColumnCount);

            var source = this.coreArray;
            var dest = buf.coreArray;

            var len = source.Length;

            for (var x = 0; x < len; x++)
            {
                var src = source[x];
                var tgt = dest[x];
                var ilen = source[x].Length;

                Array.Copy(src, tgt, ilen);
            }

            return buf;

            //buf.coreArray = (double[])this.coreArray.Clone();
            //return buf;
        }

        /// <summary>
        /// Swaps each matrix entry A[i, j] with A[j, i].
        /// </summary>
        /// <returns>Transposed of this matrix.</returns>
        public Matrix Transpose()
        {
            var buf = new Matrix(this.ColumnCount, this.RowCount);

            var newMatrix = buf.coreArray;


            for (int i = 0; i < buf.RowCount; i++)
                for (int j = 0; j < buf.ColumnCount; j++)
                    //newMatrix[column*this.RowCount + row] = this.coreArray[row*this.RowCount + column];
                    buf[i, j] = this[j, i];

            // for (int row = 0; row < this.RowCount; row++)
            // for (int column = 0; column < this.ColumnCount; column++)
            //newMatrix[column*this.RowCount + row] = this.coreArray[row*this.RowCount + column];
            // buf[column, row] = this[row, column];

            buf.coreArray = newMatrix;
            return buf;
        }

        protected bool Equals(Matrix other)
        {
            if (other.RowCount != this.RowCount)
                return false;

            if (other.ColumnCount != this.ColumnCount)
                return false;

            var n = other.rowCount;
            var m = other.columnCount;

            /*
            for (int i = 0; i < other.coreArray.Length; i++)
            {
                if (!FuzzyEquals(this.coreArray[i], other.coreArray[i]))
                    return false;
            }
            */


            for (int i = 0; i < n; i++)
                for (int j = 0; j < m; j++)
                {
                    if (!FuzzyEquals(this.coreArray[i][j], other.coreArray[i][j]))
                        return false;
                }

            return true;
        }

        public override bool Equals(object obj)
        {
            if (ReferenceEquals(null, obj)) return false;
            if (ReferenceEquals(this, obj)) return true;
            if (obj.GetType() != this.GetType()) return false;
            return Equals((Matrix)obj);
        }

        /// <summary>
        /// Checks if number of rows equals number of columns.
        /// </summary>
        /// <returns>True iff matrix is n by n.</returns>
        public bool IsSquare()
        {
            return (this.columnCount == this.rowCount);
        }

        public bool IsSymmetric()
        {
            for (int i = 0; i < this.rowCount; i++)
                for (int j = i + 1; j < this.columnCount; j++)
                    if (!FuzzyEquals(this[i, j], this[j, i]))
                    {
                        return false;
                    }

            return true;
        }

        /// <summary>
        /// Checks if A[i, j] == 0 for i < j.
        /// </summary>
        /// <returns>True iff matrix is upper trapeze.</returns>
        public bool IsUpperTrapeze()
        {
            for (int j = 1; j <= columnCount; j++)
                for (int i = j + 1; i <= rowCount; i++)
                    //if (!FuzzyEquals(this.coreArray[j * this.RowCount + i], 0))
                    if (!FuzzyEquals(this.coreArray[i][j], 0))
                        return false;

            return true;
        }

        /// <summary>
        /// Checks if A[i, j] == 0 for i > j.
        /// </summary>
        /// <returns>True iff matrix is lower trapeze.</returns>
        public bool IsLowerTrapeze()
        {
            for (int i = 1; i <= rowCount; i++)
                for (int j = i + 1; j <= columnCount; j++)
                    //if (!FuzzyEquals(this.coreArray[j * this.RowCount + i], 0.0))
                    if (!FuzzyEquals(this.coreArray[i][j], 0.0))
                        return false;

            return true;
        }

        /// <summary>
        /// Checks if matrix is lower or upper trapeze.
        /// </summary>
        /// <returns>True iff matrix is trapeze.</returns>
        public bool IsTrapeze()
        {
            return (this.IsUpperTrapeze() || this.IsLowerTrapeze());
        }

        /// <summary>
        /// Checks if matrix is trapeze and square.
        /// </summary>
        /// <returns>True iff matrix is triangular.</returns>
        public bool IsTriangular()
        {
            return (this.IsLowerTriangular() || this.IsUpperTriangular());
        }

        /// <summary>
        /// Checks if matrix is square and upper trapeze.
        /// </summary>
        /// <returns>True iff matrix is upper triangular.</returns>
        public bool IsUpperTriangular()
        {
            return (this.IsSquare() && this.IsUpperTrapeze());
        }

        /// <summary>
        /// Checks if matrix is square and lower trapeze.
        /// </summary>
        /// <returns>True iff matrix is lower triangular.</returns>
        public bool IsLowerTriangular()
        {
            return (this.IsSquare() && this.IsLowerTrapeze());
        }

        /// <summary>
        /// Swaps rows at specified indices. The latter do not have to be ordered.
        /// When equal, nothing is done.
        /// </summary>
        /// <param name="i1">One-based index of first row.</param>
        /// <param name="i2">One-based index of second row.</param>        
        public void SwapRows(int i1, int i2)
        {
            if (i1 < 0 || i1 >= rowCount || i2 < 0 || i2 >= rowCount)
                throw new ArgumentException("Indices must be positive and <= number of rows.");

            if (i1 == i2)
                return;

            var buff = this.coreArray[i1];
            this.coreArray[i1] = this.coreArray[i2];
            this.coreArray[i2] = buff;

            /*
            for (int i = 0; i < columnCount; i++)
            {
                var tmp = this[i1, i];

                this[i1, i] = this[i2, i];
                this[i2, i] = tmp;
            }
            */
        }

        /// <summary>
        /// Swaps columns at specified indices. The latter do not have to be ordered.
        /// When equal, nothing is done.
        /// </summary>
        /// <param name="j1">One-based index of first col.</param>
        /// <param name="j2">One-based index of second col.</param>       
        public void SwapColumns(int j1, int j2)
        {
            if (j1 <= 0 || j1 > columnCount || j2 <= 0 || j2 > columnCount)
                throw new ArgumentException("Indices must be positive and <= number of cols.");

            if (j1 == j2)
                return;

            for (int i = 0; i < rowCount; i++)
            {
                var tmp = this.coreArray[i][j1];

                this.coreArray[i][j1] = this.coreArray[i][j2];
                this.coreArray[i][j2] = tmp;
            }

            //var j1Col = this.ExtractColumn(j1).coreArray;
            //var j2Col = this.ExtractColumn(j2).coreArray;

            //this.SetRow(j1, j2Col);
            //this.SetRow(j2, j1Col);
        }

        /// <summary>
        /// Retrieves row vector at specfifed index and deletes it from matrix.
        /// </summary>
        /// <param name="i">One-based index at which to extract.</param>
        /// <returns>Row vector.</returns>
        public Matrix ExtractRow(int i)
        {
            if (i >= this.RowCount || i < 0)
                throw new ArgumentOutOfRangeException("i");

            var mtx = new Matrix(1, this.ColumnCount);

            Array.Copy(this.coreArray[0], mtx.coreArray[0], this.coreArray[0].Length);

            //for (int j = 0; j < this.ColumnCount; j++)
            //{
            //    mtx.coreArray[j] = this.coreArray[j * this.RowCount + i];
            //}

            return mtx;
        }

        /// <summary>
        /// Retrieves column vector at specfifed index and deletes it from matrix.
        /// </summary>
        /// <param name="j">One-based index at which to extract.</param>
        /// <returns>Row vector.</returns>
        public Matrix ExtractColumn(int j)
        {
            if (j >= this.ColumnCount || j < 0)
                throw new ArgumentOutOfRangeException("j");

            var mtx = new Matrix(this.RowCount, 1);


            for (int i = 0; i < this.RowCount; i++)
            {
                mtx.coreArray[i][0] = this.coreArray[i][j];
            }

            return mtx;
        }

        /// <summary>
        /// Gets the determinant of matrix
        /// </summary>
        /// <returns></returns>
        public double Determinant()
        {
            // Seems working good!

            if (!IsSquare())
                throw new InvalidOperationException();

            var clone = this.Clone();

            var n = this.rowCount;

            var sign = 1.0;

            var epsi1on = 1e-10 * MinAbs(clone.coreArray);

            if (epsi1on == 0)
                epsi1on = 1e-9;

            // this[row,column] = this.coreArray[column*this.rowCount + row]

            for (var i = 0; i < n - 1; i++)
            {
                if (System.Math.Abs(clone[i, i]) < epsi1on)
                {
                    var firstNonZero = -1;

                    for (var k = i + 1; k < n; k++)
                        if (System.Math.Abs(clone[k, i]) > epsi1on)
                            firstNonZero = k;

                    if (firstNonZero == -1)
                        throw new OperationCanceledException();
                    else
                    {
                        clone.SwapRows(firstNonZero, i);
                        sign = -sign;
                    }
                }



                for (var j = i + 1; j < n; j++)
                {
                    //var alfa = (clone.coreArray[j * n + i] / clone.coreArray[i * n + i]);
                    var alfa = (clone.coreArray[i][j] / clone.coreArray[i][i]);

                    for (var k = i; k < n; k++)
                    {
                        //clone.coreArray[j * n + k] -= alfa * clone.coreArray[i * n + k];
                        clone.coreArray[k][j] -= alfa * clone.coreArray[k][i];
                    }
                }
            }

            var buf = sign;

            var arr = new double[n];

            for (var i = 0; i < n; i++)
                //arr[i] = clone.coreArray[i * n + i];
                arr[i] = clone.coreArray[i][i];

            Array.Sort(arr);

            for (var i = 0; i < n; i++)
                buf = buf * arr[n - i - 1];

            return buf;
        }

        private static double MinAbs(double[] arr)
        {
            var min = double.MaxValue;

            foreach (var d in arr)
            {
                if (Math.Abs(d) < min)
                    min = Math.Abs(d);
            }

            return min;
        }

        private static double MinAbs(double[][] arr)
        {
            var min = double.MaxValue;

            foreach (var a in arr)
            {
                foreach (var d in a)
                    if (Math.Abs(d) < min)
                        min = Math.Abs(d);
            }

            return min;
        }

        /// <summary>
        /// Calculates the inverse of matrix
        /// </summary>
        /// <returns>inverse of this matrix</returns>
        public Matrix Inverse()
        {
            if (!IsSquare())
                throw new InvalidOperationException();

            //seems working good!
            var n = this.rowCount;
            var clone = this.Clone();
            var eye = Eye(n);

            var epsi1on = 1e-10 * MinAbs(clone.coreArray);

            if (epsi1on == 0)
                epsi1on = 1e-9;

            /**/

            var perm = new List<int>();

            var clonea = clone.coreArray;
            var eyea = eye.coreArray;

            for (var j = 0; j < n - 1; j++)
            {
                for (var i = j + 1; i < n; i++)
                {
                    //if (System.Math.Abs(clonea[j + j * n]) < epsi1on)
                    if (System.Math.Abs(clonea[j][j]) < epsi1on)
                    {
                        var firstNonZero = -1;

                        for (var k = j + 1; k < n; k++)
                            //if (System.Math.Abs(clonea[k + j * n]) > epsi1on)
                            if (System.Math.Abs(clonea[k][j]) > epsi1on)
                                firstNonZero = k;

                        if (firstNonZero == -1)
                            throw new OperationCanceledException();
                        else
                        {
                            clone.SwapRows(firstNonZero, j);
                            eye.SwapRows(firstNonZero, j);

                            perm.Add(j);
                            perm.Add(firstNonZero);
                        }
                    }

                    //var alfa = clonea[i + j * n] / clonea[j + j * n];
                    var alfa = clonea[i][j] / clonea[j][j];

                    for (var k = 0; k < n; k++)
                    {
                        //clonea[i + k * n] -= alfa * clonea[j + k * n];
                        clonea[i][k] -= alfa * clonea[j][k];

                        //eyea[i + k * n] -= alfa * eyea[j + k * n];
                        eyea[i][k] -= alfa * eyea[j][k];
                    }
                }
            }

            /**/

            for (var j = n - 1; j > 0; j--)
            {
                for (var i = j - 1; i >= 0; i--)
                {
                    //if (System.Math.Abs(clonea[j + j * n]) < epsi1on)
                    if (System.Math.Abs(clonea[j][j]) < epsi1on)
                    {
                        var firstNonZero = -1;

                        for (var k = j - 1; k >= 0; k--)
                            //if (System.Math.Abs(clonea[k + j * n]) > epsi1on)
                            if (System.Math.Abs(clonea[k][j]) > epsi1on)
                                firstNonZero = k;

                        if (firstNonZero == -1)
                            throw new OperationCanceledException();
                        else
                        {
                            clone.SwapRows(firstNonZero, j);
                            eye.SwapRows(firstNonZero, j);

                            perm.Add(j);
                            perm.Add(firstNonZero);
                        }
                    }

                    //var alfa = clonea[i + j * n] / clonea[j + j * n];
                    var alfa = clonea[i][j] / clonea[j][j];

                    for (var k = n - 1; k >= 0; k--)
                    {
                        //clonea[i + k * n] -= alfa * clonea[j + k * n];
                        clonea[i][k] -= alfa * clonea[j][k];

                        //eyea[i + k * n] -= alfa * eyea[j + k * n];
                        eyea[i][k] -= alfa * eyea[j][k];
                    }
                }
            }

            /**/

            for (var i = 0; i < n; i++)
            {
                //var alfa = 1 / clonea[i + i * n];
                var alfa = 1 / clonea[i][i];

                for (var j = 0; j < n; j++)
                {
                    //clonea[i + j * n] *= alfa;
                    clonea[i][j] *= alfa;

                    //eyea[i + j * n] *= alfa;
                    eyea[i][j] *= alfa;
                }
            }

            /**/

            return eye;
        }


        /// <summary>
        /// determines the equation of parameters with fuzzy approach.
        /// </summary>
        /// <remarks>
        /// Because of rounding errors in computer, this method will be used for checking equality of numbers</remarks>
        /// <param name="d1">The d1.</param>
        /// <param name="d2">The d2.</param>
        /// <returns>true, if d1 equals with fuzzy approach to d2, otherwise false</returns>
        private bool FuzzyEquals(double d1, double d2)
        {
            return Math.Abs(d1 - d2) < Epsilon;
        }

        /// <summary>
        /// The Threshold of equality of two double precision members.
        /// </summary>
        public double Epsilon = 1e-12;

        /// <summary>
        /// multiply two matrixes in a cache friendly way. for large matrixes it does have benefits
        /// </summary>
        /// <param name="mat1"></param>
        /// <param name="mat2"></param>
        /// <returns></returns>
        public static Matrix FastMultiply(Matrix mat1, Matrix mat2,bool mat2IsSymmetric = false)
        {
            //ThrowIf(mat1.RowCount != mat2.RowCount || mat1.ColumnCount != mat2.ColumnCount,
            //    "Inconsistent matrix sizes");

            if (mat1.ColumnCount != mat2.RowCount)
                throw new InvalidOperationException("No consistent dimensions");

            var buf = new Matrix(mat1.RowCount, mat2.ColumnCount);

            Matrix m2c = mat2IsSymmetric?mat2:mat2.Transpose();


            for (int i = 0; i < mat1.rowCount; i++)
                for (int j = 0; j < mat2.columnCount; j++)
                {
                    buf.coreArray[i][j] = DotProd(mat1.coreArray[i], m2c.coreArray[j]); ;
                }

            return buf;
        }

        private static double DotProd(double[] a1, double[] a2)
        {
            /**/
            var tt = 0.0;

            for (var i = 0; i < a1.Length; i++)
                tt += a1[i] * a2[i];

            return tt;
            /* */

            var n = a1.Length - a1.Length%3;

            var v1 = 0.0;
            var v2 = 0.0;
            var v3 = 0.0;

            for (var i = 0; i < n; i+=3)
            {
                v1 += a1[i] * a2[i];
                v2 += a1[i + 1] * a2[i + 1];
                v3 += a1[i + 2] * a2[i + 2];
            }

            for (var i = n; i < a1.Length; i++)
            {
                v1 += a1[i] * a2[i];
            }

            return v1 + v2 + v3;
        }

        public static Matrix FromColMajorMatrix(int rowCount,int colCount, double[] coreArr)
        {
            var buf = new Matrix(rowCount, colCount);

            for (var i = 0; i < rowCount; i++)
                for (var j = 0; j < colCount; j++)
                {
                    buf[i, j] = coreArr[j * rowCount + i];
                }

            return buf;
        }


        /// <summary>
        /// Does the transpose operation inplace
        /// </summary>
        public void InPlaceTranspose()
        {
            if (!this.IsSquare())
                throw new Exception();

            var n = this.rowCount;

            for (var i = 0; i < n; i++)
                for (var j = 0; j < i; j++)
                {
                    if (i == j)
                        continue;

                    var vij = this[i, j];
                    var vji = this[j, i];

                    this[i, j] = vji;
                    this[j, i] = vij;
                }
        }

        /// <summary>
        /// the ABS Maximum of members
        /// </summary>
        /// <returns></returns>
        public double MaxAbsMember()
        {
            var max = this.coreArray.Max(i => i.Max(j => Math.Abs(j)));

            return max;
        }

        public static Matrix FromJagged(double[][] core)
        {
            var buf = new Matrix();
            buf.rowCount = core.Length;
            buf.columnCount = core[0].Length;
            buf.coreArray = core;


            return buf;
        }

        public static Matrix From2D(double[,] core)
        {
            var rowCount = core.GetLength(0);
            var columnCount = core.GetLength(1);

            var buf = new Matrix(rowCount, columnCount);

            for (var i = 0; i < rowCount; i++)
                for (var j = 0; j < columnCount; j++)
                    buf[i, j] = core[i, j];

            return buf;
        }


        /// <summary>
        /// does the r' . X . r multiply operation with fast methods, where X is symmetric
        /// </summary>
        /// <param name="r">the r matrix</param>
        /// <param name="x">the x matrix, symmetric</param>
        /// /// <param name="x">the result</param>
        /// <returns></returns>
        public static void RtXR_Operation(Matrix r, Matrix x, Matrix result)
        {
            var n = r.rowCount;

            if (!r.IsSquare() || !x.IsSquare() || !result.IsSquare() || n != x.rowCount || n != result.rowCount)
                throw new Exception();

            var li = new double[n];


            r.InPlaceTranspose();

            for (var i = 0; i < n; i++)
            {
                #region computation of Li
                for (var j = 0; j < n; j++)
                    li[j] = DotProd(r.coreArray[i], x.coreArray[j]);
                #endregion

                for (var q = 0; q < n; q++)
                {
                    result[i, q] = DotProd(li, r.coreArray[q]);
                }
            }

            r.InPlaceTranspose();
        }
    }
}