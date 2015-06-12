# LightNetMatrix
A simple, lightweight and resonably fast dense double matrix for .NET. All codes are packed inside a single class in single .cs file which makes it simply possible to include it into your C# code.

## Usages
its very simple to use:

```c#
Matrix m = RandomMatrix(5, 6);   //creates a new 5x6 matrix ( 5 row and 6 column )
Matrix inv = m.Inverse();        // computes inverse of m into inv
double det = m.Determinant();    // computes determinant of m into det

        
```

#Features
- Reasonable Performance
- Determinant of square Matrix
- Both Inverse and Transpose of square matrix.
- Cholesky decomposition of matrix.
- Multyply, plus and minus operations overloaded.
- Compatible with .NET Portable framework.
- CLS Compliant.
- Serializable (both xml and binary)

## Performance
Because of using one dimensional array for storing data, performance is far higher than 2d array in .NET and is somehow around jagged array. Anyways time consuming parts (like inverse and determinant) are not highly optimized for performance. On a i7 pc inverse of matrix of 300x300 took ~160ms and determinant took ~50ms in release mode and no debugger attached.

## Compatibility
it is compatible with portable .Net framwork (using in Android and IOS) and .NET 2.0+ without any changes.

#Notes
If you are dealing with very large matrixes, you maybe have to use sparse matrix instead of dense matrix and also if you are trying to solve a large equation system with inversing the matrix, you maybe should use other methods like cholesky decomposition or other ones (see [this guidline!](http://www.johndcook.com/blog/2010/01/19/dont-invert-that-matrix/)). For sparse matrix library in .NET Csparse.net is a wonderfull one with cholesky decomposition!
