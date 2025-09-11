#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <random>
#include <ctime>
#include <algorithm>
#include <filesystem>
#include <sstream>
#include <regex>
#include <tuple>
#include <utility>
#include <string>
#include <cstdlib>
#include "PolyUtils.h"
#include "PolygonClipping.h"
#include "potential_cal.h"

using namespace std;

double LJPotentialBetweenTwoAtoms(double Sigma, double Epsilon, double distance) {
/* Calculates the Lennard-Jones potential energy between two atoms based on their distance.
Inputs:
  Sigma   - diameter of an atom
  Epsilon - depth of potential well
  distance - distance between atoms
Output: Potential energy as a double
*/
    double SigmaDividedByDistance = Sigma / distance;
    double Energy = 4 * Epsilon * (pow(SigmaDividedByDistance, 12) - pow(SigmaDividedByDistance, 6));
    return Energy;
}

double LJPotentialBetweenTwoParticles(double Sigma, double Epsilon, double LJCutOff,
                                   const double atomList1[][3], int size1,
                                   const double atomList2[][3], int size2, double aspectratio,
                                   double BoxLength, string shape) {
/* Calculates the total Lennard-Jones potential between two faceted particles made of multiple atoms.
Inputs:
  Sigma       - diameter of an atom
  Epsilon     - depth of potential well
  LJCutOff    - cutoff distance for interaction
  atomList1   - coordinates of atoms in particle 1
  size1       - number of atoms per side for particle 1
  atomList2   - coordinates of atoms in particle 2
  size2       - number of atoms per side for particle 2
  aspectratio - height / length of particle shape
  BoxLength   - length of the simulation box (for periodic boundary)
  shape       - shape of the particles ("cuboid" or "triangularprism")
Output: Total Lennard-Jones potential energy between the two particles as a double
*/
    double CubeCubeEnergy = 0;
    int row1, row2;
    if (shape == "cuboid") {
      row1 = size1 * size1 * size1 * aspectratio; // number of atoms/CG-beads in particle 1
      row2 = size2 * size2 * size2 * aspectratio; // number of atoms/CG-beads in particle 2
    }

    if (shape == "triangularprism") {
      row1 = size1 * size1 * (size1+1) * aspectratio / 2.0;
      row2 = size2 * size2 * (size2+1) * aspectratio / 2.0;
    }

    for (int i = 0; i < row1; ++i) {
        for (int j = 0; j < row2; ++j) {
            double Atom1X = atomList1[i][0];
            double Atom1Y = atomList1[i][1];
            double Atom1Z = atomList1[i][2];
            double Atom2X = atomList2[j][0];
            double Atom2Y = atomList2[j][1];
            double Atom2Z = atomList2[j][2];

            double dx = Atom1X - Atom2X;
            dx = dx - BoxLength * round(dx / BoxLength); // periodic boundary condition
            double dy = Atom1Y - Atom2Y;
            dy = dy - BoxLength * round(dy / BoxLength);
            double dz = Atom1Z - Atom2Z;
            dz = dz - BoxLength * round(dz / BoxLength);
            double distance = sqrt(dx * dx + dy * dy + dz * dz); // distance between two atoms

            if (distance <= LJCutOff) {
              double energy = LJPotentialBetweenTwoAtoms(Sigma, Epsilon, distance);
              CubeCubeEnergy += energy;
            }
        }
    }

    return CubeCubeEnergy;
}

typedef double Vector3[3];
typedef double Matrix3x3[3][3];

void MultiplyMatrixVector(const Matrix3x3& matrix, const Vector3 input, Vector3 result) {
/* Multiplies a 3x3 matrix by a 3D vector and stores the result.
Inputs:
  matrix - 3x3 matrix
  input  - 3D input vector
Output: result - 3D output vector (modified inside the function)
*/
    for (int i = 0; i < 3; ++i) {
        result[i] = 0;
        for (int j = 0; j < 3; ++j) {
            result[i] += input[j] * matrix[j][i];
        }
    }
}

void vectorSubtraction(const double a[3], const double b[3], double result[3]) {
/* Subtracts vector b from vector a and stores the result.
Inputs:
  a      - first 3D vector
  b      - second 3D vector to subtract from a
Output: result - resulting 3D vector after subtraction
*/
    result[0] = a[0] - b[0];
    result[1] = a[1] - b[1];
    result[2] = a[2] - b[2];
}

void normalizeVector(double vec[3]) {
/* Normalizes a 3D vector to have a magnitude of 1.
Inputs:
  vec - 3D vector to be normalized (modified in place)
Output:
  vec - normalized 3D vector
*/
    double norm = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
    vec[0] /= norm;
    vec[1] /= norm;
    vec[2] /= norm;
}

void vectorAddition(const double a[3], const double b[3], double result[3]) {
/* Adds two 3D vectors and stores the result.
Inputs:
  a      - first 3D vector
  b      - second 3D vector
Output: result - resulting 3D vector after addition
*/
    result[0] = a[0] + b[0];
    result[1] = a[1] + b[1];
    result[2] = a[2] + b[2];
}

void vectorScalarProduct(const double vec[3], double scalar, double result[3]) {
/* Multiplies a 3D vector by a scalar and stores the result.
Inputs:
  vec    - input 3D vector
  scalar - scalar value to multiply with the vector
Output: result - resulting 3D vector after multiplication
*/
    result[0] = scalar * vec[0];
    result[1] = scalar * vec[1];
    result[2] = scalar * vec[2];
}

void crossProduct(const double u[3], const double v[3], double result[3]) {
/* Computes the cross product of two 3D vectors u and v, stores in result.
Inputs:
  u      - first 3D vector
  v      - second 3D vector
Output: result - resulting 3D vector perpendicular to u and v
*/
    result[0] = u[1] * v[2] - u[2] * v[1];
    result[1] = u[2] * v[0] - u[0] * v[2];
    result[2] = u[0] * v[1] - u[1] * v[0];
}

double dotProduct(const double u[3], const double v[3]) {
/* Computes the dot product of two 3D vectors u and v.
Inputs:
  u - first 3D vector
  v - second 3D vector
Output:
  Returns the scalar dot product as a double
*/

    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
}

void computePlaneEquation(const double vertices[][3], double &A, double &B, double &C, double &D) {
/* Computes the coefficients A, B, C, D of the plane equation Ax + By + Cz + D = 0
   defined by vertices of a polygon.
Inputs:
  vertices - array of four 3D points defining the polygon
Outputs:
  A, B, C, D - coefficients of the plane equation (passed by reference)
*/
    double v1[3] = { vertices[1][0] - vertices[0][0], vertices[1][1] - vertices[0][1], vertices[1][2] - vertices[0][2] };
    double v2[3] = { vertices[2][0] - vertices[0][0], vertices[2][1] - vertices[0][1], vertices[2][2] - vertices[0][2] };

    double normal[3];
    crossProduct(v1, v2, normal);
    normalizeVector(normal);

    A = normal[0];
    B = normal[1];
    C = normal[2];

    // Ensure A is negative
    if (A > 0) {
        A = -A;
        B = -B;
        C = -C;
    }

    D = - (A * vertices[0][0] + B * vertices[0][1] + C * vertices[0][2]);
}

bool compareIndices(const int& i, const int& j, int column, const double arr[4][3]) {
/* Compares two indices based on the values in a specified column of a 2D array.
Inputs:
  i      - first index
  j      - second index
  column - column number to compare
  arr    - 2D array of doubles
Output:
  Returns true if arr[i][column] < arr[j][column], else false
*/
    return arr[i][column] < arr[j][column];
}

vector<int> sortIndicesByDesiredColumn(const double arr[4][3], int column) {
/* Sorts the indices of a 2D array based on the values in the specified column.
Inputs:
  arr    - 2D array of doubles
  column - column number to sort by
Output:
  Returns a vector of indices sorted by the specified column
*/
    // Create an array of indices
    vector<int> indices(4);
    iota(indices.begin(), indices.end(), 0);

    // Sort indices based on the second column
    sort(indices.begin(), indices.end(),
              [&arr, column](const int& i, const int& j) {
                  return compareIndices(i, j, column, arr);
              });

    return indices;
}


void ParticleRotation(const Vector3& CubeCentroid, const Vector3 VectorX, const Vector3 VectorY, const Vector3 VectorZ, const Vector3& RotationCenter,
                  const Vector3& RotationAngles, Vector3& CubeCentroidResult, Vector3 CubeVectorResult[3]) {
/* Rotates particles around the rotation center based on alpha, beta, gamma angles using rotation matrices.
Inputs:
  CubeCentroid      - original centroid of the cube (before rotation)
  VectorX, VectorY, VectorZ - axis vectors defining the orientation of the cube
  RotationCenter    - center of rotation
  RotationAngles    - angles [alpha, beta, gamma] in degrees to rotate about the X, Y, and Z axes respectively
Outputs:
  CubeCentroidResult - rotated centroid position
  CubeVectorResult   - rotated vectors for the cube's axes (X, Y, Z) after rotation
*/
    Vector3 RotVectorX, RotVectorY, RotVectorZ;

    double Alpha = RotationAngles[0] * M_PI / 180.0;
    double Beta = RotationAngles[1] * M_PI / 180.0;
    double Gamma = RotationAngles[2] * M_PI / 180.0;

    Matrix3x3 RotationX = {{1, 0, 0}, {0, cos(Alpha), -sin(Alpha)}, {0, sin(Alpha), cos(Alpha)}}; // Rx rotation matrix
    Matrix3x3 RotationY = {{cos(Beta), 0, sin(Beta)}, {0, 1, 0}, {-sin(Beta), 0, cos(Beta)}}; // Ry rotation matrix
    Matrix3x3 RotationZ = {{cos(Gamma), -sin(Gamma), 0}, {sin(Gamma), cos(Gamma), 0}, {0, 0, 1}}; // Rz rotation matrix

    Vector3 RotationCenter_To_CubeCentroid;
    for (int i = 0; i < 3; ++i) {
        RotationCenter_To_CubeCentroid[i] = CubeCentroid[i] - RotationCenter[i]; // Find the displacement vector from rotation center to particle center
    }

    Vector3 tempresult1;
    Vector3 tempresult2;
    MultiplyMatrixVector(RotationX, VectorX, tempresult1);
    MultiplyMatrixVector(RotationY, tempresult1, tempresult2);
    MultiplyMatrixVector(RotationZ, tempresult2, RotVectorX); // RotVectorX is the rotated VectorX = RzRyRxVectorX

    MultiplyMatrixVector(RotationX, VectorY, tempresult1);
    MultiplyMatrixVector(RotationY, tempresult1, tempresult2);
    MultiplyMatrixVector(RotationZ, tempresult2, RotVectorY); // RotVectorY is the rotated VectorY = RzRyRxVectorY

    MultiplyMatrixVector(RotationX, VectorZ, tempresult1);
    MultiplyMatrixVector(RotationY, tempresult1, tempresult2);
    MultiplyMatrixVector(RotationZ, tempresult2, RotVectorZ); // RotVectorZ is the rotated VectorZ = RzRyRxVectorZ

    MultiplyMatrixVector(RotationX, RotationCenter_To_CubeCentroid, tempresult1);
    MultiplyMatrixVector(RotationY, tempresult1, tempresult2);
    MultiplyMatrixVector(RotationZ, tempresult2, RotationCenter_To_CubeCentroid); // New rotated particle center

    double lengthX = sqrt(RotVectorX[0] * RotVectorX[0] + RotVectorX[1] * RotVectorX[1] + RotVectorX[2] * RotVectorX[2]);
    double lengthY = sqrt(RotVectorY[0] * RotVectorY[0] + RotVectorY[1] * RotVectorY[1] + RotVectorY[2] * RotVectorY[2]);
    double lengthZ = sqrt(RotVectorZ[0] * RotVectorZ[0] + RotVectorZ[1] * RotVectorZ[1] + RotVectorZ[2] * RotVectorZ[2]);

    RotVectorX[0] /= lengthX; // normalization of vectors
    RotVectorX[1] /= lengthX;
    RotVectorX[2] /= lengthX;

    RotVectorY[0] /= lengthY;
    RotVectorY[1] /= lengthY;
    RotVectorY[2] /= lengthY;

    RotVectorZ[0] /= lengthZ;
    RotVectorZ[1] /= lengthZ;
    RotVectorZ[2] /= lengthZ;

    for (int i = 0; i < 3; ++i) {
        CubeCentroidResult[i] = RotationCenter_To_CubeCentroid[i] + RotationCenter[i];
        CubeVectorResult[0][i] = RotVectorX[i];
        CubeVectorResult[1][i] = RotVectorY[i];
        CubeVectorResult[2][i] = RotVectorZ[i];
    }
}

void matrixMultiply(const double A[3][3], const double B[3][3], double result[3][3]) {
/* Multiplies two 3x3 matrices A and B, storing the result in the result matrix.
Inputs:
  A - first 3x3 matrix
  B - second 3x3 matrix
Output:
  result - 3x3 matrix containing the product of A and B
*/
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result[i][j] = 0;
            for (int k = 0; k < 3; ++k) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void findVerticesOfCube(double ParticleSideLength, double Height, double COMX, double COMY, double COMZ, const Vector3 VectorX, const Vector3 VectorY, const Vector3 VectorZ, double Vertex[][3]){
/* Returns positions of vertices based on the center of mass and orientation of the cube.
Inputs:
  ParticleSideLength - length of the side of the cube
  Height             - height of the cube (aspect ratio consideration)
  COMX, COMY, COMZ   - center of mass coordinates
  VectorX, VectorY, VectorZ - orientation vectors for cube
Output:
  Vertex - 2D array containing the positions of the 8 vertices of the cube
*/
    double aspectratio = Height/ParticleSideLength;
    double matrix[8][3] = {
        {-1, +1, +aspectratio},
        {-1, -1, +aspectratio},
        {-1, -1, -aspectratio},
        {-1, +1, -aspectratio},
        {+1, +1, +aspectratio},
        {+1, -1, +aspectratio},
        {+1, -1, -aspectratio},
        {+1, +1, -aspectratio}
    };

    double x,y,z;
    for (int i = 0; i < 8; ++i) {
        x = matrix[i][0];
        y = matrix[i][1];
        z = matrix[i][2];
        Vertex[i][0] = COMX + (ParticleSideLength / 2.0) * (x * VectorX[0] + y * VectorY[0] + z * VectorZ[0]);
        Vertex[i][1] = COMY + (ParticleSideLength / 2.0) * (x * VectorX[1] + y * VectorY[1] + z * VectorZ[1]);
        Vertex[i][2] = COMZ + (ParticleSideLength / 2.0) * (x * VectorX[2] + y * VectorY[2] + z * VectorZ[2]);
    }
}

void findVerticesOfTriangularPrism(double CubeSideLength, double Height, double AASigma, double COMX, double COMY, double COMZ, const Vector3 VectorX, const Vector3 VectorY, const Vector3 VectorZ, double Vertex[][3]){
/* Returns positions of vertices based on the center of mass and orientation of a triangular prism.
Inputs:
  CubeSideLength - length of the side of the triangular prism base
  Height          - height of the triangular prism
  AASigma         - diameter of atom
  COMX, COMY, COMZ - center of mass coordinates
  VectorX, VectorY, VectorZ - orientation vectors for triangular prism
Output:
  Vertex - 2D array containing the positions of the 6 vertices of the triangular prism
*/

    double aspectratio = Height/CubeSideLength;
    double side = CubeSideLength/AASigma;
    double xright = (side-1)*sqrt(3)/(3.0*side)+1/side;
    double xleft = -(side-1)*sqrt(3)*2.0/(3.0*side)-1/side;
    double matrix[6][3] = {
        {xright, -1, +aspectratio},
        {xright, +1, +aspectratio},
        {xright, +1, -aspectratio},
        {xright, -1, -aspectratio},
        {xleft, 0, +aspectratio},
        {xleft, 0, -aspectratio}
    };

    double x,y,z;
    for (int i = 0; i < 6; ++i) {
        x = matrix[i][0];
        y = matrix[i][1];
        z = matrix[i][2];
        Vertex[i][0] = COMX + (CubeSideLength / 2.0) * (x * VectorX[0] + y * VectorY[0] + z * VectorZ[0]);
        Vertex[i][1] = COMY + (CubeSideLength / 2.0) * (x * VectorX[1] + y * VectorY[1] + z * VectorZ[1]);
        Vertex[i][2] = COMZ + (CubeSideLength / 2.0) * (x * VectorX[2] + y * VectorY[2] + z * VectorZ[2]);
    }
}

void Reorientation_cuboid(int i, int j, double Cube2SideLength, double Cube2Height, double BoxLength,
                   Vector3 Particle1Centroid, Vector3 Particle2Centroid,
                   Vector3 Particle1VectorsX, Vector3 Particle1VectorsY, Vector3 Particle1VectorsZ,
                   Vector3 Particle2VectorsX, Vector3 Particle2VectorsY, Vector3 Particle2VectorsZ, double Vertex[][3],
                   double y_coord11[], double z_coord11[], double &ds, double lc) {
/* Reorients two cubes to determine the reference particle and align them such that the reference particle is aligned along coordinate axes.
Inputs:
  i, j                   - indices of the two particles (cube 1 and cube 2)
  Cube2SideLength        - side length of the second cube
  Cube2Height            - height of the second cube
  BoxLength              - the length of the periodic box
  Particle1Centroid      - centroid of the first particle
  Particle2Centroid      - centroid of the second particle
  Particle1VectorsX, Particle1VectorsY, Particle1VectorsZ - orientation vectors for the first cube
  Particle2VectorsX, Particle2VectorsY, Particle2VectorsZ - orientation vectors for the second cube
Output:
  Vertex                 - positions of the vertices of the second cube after reorientation
  y_coord11, z_coord11   - y-z coordinates of interacting face of reference particle
  ds                     - x-distance of center of interacting face of reference particle from its center
  lc                     - half of side length of reference particle in the units of AASigma
*/


    // Create temporary variables to store intermediate values
    double Size1Temp, Size2Temp;
    Vector3 Particle1CentroidTemp, Particle2CentroidTemp;
    Vector3 Particle1VectorsXTemp, Particle1VectorsYTemp, Particle1VectorsZTemp;
    Vector3 Particle2VectorsXTemp, Particle2VectorsYTemp, Particle2VectorsZTemp;

    Vector3 P1P2;
    vectorSubtraction(Particle2Centroid, Particle1Centroid, P1P2); // vector between centers of particle 1 and 2
    for (int k = 0; k < 3; ++k) {
       P1P2[k] = P1P2[k]-BoxLength*round(P1P2[k]/BoxLength); // periodic boundary condition
    }

    double P1P2ScaledX = P1P2[0]/Cube2SideLength;
    double P1P2ScaledY = P1P2[1]/Cube2SideLength;
    double P1P2ScaledZ = P1P2[2]/Cube2Height;
    double P1P2NormScaled = sqrt(P1P2ScaledX * P1P2ScaledX + P1P2ScaledY * P1P2ScaledY + P1P2ScaledZ * P1P2ScaledZ);
    //cout << "P1P2 " << P1P2[0] << " " << P1P2[1] << " " << P1P2[2] << endl;

    //Find angles between P1P2 and coordinate axes vectors
    double AngleBetweenP1P2AndX = std::acos(std::max(-1.0, std::min(1.0, P1P2ScaledX / P1P2NormScaled)));
    double AngleBetweenP1P2AndXNegative = std::acos(std::max(-1.0, std::min(1.0, -P1P2ScaledX / P1P2NormScaled)));
    double AngleBetweenP1P2AndY = std::acos(std::max(-1.0, std::min(1.0, P1P2ScaledY / P1P2NormScaled)));
    double AngleBetweenP1P2AndYNegative = std::acos(std::max(-1.0, std::min(1.0, -P1P2ScaledY / P1P2NormScaled)));
    double AngleBetweenP1P2AndZ = std::acos(std::max(-1.0, std::min(1.0, P1P2ScaledZ / P1P2NormScaled)));
    double AngleBetweenP1P2AndZNegative = std::acos(std::max(-1.0, std::min(1.0, -P1P2ScaledZ / P1P2NormScaled)));

    // Find out which axis (and direction) has the smallest angle with the P1P2 vector
    double MinimumAngle = min({AngleBetweenP1P2AndX, AngleBetweenP1P2AndXNegative, AngleBetweenP1P2AndY,
                                    AngleBetweenP1P2AndYNegative, AngleBetweenP1P2AndZ, AngleBetweenP1P2AndZNegative});

    //cout << "Minimum angle " << MinimumAngle << endl;
    //cout << AngleBetweenP1P2AndX << " " << AngleBetweenP1P2AndXNegative << " " << AngleBetweenP1P2AndY << " " << AngleBetweenP1P2AndYNegative << " " << AngleBetweenP1P2AndZ << " " << AngleBetweenP1P2AndZNegative << endl;

    double temp1, temp2;
    bool swapx = false;
    bool swapy = false;
    bool swapz = false;

    if (AngleBetweenP1P2AndXNegative == MinimumAngle) { // Particle 1 and 2 assignment is swapped
       swapx = true;
       for (int k = 0; k < 3; ++k) {
           temp1 = -P1P2[k];
           P1P2[k] = temp1;
           temp1 = Particle1VectorsX[k];
           temp2 = Particle2VectorsX[k];
           Particle1VectorsX[k] = temp2;
           Particle2VectorsX[k] = temp1;
           temp1 = Particle1VectorsY[k];
           temp2 = Particle2VectorsY[k];
           Particle1VectorsY[k] = temp2;
           Particle2VectorsY[k] = temp1;
           temp1 = Particle1VectorsZ[k];
           temp2 = Particle2VectorsZ[k];
           Particle1VectorsZ[k] = temp2;
           Particle2VectorsZ[k] = temp1;
        }
    } else if (AngleBetweenP1P2AndYNegative == MinimumAngle) { // Particle 1 and 2 assignment is swapped
       swapy = true;
       for (int k = 0; k < 3; ++k) {
           temp1 = -P1P2[k];
           P1P2[k] = temp1;
           temp1 = Particle1VectorsX[k];
           temp2 = Particle2VectorsX[k];
           Particle1VectorsX[k] = temp2;
           Particle2VectorsX[k] = temp1;
           temp1 = Particle1VectorsY[k];
           temp2 = Particle2VectorsY[k];
           Particle1VectorsY[k] = temp2;
           Particle2VectorsY[k] = temp1;
           temp1 = Particle1VectorsZ[k];
           temp2 = Particle2VectorsZ[k];
           Particle1VectorsZ[k] = temp2;
           Particle2VectorsZ[k] = temp1;
        }
    } else if (AngleBetweenP1P2AndZNegative == MinimumAngle) { // Particle 1 and 2 assignment is swapped
       swapz = true;
       for (int k = 0; k < 3; ++k) {
           temp1 = -P1P2[k];
           P1P2[k] = temp1;
           temp1 = Particle1VectorsX[k];
           temp2 = Particle2VectorsX[k];
           Particle1VectorsX[k] = temp2;
           Particle2VectorsX[k] = temp1;
           temp1 = Particle1VectorsY[k];
           temp2 = Particle2VectorsY[k];
           Particle1VectorsY[k] = temp2;
           Particle2VectorsY[k] = temp1;
           temp1 = Particle1VectorsZ[k];
           temp2 = Particle2VectorsZ[k];
           Particle1VectorsZ[k] = temp2;
           Particle2VectorsZ[k] = temp1;
        }
    }

    // Define the angles to rotate cube 2 about cube 1
    double VirtualMoveAlpha, VirtualMoveBeta, VirtualMoveGamma;

    /*double sinBeta1 = Particle1VectorsX[2];
    double cosBeta1 = sqrt(1-sinBeta1*sinBeta1);
    double sinAlpha1 = -Particle1VectorsY[2]/cosBeta1;
    double cosAlpha1 = Particle1VectorsZ[2]/cosBeta1;
    double sinGamma1 = -Particle1VectorsX[1]/cosBeta1;
    double cosGamma1 = Particle1VectorsX[0]/cosBeta1;

    double sinBeta2 = Particle2VectorsX[2];
    double cosBeta2 = sqrt(1-sinBeta2*sinBeta2);
    double sinAlpha2 = -Particle2VectorsY[2]/cosBeta2;
    double cosAlpha2 = Particle2VectorsZ[2]/cosBeta2;
    double sinGamma2 = -Particle2VectorsX[1]/cosBeta2;
    double cosGamma2 = Particle2VectorsX[0]/cosBeta2;
    */

    // Generate rotation matrix for particle 1
    double R1[3][3] = {{Particle1VectorsX[0], Particle1VectorsX[1], Particle1VectorsX[2]},
                       {Particle1VectorsY[0], Particle1VectorsY[1], Particle1VectorsY[2]},
                       {Particle1VectorsZ[0], Particle1VectorsZ[1], Particle1VectorsZ[2]}};

    // Calculate the inverse matrix of cube 1, rotation matrix is orthogonal matrix, so inverse of it is equal to its transpose
    double inv_1[3][3];
    //cout << "Inverse rotation matrix" << endl;
    for (int k = 0; k < 3; ++k) {
        for (int l = 0; l < 3; ++l) {
            inv_1[k][l] = R1[l][k];
            //cout << inv_1[k][l] << " " ;
        }
        //cout << endl;
    }

    Vector3 RotatedParticle2Centroid, RotatedParticle2VectorsX, RotatedParticle2VectorsY, RotatedParticle2VectorsZ;
    //cout << "Particle2VectorsX[0] " << Particle2VectorsX[0] << " Particle2VectorsX[1] " << Particle2VectorsX[1] << " Particle2VectorsX[2] " << Particle2VectorsX[2] << endl;
    MultiplyMatrixVector(inv_1, P1P2, RotatedParticle2Centroid);
    MultiplyMatrixVector(inv_1, Particle2VectorsX, RotatedParticle2VectorsX);
    //cout << "RotatedParticle2VectorsX[0] " << RotatedParticle2VectorsX[0] << " RotatedParticle2VectorsX[1] " << RotatedParticle2VectorsX[1] << " RotatedParticle2VectorsX[2] " << RotatedParticle2VectorsX[2] << endl;
    MultiplyMatrixVector(inv_1, Particle2VectorsY, RotatedParticle2VectorsY);
    //cout << "RotatedParticle2VectorsY[0] " << RotatedParticle2VectorsY[0] << " RotatedParticle2VectorsY[1] " << RotatedParticle2VectorsY[1] << " RotatedParticle2VectorsY[2] " << RotatedParticle2VectorsY[2] << endl;
    MultiplyMatrixVector(inv_1, Particle2VectorsZ, RotatedParticle2VectorsZ);
    //cout << "RotatedParticle2VectorsZ[0] " << RotatedParticle2VectorsZ[0] << " RotatedParticle2VectorsZ[1] " << RotatedParticle2VectorsZ[1] << " RotatedParticle2VectorsZ[2] " << RotatedParticle2VectorsZ[2] << endl;

    double CubeVertex[8][3];
    // Find the vertices of rotated particle 2
    findVerticesOfCube(Cube2SideLength,Cube2Height,RotatedParticle2Centroid[0],RotatedParticle2Centroid[1],RotatedParticle2Centroid[2],RotatedParticle2VectorsX,RotatedParticle2VectorsY,RotatedParticle2VectorsZ,CubeVertex);

    double RotatedParticleScaledX = RotatedParticle2Centroid[0]/Cube2SideLength;
    double RotatedParticleScaledY = RotatedParticle2Centroid[1]/Cube2SideLength;
    double RotatedParticleScaledZ = RotatedParticle2Centroid[2]/Cube2Height;
    double Norm = sqrt(RotatedParticleScaledX*RotatedParticleScaledX+RotatedParticleScaledY*RotatedParticleScaledY+RotatedParticleScaledZ*RotatedParticleScaledZ);
    // Now again compute the angles between rotated particle P1P2 vector and coordinate axes
    AngleBetweenP1P2AndX = std::acos(std::max(-1.0, std::min(1.0, RotatedParticleScaledX / Norm)));
    AngleBetweenP1P2AndXNegative = std::acos(std::max(-1.0, std::min(1.0, -RotatedParticleScaledX / Norm)));
    AngleBetweenP1P2AndY = std::acos(std::max(-1.0, std::min(1.0, RotatedParticleScaledY / Norm)));
    AngleBetweenP1P2AndYNegative = std::acos(std::max(-1.0, std::min(1.0, -RotatedParticleScaledY / Norm)));
    AngleBetweenP1P2AndZ = std::acos(std::max(-1.0, std::min(1.0, RotatedParticleScaledZ / Norm)));
    AngleBetweenP1P2AndZNegative = std::acos(std::max(-1.0, std::min(1.0, -RotatedParticleScaledZ / Norm)));
    // Find out which axis (and direction) has the smallest angle with the rotated P1P2 vector
    MinimumAngle = min({AngleBetweenP1P2AndX, AngleBetweenP1P2AndXNegative, AngleBetweenP1P2AndY,
                                    AngleBetweenP1P2AndYNegative, AngleBetweenP1P2AndZ, AngleBetweenP1P2AndZNegative});

    //cout << "Minimum angle " << MinimumAngle << endl;
    //cout << AngleBetweenP1P2AndX << " " << AngleBetweenP1P2AndXNegative << " " << AngleBetweenP1P2AndY << " " << AngleBetweenP1P2AndYNegative << " " << AngleBetweenP1P2AndZ << " " << AngleBetweenP1P2AndZNegative << endl;

    double aspectratio = Cube2Height/Cube2SideLength;
    // Now, you can determine the interacting facet y and z coordinates and final vertices of second particle
    if (AngleBetweenP1P2AndX == MinimumAngle) {
        Particle2Centroid[0] = RotatedParticle2Centroid[0];
        Particle2Centroid[1] = RotatedParticle2Centroid[1];
        Particle2Centroid[2] = RotatedParticle2Centroid[2];
        y_coord11[0] = -lc; //{-lc,-lc,lc,lc};
        y_coord11[1] = lc;
        y_coord11[2] = lc;
        y_coord11[3] = -lc;
        z_coord11[0] = lc*aspectratio; //{-lc*aspectratio,lc*aspectratio,lc*aspectratio,-lc*aspectratio};
        z_coord11[1] = lc*aspectratio;
        z_coord11[2] = -lc*aspectratio;
        z_coord11[3] = -lc*aspectratio;
        ds = lc;
        for (int k = 0; k < 8; ++k) {
           Vertex[k][0] = CubeVertex[k][0];
           Vertex[k][1] = CubeVertex[k][1];
           Vertex[k][2] = CubeVertex[k][2];
        }
    } else if (AngleBetweenP1P2AndY == MinimumAngle) {
        Particle2Centroid[0] = RotatedParticle2Centroid[1];
        Particle2Centroid[1] = -RotatedParticle2Centroid[0];
        Particle2Centroid[2] = RotatedParticle2Centroid[2];
        y_coord11[0] = -lc; //{-lc,-lc,lc,lc};
        y_coord11[1] = lc;
        y_coord11[2] = lc;
        y_coord11[3] = -lc;
        z_coord11[0] = lc*aspectratio; //{-lc*aspectratio,lc*aspectratio,lc*aspectratio,-lc*aspectratio};
        z_coord11[1] = lc*aspectratio;
        z_coord11[2] = -lc*aspectratio;
        z_coord11[3] = -lc*aspectratio;
        ds = lc;
        for (int k = 0; k < 8; ++k) {
           Vertex[k][0] = CubeVertex[k][1];
           Vertex[k][1] = -CubeVertex[k][0];
           Vertex[k][2] = CubeVertex[k][2];
        }
    } else if (AngleBetweenP1P2AndZ == MinimumAngle) {
        Particle2Centroid[0] = RotatedParticle2Centroid[2];
        Particle2Centroid[1] = RotatedParticle2Centroid[1];
        Particle2Centroid[2] = -RotatedParticle2Centroid[0];
        y_coord11[0] = -lc; //{-lc,-lc,lc,lc};
        y_coord11[1] = lc;
        y_coord11[2] = lc;
        y_coord11[3] = -lc;
        z_coord11[0] = lc; //{-lc,lc,lc,-lc};
        z_coord11[1] = lc;
        z_coord11[2] = -lc;
        z_coord11[3] = -lc;
        ds = lc*aspectratio;
        for (int k = 0; k < 8; ++k) {
           Vertex[k][0] = CubeVertex[k][2];
           Vertex[k][1] = CubeVertex[k][1];
           Vertex[k][2] = -CubeVertex[k][0];
        }
    } else if (AngleBetweenP1P2AndXNegative == MinimumAngle) { // in case reorienting centroid again results in -x axis
        Particle2Centroid[0] = -RotatedParticle2Centroid[0];
        Particle2Centroid[1] = -RotatedParticle2Centroid[1];
        Particle2Centroid[2] = -RotatedParticle2Centroid[2];
        y_coord11[0] = -lc; //{-lc,-lc,lc,lc};
        y_coord11[1] = lc;
        y_coord11[2] = lc;
        y_coord11[3] = -lc;
        z_coord11[0] = lc*aspectratio; //{-lc*aspectratio,lc*aspectratio,lc*aspectratio,-lc*aspectratio};
        z_coord11[1] = lc*aspectratio;
        z_coord11[2] = -lc*aspectratio;
        z_coord11[3] = -lc*aspectratio;
        ds = lc;
        for (int k = 0; k < 8; ++k) {
           Vertex[k][0] = -CubeVertex[k][0];
           Vertex[k][1] = -CubeVertex[k][1];
           Vertex[k][2] = -CubeVertex[k][2];
        }
    } else if (AngleBetweenP1P2AndYNegative == MinimumAngle) {
        Particle2Centroid[0] = -RotatedParticle2Centroid[1];
        Particle2Centroid[1] = RotatedParticle2Centroid[0];
        Particle2Centroid[2] = -RotatedParticle2Centroid[2];
        y_coord11[0] = -lc; //{-lc,-lc,lc,lc};
        y_coord11[1] = lc;
        y_coord11[2] = lc;
        y_coord11[3] = -lc;
        z_coord11[0] = lc*aspectratio; //{-lc*aspectratio,lc*aspectratio,lc*aspectratio,-lc*aspectratio};
        z_coord11[1] = lc*aspectratio;
        z_coord11[2] = -lc*aspectratio;
        z_coord11[3] = -lc*aspectratio;
        ds = lc;
        for (int k = 0; k < 8; ++k) {
           Vertex[k][0] = -CubeVertex[k][1];
           Vertex[k][1] = CubeVertex[k][0];
           Vertex[k][2] = -CubeVertex[k][2];
        }
    } else if (AngleBetweenP1P2AndZNegative == MinimumAngle) {
        Particle2Centroid[0] = -RotatedParticle2Centroid[2];
        Particle2Centroid[1] = -RotatedParticle2Centroid[1];
        Particle2Centroid[2] = RotatedParticle2Centroid[0];
        y_coord11[0] = -lc; //{-lc,-lc,lc,lc};
        y_coord11[1] = lc;
        y_coord11[2] = lc;
        y_coord11[3] = -lc;
        z_coord11[0] = lc; //{-lc,lc,lc,-lc};
        z_coord11[1] = lc;
        z_coord11[2] = -lc;
        z_coord11[3] = -lc;
        ds = lc*aspectratio;
        for (int k = 0; k < 8; ++k) {
           Vertex[k][0] = -CubeVertex[k][2];
           Vertex[k][1] = -CubeVertex[k][1];
           Vertex[k][2] = CubeVertex[k][0];
        }
    }

    // Reset cube 1 orientation to standard basis vectors
    Particle1Centroid[0] = 0.0;
    Particle1Centroid[1] = 0.0;
    Particle1Centroid[2] = 0.0;
    Particle1VectorsX[0] = 1.0;
    Particle1VectorsX[1] = 0.0;
    Particle1VectorsX[2] = 0.0;
    Particle1VectorsY[0] = 0.0;
    Particle1VectorsY[1] = 1.0;
    Particle1VectorsY[2] = 0.0;
    Particle1VectorsZ[0] = 0.0;
    Particle1VectorsZ[1] = 0.0;
    Particle1VectorsZ[2] = 1.0;

   for (int k = 0; k < 3; ++k) {
      Particle2VectorsX[k] = RotatedParticle2VectorsX[k];
      Particle2VectorsY[k] = RotatedParticle2VectorsY[k];
      Particle2VectorsZ[k] = RotatedParticle2VectorsZ[k];
   }


   // Find the minimum value along the first column
   auto minRow = min_element(Vertex, Vertex+8,
        [](const auto& a, const auto& b) {
            return a[0] < b[0];
        });

    // Print the Vertex array
    /*std::cout << "Vertex array:" << std::endl;
    for (int k = 0; k < 8; ++k) {
        std::cout << "Vertex[" << k << "]: ("
                  << Vertex[k][0] << ", "
                  << Vertex[k][1] << ", "
                  << Vertex[k][2] << ")" << std::endl;
    }*/


    int minindex = distance(Vertex, minRow); // this gives the index of vertex that has lowest x value (closest vertex of particle 2 in x-direction)
    double minx = Vertex[minindex][0];
    double miny = Vertex[minindex][1];
    double minz = Vertex[minindex][2];

    // Step 2: Collect indices with minimum x
    vector<int> minXIndices;
    for (int i = 0; i < 8; ++i) {
       if (abs(Vertex[i][0] - minx) < 0.01) {
          minXIndices.push_back(i);
       }
    }

    // Step 3: Find the one among those with smallest distance to the origin of interacting facet of particle 1
    int finalindex = minXIndices[0]; // Initialize with the first one
    for (int idx : minXIndices) {
       double dist1 = Vertex[idx][1] * Vertex[idx][1] + Vertex[idx][2] * Vertex[idx][2];
       double dist2 = Vertex[finalindex][1] * Vertex[finalindex][1] + Vertex[finalindex][2] * Vertex[finalindex][2];
       if (dist1 < dist2) {
          finalindex = idx;
       }
    }

    miny = Vertex[finalindex][1];
    minz = Vertex[finalindex][2];

    //cout << "Closest vertex " << minx << " " << miny << " " << minz << endl;

    bool swapref = false;
    // Using the closest vertex coordinates (minx,miny,minz), we can check whether this vertex is within the boundaries of interacting face of particle 1
    if (-Cube2Height/2.0 <= minz && minz <= Cube2Height/2.0 && -Cube2SideLength/2.0 <= miny && miny <= Cube2SideLength/2.0) { //if yes, no need for swap
        swapref = false;
    } else {
        swapref = true; // if no, reference particle needs to changed by swapping reference assignment
    }

    //cout << "swapref " << swapref << endl; cout << "P1P2 " << P1P2[0] << " " << P1P2[1] << " " << P1P2[2] << endl;

    if (swapref) { // reference particle is changed
        for (int k = 0; k < 3; ++k) {
           temp1 = -P1P2[k];
           P1P2[k] = temp1;
           temp1 = Particle1VectorsX[k];
           temp2 = Particle2VectorsX[k];
           Particle1VectorsX[k] = temp2;
           Particle2VectorsX[k] = temp1;
           temp1 = Particle1VectorsY[k];
           temp2 = Particle2VectorsY[k];
           Particle1VectorsY[k] = temp2;
           Particle2VectorsY[k] = temp1;
           temp1 = Particle1VectorsZ[k];
           temp2 = Particle2VectorsZ[k];
           Particle1VectorsZ[k] = temp2;
           Particle2VectorsZ[k] = temp1;
        }
        // Now, you do the same operation as done above
        // Generate rotation matrix for new particle 1 (reference particle)
        double R1[3][3] = {{Particle1VectorsX[0], Particle1VectorsX[1], Particle1VectorsX[2]},
                           {Particle1VectorsY[0], Particle1VectorsY[1], Particle1VectorsY[2]},
                           {Particle1VectorsZ[0], Particle1VectorsZ[1], Particle1VectorsZ[2]}};

        // Calculate the inverse matrix of cube 1, rotation matrix is orthogonal matrix, so inverse of it is equal to its transpose
        double inv_1[3][3];
        //cout << "Inverse rotation matrix" << endl;
        for (int k = 0; k < 3; ++k) {
           for (int l = 0; l < 3; ++l) {
              inv_1[k][l] = R1[l][k];
           //cout << inv_1[k][l] << " " ;
           }
           //cout << endl;
        }

        MultiplyMatrixVector(inv_1, P1P2, RotatedParticle2Centroid);
        MultiplyMatrixVector(inv_1, Particle2VectorsX, RotatedParticle2VectorsX);
        //cout << "RotatedParticle2VectorsX[0] " << RotatedParticle2VectorsX[0] << " RotatedParticle2VectorsX[1] " << RotatedParticle2VectorsX[1] << " RotatedParticle2VectorsX[2] " << RotatedParticle2VectorsX[2] << endl;
        MultiplyMatrixVector(inv_1, Particle2VectorsY, RotatedParticle2VectorsY);
        //cout << "RotatedParticle2VectorsY[0] " << RotatedParticle2VectorsY[0] << " RotatedParticle2VectorsY[1] " << RotatedParticle2VectorsY[1] << " RotatedParticle2VectorsY[2] " << RotatedParticle2VectorsY[2] << endl;
        MultiplyMatrixVector(inv_1, Particle2VectorsZ, RotatedParticle2VectorsZ);
        //cout << "RotatedParticle2VectorsZ[0] " << RotatedParticle2VectorsZ[0] << " RotatedParticle2VectorsZ[1] " << RotatedParticle2VectorsZ[1] << " RotatedParticle2VectorsZ[2] " << RotatedParticle2VectorsZ[2] << endl;

        double CubeVertex[8][3];
        findVerticesOfCube(Cube2SideLength,Cube2Height,RotatedParticle2Centroid[0],RotatedParticle2Centroid[1],RotatedParticle2Centroid[2],RotatedParticle2VectorsX,RotatedParticle2VectorsY,RotatedParticle2VectorsZ,CubeVertex);

        double RotatedParticleScaledX = RotatedParticle2Centroid[0]/Cube2SideLength;
        double RotatedParticleScaledY = RotatedParticle2Centroid[1]/Cube2SideLength;
        double RotatedParticleScaledZ = RotatedParticle2Centroid[2]/Cube2Height;
        double Norm = sqrt(RotatedParticleScaledX*RotatedParticleScaledX+RotatedParticleScaledY*RotatedParticleScaledY+RotatedParticleScaledZ*RotatedParticleScaledZ);
        AngleBetweenP1P2AndX = std::acos(std::max(-1.0, std::min(1.0, RotatedParticleScaledX / Norm)));
        AngleBetweenP1P2AndXNegative = std::acos(std::max(-1.0, std::min(1.0, -RotatedParticleScaledX / Norm)));
        AngleBetweenP1P2AndY = std::acos(std::max(-1.0, std::min(1.0, RotatedParticleScaledY / Norm)));
        AngleBetweenP1P2AndYNegative = std::acos(std::max(-1.0, std::min(1.0, -RotatedParticleScaledY / Norm)));
        AngleBetweenP1P2AndZ = std::acos(std::max(-1.0, std::min(1.0, RotatedParticleScaledZ / Norm)));
        AngleBetweenP1P2AndZNegative = std::acos(std::max(-1.0, std::min(1.0, -RotatedParticleScaledZ / Norm)));
        // Find out which axis (and direction) has the smallest angle with the P1P2 vector
        MinimumAngle = min({AngleBetweenP1P2AndX, AngleBetweenP1P2AndXNegative, AngleBetweenP1P2AndY,
                                    AngleBetweenP1P2AndYNegative, AngleBetweenP1P2AndZ, AngleBetweenP1P2AndZNegative});

        //cout << "Minimum angle " << MinimumAngle << endl;
        //cout << AngleBetweenP1P2AndX << " " << AngleBetweenP1P2AndXNegative << " " << AngleBetweenP1P2AndY << " " << AngleBetweenP1P2AndYNegative << " " << AngleBetweenP1P2AndZ << " " << AngleBetweenP1P2AndZNegative << endl;

        double aspectratio = Cube2Height/Cube2SideLength;

        if (AngleBetweenP1P2AndX == MinimumAngle) {
            Particle2Centroid[0] = RotatedParticle2Centroid[0];
            Particle2Centroid[1] = RotatedParticle2Centroid[1];
            Particle2Centroid[2] = RotatedParticle2Centroid[2];
            y_coord11[0] = -lc; //{-lc,lc,lc,-lc};
            y_coord11[1] = lc;
            y_coord11[2] = lc;
            y_coord11[3] = -lc;
            z_coord11[0] = lc*aspectratio; //lc*aspectratio,lc*aspectratio,-lc*aspectratio,-lc*aspectratio};
            z_coord11[1] = lc*aspectratio;
            z_coord11[2] = -lc*aspectratio;
            z_coord11[3] = -lc*aspectratio;
            ds = lc;
            for (int k = 0; k < 8; ++k) {
               Vertex[k][0] = CubeVertex[k][0];
               Vertex[k][1] = CubeVertex[k][1];
               Vertex[k][2] = CubeVertex[k][2];
            }
        } else if (AngleBetweenP1P2AndY == MinimumAngle) {
            Particle2Centroid[0] = RotatedParticle2Centroid[1];
            Particle2Centroid[1] = -RotatedParticle2Centroid[0];
            Particle2Centroid[2] = RotatedParticle2Centroid[2];
            y_coord11[0] = -lc; //{-lc,lc,lc,-lc};
            y_coord11[1] = lc;
            y_coord11[2] = lc;
            y_coord11[3] = -lc;
            z_coord11[0] = lc*aspectratio; //{lc*aspectratio,lc*aspectratio,-lc*aspectratio,-lc*aspectratio};
            z_coord11[1] = lc*aspectratio;
            z_coord11[2] = -lc*aspectratio;
            z_coord11[3] = -lc*aspectratio;
            ds = lc;
            for (int k = 0; k < 8; ++k) {
               Vertex[k][0] = CubeVertex[k][1];
               Vertex[k][1] = -CubeVertex[k][0];
               Vertex[k][2] = CubeVertex[k][2];
            }
        } else if (AngleBetweenP1P2AndZ == MinimumAngle) {
            Particle2Centroid[0] = RotatedParticle2Centroid[2];
            Particle2Centroid[1] = RotatedParticle2Centroid[1];
            Particle2Centroid[2] = -RotatedParticle2Centroid[0];
            y_coord11[0] = -lc; //{-lc,lc,lc,-lc};
            y_coord11[1] = lc;
            y_coord11[2] = lc;
            y_coord11[3] = -lc;
            z_coord11[0] = lc; //{lc,lc,-lc,-lc};
            z_coord11[1] = lc;
            z_coord11[2] = -lc;
            z_coord11[3] = -lc;
            ds = lc*aspectratio;
            for (int k = 0; k < 8; ++k) {
               Vertex[k][0] = CubeVertex[k][2];
               Vertex[k][1] = CubeVertex[k][1];
               Vertex[k][2] = -CubeVertex[k][0];
            }
        } else if (AngleBetweenP1P2AndXNegative == MinimumAngle) { // in case reorienting centroid again results in -x axis
            Particle2Centroid[0] = -RotatedParticle2Centroid[0];
            Particle2Centroid[1] = -RotatedParticle2Centroid[1];
            Particle2Centroid[2] = -RotatedParticle2Centroid[2];
            y_coord11[0] = -lc; //{-lc,lc,lc,-lc};
            y_coord11[1] = lc;
            y_coord11[2] = lc;
            y_coord11[3] = -lc;
            z_coord11[0] = lc*aspectratio; //{lc*aspectratio,lc*aspectratio,-lc*aspectratio,-lc*aspectratio};
            z_coord11[1] = lc*aspectratio;
            z_coord11[2] = -lc*aspectratio;
            z_coord11[3] = -lc*aspectratio;
            ds = lc;
            for (int k = 0; k < 8; ++k) {
               Vertex[k][0] = -CubeVertex[k][0];
               Vertex[k][1] = -CubeVertex[k][1];
               Vertex[k][2] = -CubeVertex[k][2];
            }
        } else if (AngleBetweenP1P2AndYNegative == MinimumAngle) {
            Particle2Centroid[0] = -RotatedParticle2Centroid[1];
            Particle2Centroid[1] = RotatedParticle2Centroid[0];
            Particle2Centroid[2] = -RotatedParticle2Centroid[2];
            y_coord11[0] = -lc; //{-lc,lc,lc,-lc};
            y_coord11[1] = lc;
            y_coord11[2] = lc;
            y_coord11[3] = -lc;
            z_coord11[0] = lc*aspectratio; //{lc*aspectratio,lc*aspectratio,-lc*aspectratio, -lc*aspectratio};
            z_coord11[1] = lc*aspectratio;
            z_coord11[2] = -lc*aspectratio;
            z_coord11[3] = -lc*aspectratio;
            ds = lc;
            for (int k = 0; k < 8; ++k) {
               Vertex[k][0] = -CubeVertex[k][1];
               Vertex[k][1] = CubeVertex[k][0];
               Vertex[k][2] = -CubeVertex[k][2];
            }
        } else if (AngleBetweenP1P2AndZNegative == MinimumAngle) {
            Particle2Centroid[0] = -RotatedParticle2Centroid[2];
            Particle2Centroid[1] = -RotatedParticle2Centroid[1];
            Particle2Centroid[2] = RotatedParticle2Centroid[0];
            y_coord11[0] = -lc; //-{lc,lc,lc,-lc};
            y_coord11[1] = lc;
            y_coord11[2] = lc;
            y_coord11[3] = -lc;
            z_coord11[0] = lc; //{lc,lc,-lc,-lc};
            z_coord11[1] = lc;
            z_coord11[2] = -lc;
            z_coord11[3] = -lc;
            ds = lc*aspectratio;
            for (int k = 0; k < 8; ++k) {
               Vertex[k][0] = -CubeVertex[k][2];
               Vertex[k][1] = -CubeVertex[k][1];
               Vertex[k][2] = CubeVertex[k][0];
            }
        }

        //cout << "y_cooord11 " << y_coord11[0] << " " << y_coord11[1] << " " << y_coord11[2] << " " << y_coord11[3] << endl;
        //cout << "z_cooord11 " << z_coord11[0] << " " << z_coord11[1] << " " << z_coord11[2] << " " << z_coord11[3] << endl;

        // Reset cube 1 orientation to standard basis vectors
        Particle1Centroid[0] = 0.0;
        Particle1Centroid[1] = 0.0;
        Particle1Centroid[2] = 0.0;
        Particle1VectorsX[0] = 1.0;
        Particle1VectorsX[1] = 0.0;
        Particle1VectorsX[2] = 0.0;
        Particle1VectorsY[0] = 0.0;
        Particle1VectorsY[1] = 1.0;
        Particle1VectorsY[2] = 0.0;
        Particle1VectorsZ[0] = 0.0;
        Particle1VectorsZ[1] = 0.0;
        Particle1VectorsZ[2] = 1.0;

        for (int k = 0; k < 3; ++k) {
          Particle2VectorsX[k] = RotatedParticle2VectorsX[k];
          Particle2VectorsY[k] = RotatedParticle2VectorsY[k];
          Particle2VectorsZ[k] = RotatedParticle2VectorsZ[k];
        }
   }
}


void Reorientation_triangularprism(int i, int j, double Cube2SideLength, double Cube2Height, double AASigma, double BoxLength,
                   Vector3 Particle1Centroid, Vector3 Particle2Centroid,
                   Vector3 Particle1VectorsX, Vector3 Particle1VectorsY, Vector3 Particle1VectorsZ,
                   Vector3 Particle2VectorsX, Vector3 Particle2VectorsY, Vector3 Particle2VectorsZ, double Vertex[][3],
                   vector<double>& y_coord11, vector<double>& z_coord11, double &ds, double lc) {
/* Reorients two triangular prisms to determine the reference particle and align them such that the reference particle is aligned along the coordinate axes.
Inputs:
  i, j                   - indices of the two particles (triangular prism 1 and triangular prism 2)
  Cube2SideLength        - side length of the second triangular prism
  Cube2Height            - height of the second triangular prism
  AASigma                - scaling factor for the particle geometry
  BoxLength              - the length of the periodic box
  Particle1Centroid      - centroid of the first particle (triangular prism 1)
  Particle2Centroid      - centroid of the second particle (triangular prism 2)
  Particle1VectorsX, Particle1VectorsY, Particle1VectorsZ - orientation vectors for the first triangular prism
  Particle2VectorsX, Particle2VectorsY, Particle2VectorsZ - orientation vectors for the second triangular prism
Output:
  Vertex                 - positions of the vertices of the second triangular prism after reorientation
  y_coord11, z_coord11   - y-z coordinates of the interacting face of reference particle
  ds                     - x-distance of center of interacting face of reference particle from its center
  lc                     - half of side length of reference particle in the units of AASigma
*/

    // Create temporary variables to store intermediate values
    double Size1Temp, Size2Temp;
    Vector3 Particle1CentroidTemp, Particle2CentroidTemp;
    Vector3 Particle1VectorsXTemp, Particle1VectorsYTemp, Particle1VectorsZTemp;
    Vector3 Particle2VectorsXTemp, Particle2VectorsYTemp, Particle2VectorsZTemp;

    Vector3 P1P2;
    vectorSubtraction(Particle2Centroid, Particle1Centroid, P1P2);
    for (int k = 0; k < 3; ++k) {
       P1P2[k] = P1P2[k]-BoxLength*round(P1P2[k]/BoxLength);
    }

    double P1P2ScaledX = P1P2[0]/Cube2SideLength;
    double P1P2ScaledY = P1P2[1]/Cube2SideLength;
    double P1P2ScaledZ = P1P2[2]/Cube2Height;
    double P1P2NormScaled = sqrt(P1P2ScaledX * P1P2ScaledX + P1P2ScaledY * P1P2ScaledY + P1P2ScaledZ * P1P2ScaledZ);
    //cout << "P1P2 " << P1P2[0] << " " << P1P2[1] << " " << P1P2[2] << endl;

    double AngleBetweenP1P2AndX = std::acos(std::max(-1.0, std::min(1.0, P1P2ScaledX / P1P2NormScaled)));
    double AngleBetweenP1P2AndXNegative = std::acos(std::max(-1.0, std::min(1.0, -P1P2ScaledX / P1P2NormScaled)));
    double AngleBetweenP1P2AndY = std::acos(std::max(-1.0, std::min(1.0, P1P2ScaledY / P1P2NormScaled)));
    double AngleBetweenP1P2AndYNegative = std::acos(std::max(-1.0, std::min(1.0, -P1P2ScaledY / P1P2NormScaled)));
    double AngleBetweenP1P2AndZ = std::acos(std::max(-1.0, std::min(1.0, P1P2ScaledZ / P1P2NormScaled)));
    double AngleBetweenP1P2AndZNegative = std::acos(std::max(-1.0, std::min(1.0, -P1P2ScaledZ / P1P2NormScaled)));

    // Find out which axis (and direction) has the smallest angle with the P1P2 vector
    double MinimumAngle = min({AngleBetweenP1P2AndX, AngleBetweenP1P2AndXNegative, AngleBetweenP1P2AndY,
                                    AngleBetweenP1P2AndYNegative, AngleBetweenP1P2AndZ, AngleBetweenP1P2AndZNegative});

    //cout << "Minimum angle " << MinimumAngle << endl;
    //cout << AngleBetweenP1P2AndX << " " << AngleBetweenP1P2AndXNegative << " " << AngleBetweenP1P2AndY << " " << AngleBetweenP1P2AndYNegative << " " << AngleBetweenP1P2AndZ << " " << AngleBetweenP1P2AndZNegative << endl;

    double temp1, temp2;
    bool swapx = false;
    bool swapy = false;
    bool swapz = false;

    if (AngleBetweenP1P2AndXNegative == MinimumAngle) {
       swapx = true;
       for (int k = 0; k < 3; ++k) {
           temp1 = -P1P2[k];
           P1P2[k] = temp1;
           temp1 = Particle1VectorsX[k];
           temp2 = Particle2VectorsX[k];
           Particle1VectorsX[k] = temp2;
           Particle2VectorsX[k] = temp1;
           temp1 = Particle1VectorsY[k];
           temp2 = Particle2VectorsY[k];
           Particle1VectorsY[k] = temp2;
           Particle2VectorsY[k] = temp1;
           temp1 = Particle1VectorsZ[k];
           temp2 = Particle2VectorsZ[k];
           Particle1VectorsZ[k] = temp2;
           Particle2VectorsZ[k] = temp1;
        }
    } else if (AngleBetweenP1P2AndYNegative == MinimumAngle) {
       swapy = true;
       for (int k = 0; k < 3; ++k) {
           temp1 = -P1P2[k];
           P1P2[k] = temp1;
           temp1 = Particle1VectorsX[k];
           temp2 = Particle2VectorsX[k];
           Particle1VectorsX[k] = temp2;
           Particle2VectorsX[k] = temp1;
           temp1 = Particle1VectorsY[k];
           temp2 = Particle2VectorsY[k];
           Particle1VectorsY[k] = temp2;
           Particle2VectorsY[k] = temp1;
           temp1 = Particle1VectorsZ[k];
           temp2 = Particle2VectorsZ[k];
           Particle1VectorsZ[k] = temp2;
           Particle2VectorsZ[k] = temp1;
        }
    } else if (AngleBetweenP1P2AndZNegative == MinimumAngle) {
       swapz = true;
       for (int k = 0; k < 3; ++k) {
           temp1 = -P1P2[k];
           P1P2[k] = temp1;
           temp1 = Particle1VectorsX[k];
           temp2 = Particle2VectorsX[k];
           Particle1VectorsX[k] = temp2;
           Particle2VectorsX[k] = temp1;
           temp1 = Particle1VectorsY[k];
           temp2 = Particle2VectorsY[k];
           Particle1VectorsY[k] = temp2;
           Particle2VectorsY[k] = temp1;
           temp1 = Particle1VectorsZ[k];
           temp2 = Particle2VectorsZ[k];
           Particle1VectorsZ[k] = temp2;
           Particle2VectorsZ[k] = temp1;
        }
    }

    // Define the angles to rotate cube 2 about cube 1
    double VirtualMoveAlpha, VirtualMoveBeta, VirtualMoveGamma;

    /*double sinBeta1 = Particle1VectorsX[2];
    double cosBeta1 = sqrt(1-sinBeta1*sinBeta1);
    double sinAlpha1 = -Particle1VectorsY[2]/cosBeta1;
    double cosAlpha1 = Particle1VectorsZ[2]/cosBeta1;
    double sinGamma1 = -Particle1VectorsX[1]/cosBeta1;
    double cosGamma1 = Particle1VectorsX[0]/cosBeta1;

    double sinBeta2 = Particle2VectorsX[2];
    double cosBeta2 = sqrt(1-sinBeta2*sinBeta2);
    double sinAlpha2 = -Particle2VectorsY[2]/cosBeta2;
    double cosAlpha2 = Particle2VectorsZ[2]/cosBeta2;
    double sinGamma2 = -Particle2VectorsX[1]/cosBeta2;
    double cosGamma2 = Particle2VectorsX[0]/cosBeta2;
    */
    /*if (abs(Particle1VectorsY[0]) > abs(Particle1VectorsX[0])) {
       for (int k = 0; k < 3; ++k) {
          temp1 = Particle1VectorsX[k];
          temp2 = Particle1VectorsY[k];
          Particle1VectorsX[k] = temp2;
          Particle1VectorsY[k] = temp1;
       }
    } else if (abs(Particle1VectorsZ[0]) > abs(Particle1VectorsX[0])) {
       for (int k = 0; k < 3; ++k) {
          temp1 = Particle1VectorsX[k];
          temp2 = Particle1VectorsZ[k];
          Particle1VectorsX[k] = temp2;
          Particle1VectorsZ[k] = temp1;
       }
    }

    if (abs(Particle1VectorsX[1]) > abs(Particle1VectorsY[1])) {
       for (int k = 0; k < 3; ++k) {
          temp1 = Particle1VectorsX[k];
          temp2 = Particle1VectorsY[k];
          Particle1VectorsX[k] = temp2;
          Particle1VectorsY[k] = temp1;
       }
    } else if (abs(Particle1VectorsZ[1]) > abs(Particle1VectorsY[1])) {
       for (int k = 0; k < 3; ++k) {
          temp1 = Particle1VectorsY[k];
          temp2 = Particle1VectorsZ[k];
          Particle1VectorsY[k] = temp2;
          Particle1VectorsZ[k] = temp1;
       }
    }

    if (abs(Particle1VectorsX[2]) > abs(Particle1VectorsZ[2])) {
       for (int k = 0; k < 3; ++k) {
          temp1 = Particle1VectorsX[k];
          temp2 = Particle1VectorsZ[k];
          Particle1VectorsX[k] = temp2;
          Particle1VectorsZ[k] = temp1;
       }
    } else if (abs(Particle1VectorsY[2]) > abs(Particle1VectorsZ[2])) {
       for (int k = 0; k < 3; ++k) {
          temp1 = Particle1VectorsY[k];
          temp2 = Particle1VectorsZ[k];
          Particle1VectorsY[k] = temp2;
          Particle1VectorsZ[k] = temp1;
       }
    }

   if (abs(Particle2VectorsY[0]) > abs(Particle2VectorsX[0])) {
       for (int k = 0; k < 3; ++k) {
          temp1 = Particle2VectorsX[k];
          temp2 = Particle2VectorsY[k];
          Particle2VectorsX[k] = temp2;
          Particle2VectorsY[k] = temp1;
       }
    } else if (abs(Particle2VectorsZ[0]) > abs(Particle2VectorsX[0])) {
       for (int k = 0; k < 3; ++k) {
          temp1 = Particle2VectorsX[k];
          temp2 = Particle2VectorsZ[k];
          Particle2VectorsX[k] = temp2;
          Particle2VectorsZ[k] = temp1;
       }
    }

    if (abs(Particle2VectorsX[1]) > abs(Particle2VectorsY[1])) {
       for (int k = 0; k < 3; ++k) {
          temp1 = Particle2VectorsX[k];
          temp2 = Particle2VectorsY[k];
          Particle2VectorsX[k] = temp2;
          Particle2VectorsY[k] = temp1;
       }
    } else if (abs(Particle2VectorsZ[1]) > abs(Particle2VectorsY[1])) {
       for (int k = 0; k < 3; ++k) {
          temp1 = Particle2VectorsY[k];
          temp2 = Particle2VectorsZ[k];
          Particle2VectorsY[k] = temp2;
          Particle2VectorsZ[k] = temp1;
       }
    }

    if (abs(Particle2VectorsX[2]) > abs(Particle2VectorsZ[2])) {
       for (int k = 0; k < 3; ++k) {
          temp1 = Particle2VectorsX[k];
          temp2 = Particle2VectorsZ[k];
          Particle2VectorsX[k] = temp2;
          Particle2VectorsZ[k] = temp1;
       }
    } else if (abs(Particle2VectorsY[2]) > abs(Particle2VectorsZ[2])) {
       for (int k = 0; k < 3; ++k) {
          temp1 = Particle2VectorsY[k];
          temp2 = Particle2VectorsZ[k];
          Particle2VectorsY[k] = temp2;
          Particle2VectorsZ[k] = temp1;
       }
    }*/

    // Generate rotation matrix for cube 1
    double R1[3][3] = {{Particle1VectorsX[0], Particle1VectorsX[1], Particle1VectorsX[2]},
                       {Particle1VectorsY[0], Particle1VectorsY[1], Particle1VectorsY[2]},
                       {Particle1VectorsZ[0], Particle1VectorsZ[1], Particle1VectorsZ[2]}};

    /*double inv_1x[3][3] = {{1.0, 0.0, 0.0},
                          {0.0, cosAlpha1, sinAlpha1},
                          {0.0, -sinAlpha1, cosAlpha1}};

    double inv_1y[3][3] = {{cosBeta1, 0.0, -sinBeta1},
                           {0.0, 1.0, 0.0},
                           {sinBeta1, 0.0, cosBeta1}};

    double inv_1z[3][3] = {{cosGamma1, sinGamma1, 0.0},
                           {-sinGamma1, cosGamma1, 0.0},
                           {0.0, 0.0, 1.0}};

    cout << "Inverse rotation x matrix" << endl;
    for (int k = 0; k < 3; ++k) {
        for (int l = 0; l < 3; ++l) {
            cout << inv_1x[k][l] << " " ;
        }
        cout << endl;
    }

    cout << "Inverse rotation y matrix" << endl;
    for (int k = 0; k < 3; ++k) {
        for (int l = 0; l < 3; ++l) {
            cout << inv_1y[k][l] << " " ;
        }
        cout << endl;
    }
    cout << "Inverse rotation z matrix" << endl;
    for (int k = 0; k < 3; ++k) {
        for (int l = 0; l < 3; ++l) {
            cout << inv_1z[k][l] << " " ;
        }
        cout << endl;
    }

    double R2_x[3][3] = {{1.0, 0.0, 0.0},
                          {0.0, cosAlpha2, -sinAlpha2},
                          {0.0, sinAlpha2, cosAlpha2}};

    double R2_y[3][3] = {{cosBeta2, 0.0, sinBeta2},
                           {0.0, 1.0, 0.0},
                           {-sinBeta2, 0.0, cosBeta2}};

    double R2_z[3][3] = {{cosGamma2, -sinGamma2, 0.0},
                           {sinGamma2, cosGamma2, 0.0},
                           {0.0, 0.0, 1.0}};

    double resultx[3][3], resulty[3][3], resultz[3][3];
    matrixMultiply(R2_x, inv_1x, resultx);
    matrixMultiply(R2_y, inv_1y, resulty);
    matrixMultiply(R2_z, inv_1z, resultz);

    cout << "Resultx" << endl;
    for (int k = 0; k < 3; ++k) {
        for (int l = 0; l < 3; ++l) {
            cout << resultx[k][l] << " " ;
        }
        cout << endl;
    }

    cout << "Resulty" << endl;
    for (int k = 0; k < 3; ++k) {
        for (int l = 0; l < 3; ++l) {
            cout << resulty[k][l] << " " ;
        }
        cout << endl;
    }

    cout << "Resultz" << endl;
    for (int k = 0; k < 3; ++k) {
        for (int l = 0; l < 3; ++l) {
            cout << resultz[k][l] << " " ;
        }
        cout << endl;
    }

    double ryz[3][3], inv_1[3][3];
    matrixMultiply(resulty, resultz, ryz);
    matrixMultiply(resultx, ryz, inv_1);
    */
    // Calculate the inverse matrix of cube 1, rotation matrix is orthogonal matrix, so inverse of it is equal to its transpose
    double inv_1[3][3];
    //cout << "Inverse rotation matrix" << endl;
    for (int k = 0; k < 3; ++k) {
        for (int l = 0; l < 3; ++l) {
            inv_1[k][l] = R1[l][k];
            //cout << inv_1[k][l] << " " ;
        }
        //cout << endl;
    }

    Vector3 RotatedParticle2Centroid, RotatedParticle2VectorsX, RotatedParticle2VectorsY, RotatedParticle2VectorsZ;
    //cout << "Particle2VectorsX[0] " << Particle2VectorsX[0] << " Particle2VectorsX[1] " << Particle2VectorsX[1] << " Particle2VectorsX[2] " << Particle2VectorsX[2] << endl;
    MultiplyMatrixVector(inv_1, P1P2, RotatedParticle2Centroid);
    MultiplyMatrixVector(inv_1, Particle2VectorsX, RotatedParticle2VectorsX);
    //cout << "RotatedParticle2VectorsX[0] " << RotatedParticle2VectorsX[0] << " RotatedParticle2VectorsX[1] " << RotatedParticle2VectorsX[1] << " RotatedParticle2VectorsX[2] " << RotatedParticle2VectorsX[2] << endl;
    MultiplyMatrixVector(inv_1, Particle2VectorsY, RotatedParticle2VectorsY);
    //cout << "RotatedParticle2VectorsY[0] " << RotatedParticle2VectorsY[0] << " RotatedParticle2VectorsY[1] " << RotatedParticle2VectorsY[1] << " RotatedParticle2VectorsY[2] " << RotatedParticle2VectorsY[2] << endl;
    MultiplyMatrixVector(inv_1, Particle2VectorsZ, RotatedParticle2VectorsZ);
    //cout << "RotatedParticle2VectorsZ[0] " << RotatedParticle2VectorsZ[0] << " RotatedParticle2VectorsZ[1] " << RotatedParticle2VectorsZ[1] << " RotatedParticle2VectorsZ[2] " << RotatedParticle2VectorsZ[2] << endl;

    double CubeVertex[6][3];
    findVerticesOfTriangularPrism(Cube2SideLength,Cube2Height,AASigma,RotatedParticle2Centroid[0],RotatedParticle2Centroid[1],RotatedParticle2Centroid[2],RotatedParticle2VectorsX,RotatedParticle2VectorsY,RotatedParticle2VectorsZ,CubeVertex);

    double RotatedParticleScaledX = RotatedParticle2Centroid[0]/Cube2SideLength;
    double RotatedParticleScaledY = RotatedParticle2Centroid[1]/Cube2SideLength;
    double RotatedParticleScaledZ = RotatedParticle2Centroid[2]/Cube2Height;
    double Norm = sqrt(RotatedParticleScaledX*RotatedParticleScaledX+RotatedParticleScaledY*RotatedParticleScaledY+RotatedParticleScaledZ*RotatedParticleScaledZ);
    AngleBetweenP1P2AndX = std::acos(std::max(-1.0, std::min(1.0, RotatedParticleScaledX / Norm)));
    AngleBetweenP1P2AndXNegative = std::acos(std::max(-1.0, std::min(1.0, -RotatedParticleScaledX / Norm)));
    AngleBetweenP1P2AndY = std::acos(std::max(-1.0, std::min(1.0, RotatedParticleScaledY / Norm)));
    AngleBetweenP1P2AndYNegative = std::acos(std::max(-1.0, std::min(1.0, -RotatedParticleScaledY / Norm)));
    AngleBetweenP1P2AndZ = std::acos(std::max(-1.0, std::min(1.0, RotatedParticleScaledZ / Norm)));
    AngleBetweenP1P2AndZNegative = std::acos(std::max(-1.0, std::min(1.0, -RotatedParticleScaledZ / Norm)));
    // Find out which axis (and direction) has the smallest angle with the P1P2 vector
    MinimumAngle = min({AngleBetweenP1P2AndX, AngleBetweenP1P2AndXNegative, AngleBetweenP1P2AndY,
                                    AngleBetweenP1P2AndYNegative, AngleBetweenP1P2AndZ, AngleBetweenP1P2AndZNegative});

    //cout << "Minimum angle " << MinimumAngle << endl;
    //cout << AngleBetweenP1P2AndX << " " << AngleBetweenP1P2AndXNegative << " " << AngleBetweenP1P2AndY << " " << AngleBetweenP1P2AndYNegative << " " << AngleBetweenP1P2AndZ << " " << AngleBetweenP1P2AndZNegative << endl;

    double aspectratio = Cube2Height/Cube2SideLength;

    if (AngleBetweenP1P2AndX == MinimumAngle) {
        Particle2Centroid[0] = RotatedParticle2Centroid[0];
        Particle2Centroid[1] = RotatedParticle2Centroid[1];
        Particle2Centroid[2] = RotatedParticle2Centroid[2];
        y_coord11.push_back(-lc); //{-lc,-lc,lc,lc};
        y_coord11.push_back(-lc);
        y_coord11.push_back(lc);
        y_coord11.push_back(lc);
        z_coord11.push_back(-lc*aspectratio); //{-lc*aspectratio,lc*aspectratio,lc*aspectratio,-lc*aspectratio};
        z_coord11.push_back(lc*aspectratio);
        z_coord11.push_back(lc*aspectratio);
        z_coord11.push_back(-lc*aspectratio);
        ds = (lc-0.5)/sqrt(3)+0.5;
        for (int k = 0; k < 6; ++k) {
           Vertex[k][0] = CubeVertex[k][0];
           Vertex[k][1] = CubeVertex[k][1];
           Vertex[k][2] = CubeVertex[k][2];
        }
    } else if (AngleBetweenP1P2AndY == MinimumAngle) {
        Particle2Centroid[0] = -0.5*RotatedParticle2Centroid[0]+0.5*RotatedParticle2Centroid[1]*sqrt(3);
        Particle2Centroid[1] = -0.5*RotatedParticle2Centroid[0]*sqrt(3)-0.5*RotatedParticle2Centroid[1];
        Particle2Centroid[2] = RotatedParticle2Centroid[2];
        y_coord11.push_back(-lc); //{-lc,-lc,lc,lc};
        y_coord11.push_back(-lc);
        y_coord11.push_back(lc);
        y_coord11.push_back(lc);
        z_coord11.push_back(-lc*aspectratio); //{-lc*aspectratio,lc*aspectratio,lc*aspectratio,-lc*aspectratio};
        z_coord11.push_back(lc*aspectratio);
        z_coord11.push_back(lc*aspectratio);
        z_coord11.push_back(-lc*aspectratio);
        ds = (lc-0.5)/sqrt(3)+0.5;
        for (int k = 0; k < 6; ++k) {
           Vertex[k][0] = -0.5*CubeVertex[k][0]+0.5*CubeVertex[k][1]*sqrt(3);
           Vertex[k][1] = -0.5*CubeVertex[k][0]*sqrt(3)-0.5*CubeVertex[k][1];
           Vertex[k][2] = CubeVertex[k][2];
        }
    } else if (AngleBetweenP1P2AndZ == MinimumAngle) {
        Particle2Centroid[0] = RotatedParticle2Centroid[2];
        Particle2Centroid[1] = RotatedParticle2Centroid[1];
        Particle2Centroid[2] = -RotatedParticle2Centroid[0];
        y_coord11.push_back(-lc);
        y_coord11.push_back(lc);
        y_coord11.push_back(0.0);
        z_coord11.push_back(-(lc-0.5)/sqrt(3));
        z_coord11.push_back(-(lc-0.5)/sqrt(3));
        z_coord11.push_back((lc-0.5)*2.0/sqrt(3));
        ds = lc*aspectratio;
        for (int k = 0; k < 6; ++k) {
           Vertex[k][0] = CubeVertex[k][2];
           Vertex[k][1] = CubeVertex[k][1];
           Vertex[k][2] = -CubeVertex[k][0];
        }
    } else if (AngleBetweenP1P2AndXNegative == MinimumAngle) { // in case reorienting centroid again results in -x axis
        if (AngleBetweenP1P2AndY <= AngleBetweenP1P2AndYNegative) {
           Particle2Centroid[0] = -0.5*RotatedParticle2Centroid[0]+0.5*RotatedParticle2Centroid[1]*sqrt(3);
           Particle2Centroid[1] = -0.5*RotatedParticle2Centroid[0]*sqrt(3)-0.5*RotatedParticle2Centroid[1];
           Particle2Centroid[2] = RotatedParticle2Centroid[2];
           y_coord11.push_back(-lc); //{-lc,-lc,lc,lc};
           y_coord11.push_back(-lc);
           y_coord11.push_back(lc);
           y_coord11.push_back(lc);
           z_coord11.push_back(-lc*aspectratio); //{-lc*aspectratio,lc*aspectratio,lc*aspectratio,-lc*aspectratio};
           z_coord11.push_back(lc*aspectratio);
           z_coord11.push_back(lc*aspectratio);
           z_coord11.push_back(-lc*aspectratio);
           ds = (lc-0.5)/sqrt(3)+0.5;
           for (int k = 0; k < 6; ++k) {
              Vertex[k][0] = -0.5*CubeVertex[k][0]+0.5*CubeVertex[k][1]*sqrt(3);
              Vertex[k][1] = -0.5*CubeVertex[k][0]*sqrt(3)-0.5*CubeVertex[k][1];
              Vertex[k][2] = CubeVertex[k][2];
           }
        } else {
           Particle2Centroid[0] = -0.5*RotatedParticle2Centroid[0]-0.5*RotatedParticle2Centroid[1]*sqrt(3);
           Particle2Centroid[1] = 0.5*RotatedParticle2Centroid[0]*sqrt(3)-0.5*RotatedParticle2Centroid[1];
           Particle2Centroid[2] = RotatedParticle2Centroid[2];
           y_coord11.push_back(-lc); //{-lc,-lc,lc,lc};
           y_coord11.push_back(-lc);
           y_coord11.push_back(lc);
           y_coord11.push_back(lc);
           z_coord11.push_back(-lc*aspectratio); //{-lc*aspectratio,lc*aspectratio,lc*aspectratio,-lc*aspectratio};
           z_coord11.push_back(lc*aspectratio);
           z_coord11.push_back(lc*aspectratio);
           z_coord11.push_back(-lc*aspectratio);
           ds = (lc-0.5)/sqrt(3)+0.5;
           for (int k = 0; k < 6; ++k) {
              Vertex[k][0] = -0.5*CubeVertex[k][0]-0.5*CubeVertex[k][1]*sqrt(3);
              Vertex[k][1] = 0.5*CubeVertex[k][0]*sqrt(3)-0.5*CubeVertex[k][1];
              Vertex[k][2] = CubeVertex[k][2];
           }
        }
    } else if (AngleBetweenP1P2AndYNegative == MinimumAngle) {
        Particle2Centroid[0] = -0.5*RotatedParticle2Centroid[0]-0.5*RotatedParticle2Centroid[1]*sqrt(3);
        Particle2Centroid[1] = 0.5*RotatedParticle2Centroid[0]*sqrt(3)-0.5*RotatedParticle2Centroid[1];
        Particle2Centroid[2] = RotatedParticle2Centroid[2];
        y_coord11.push_back(-lc); //{-lc,-lc,lc,lc};
        y_coord11.push_back(-lc);
        y_coord11.push_back(lc);
        y_coord11.push_back(lc);
        z_coord11.push_back(-lc*aspectratio); //{-lc*aspectratio,lc*aspectratio,lc*aspectratio,-lc*aspectratio};
        z_coord11.push_back(lc*aspectratio);
        z_coord11.push_back(lc*aspectratio);
        z_coord11.push_back(-lc*aspectratio);
        ds = (lc-0.5)/sqrt(3)+0.5;
        for (int k = 0; k < 6; ++k) {
           Vertex[k][0] = -0.5*CubeVertex[k][0]-0.5*CubeVertex[k][1]*sqrt(3);
           Vertex[k][1] = 0.5*CubeVertex[k][0]*sqrt(3)-0.5*CubeVertex[k][1];
           Vertex[k][2] = CubeVertex[k][2];
        }
    } else if (AngleBetweenP1P2AndZNegative == MinimumAngle) {
        Particle2Centroid[0] = -RotatedParticle2Centroid[2];
        Particle2Centroid[1] = RotatedParticle2Centroid[1];
        Particle2Centroid[2] = -RotatedParticle2Centroid[0];
        y_coord11.push_back(-lc);
        y_coord11.push_back(lc);
        y_coord11.push_back(0.0);
        z_coord11.push_back(-(lc-0.5)/sqrt(3));
        z_coord11.push_back(-(lc-0.5)/sqrt(3));
        z_coord11.push_back((lc-0.5)*2.0/sqrt(3));
        ds = lc*aspectratio;
        for (int k = 0; k < 6; ++k) {
           Vertex[k][0] = -CubeVertex[k][2];
           Vertex[k][1] = CubeVertex[k][1];
           Vertex[k][2] = -CubeVertex[k][0];
        }
    }
    // Reset cube 1 orientation to standard basis vectors
    Particle1Centroid[0] = 0.0;
    Particle1Centroid[1] = 0.0;
    Particle1Centroid[2] = 0.0;
    Particle1VectorsX[0] = 1.0;
    Particle1VectorsX[1] = 0.0;
    Particle1VectorsX[2] = 0.0;
    Particle1VectorsY[0] = 0.0;
    Particle1VectorsY[1] = 1.0;
    Particle1VectorsY[2] = 0.0;
    Particle1VectorsZ[0] = 0.0;
    Particle1VectorsZ[1] = 0.0;
    Particle1VectorsZ[2] = 1.0;

   for (int k = 0; k < 3; ++k) {
      Particle2VectorsX[k] = RotatedParticle2VectorsX[k];
      Particle2VectorsY[k] = RotatedParticle2VectorsY[k];
      Particle2VectorsZ[k] = RotatedParticle2VectorsZ[k];
   }

}

double calc_ver_cuboid(double Particle2VerticesEdge[][3],double AASigma, double catr1, double atre1, double catr2, double atre2, double crep, double repe, double Epsilon, double lc, double aspectratio, double y_coord11[], double z_coord11[], double ds, bool error_flag) {
/* Calculates the interaction energy between two particles (cuboid) based on the vertex coordinates of the second particle.
The function computes the potential energy by comparing the positions of the vertices and determining the configuration of particles,
such as face-face configuration or other types of interactions (edge-face, corner-face, etc.).

Inputs:
  Particle2VerticesEdge  - Vertex coordinates of the second particle (cuboid) after reorientation
  AASigma                - Scaling factor for particle geometry (typically related to particle size)
  catr1, atre1, catr2, atre2 - Constants for energy calculation (specific to interaction model, refer to original Brian Lee paper for more detail)
  crep                   - Repulsive potential constant
  repe                   - Repulsive potential energy parameter
  Epsilon                - Energy parameter for calculating the interaction strength
  lc                     - Half of the side length of the reference particle, scaled by AASigma
  aspectratio            - Aspect ratio of the cuboid (height/side length)
  y_coord11, z_coord11   - Coordinates of the interacting face of the reference particle
  ds                     - x-distance between the centers of interacting face of reference particle from its own center
  error_flag             - Flag to track if an error occurred during calculation

Output:
  Returns the calculated total energy (utot) based on the interaction of the two particles.
*/
    vector<vector<int>> facever = {{0, 3, 2, 1}, {2, 3, 7, 6}, {4, 7, 6, 5}, {0, 4, 5, 1}, {1, 2, 6, 5}, {0, 3, 7, 4}};
    vector<vector<int>> fid2d = {{0, 3}, {0, 5}, {3, 5}, {0, 4}, {3, 4}, {0, 1}, {1, 4}, {1, 5}, {2, 3}, {2, 5}, {2, 4}, {1, 2}};
    vector<vector<int>> fid3d = {{0, 3, 5}, {0, 3, 4}, {0, 1, 4}, {0, 1, 5}, {2, 3, 5}, {2, 3, 4}, {1, 2, 4}, {1, 2, 5}};

    // Find the minimum value along the first column (corresponding to minimum value along x-direction)
    auto minRow = min_element(Particle2VerticesEdge, Particle2VerticesEdge+8,
        [](const auto& a, const auto& b) {
            return a[0] < b[0];
        });

    // Get the row index of the minimum value
    int minindex = distance(Particle2VerticesEdge, minRow);
    //cout << "minindex " << minindex << endl;
    int mincount = 0;
    double minval = Particle2VerticesEdge[minindex][0];
    //cout << "minval " << minval << endl;
    double minedgeY = Particle2VerticesEdge[minindex][1];
    double minedgeZ = Particle2VerticesEdge[minindex][2];
    vector<int> minind;

    // Calculate number of indices with smallest x to determine if the particle is in parallel, edge-face, and other configurations
    for (int i = 0; i < 8; ++i) {
       if (abs(Particle2VerticesEdge[i][0]-minval) < 0.01) {
          mincount++;
          minind.push_back(i);
       }
    }

    //cout << "comx" << comx << endl;
    //cout << "mincount " << mincount << endl;
    //cout << "mindind" << endl;
    //for (int i = 0; i < minind.size(); ++i) {
    //    std::cout << minind[i] << " ";
    //}
    //cout << endl;
    int mem_size=256;
    initial_memory(mem_size);

    double comx = 0.0;
    double comy = 0.0;
    double comz = 0.0;

    double ver[8][3];
    for (int i = 0; i < 8; ++i) {
       for (int j = 0; j < 3; ++j) {
          ver[i][j] = Particle2VerticesEdge[i][j] / AASigma; // store vertex coordinates in terms of AASigma
          if (j == 0) {
             ver[i][j] = ver[i][j] - minval/AASigma+lc; // you will add this subtraction back in the following code, but this brings closest vertex to x value of particle 1 interacting face
          }
       }
    }

    double utot, uval;
    if (mincount == 4) { // face-face configuration
       int rcount = 0;
       double tempv[4][3];

       for (int i = 0; i < mincount; ++i) {
          copy(ver[minind[i]], ver[minind[i]] + 3, tempv[i]);
       }

       /*cout <<"tempv" << endl;
       for (int i = 0; i < 4; ++i) {
          for (int j = 0; j < 3; ++j) {
             std::cout << tempv[i][j] << " ";
          }
          std::cout << std::endl;
       }*/

       double facevs[4][3];
       int column = 2; // third column
       vector<int> sortedIndices = sortIndicesByDesiredColumn(tempv,column); //sort vertices based on their z-values
       // Assign rows based on the sorted indices
       for (int i = 0; i < mincount; ++i) {
          int originalIndex = sortedIndices[i];
          copy(begin(tempv[originalIndex]), end(tempv[originalIndex]), begin(facevs[i]));
       }
       // Assign values to facev based on the sorted tempv
       double facev[4][3];

       if (facevs[2][1] > facevs[3][1]) {
          copy(facevs[3], facevs[3] + 3, facev[0]);
          copy(facevs[2], facevs[2] + 3, facev[1]);
       } else {
          copy(facevs[2], facevs[2] + 3, facev[0]);
          copy(facevs[3], facevs[3] + 3, facev[1]);
       }
       if (facevs[0][1] > facevs[1][1]) {
          copy(facevs[0], facevs[0] + 3, facev[2]);
          copy(facevs[1], facevs[1] + 3, facev[3]);
       } else {
          copy(facevs[1], facevs[1] + 3, facev[2]);
          copy(facevs[0], facevs[0] + 3, facev[3]);
       }

       /*cout <<"facev" << endl;
       for (int i = 0; i < 4; ++i) {
          for (int j = 0; j < 3; ++j) {
             std::cout << facev[i][j] << " ";
          }
          std::cout << std::endl;
       }*/

       utot = 0;
       double tempface[4][3];

       for (int j = 0; j < 4; ++j) {
          for (int k = 0; k < 3; ++k) {
             tempface[j][k] = facev[j][k];
             if (k == 0) {
                tempface[j][k] += 1 + minval/AASigma - lc - ds; // this additional plus 1 is because two touching atoms are at sigma distance
             }
          }
       }


       /*std::cout << "tempface:" << std::endl;
       for (int i = 0; i < 4; ++i) {
          for (int j = 0; j < 3; ++j) {
             std::cout << tempface[i][j] << " ";
          }
          std::cout << std::endl;
       }*/

       int N1 = 4; // number of vertices of interacting face 1
       int N2 = 4; // number of vertices of interacting face 2

       double y_coord22[]={tempface[0][1],tempface[1][1],tempface[2][1],tempface[3][1]}; // y coordinates of 4 vertices of particle 2 interacting face
       double z_coord22[]={tempface[0][2],tempface[1][2],tempface[2][2],tempface[3][2]}; // z coordinates of 4 vertices of particle 2 interacting face

       // Surface norm is the unit normal vector of Ax+By+Cz+D plane formed by these 4 vertices of interacting face of particle 2 such that {A,B,C,D}
       double A, B, C, D;
       computePlaneEquation(tempface, A, B, C, D);
       //cout << "A " << A << " B " << B << " C " << C << " D " << D <<endl;
       double surface_norm[] = {A,B,C,D};
       //cout << AASigma << " " << lc << " " << atre1 << " " << atre2 << " " << repe << " " << catr1 << " " << catr2 << " " << crep << " " << Epsilon << endl;
       int i;
       error_flag = false;
       i=get_regions(y_coord11, z_coord11, y_coord22, z_coord22, N1, N2, surface_norm, error_flag);
       //cout << "error_flag_main " << error_flag << endl;
       uval = surface_energy(AASigma,lc,atre1,atre2,repe,catr1,catr2,crep,Epsilon,y_coord11,z_coord11,y_coord22,z_coord22,N1,N2,surface_norm);
       //cout << "uval " << uval << endl;
       utot += uval;
    } else if (mincount == 2) { // face-edge configuration
       double rcount = 0;
       utot = 0;
       int fid;
       // based on which 2 vertices are the closest ones to particle 1 interacting face, the interacting face of particle 2 is determined from the list.
       if (minind[0] == 0) {
          if (minind[1] == 1) {
             fid = 0;
          } else if (minind[1] == 3) {
             fid = 1;
          } else if (minind[1] == 4) {
             fid = 2;
          }
       } else if (minind[0] == 1) {
          if (minind[1] == 2) {
             fid = 3;
          } else if (minind[1] == 5) {
             fid = 4;
          }
       } else if (minind[0] == 2) {
          if (minind[1] == 3) {
             fid = 5;
          } else if (minind[1] == 6) {
             fid = 6;
          }
       } else if (minind[0] == 3) {
          if (minind[1] == 7) {
             fid = 7;
          }
       } else if (minind[0] == 4) {
          if (minind[1] == 5) {
             fid = 8;
          } else if (minind[1] == 7) {
             fid = 9;
          }
       } else if (minind[0] == 5) {
          if (minind[1] == 6) {
             fid = 10;
          }
       } else if (minind[0] == 6) {
          if (minind[1] == 7) {
             fid = 11;
          }
       }
       double tempv[4][3];
       for (int i1 = 0; i1 < 2; ++i1) {
          for (int i2 = 0; i2 < 4; ++i2) {
             int row = fid2d[fid][i1];
             int faceIndex = facever[row][i2]; // the interacting face determined based on the information of closest 2 vertices
             for (int i3 = 0; i3 < 3; ++i3) {
                tempv[i2][i3] = ver[faceIndex][i3];
             }
          }
          /*cout <<"tempv" << endl;
          for (int i = 0; i < 4; ++i) {
             for (int j = 0; j < 3; ++j) {
                std::cout << tempv[i][j] << " ";
             }
             std::cout << std::endl;
          }*/
          double facevs[4][3];
          int column = 2; // third column
          vector<int> sortedIndices = sortIndicesByDesiredColumn(tempv,column); // sort them based on their z-values
          for (int i = 0; i < 4; ++i) {
             int originalIndex = sortedIndices[i];
             copy(begin(tempv[originalIndex]), end(tempv[originalIndex]), begin(facevs[i]));
          }

          double facev[4][3];
          if (facevs[2][1] > facevs[3][1]) {
             facev[0][0] = facevs[3][0];
             facev[0][1] = facevs[3][1];
             facev[0][2] = facevs[3][2];

             facev[1][0] = facevs[2][0];
             facev[1][1] = facevs[2][1];
             facev[1][2] = facevs[2][2];
          } else {
             facev[0][0] = facevs[2][0];
             facev[0][1] = facevs[2][1];
             facev[0][2] = facevs[2][2];

             facev[1][0] = facevs[3][0];
             facev[1][1] = facevs[3][1];
             facev[1][2] = facevs[3][2];
          }

          if (facevs[0][1] > facevs[1][1]) {
             facev[2][0] = facevs[0][0];
             facev[2][1] = facevs[0][1];
             facev[2][2] = facevs[0][2];

             facev[3][0] = facevs[1][0];
             facev[3][1] = facevs[1][1];
             facev[3][2] = facevs[1][2];
          } else {
             facev[2][0] = facevs[1][0];
             facev[2][1] = facevs[1][1];
             facev[2][2] = facevs[1][2];

             facev[3][0] = facevs[0][0];
             facev[3][1] = facevs[0][1];
             facev[3][2] = facevs[0][2];
          }
          /*cout <<"facev" << endl;
          for (int i = 0; i < 4; ++i) {
             for (int j = 0; j < 3; ++j) {
                std::cout << facev[i][j] << " ";
             }
             std::cout << std::endl;
          }*/

          double tempface[4][3];

          for (int j = 0; j < 4; ++j) {
             for (int k = 0; k < 3; ++k) {
                tempface[j][k] = facev[j][k];
                if (k == 0) {
                   tempface[j][k] += 1 + minval/AASigma - lc - ds; // this additional plus 1 is because two touching atoms are at 1 sigma distance
                }
             }
          }

             /*std::cout << "tempface:" << std::endl;
             for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 3; ++j) {
                   std::cout << tempface[i][j] << " ";
                }
                std::cout << std::endl;
             }*/

          // Printing tempregion
             /*std::cout << "tempregion:" << std::endl;
             for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 2; ++j) {
                   std::cout << tempregion[i][j] << " ";
                }
                std::cout << std::endl;
             }*/

          int N1 = 4; // numbers of vertices of particle 1 interacting face
          int N2 = 4; // numbers of vertices of particle 2 interacting face
          double y_coord22[]={tempface[0][1],tempface[1][1],tempface[2][1],tempface[3][1]}; // y coordinates of 4 vertices of particle 2 interacting face
          double z_coord22[]={tempface[0][2],tempface[1][2],tempface[2][2],tempface[3][2]}; // z coordinates of 4 vertices of particle 2 interacting face

          // Surface norm is the unit normal vector of Ax+By+Cz+D plane formed by these 4 vertices such that {A,B,C,D}
          double A, B, C, D;
          computePlaneEquation(tempface, A, B, C, D);
          //cout << "A " << A << " B " << B << " C " << C << " D " << D <<endl;
          double surface_norm[]={A,B,C,D};
          int i;
          error_flag = false;
          i=get_regions(y_coord11, z_coord11, y_coord22, z_coord22, N1, N2, surface_norm, error_flag);
          //cout << "error_flag_main " << error_flag << endl;
          //cout << "surface_norm " << surface_norm[0] << " " << surface_norm[1] << " " << surface_norm[2] << " " << surface_norm[3] << endl;
          //cout << "y_coord11 " << y_coord11[0] << " " << y_coord11[1] << " " << y_coord11[2] << " " << y_coord11[3] << endl;
          //cout << "z_coord11 " << z_coord11[0] << " " << z_coord11[1] << " " << z_coord11[2] << " " << z_coord11[3] << endl;
          //cout << "y_coord22 " << y_coord22[0] << " " << y_coord22[1] << " " << y_coord22[2] << " " << y_coord22[3] << endl;
          //cout << "z_coord22 " << z_coord22[0] << " " << z_coord22[1] << " " << z_coord22[2] << " " << z_coord22[3] << endl;
          uval = surface_energy(AASigma,lc,atre1,atre2,repe,catr1,catr2,crep,Epsilon,y_coord11,z_coord11,y_coord22,z_coord22,N1,N2,surface_norm);
          //cout << "uval " << uval << endl;
          utot += uval;
       }
    }  else if (mincount == 1) { // face-vertex configuration
       double rcount = 0;
       utot = 0;
       double tempv[4][3];

       for (int i1 = 0; i1 < 3; ++i1) {
          for (int i2 = 0; i2 < 4; ++i2) {
             for (int j = 0; j < 3; ++j) {
                int row = fid3d[minindex][i1];
                int faceindex = facever[row][i2]; // interacting face determined based on the information of closest vertex
                tempv[i2][j] = ver[faceindex][j];
             }
          }
          /*cout <<"tempv" << endl;
          for (int i = 0; i < 4; ++i) {
             for (int j = 0; j < 3; ++j) {
                std::cout << tempv[i][j] << " ";
             }
             std::cout << std::endl;
          }*/
          double facevs[4][3];
          int column = 2; // third column
          vector<int> sortedIndices = sortIndicesByDesiredColumn(tempv,column); // sort them based on their z-values
          for (int i = 0; i < 4; ++i) {
             int originalIndex = sortedIndices[i];
             copy(begin(tempv[originalIndex]), end(tempv[originalIndex]), begin(facevs[i]));
          }

          double facev[4][3];
          if (facevs[2][1] > facevs[3][1]) {
             facev[0][0] = facevs[3][0];
             facev[0][1] = facevs[3][1];
             facev[0][2] = facevs[3][2];

             facev[1][0] = facevs[2][0];
             facev[1][1] = facevs[2][1];
             facev[1][2] = facevs[2][2];
          } else {
             facev[0][0] = facevs[2][0];
             facev[0][1] = facevs[2][1];
             facev[0][2] = facevs[2][2];

             facev[1][0] = facevs[3][0];
             facev[1][1] = facevs[3][1];
             facev[1][2] = facevs[3][2];
          }

          if (facevs[0][1] > facevs[1][1]) {
             facev[2][0] = facevs[0][0];
             facev[2][1] = facevs[0][1];
             facev[2][2] = facevs[0][2];

             facev[3][0] = facevs[1][0];
             facev[3][1] = facevs[1][1];
             facev[3][2] = facevs[1][2];
          } else {
             facev[2][0] = facevs[1][0];
             facev[2][1] = facevs[1][1];
             facev[2][2] = facevs[1][2];

             facev[3][0] = facevs[0][0];
             facev[3][1] = facevs[0][1];
             facev[3][2] = facevs[0][2];
          }
          /*cout <<"facev" << endl;
          for (int i = 0; i < 4; ++i) {
             for (int j = 0; j < 3; ++j) {
                std::cout << facev[i][j] << " ";
             }
             std::cout << std::endl;
          }*/
          double tempface[4][3];

          for (int j = 0; j < 4; ++j) {
             for (int k = 0; k < 3; ++k) {
                tempface[j][k] = facev[j][k];
                if (k == 0) {
                   tempface[j][k] += 1 + minval/AASigma - lc - ds; // this additional plus 1 is because two touching atoms are at 1 sigma distance
                }
             }
          }
             /*std::cout << "tempface:" << std::endl;
             for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 3; ++j) {
                   std::cout << tempface[i][j] << " ";
                }
                std::cout << std::endl;
             }

             // Printing tempregion
             std::cout << "tempregion:" << std::endl;
             for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 2; ++j) {
                   std::cout << tempregion[i][j] << " ";
                }
                std::cout << std::endl;
             }*/

          int N1 = 4; // number of vertices of particle 1 interacting face
          int N2 = 4; // number of vertices of particle 2 interacting face
          double y_coord22[]={tempface[0][1],tempface[1][1],tempface[2][1],tempface[3][1]}; // y coordinates of 4 vertices of particle 2 interacting face
          double z_coord22[]={tempface[0][2],tempface[1][2],tempface[2][2],tempface[3][2]}; // z coordinates of 4 vertices of particle 2 interacting face

          // Surface norm is the unit normal vector of Ax+By+Cz+D plane formed by these 4 vertices such that {A,B,C,D}
          double A, B, C, D;
          computePlaneEquation(tempface, A, B, C, D);
          //cout << "A " << A << " B " << B << " C " << C << " D " << D <<endl;
          double surface_norm[]={A,B,C,D};
          int i;
          error_flag = false;
          i=get_regions(y_coord11, z_coord11, y_coord22, z_coord22, N1, N2, surface_norm, error_flag);
          //cout << "error_flag_main " << error_flag << endl;
          uval = surface_energy(AASigma,lc,atre1,atre2,repe,catr1,catr2,crep,Epsilon,y_coord11,z_coord11,y_coord22,z_coord22,N1,N2,surface_norm);
          //cout << "uval " << uval << endl;
          utot += uval;
       }
    } else {
       utot = 0;
    }

    return utot;
}

double calc_ver_triangularprism(double Particle2VerticesEdge[][3],double AASigma, double catr1, double atre1, double catr2, double atre2, double crep, double repe, double Epsilon, double lc, double aspectratio, vector<double>& y_coord11, vector<double>& z_coord11, double ds, bool error_flag) {
/*
Calculates the interaction energy between two triangular prism particles based on the vertex coordinates of the second particle.
The function computes the potential energy by comparing the positions of the vertices and determining the configuration of particles,
such as face-face, edge-face, or vertex-face interactions.

Inputs:
  Particle2VerticesEdge  - Vertex coordinates of the second particle (triangular prism) after reorientation
  AASigma                - Atom diameter
  catr1, atre1, catr2, atre2 - Constants for energy calculation specific to the interaction model (refer to the original Brian Lee paper for details)
  crep                   - Repulsive potential constant
  repe                   - Repulsive potential energy parameter
  Epsilon                - Energy parameter used for calculating the interaction strength
  lc                     - Half of the side length of the reference particle, scaled by AASigma
  aspectratio            - Aspect ratio of the triangular prism (height/side length)
  y_coord11, z_coord11   - Coordinates of the interacting face of the reference particle
  ds                     - x-distance between the centers of interacting face of the reference particle from its own center
  error_flag             - Flag to indicate if an error occurred during the calculation

Output:
  Returns the calculated total energy (utot) based on the interaction between the two particles.
*/

    vector<vector<int>> facever = {{3, 0, 1, 2}, {2, 1, 4, 5}, {5, 4, 0, 3}, {0, 4, 1}, {3, 5, 2}};
    vector<vector<int>> fid2d = {{0, 3}, {0, 2}, {2, 3}, {0, 1}, {1, 3}, {0, 4}, {1, 4}, {2, 4}, {1, 2}};
    vector<vector<int>> fid3d = {{0, 2, 3}, {0, 1, 3}, {0, 1, 4}, {0, 2, 4}, {1, 2, 3}, {1, 2, 4}};

    // Find the minimum value along the first column
    auto minRow = min_element(Particle2VerticesEdge, Particle2VerticesEdge+6,
        [](const auto& a, const auto& b) {
            return a[0] < b[0];
        });

    // Get the row index of the minimum value
    int minindex = distance(Particle2VerticesEdge, minRow);
    //cout << "minindex " << minindex << endl;
    int mincount = 0;
    double minval = Particle2VerticesEdge[minindex][0];
    //cout << "minval " << minval << endl;
    double minedgeY = Particle2VerticesEdge[minindex][1];
    double minedgeZ = Particle2VerticesEdge[minindex][2];
    vector<int> minind;

    for (int i = 0; i < 6; ++i) {
       if (abs(Particle2VerticesEdge[i][0]-minval) < 0.01) {
          mincount++;
          minind.push_back(i);
       }
    }

    //cout << "comx" << comx << endl;
    //cout << "mincount " << mincount << endl;
    //cout << "mindind" << endl;
    //for (int i = 0; i < minind.size(); ++i) {
    //    std::cout << minind[i] << " ";
    //}
    //cout << endl;
    int mem_size=256;
    initial_memory(mem_size);

    double comx = 0.0;
    double comy = 0.0;
    double comz = 0.0;

    double ver[6][3];
    for (int i = 0; i < 6; ++i) {
       for (int j = 0; j < 3; ++j) {
          ver[i][j] = Particle2VerticesEdge[i][j] / AASigma;
          if (j == 0) {
             ver[i][j] = ver[i][j] - minval/AASigma+lc;
          }
       }
    }

    double utot, uval;
    if (mincount == 4) {
       int rcount = 0;
       double tempv[4][3];

       for (int i = 0; i < mincount; ++i) {
          copy(ver[minind[i]], ver[minind[i]] + 3, tempv[i]);
       }

       /*cout <<"tempv" << endl;
       for (int i = 0; i < 4; ++i) {
          for (int j = 0; j < 3; ++j) {
             std::cout << tempv[i][j] << " ";
          }
          std::cout << std::endl;
       }*/

       double facevs[4][3];
       int column = 2; // third column
       vector<int> sortedIndices = sortIndicesByDesiredColumn(tempv,column);
       // Assign rows based on the sorted indices
       for (int i = 0; i < mincount; ++i) {
          int originalIndex = sortedIndices[i];
          copy(begin(tempv[originalIndex]), end(tempv[originalIndex]), begin(facevs[i]));
       }
       // Assign values to facev based on the sorted tempv
       double facev[4][3];

       if (facevs[2][1] > facevs[3][1]) {
          copy(facevs[3], facevs[3] + 3, facev[0]);
          copy(facevs[2], facevs[2] + 3, facev[1]);
       } else {
          copy(facevs[2], facevs[2] + 3, facev[0]);
          copy(facevs[3], facevs[3] + 3, facev[1]);
       }
       if (facevs[0][1] > facevs[1][1]) {
          copy(facevs[0], facevs[0] + 3, facev[2]);
          copy(facevs[1], facevs[1] + 3, facev[3]);
       } else {
          copy(facevs[1], facevs[1] + 3, facev[2]);
          copy(facevs[0], facevs[0] + 3, facev[3]);
       }

       /*cout <<"facev" << endl;
       for (int i = 0; i < 4; ++i) {
          for (int j = 0; j < 3; ++j) {
             std::cout << facev[i][j] << " ";
          }
          std::cout << std::endl;
       }*/

       utot = 0;
       double tempface[4][3];

       for (int j = 0; j < 4; ++j) {
          for (int k = 0; k < 3; ++k) {
             tempface[j][k] = facev[j][k];
             if (k == 0) {
                tempface[j][k] += 1 + minval/AASigma - lc - ds;
             }
          }
       }


       /*std::cout << "tempface:" << std::endl;
       for (int i = 0; i < 4; ++i) {
          for (int j = 0; j < 3; ++j) {
             std::cout << tempface[i][j] << " ";
          }
          std::cout << std::endl;
       }*/

       int N1 = y_coord11.size();
       int N2 = 4;

       double y_coord22[]={tempface[0][1],tempface[1][1],tempface[2][1],tempface[3][1]};
       double z_coord22[]={tempface[0][2],tempface[1][2],tempface[2][2],tempface[3][2]};

       // Surface norm is the unit normal vector of Ax+By+Cz+D plane formed by these 4 vertices such that {A,B,C,D}
       double A, B, C, D;
       computePlaneEquation(tempface, A, B, C, D);
       //cout << "A " << A << " B " << B << " C " << C << " D " << D <<endl;
       double surface_norm[] = {A,B,C,D};
       //cout << AASigma << " " << lc << " " << atre1 << " " << atre2 << " " << repe << " " << catr1 << " " << catr2 << " " << crep << " " << Epsilon << endl;
       int i;
       error_flag = false;
       i=get_regions(y_coord11.data(), z_coord11.data(), y_coord22, z_coord22, N1, N2, surface_norm, error_flag);
       //cout << "error_flag_main " << error_flag << endl;
       uval = surface_energy(AASigma,lc,atre1,atre2,repe,catr1,catr2,crep,Epsilon,y_coord11.data(),z_coord11.data(),y_coord22,z_coord22,N1,N2,surface_norm);
       //cout << "uval " << uval << endl;
       utot += uval;
    } else if (mincount == 2) {
       double rcount = 0;
       utot = 0;
       int fid;
       if (minind[0] == 0) {
          if (minind[1] == 1) {
             fid = 0;
          } else if (minind[1] == 3) {
             fid = 1;
          } else if (minind[1] == 4) {
             fid = 2;
          }
       } else if (minind[0] == 1) {
          if (minind[1] == 2) {
             fid = 3;
          } else if (minind[1] == 4) {
             fid = 4;
          }
       } else if (minind[0] == 2) {
          if (minind[1] == 3) {
             fid = 5;
          } else if (minind[1] == 5) {
             fid = 6;
          }
       } else if (minind[0] == 3) {
          if (minind[1] == 5) {
             fid = 7;
          }
       } else if (minind[0] == 4) {
          if (minind[1] == 5) {
             fid = 8;
          }
       }

       for (int i1 = 0; i1 < 2; ++i1) {
          int row = fid2d[fid][i1];
          if (facever[row].size() == 4) {
            double tempv[4][3];
            for (int i2 = 0; i2 < 4; ++i2) {
               int faceIndex = facever[row][i2];
               for (int i3 = 0; i3 < 3; ++i3) {
                  tempv[i2][i3] = ver[faceIndex][i3];
               }
            }
            /*cout <<"tempv" << endl;
            for (int i = 0; i < 4; ++i) {
               for (int j = 0; j < 3; ++j) {
                  std::cout << tempv[i][j] << " ";
               }
               std::cout << std::endl;
            }*/
            double facevs[4][3];
            int column = 2; // third column
            vector<int> sortedIndices = sortIndicesByDesiredColumn(tempv,column);
            for (int i = 0; i < 4; ++i) {
               int originalIndex = sortedIndices[i];
               copy(begin(tempv[originalIndex]), end(tempv[originalIndex]), begin(facevs[i]));
            }

            double facev[4][3];
            if (facevs[2][1] > facevs[3][1]) {
               facev[0][0] = facevs[3][0];
               facev[0][1] = facevs[3][1];
               facev[0][2] = facevs[3][2];

               facev[1][0] = facevs[2][0];
               facev[1][1] = facevs[2][1];
               facev[1][2] = facevs[2][2];
            } else {
               facev[0][0] = facevs[2][0];
               facev[0][1] = facevs[2][1];
               facev[0][2] = facevs[2][2];

               facev[1][0] = facevs[3][0];
               facev[1][1] = facevs[3][1];
               facev[1][2] = facevs[3][2];
            }

            if (facevs[0][1] > facevs[1][1]) {
               facev[2][0] = facevs[0][0];
               facev[2][1] = facevs[0][1];
               facev[2][2] = facevs[0][2];

               facev[3][0] = facevs[1][0];
               facev[3][1] = facevs[1][1];
               facev[3][2] = facevs[1][2];
            } else {
               facev[2][0] = facevs[1][0];
               facev[2][1] = facevs[1][1];
               facev[2][2] = facevs[1][2];

               facev[3][0] = facevs[0][0];
               facev[3][1] = facevs[0][1];
               facev[3][2] = facevs[0][2];
            }
            /*cout <<"facev" << endl;
            for (int i = 0; i < 4; ++i) {
               for (int j = 0; j < 3; ++j) {
                  std::cout << facev[i][j] << " ";
               }
               std::cout << std::endl;
            }*/

            double tempface[4][3];

            for (int j = 0; j < 4; ++j) {
               for (int k = 0; k < 3; ++k) {
                  tempface[j][k] = facev[j][k];
                  if (k == 0) {
                     tempface[j][k] += 1 + minval/AASigma - lc - ds;
                  }
               }
            }

             /*std::cout << "tempface:" << std::endl;
             for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 3; ++j) {
                   std::cout << tempface[i][j] << " ";
                }
                std::cout << std::endl;
             }*/

            // Printing tempregion
             /*std::cout << "tempregion:" << std::endl;
             for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 2; ++j) {
                   std::cout << tempregion[i][j] << " ";
                }
                std::cout << std::endl;
             }*/

            int N1 = y_coord11.size();
            int N2 = 4;
            double y_coord22[]={tempface[0][1],tempface[1][1],tempface[2][1],tempface[3][1]};
            double z_coord22[]={tempface[0][2],tempface[1][2],tempface[2][2],tempface[3][2]};

            // Surface norm is the unit normal vector of Ax+By+Cz+D plane formed by these 4 vertices such that {A,B,C,D}
            double A, B, C, D;
            computePlaneEquation(tempface, A, B, C, D);
            //cout << "A " << A << " B " << B << " C " << C << " D " << D <<endl;
            double surface_norm[]={A,B,C,D};
            int i;
            error_flag = false;
            i=get_regions(y_coord11.data(), z_coord11.data(), y_coord22, z_coord22, N1, N2, surface_norm, error_flag);
            //cout << "error_flag_main " << error_flag << endl;
            //cout << "surface_norm " << surface_norm[0] << " " << surface_norm[1] << " " << surface_norm[2] << " " << surface_norm[3] << endl;
            //cout << "y_coord11 " << y_coord11[0] << " " << y_coord11[1] << " " << y_coord11[2] << " " << y_coord11[3] << endl;
            //cout << "z_coord11 " << z_coord11[0] << " " << z_coord11[1] << " " << z_coord11[2] << " " << z_coord11[3] << endl;
            //cout << "y_coord22 " << y_coord22[0] << " " << y_coord22[1] << " " << y_coord22[2] << " " << y_coord22[3] << endl;
            //cout << "z_coord22 " << z_coord22[0] << " " << z_coord22[1] << " " << z_coord22[2] << " " << z_coord22[3] << endl;
            uval = surface_energy(AASigma,lc,atre1,atre2,repe,catr1,catr2,crep,Epsilon,y_coord11.data(),z_coord11.data(),y_coord22,z_coord22,N1,N2,surface_norm);
            //cout << "uval " << uval << endl;
          } else {
            double tempv[3][3];
            for (int i2 = 0; i2 < 3; ++i2) {
               int faceIndex = facever[row][i2];
               for (int i3 = 0; i3 < 3; ++i3) {
                  tempv[i2][i3] = ver[faceIndex][i3];
               }
            }
            double tempface[3][3];

            for (int j = 0; j < 3; ++j) {
               for (int k = 0; k < 3; ++k) {
                  tempface[j][k] = tempv[j][k];
                  if (k == 0) {
                     tempface[j][k] += 1 + minval/AASigma - lc - ds;
                  }
               }
            }

             /*std::cout << "tempface:" << std::endl;
             for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                   std::cout << tempface[i][j] << " ";
                }
                std::cout << std::endl;
             }*/

            // Printing tempregion
             /*std::cout << "tempregion:" << std::endl;
             for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 2; ++j) {
                   std::cout << tempregion[i][j] << " ";
                }
                std::cout << std::endl;
             }*/

            int N1 = y_coord11.size();
            int N2 = 3;
            double y_coord22[]={tempface[0][1],tempface[1][1],tempface[2][1]};
            double z_coord22[]={tempface[0][2],tempface[1][2],tempface[2][2]};

            // Surface norm is the unit normal vector of Ax+By+Cz+D plane formed by these 4 vertices such that {A,B,C,D}
            double A, B, C, D;
            computePlaneEquation(tempface, A, B, C, D);
            //cout << "A " << A << " B " << B << " C " << C << " D " << D <<endl;
            double surface_norm[]={A,B,C,D};
            int i;
            error_flag = false;
            i=get_regions(y_coord11.data(), z_coord11.data(), y_coord22, z_coord22, N1, N2, surface_norm, error_flag);
            //cout << "error_flag_main " << error_flag << endl;
            //cout << "surface_norm " << surface_norm[0] << " " << surface_norm[1] << " " << surface_norm[2] << " " << surface_norm[3] << endl;
            //cout << "y_coord11 " << y_coord11[0] << " " << y_coord11[1] << " " << y_coord11[2] << " " << y_coord11[3] << endl;
            //cout << "z_coord11 " << z_coord11[0] << " " << z_coord11[1] << " " << z_coord11[2] << " " << z_coord11[3] << endl;
            //cout << "y_coord22 " << y_coord22[0] << " " << y_coord22[1] << " " << y_coord22[2] << " " << y_coord22[3] << endl;
            //cout << "z_coord22 " << z_coord22[0] << " " << z_coord22[1] << " " << z_coord22[2] << " " << z_coord22[3] << endl;
            uval = surface_energy(AASigma,lc,atre1,atre2,repe,catr1,catr2,crep,Epsilon,y_coord11.data(),z_coord11.data(),y_coord22,z_coord22,N1,N2,surface_norm);
            //cout << "uval " << uval << endl;
          }
          utot += uval;
       }
    }  else if (mincount == 1) {
       double rcount = 0;
       utot = 0;
       for (int i1 = 0; i1 < 3; ++i1) {
          int row = fid3d[minindex][i1];
          if (facever[row].size() == 4) {
            double tempv[4][3];
            for (int i2 = 0; i2 < 4; ++i2) {
               for (int j = 0; j < 3; ++j) {
                  int faceindex = facever[row][i2];
                  tempv[i2][j] = ver[faceindex][j];
               }
            }
            /*cout <<"tempv" << endl;
            for (int i = 0; i < 4; ++i) {
               for (int j = 0; j < 3; ++j) {
                  std::cout << tempv[i][j] << " ";
               }
               std::cout << std::endl;
            }*/
            double facevs[4][3];
            int column = 2; // third column
            vector<int> sortedIndices = sortIndicesByDesiredColumn(tempv,column);
            for (int i = 0; i < 4; ++i) {
               int originalIndex = sortedIndices[i];
               copy(begin(tempv[originalIndex]), end(tempv[originalIndex]), begin(facevs[i]));
            }

            double facev[4][3];
            if (facevs[2][1] > facevs[3][1]) {
               facev[0][0] = facevs[3][0];
               facev[0][1] = facevs[3][1];
               facev[0][2] = facevs[3][2];

               facev[1][0] = facevs[2][0];
               facev[1][1] = facevs[2][1];
               facev[1][2] = facevs[2][2];
            } else {
               facev[0][0] = facevs[2][0];
               facev[0][1] = facevs[2][1];
               facev[0][2] = facevs[2][2];

               facev[1][0] = facevs[3][0];
               facev[1][1] = facevs[3][1];
               facev[1][2] = facevs[3][2];
            }

            if (facevs[0][1] > facevs[1][1]) {
               facev[2][0] = facevs[0][0];
               facev[2][1] = facevs[0][1];
               facev[2][2] = facevs[0][2];

               facev[3][0] = facevs[1][0];
               facev[3][1] = facevs[1][1];
               facev[3][2] = facevs[1][2];
            } else {
               facev[2][0] = facevs[1][0];
               facev[2][1] = facevs[1][1];
               facev[2][2] = facevs[1][2];

               facev[3][0] = facevs[0][0];
               facev[3][1] = facevs[0][1];
               facev[3][2] = facevs[0][2];
            }
            /*cout <<"facev" << endl;
            for (int i = 0; i < 4; ++i) {
               for (int j = 0; j < 3; ++j) {
                std::cout << facev[i][j] << " ";
               }
               std::cout << std::endl;
            }*/
            double tempface[4][3];

            for (int j = 0; j < 4; ++j) {
               for (int k = 0; k < 3; ++k) {
                  tempface[j][k] = facev[j][k];
                  if (k == 0) {
                     tempface[j][k] += 1 + minval/AASigma - lc - ds;
                  }
               }
            }
             /*std::cout << "tempface:" << std::endl;
             for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 3; ++j) {
                   std::cout << tempface[i][j] << " ";
                }
                std::cout << std::endl;
             }

             // Printing tempregion
             std::cout << "tempregion:" << std::endl;
             for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 2; ++j) {
                   std::cout << tempregion[i][j] << " ";
                }
                std::cout << std::endl;
             }*/

            int N1 = y_coord11.size();
            int N2 = 4;
            double y_coord22[]={tempface[0][1],tempface[1][1],tempface[2][1],tempface[3][1]};
            double z_coord22[]={tempface[0][2],tempface[1][2],tempface[2][2],tempface[3][2]};

            // Surface norm is the unit normal vector of Ax+By+Cz+D plane formed by these 4 vertices such that {A,B,C,D}
            double A, B, C, D;
            computePlaneEquation(tempface, A, B, C, D);
            //cout << "A " << A << " B " << B << " C " << C << " D " << D <<endl;
            double surface_norm[]={A,B,C,D};
            int i;
            error_flag = false;
            i=get_regions(y_coord11.data(), z_coord11.data(), y_coord22, z_coord22, N1, N2, surface_norm, error_flag);
            //cout << "error_flag_main " << error_flag << endl;
            uval = surface_energy(AASigma,lc,atre1,atre2,repe,catr1,catr2,crep,Epsilon,y_coord11.data(),z_coord11.data(),y_coord22,z_coord22,N1,N2,surface_norm);
            //cout << "uval " << uval << endl;
          } else {
            double tempv[3][3];
            for (int i2 = 0; i2 < 3; ++i2) {
               for (int j = 0; j < 3; ++j) {
                  int faceindex = facever[row][i2];
                  tempv[i2][j] = ver[faceindex][j];
               }
            }
            /*cout <<"tempv" << endl;
            for (int i = 0; i < 3; ++i) {
               for (int j = 0; j < 3; ++j) {
                  std::cout << tempv[i][j] << " ";
               }
               std::cout << std::endl;
            }*/

            double tempface[3][3];

            for (int j = 0; j < 3; ++j) {
               for (int k = 0; k < 3; ++k) {
                  tempface[j][k] = tempv[j][k];
                  if (k == 0) {
                     tempface[j][k] += 1 + minval/AASigma - lc - ds;
                  }
               }
            }
            int N1 = y_coord11.size();
            int N2 = 3;
            double y_coord22[]={tempface[0][1],tempface[1][1],tempface[2][1]};
            double z_coord22[]={tempface[0][2],tempface[1][2],tempface[2][2]};

            // Surface norm is the unit normal vector of Ax+By+Cz+D plane formed by these 4 vertices such that {A,B,C,D}
            double A, B, C, D;
            computePlaneEquation(tempface, A, B, C, D);
            //cout << "A " << A << " B " << B << " C " << C << " D " << D <<endl;
            double surface_norm[]={A,B,C,D};
            int i;
            error_flag = false;
            i=get_regions(y_coord11.data(), z_coord11.data(), y_coord22, z_coord22, N1, N2, surface_norm, error_flag);
            //cout << "error_flag_main " << error_flag << endl;
            uval = surface_energy(AASigma,lc,atre1,atre2,repe,catr1,catr2,crep,Epsilon,y_coord11.data(),z_coord11.data(),y_coord22,z_coord22,N1,N2,surface_norm);
            //cout << "uval " << uval << endl;
          }
          utot += uval;
       }
    } else {
       utot = 0;
    }

    return utot;
}

// This function checks if two interacting particles overlap with each other or not
bool CheckOverlapForTwoParticles(const int i, const int j, double BoxLength, double AASigma, double ParticleSideLength1, double ParticleSideLength2,
                              double ParticleSideHeight1, double ParticleSideHeight2,
                              const Vector3 Particle1Centroid, const Vector3 Particle2Centroid,
                              const Vector3 Particle1VectorX, const Vector3 Particle1VectorY, const Vector3 Particle1VectorZ,
                              const Vector3 Particle2VectorX, const Vector3 Particle2VectorY, const Vector3 Particle2VectorZ, string shape) {

/* Checks whether two interacting particles overlap in a 3D space or not.

The function checks if two particles overlap by calculating and comparing their projected distances on the x, y, and z axes. It supports different particle
shapes (e.g., cuboids) and verifies their positions, distances, and configurations to determine whether they overlap or not.

Inputs:
  i, j                   - Particle indices for the two interacting particles
  BoxLength              - Length of the periodic box in each direction
  AASigma                - Atom diameter
  ParticleSideLength1    - Side length of the first particle
  ParticleSideLength2    - Side length of the second particle
  ParticleSideHeight1    - Height of the first particle
  ParticleSideHeight2    - Height of the second particle
  Particle1Centroid      - Centroid coordinates of the first particle
  Particle2Centroid      - Centroid coordinates of the second particle
  Particle1VectorX, Particle1VectorY, Particle1VectorZ - orientation vectors for the first particle
  Particle2VectorX, Particle2VectorY, Particle2VectorZ - orientation vectors for the second particle
  shape                  - Type of particle shape (e.g., "cuboid")

Output:
  Returns a boolean value indicating whether the two particles overlap:
    - `true`  if the particles overlap
    - `false` if the particles do not overlap

The function works as follows:
1. It computes the positions and orientations of the two particles based on their centroids and orientation vectors.
2. The vertices of the selected shape can be computed.
3. It then checks if the particles are far enough apart such that they cannot possibly overlap. If the particles are too far apart, the function returns `false` immediately.
4. If the particles are close enough, it checks if any vertex of one particle is inside the other or if any edge of one particle intersects with the other. This is done through geometric comparisons using their vertices and edges.
5. Lastly, the function projects both particles along the x, y, and z axes. If there is no overlap along any of these axes, it returns `false`. If overlap occurs along all axes, the function returns `true`.

Notes:
  - The function supports periodic boundary conditions, as indicated by the `BoxLength` and modulo operations when calculating vertex positions.
  - Currently, it focuses on cuboids and triangular prism, but the framework can be extended to other shapes by adding additional checks and calculations for different geometries.
  - Numeric precision issues are accounted for with small thresholds (e.g., `1e-5`).
*/

    double Size1 = ParticleSideLength1 / AASigma;
    double Size2 = ParticleSideLength2 / AASigma;
    Vector3 ParticleC1, ParticleC2, VectorX_1, VectorY_1, VectorZ_1, VectorX_2, VectorY_2, VectorZ_2;
    for (int k = 0; k < 3; ++k) {
       ParticleC1[k] = Particle1Centroid[k];
       ParticleC2[k] = Particle2Centroid[k];
       VectorX_1[k] = Particle1VectorX[k];
       VectorY_1[k] = Particle1VectorY[k];
       VectorZ_1[k] = Particle1VectorZ[k];
       VectorX_2[k] = Particle2VectorX[k];
       VectorY_2[k] = Particle2VectorY[k];
       VectorZ_2[k] = Particle2VectorZ[k];
    }

    if (shape == "cuboid") {
      double cube1initial[8][3];
      // vertices of particle 1
      findVerticesOfCube(ParticleSideLength1,ParticleSideHeight1,ParticleC1[0],ParticleC1[1],ParticleC1[2],VectorX_1,VectorY_1,VectorZ_1,cube1initial);
      double cube2initial[8][3];
      // vertices of particle 2
      findVerticesOfCube(ParticleSideLength2,ParticleSideHeight2,ParticleC2[0],ParticleC2[1],ParticleC2[2],VectorX_2,VectorY_2,VectorZ_2,cube2initial);

    /*cout << "Vertex 1 initial" << endl;
    for (int i = 0; i < 8; ++i) {
       for (int j = 0; j < 3; ++j) {
          std::cout << cube1initial[i][j] << " ";
       }
       std::cout << std::endl;
    }
    cout << "Vertex 2 initial" << endl;
    for (int i = 0; i < 8; ++i) {
       for (int j = 0; j < 3; ++j) {
          std::cout << cube2initial[i][j] << " ";
       }
       std::cout << std::endl;
    }*/
    //cout << "ParticleC1 " << ParticleC1[0] << " " << ParticleC1[1] << " " << ParticleC1[2] << endl;
    //cout << "ParticleC2 " << ParticleC2[0] << " " << ParticleC2[1] << " " << ParticleC2[2] << endl;
      double particle2Vertex[8][3];
      double y_coord11[4]; // y coordinates of interacting face of particle 1
      double z_coord11[4]; // z coordinates of interacting face of particle 1
      double ds;
      double lc = Size1/2.0; // half of particle 1 edge length
      // reorient them such that particle 1 aligns along coordinate axes
      Reorientation_cuboid(i, j, ParticleSideLength2, ParticleSideHeight2, BoxLength, ParticleC1, ParticleC2, VectorX_1, VectorY_1, VectorZ_1, VectorX_2, VectorY_2, VectorZ_2, particle2Vertex, y_coord11, z_coord11, ds, lc);

      double particle1Vertex[8][3]; // vertices of particle 1
      for (int i = 0; i < 8; ++i) {
         if(i < 4) {
              particle1Vertex[i][0] = -ds * AASigma;
              particle1Vertex[i][1] = y_coord11[i] * AASigma;
              particle1Vertex[i][2] = z_coord11[i] * AASigma;
         }
         else {
              particle1Vertex[i][0] = ds * AASigma;
              particle1Vertex[i][1] = y_coord11[i-4] * AASigma;
              particle1Vertex[i][2] = z_coord11[i-4] * AASigma;
         }
      }

      double Cube1X = ParticleC1[0];
      double Cube1Y = ParticleC1[1];
      double Cube1Z = ParticleC1[2];
      double Cube2X = ParticleC2[0];
      double Cube2Y = ParticleC2[1];
      double Cube2Z = ParticleC2[2];
      //cout << "Cube1 " << Cube1X << " " << Cube1Y << " " << Cube1Z << endl;
      //cout << "Cube2 " << Cube2X << " " << Cube2Y << " " << Cube2Z << endl;

    /*cout << "cube1Vertex " << endl;
    for (int i = 0; i < 8; ++i) {
       for (int j = 0; j < 3; ++j) {
          std::cout << cube1Vertex[i][j] << " ";
       }
       std::cout << std::endl;
    }

    cout << "cube2Vertex " << endl;
    for (int i = 0; i < 8; ++i) {
       for (int j = 0; j < 3; ++j) {
          std::cout << cube2Vertex[i][j] << " ";
       }
       std::cout << std::endl;
    }*/

      if (max(Size1, Size2) == 1) { // if spherical single particles
          // Single particle
          double center_distance_sq = pow((Cube1X - Cube2X), 2) + pow((Cube1Y - Cube2Y), 2) + pow((Cube1Z - Cube2Z), 2);
          double max_val_sq = pow((ParticleSideLength1 + ParticleSideLength2), 2) / 4.0;

          if (center_distance_sq > max_val_sq) {
              return false;
          }
      } else { // if faceted particles
          double center_distance_sq = pow((Cube1X - Cube2X), 2) + pow((Cube1Y - Cube2Y), 2) + pow((Cube1Z - Cube2Z), 2);
          double mult = max(1.0,ParticleSideHeight2/ParticleSideLength2);
          double max_val_sq = 2 * (pow((ParticleSideLength1 + ParticleSideLength2), 2)) / 4.0 + pow(ParticleSideHeight2, 2);

          if (center_distance_sq > max_val_sq) { // if the distance is larger than that, there is no chance for overlap, so return false immediately
              //cout << "Here1" << endl;
              return false;
          }

          // Check if any vertex is inside the other cube
          double ParticleSideLengthHalf1 = ParticleSideLength1/2.0;
          double CubeSideHeightHalf1 = ParticleSideHeight1/2.0;
          for (int i = 0; i < 8; ++i) {
             double vertex_x = particle2Vertex[i][0];
             double vertex_y = particle2Vertex[i][1];
             double vertex_z = particle2Vertex[i][2];
             vertex_x = vertex_x-BoxLength*round(vertex_x/BoxLength);
             vertex_y = vertex_y-BoxLength*round(vertex_y/BoxLength);
             vertex_z = vertex_z-BoxLength*round(vertex_z/BoxLength);
             //cout << "vertex " << vertex_x << " " << vertex_y << " " << vertex_z << endl;
             bool check1 = -ds * AASigma - vertex_x < -1e-5; //-CubeSideLengthHalf1 < vertex_x, -1e-5 for numeric error precision
             bool check2 = vertex_x - ds * AASigma <-1e-5; // vertex_x < CubeSideLengthHalf1, -1e-5 for numeric error precision
             bool check3 = y_coord11[0] * AASigma - vertex_y < -1e-5; //-CubeSideLengthHalf1 < vertex_y, -1e-5 for numeric error precision
             bool check4 = vertex_y - y_coord11[2] * AASigma < -1e-5; //vertex_y < CubeSideLengthHalf1, -1e-5 for numeric error precision
             bool check5 = z_coord11[0] * AASigma - vertex_z < -1e-5; //-CubeSideHeightHalf1 < vertex_z, -1e-5 for numeric error precision
             bool check6 = vertex_z - z_coord11[2] * AASigma < -1e-5; // vertex_z < CubeSideHeightHalf1, -1e-5 for numeric error precision
             if (check1 && check2 && check3 && check4 && check5 && check6) {
                //cout << "Vertex inside overlap" << endl;
                return true;
             }
          }
          // Check if any edge intersects with other cube interior
          vector<int> sortedIndices = sortIndicesByDesiredColumn(particle2Vertex,0);
          for (int i = 1; i < 4; ++i) {
             int min_ind = sortedIndices[0];
             int next_ind = sortedIndices[i];
             double slope_x = (particle2Vertex[next_ind][0]-particle2Vertex[min_ind][0])/Size2;
             double slope_y = (particle2Vertex[next_ind][1]-particle2Vertex[min_ind][1])/Size2;
             double slope_z = (particle2Vertex[next_ind][2]-particle2Vertex[min_ind][2])/Size2;
             for (int j = 0; j < Size2; ++j) {
                double line_x = slope_x*j+particle2Vertex[min_ind][0];
                double line_y = slope_y*j+particle2Vertex[min_ind][1];
                double line_z = slope_z*j+particle2Vertex[min_ind][2];
                bool check1 = -ds * AASigma - line_x < -1e-5; //-CubeSideLengthHalf1 < line_x, -1e-5 for numeric error precision
                bool check2 = line_x - ds * AASigma < -1e-5; // line_x < CubeSideLengthHalf1, -1e-5 for numeric error precision
                bool check3 = y_coord11[0] * AASigma - line_y < -1e-5; //-CubeSideLengthHalf1 < line_y, -1e-5 for numeric error precision
                bool check4 = line_y - y_coord11[2] * AASigma < -1e-5; //line_y < CubeSideLengthHalf1, -1e-5 for numeric error precision
                bool check5 = z_coord11[0] * AASigma - line_z < -1e-5; //-CubeSideHeightHalf1 < line_z, -1e-5 for numeric error precision
                bool check6 = line_z - z_coord11[2] * AASigma < -1e-5; // line_z < CubeSideHeightHalf1, -1e-5 for numeric error precision
                if (check1 && check2 && check3 && check4 && check5 && check6) {
                   //cout << "Edge intersect overlap" << endl;
                   return true;
                }
             }
          }


        /*cout << "Vertex 1 final" << endl;
        for (int i = 0; i < 8; ++i) {
           for (int j = 0; j < 3; ++j) {
              std::cout << cube1Vertex[i][j] << " ";
           }
           std::cout << std::endl;
        }
        cout << "Vertex 2 final" << endl;
        for (int i = 0; i < 8; ++i) {
           for (int j = 0; j < 3; ++j) {
              std::cout << cube2Vertex[i][j] << " ";
           }
           std::cout << std::endl;
        }*/

          // If The algorithm projects both particles on x, y, z directions; if overlap detected on all these directions, there is overlap; otherwise not.
          for (int i = 0; i < 3; ++i) {
              Vector3 axis_vec = {0.0, 0.0, 0.0};
              axis_vec[i] = 1.0;
              double max_proj1 = 0.0;
              double max_proj2 = 0.0;
              double min_proj1 = 1000000000.0; //very big number to initialize
              double min_proj2 = 1000000000.0; //very big number to initialize
              for (int j = 0; j < 8; ++j) {
                  double proj1 = particle1Vertex[j][0]*axis_vec[0]+particle1Vertex[j][1]*axis_vec[1]+particle1Vertex[j][2]*axis_vec[2];
                  double proj2 = particle2Vertex[j][0]*axis_vec[0]+particle2Vertex[j][1]*axis_vec[1]+particle2Vertex[j][2]*axis_vec[2];
                  if (proj1 < min_proj1) {
                     min_proj1 = proj1;
                  }
                  if (proj2 < min_proj2) {
                     min_proj2 = proj2;
                  }
                  if (proj1 > max_proj1) {
                     max_proj1 = proj1;
                  }
                  if (proj2 > max_proj2) {
                     max_proj2 = proj2;
                  }
              }
              //cout << axis_vec[0] << " " << axis_vec[1] << " " << axis_vec[2] << endl;
              //cout << min_proj1 << " " << min_proj2 << " " << max_proj1 << " " << max_proj2 << endl;
              bool check1 = max_proj1 - min_proj2 < 1e-5; // max_proj1 < min_proj2, 1e-5 for numeric error precision
              bool check2 = max_proj2 - min_proj1 < 1e-5; // max_proj2 < min_proj1, 1e-5 for numeric error precision
              if (check1 || check2) {
                 //cout << "No Projection overlap"<<endl;
                 return false;
              }
          }

          // If The algorithm projects both particles onto the xy, xz, yz planes; if overlap detected on all these planes, there is overlap; otherwise not.
          for (int i = 0; i < 3; ++i) {
              for (int j = 0; j < 3; ++j) {
                  if (i != j) {
                     Vector3 axis_vec = {0.0, 0.0, 0.0};
                     axis_vec[i] = 1.0;
                     axis_vec[j] = 1.0;
                     double max_proj1 = 0.0;
                     double max_proj2 = 0.0;
                     double min_proj1 = 1000000000.0; //very big number to initialize
                     double min_proj2 = 1000000000.0; //very big number to initialize
                     for (int j = 0; j < 8; ++j) {
                         double proj1 = particle1Vertex[j][0]*axis_vec[0]+particle1Vertex[j][1]*axis_vec[1]+particle1Vertex[j][2]*axis_vec[2];
                         double proj2 = particle2Vertex[j][0]*axis_vec[0]+particle2Vertex[j][1]*axis_vec[1]+particle2Vertex[j][2]*axis_vec[2];
                         if (proj1 < min_proj1) {
                            min_proj1 = proj1;
                         }
                         if (proj2 < min_proj2) {
                            min_proj2 = proj2;
                         }
                         if (proj1 > max_proj1) {
                            max_proj1 = proj1;
                         }
                         if (proj2 > max_proj2) {
                            max_proj2 = proj2;
                         }
                     }
                     //cout << axis_vec[0] << " " << axis_vec[1] << " " << axis_vec[2] << endl;
                     //cout << min_proj1 << " " << min_proj2 << " " << max_proj1 << " " << max_proj2 << endl;
                     bool check1 = max_proj1 - min_proj2 < 1e-5; // max_proj1 < min_proj2, 1e-5 for numeric error precision
                     bool check2 = max_proj2 - min_proj1 < 1e-5; // max_proj2 < min_proj1, 1e-5 for numeric error precision
                     if (check1 || check2) {
                        //cout << "No Projection 2 overlap"<<endl;
                        return false;
                     }
                  }
              }
          }
      }
      return true;
    } else if (shape == "triangularprism") {
      double cube1initial[6][3];
      findVerticesOfTriangularPrism(ParticleSideLength1,ParticleSideHeight1,AASigma,ParticleC1[0],ParticleC1[1],ParticleC1[2],VectorX_1,VectorY_1,VectorZ_1,cube1initial);
      double cube2initial[6][3];
      findVerticesOfTriangularPrism(ParticleSideLength2,ParticleSideHeight2,AASigma,ParticleC2[0],ParticleC2[1],ParticleC2[2],VectorX_2,VectorY_2,VectorZ_2,cube2initial);

    /*cout << "Vertex 1 initial" << endl;
    for (int i = 0; i < 8; ++i) {
       for (int j = 0; j < 3; ++j) {
          std::cout << cube1initial[i][j] << " ";
       }
       std::cout << std::endl;
    }
    cout << "Vertex 2 initial" << endl;
    for (int i = 0; i < 8; ++i) {
       for (int j = 0; j < 3; ++j) {
          std::cout << cube2initial[i][j] << " ";
       }
       std::cout << std::endl;
    }*/
    //cout << "ParticleC1 " << ParticleC1[0] << " " << ParticleC1[1] << " " << ParticleC1[2] << endl;
    //cout << "ParticleC2 " << ParticleC2[0] << " " << ParticleC2[1] << " " << ParticleC2[2] << endl;
      double cube2Vertex[6][3];
      vector<double> y_coord11;
      vector<double> z_coord11;
      double ds;
      double lc = Size1/2.0;
      Reorientation_triangularprism(i, j, ParticleSideLength2, ParticleSideHeight2, AASigma, BoxLength, ParticleC1, ParticleC2, VectorX_1, VectorY_1, VectorZ_1, VectorX_2, VectorY_2, VectorZ_2, cube2Vertex, y_coord11, z_coord11, ds, lc);
      int size = y_coord11.size();
      //cout << "size " << size << endl;

      double cube1Vertex[6][3];
      for (int i = 0; i < 6; ++i) {
         if(i < size) {
              cube1Vertex[i][0] = ds * AASigma;
              cube1Vertex[i][1] = y_coord11[i] * AASigma;
              cube1Vertex[i][2] = z_coord11[i] * AASigma;
         }
         else {
              cube1Vertex[i][0] = -ds * AASigma;
              cube1Vertex[i][1] = y_coord11[i-size] * AASigma;
              cube1Vertex[i][2] = z_coord11[i-size] * AASigma;
         }
      }

      double Cube1X = ParticleC1[0];
      double Cube1Y = ParticleC1[1];
      double Cube1Z = ParticleC1[2];
      double Cube2X = ParticleC2[0];
      double Cube2Y = ParticleC2[1];
      double Cube2Z = ParticleC2[2];
      //cout << "Cube1 " << Cube1X << " " << Cube1Y << " " << Cube1Z << endl;
      //cout << "Cube2 " << Cube2X << " " << Cube2Y << " " << Cube2Z << endl;

    /*cout << "cube1Vertex " << endl;
    for (int i = 0; i < 6; ++i) {
       for (int j = 0; j < 3; ++j) {
          std::cout << cube1Vertex[i][j] << " ";
       }
       std::cout << std::endl;
    }

    cout << "cube2Vertex " << endl;
    for (int i = 0; i < 6; ++i) {
       for (int j = 0; j < 3; ++j) {
          std::cout << cube2Vertex[i][j] << " ";
       }
       std::cout << std::endl;
    }*/

      if (max(Size1, Size2) == 1) {
          // Single particle
          double center_distance_sq = pow((Cube1X - Cube2X), 2) + pow((Cube1Y - Cube2Y), 2) + pow((Cube1Z - Cube2Z), 2);
          double max_val_sq = pow((ParticleSideLength1 + ParticleSideLength2), 2) / 4.0;

          if (center_distance_sq > max_val_sq) {
              return false;
          }
      } else {
          double center_distance_sq = pow((Cube1X - Cube2X), 2) + pow((Cube1Y - Cube2Y), 2) + pow((Cube1Z - Cube2Z), 2);
          double mult = max(1.0,ParticleSideHeight2/ParticleSideLength2);
          double max_val_sq = 2 * (pow((ParticleSideLength1 + ParticleSideLength2), 2)) / 4.0 + pow(ParticleSideHeight2, 2);

          if (center_distance_sq > max_val_sq) {
              //cout << "Here1" << endl;
              return false;
          }

          // Check if any vertex is inside the other cube
          double CubeSideLengthHalf1 = ParticleSideLength1/2.0;
          double CubeSideHeightHalf1 = ParticleSideHeight1/2.0;
          for (int i = 0; i < 6; ++i) {
             double vertex_x = cube2Vertex[i][0];
             double vertex_y = cube2Vertex[i][1];
             double vertex_z = cube2Vertex[i][2];
             vertex_x = vertex_x-BoxLength*round(vertex_x/BoxLength);
             vertex_y = vertex_y-BoxLength*round(vertex_y/BoxLength);
             vertex_z = vertex_z-BoxLength*round(vertex_z/BoxLength);
             //cout << "vertex " << vertex_x << " " << vertex_y << " " << vertex_z << endl;
             bool check1 = -ds * AASigma - vertex_x < -1e-5; //-x_cube1 < vertex_x, -1e-5 for numeric error precision
             bool check2 = vertex_x - ds * AASigma <-1e-5; // vertex_x < x_cube1, -1e-5 for numeric error precision
             bool check3 = y_coord11[0] * AASigma - vertex_y < -1e-5; //-y_cube1 < vertex_y, -1e-5 for numeric error precision
             bool check4 = vertex_y - y_coord11[2] * AASigma < -1e-5; //vertex_y < y_cube1, -1e-5 for numeric error precision
             bool check5 = z_coord11[0] * AASigma - vertex_z < -1e-5; //-z_cube1 < vertex_z, -1e-5 for numeric error precision
             bool check6 = vertex_z - z_coord11[2] * AASigma < -1e-5; // z_cube1 < CubeSideHeightHalf1, -1e-5 for numeric error precision
             if (size == 4) {
                check1 = -ds * AASigma - vertex_x < -1e-5; //-x_cube1 < vertex_x, -1e-5 for numeric error precision
                check2 = vertex_x - ds * AASigma <-1e-5; // vertex_x < x_cube1, -1e-5 for numeric error precision
                check3 = y_coord11[0] * AASigma + (ds * AASigma - vertex_x) / sqrt(3) - vertex_y < -1e-5; //-y_cube1 < vertex_y
                check4 = vertex_y - y_coord11[2] * AASigma + (ds * AASigma - vertex_x) / sqrt(3) < -1e-5; //vertex_y < y_cube1
                check5 = z_coord11[0] * AASigma - vertex_z < -1e-5; //-z_cube1 < vertex_z, -1e-5 for numeric error precision
                check6 = vertex_z - z_coord11[2] * AASigma < -1e-5; // z_cube1 < CubeSideHeightHalf1, -1e-5 for numeric error precision
             }
             if (check1 && check2 && check3 && check4 && check5 && check6) {
                //cout << "Vertex inside overlap" << endl;
                return true;
             }
          }
          // Check if any edge intersects with other cube interior
          vector<int> sortedIndices = sortIndicesByDesiredColumn(cube2Vertex,0);
          for (int i = 1; i < 4; ++i) {
             int min_ind = sortedIndices[0];
             int next_ind = sortedIndices[i];
             double slope_x = (cube2Vertex[next_ind][0]-cube2Vertex[min_ind][0])/Size2;
             double slope_y = (cube2Vertex[next_ind][1]-cube2Vertex[min_ind][1])/Size2;
             double slope_z = (cube2Vertex[next_ind][2]-cube2Vertex[min_ind][2])/Size2;
             for (int j = 0; j < Size2; ++j) {
                double line_x = slope_x*j+cube2Vertex[min_ind][0];
                double line_y = slope_y*j+cube2Vertex[min_ind][1];
                double line_z = slope_z*j+cube2Vertex[min_ind][2];
                bool check1 = -ds * AASigma - line_x < -1e-5; //-x_cube1 < line_x, -1e-5 for numeric error precision
                bool check2 = line_x - ds * AASigma < -1e-5; // line_x < x_cube1, -1e-5 for numeric error precision
                bool check3 = y_coord11[0] * AASigma - line_y < -1e-5; //-y_cube1 < line_y, -1e-5 for numeric error precision
                bool check4 = line_y - y_coord11[2] * AASigma < -1e-5; //line_y < y_cube1, -1e-5 for numeric error precision
                bool check5 = z_coord11[0] * AASigma - line_z < -1e-5; //-z_cube1 < line_z, -1e-5 for numeric error precision
                bool check6 = line_z - z_coord11[2] * AASigma < -1e-5; // line_z < z_cube1, -1e-5 for numeric error precision
                if (size == 4) {
                   check1 = -ds * AASigma - line_x < -1e-5; //-x_cube1 < line_x, -1e-5 for numeric error precision
                   check2 = line_x - ds * AASigma < -1e-5; // line_x < x_cube1, -1e-5 for numeric error precision
                   check3 = y_coord11[0] * AASigma + (ds * AASigma - line_x) / sqrt(3) - line_y < -1e-5; //-y_cube1 < line_y
                   check4 = line_y - y_coord11[2] * AASigma + (ds * AASigma - line_x) / sqrt(3) < -1e-5; //line_y < y_cube1
                   check5 = z_coord11[0] * AASigma - line_z < -1e-5; //-z_cube1 < line_z, -1e-5 for numeric error precision
                   check6 = line_z - z_coord11[2] * AASigma < -1e-5; // line_z < z_cube1, -1e-5 for numeric error precision
                }
                if (check1 && check2 && check3 && check4 && check5 && check6) {
                   //cout << "Edge intersect overlap" << endl;
                   return true;
                }
             }
          }


        /*cout << "Vertex 1 final" << endl;
        for (int i = 0; i < 8; ++i) {
           for (int j = 0; j < 3; ++j) {
              std::cout << cube1Vertex[i][j] << " ";
           }
           std::cout << std::endl;
        }
        cout << "Vertex 2 final" << endl;
        for (int i = 0; i < 8; ++i) {
           for (int j = 0; j < 3; ++j) {
              std::cout << cube2Vertex[i][j] << " ";
           }
           std::cout << std::endl;
        }*/


          for (int i = 0; i < 3; ++i) {
              Vector3 axis_vec = {0.0, 0.0, 0.0};
              axis_vec[i] = 1.0;
              double max_proj1 = 0.0;
              double max_proj2 = 0.0;
              double min_proj1 = 1000000000.0; //very big number to initialize
              double min_proj2 = 1000000000.0; //very big number to initialize
              for (int j = 0; j < 6; ++j) {
                  double proj1 = cube1Vertex[j][0]*axis_vec[0]+cube1Vertex[j][1]*axis_vec[1]+cube1Vertex[j][2]*axis_vec[2];
                  double proj2 = cube2Vertex[j][0]*axis_vec[0]+cube2Vertex[j][1]*axis_vec[1]+cube2Vertex[j][2]*axis_vec[2];
                  if (proj1 < min_proj1) {
                     min_proj1 = proj1;
                  }
                  if (proj2 < min_proj2) {
                     min_proj2 = proj2;
                  }
                  if (proj1 > max_proj1) {
                     max_proj1 = proj1;
                  }
                  if (proj2 > max_proj2) {
                     max_proj2 = proj2;
                  }
              }
              //cout << axis_vec[0] << " " << axis_vec[1] << " " << axis_vec[2] << endl;
              //cout << min_proj1 << " " << min_proj2 << " " << max_proj1 << " " << max_proj2 << endl;
              bool check1 = max_proj1 - min_proj2 < 1e-5; // max_proj1 < min_proj2, 1e-5 for numeric error precision
              bool check2 = max_proj2 - min_proj1 < 1e-5; // max_proj2 < min_proj1, 1e-5 for numeric error precision
              if (check1 || check2) {
                 //cout << "No Projection overlap"<<endl;
                 return false;
              }
          }

          for (int i = 0; i < 3; ++i) {
              for (int j = 0; j < 3; ++j) {
                  if (i != j) {
                     Vector3 axis_vec = {0.0, 0.0, 0.0};
                     axis_vec[i] = 1.0;
                     axis_vec[j] = 1.0;
                     double max_proj1 = 0.0;
                     double max_proj2 = 0.0;
                     double min_proj1 = 1000000000.0; //very big number to initialize
                     double min_proj2 = 1000000000.0; //very big number to initialize
                     for (int j = 0; j < 6; ++j) {
                         double proj1 = cube1Vertex[j][0]*axis_vec[0]+cube1Vertex[j][1]*axis_vec[1]+cube1Vertex[j][2]*axis_vec[2];
                         double proj2 = cube2Vertex[j][0]*axis_vec[0]+cube2Vertex[j][1]*axis_vec[1]+cube2Vertex[j][2]*axis_vec[2];
                         if (proj1 < min_proj1) {
                            min_proj1 = proj1;
                         }
                         if (proj2 < min_proj2) {
                            min_proj2 = proj2;
                         }
                         if (proj1 > max_proj1) {
                            max_proj1 = proj1;
                         }
                         if (proj2 > max_proj2) {
                            max_proj2 = proj2;
                         }
                     }
                     //cout << axis_vec[0] << " " << axis_vec[1] << " " << axis_vec[2] << endl;
                     //cout << min_proj1 << " " << min_proj2 << " " << max_proj1 << " " << max_proj2 << endl;
                     bool check1 = max_proj1 - min_proj2 < 1e-5; // max_proj1 < min_proj2, 1e-5 for numeric error precision
                     bool check2 = max_proj2 - min_proj1 < 1e-5; // max_proj2 < min_proj1, 1e-5 for numeric error precision
                     if (check1 || check2) {
                        //cout << "No Projection 2 overlap"<<endl;
                        return false;
                     }
                  }
              }
          }
      }
      return false;
    }

    return true;
}

// Initial Configuration of Particles starting from crystalline phase
void InitialConfigurationCuboid(int NumberOfParticles, double BoxLength, double AASigma, double aspectratio, string shape,
                          unordered_map<double, double>& SizeRatioDictionary,
                          vector<double>& Sizes, vector<double>& Rx, vector<double>& Ry, vector<double>& Rz,
                          double (*VectorX)[3], double (*VectorY)[3], double (*VectorZ)[3],
                          int& restartStep, ofstream& LAMMPSTrajectoryFile, ofstream& EnergyFile, ofstream& TimeFile) {
/*
    Initializes the configuration of particles starting from a crystalline phase in a cuboid-shaped box.

    The function initializes the particle system by distributing particles of varying sizes based on a given size ratio dictionary.
    The positions of particles are set in a regular grid structure, ensuring that no overlaps occur initially. The particles are
    oriented to simulate the crystalline arrangement.

    Inputs:
      NumberOfParticles       - Total number of particles in the system
      BoxLength               - Length of the periodic box in each direction
      AASigma                 - Atom diameter
      aspectratio             - Aspect ratio for the particle shapes (height vs side length)
      shape                   - Type of particle shape (e.g., "cuboid")
      SizeRatioDictionary     - Dictionary containing the particle size ratios and their relative frequencies
      Sizes                   - A vector storing the sizes of the particles
      Rx, Ry, Rz              - Vectors for storing the x, y, and z positions of the particles
      VectorX, VectorY, VectorZ -  Orientation vectors representing the orientations of the particles
      restartStep             - Step number for restarting the simulation (if applicable)
      LAMMPSTrajectoryFile    - Output file for LAMMPS trajectory data
      EnergyFile              - Output file for energy data
      TimeFile                - Output file for time data

    Output:
      None (the function modifies the `Sizes`, `Rx`, `Ry`, `Rz`, `VectorX`, `VectorY`, and `VectorZ` vectors in place)*/

    srand(2);
    vector<double> keys_list;
    vector<double> values_list;

    for (const auto& entry : SizeRatioDictionary) {
        keys_list.push_back(entry.first);
        values_list.push_back(entry.second);
    }

    double valuesSum = accumulate(values_list.begin(), values_list.end(), 0.0);
    int ParticleID = 0;
    vector<int> CumulativeNumberOfParticlesList;

    for (double size : keys_list) {
        int NumberOfCurrentSize = round(SizeRatioDictionary[size] / valuesSum * NumberOfParticles);
        CumulativeNumberOfParticlesList.push_back(NumberOfCurrentSize);
    }

    for (size_t i = 0; i < CumulativeNumberOfParticlesList.size(); ++i) {
        if (ParticleID < NumberOfParticles) {
            int NumberOfCurrentSize = CumulativeNumberOfParticlesList[i];
            for (int j = 0; j < NumberOfCurrentSize; ++j) {
                Sizes[ParticleID] = keys_list[i];
                ParticleID++;
            }
        } else {
            break; // Exit the loop if ParticleID exceeds NumberOfParticles
        }
    }

    double BoxLengthHalf = BoxLength / 2.0;
    double Length = Sizes[0] * AASigma;
    double Height = Length * aspectratio;

    double FirstPositionX = 0.5*Length/BoxLength;
    double FirstPositionY = 0.5*Length/BoxLength;
    double FirstPositionZ = 0.5*Height/BoxLength;
    double FirstCubeAlpha = 0.0;
    double FirstCubeBeta = 0.0;
    double FirstCubeGamma = 0.0;

    for (int i = 0; i < NumberOfParticles; ++i) {
        VectorX[i][0] = 1.0;
        VectorX[i][1] = 0.0;
        VectorX[i][2] = 0.0;
        VectorY[i][0] = 0.0;
        VectorY[i][1] = 1.0;
        VectorY[i][2] = 0.0;
        VectorZ[i][0] = 0.0;
        VectorZ[i][1] = 0.0;
        VectorZ[i][2] = 1.0;
    }

    Vector3 RotationCenter_To_CubeCentroid;
    Vector3 FirstCubeVectors[3];
    Vector3 CubeVectors[3];

    ParticleRotation({FirstPositionX, FirstPositionY, FirstPositionZ},
             VectorX[0],VectorY[0],VectorZ[0],
             {FirstPositionX, FirstPositionY, FirstPositionZ},
             {FirstCubeAlpha, FirstCubeBeta, FirstCubeGamma},
             RotationCenter_To_CubeCentroid,
             FirstCubeVectors);

    for (int i = 0; i < 3; ++i) {
        VectorX[0][i] = FirstCubeVectors[0][i];
        VectorY[0][i] = FirstCubeVectors[1][i];
        VectorZ[0][i] = FirstCubeVectors[2][i];
    }

    Rx[0] = BoxLength * FirstPositionX - BoxLengthHalf;
    Ry[0] = BoxLength * FirstPositionY - BoxLengthHalf;
    Rz[0] = BoxLength * FirstPositionZ - BoxLengthHalf;

    /*int layer = floor(BoxLength / Height);
    cout << "NumberOfParticles " << NumberOfParticles << endl;
    cout <<"layer " << layer << endl;
    int avg_par_layer = ceil(static_cast<double>(NumberOfParticles) / layer);
    cout << "avg_par_layer " << avg_par_layer << endl;
    double z_gap = (BoxLength - layer * Height) / layer;
    cout << "z_gap " << z_gap << endl;*/


    for (int i = 1; i < NumberOfParticles; ++i) {
        std::cout << "Insertion Of Molecule " << i + 1 << " Successful\n";
        double ParticleSideLength1 = Sizes[i] * AASigma;
        int numy = floor(BoxLength/Length);
        int numysq = numy * numy;
        double posxy = BoxLength / numy;
        double PositionX = Rx[0] + posxy * (i % numy);
        double PositionY = Ry[0] + posxy * floor((i % numysq) / numy);
        double PositionZ = Rz[0] + ParticleSideLength1 * aspectratio * floor(i / numysq);

        bool Repeat = true;
        while (Repeat) {
            Repeat = false;
            double Alpha = 0.0;
            double Beta = 0.0;
            double Gamma = 0.0;

            ParticleRotation({PositionX, PositionY, PositionZ},
                 VectorX[i], VectorY[i], VectorZ[i],
                 {PositionX, PositionY, PositionZ},
                 {Alpha, Beta, Gamma},
                 RotationCenter_To_CubeCentroid,
                 CubeVectors);

            for (int K = 0; K < i; ++K) {
                //cout <<"K " << K << endl;
                //cout << "Rx " << Rx[K] << " " << Ry[K] << " " << Rz[K] << endl;
                //cout << "Position " << PositionX << " " << PositionY << " " << PositionZ << endl;
                //cout << "index " << index << endl;
                //cout << "selectedPositions[index] " << selectedPositions[index].first << " " << selectedPositions[index].second << endl;
                //cout << "selectedPositions[K] " << selectedPositions[K].first << " " << selectedPositions[K].second << endl;
                double ParticleSideLength2 = Sizes[K] * AASigma;
                double DistanceFromAnOldAtomX = PositionX - Rx[K];
                double DistanceFromAnOldAtomY = PositionY - Ry[K];
                double DistanceFromAnOldAtomZ = PositionZ - Rz[K];

                DistanceFromAnOldAtomX -= BoxLength * round(DistanceFromAnOldAtomX / BoxLength);
                DistanceFromAnOldAtomY -= BoxLength * round(DistanceFromAnOldAtomY / BoxLength);
                DistanceFromAnOldAtomZ -= BoxLength * round(DistanceFromAnOldAtomZ / BoxLength);

                Vector3 ParticleiCentroid = {0.0, 0.0, 0.0};
                Vector3 ParticlejCentroid = {DistanceFromAnOldAtomX, DistanceFromAnOldAtomY, DistanceFromAnOldAtomZ};

                double DistanceSquare = pow(DistanceFromAnOldAtomX, 2)+pow(DistanceFromAnOldAtomY, 2)+pow(DistanceFromAnOldAtomZ, 2);
                double maxDistanceSquare = 2.0*pow(ParticleSideLength2,2)+pow(ParticleSideLength2*aspectratio,2);
                if (DistanceSquare > maxDistanceSquare) {
                    Repeat = false;
                }
                else if (CheckOverlapForTwoParticles(0,1,BoxLength,AASigma,ParticleSideLength1,ParticleSideLength2,ParticleSideLength1*aspectratio,ParticleSideLength2*aspectratio,ParticleiCentroid,ParticlejCentroid,VectorX[K],VectorY[K],VectorZ[K],CubeVectors[0],CubeVectors[1],CubeVectors[2],shape)) {
                    Repeat = true;
                    break;
                }
            }
            //cout << "Repeat " << Repeat << endl;

            for (int j = 0; j < 3; ++j) {
                VectorX[i][j] = CubeVectors[0][j];
                VectorY[i][j] = CubeVectors[1][j];
                VectorZ[i][j] = CubeVectors[2][j];
            }

            Rx[i] = PositionX;
            Ry[i] = PositionY;
            Rz[i] = PositionZ;
        }
    }
    std::cout << "Start from random configuration\n";
}

void InitialConfigurationTriangularPrism(int NumberOfParticles, double BoxLength, double AASigma, double aspectratio, string shape,
                          unordered_map<double, double>& SizeRatioDictionary,
                          vector<double>& Sizes, vector<double>& Rx, vector<double>& Ry, vector<double>& Rz,
                          double (*VectorX)[3], double (*VectorY)[3], double (*VectorZ)[3],
                          int& restartStep, ofstream& LAMMPSTrajectoryFile, ofstream& EnergyFile, ofstream& TimeFile) {
/*
    Initializes the configuration of particles starting from a crystalline phase in a triangular prism-shaped box.

    The function initializes the particle system by distributing particles of varying sizes based on a given size ratio dictionary. 
    The positions of particles are set in a regular grid structure, ensuring that no overlaps occur initially. The particles are 
    oriented to simulate the crystalline arrangement within a triangular prism geometry.

    Inputs:
      NumberOfParticles       - Total number of particles in the system
      BoxLength               - Length of the periodic box in each direction
      AASigma                 - Atom diameter
      aspectratio             - Aspect ratio for the particle shapes (height vs side length)
      shape                   - Type of particle shape (e.g., "triangular prism")
      SizeRatioDictionary     - Dictionary containing the particle size ratios and their relative frequencies
      Sizes                   - A vector storing the sizes of the particles
      Rx, Ry, Rz              - Vectors for storing the x, y, and z positions of the particles
      VectorX, VectorY, VectorZ - Orientation vectors representing the orientations of the particles
      restartStep             - Step number for restarting the simulation (if applicable)
      LAMMPSTrajectoryFile    - Output file for LAMMPS trajectory data
      EnergyFile              - Output file for energy data
      TimeFile                - Output file for time data

    Output:
      None (the function modifies the `Sizes`, `Rx`, `Ry`, `Rz`, `VectorX`, `VectorY`, and `VectorZ` vectors in place) */

    srand(2);
    vector<double> keys_list;
    vector<double> values_list;

    for (const auto& entry : SizeRatioDictionary) {
        keys_list.push_back(entry.first);
        values_list.push_back(entry.second);
    }

    double valuesSum = accumulate(values_list.begin(), values_list.end(), 0.0);
    int ParticleID = 0;
    vector<int> CumulativeNumberOfParticlesList;

    for (double size : keys_list) {
        int NumberOfCurrentSize = round(SizeRatioDictionary[size] / valuesSum * NumberOfParticles);
        CumulativeNumberOfParticlesList.push_back(NumberOfCurrentSize);
    }

    for (size_t i = 0; i < CumulativeNumberOfParticlesList.size(); ++i) {
        if (ParticleID < NumberOfParticles) {
            int NumberOfCurrentSize = CumulativeNumberOfParticlesList[i];
            for (int j = 0; j < NumberOfCurrentSize; ++j) {
                Sizes[ParticleID] = keys_list[i];
                ParticleID++;
            }
        } else {
            break; // Exit the loop if ParticleID exceeds NumberOfParticles
        }
    }

    double BoxLengthHalf = BoxLength / 2.0;
    double Length = Sizes[0] * AASigma;
    double Height = Length * aspectratio;

    double FirstPositionX = (BoxLength/12.0+AASigma/4.0+(Length-AASigma)*sqrt(3)/12.0)/BoxLength;
    double FirstPositionY = 0.0;
    double FirstPositionZ = 0.5*Height/BoxLength;
    double FirstCubeAlpha = 0.0;
    double FirstCubeBeta = 0.0;
    double FirstCubeGamma = 0.0;

    for (int i = 0; i < NumberOfParticles; ++i) {
        VectorX[i][0] = 1.0;
        VectorX[i][1] = 0.0;
        VectorX[i][2] = 0.0;
        VectorY[i][0] = 0.0;
        VectorY[i][1] = 1.0;
        VectorY[i][2] = 0.0;
        VectorZ[i][0] = 0.0;
        VectorZ[i][1] = 0.0;
        VectorZ[i][2] = 1.0;
    }

    Vector3 RotationCenter_To_CubeCentroid;
    Vector3 FirstCubeVectors[3];
    Vector3 CubeVectors[3];

    ParticleRotation({FirstPositionX, FirstPositionY, FirstPositionZ},
             VectorX[0],VectorY[0],VectorZ[0],
             {FirstPositionX, FirstPositionY, FirstPositionZ},
             {FirstCubeAlpha, FirstCubeBeta, FirstCubeGamma},
             RotationCenter_To_CubeCentroid,
             FirstCubeVectors);

    for (int i = 0; i < 3; ++i) {
        VectorX[0][i] = FirstCubeVectors[0][i];
        VectorY[0][i] = FirstCubeVectors[1][i];
        VectorZ[0][i] = FirstCubeVectors[2][i];
    }

    Rx[0] = BoxLength * FirstPositionX - BoxLengthHalf;
    Ry[0] = BoxLength * FirstPositionY - BoxLengthHalf;
    Rz[0] = BoxLength * FirstPositionZ - BoxLengthHalf;

    /*int layer = floor(BoxLength / Height);
    cout << "NumberOfParticles " << NumberOfParticles << endl;
    cout <<"layer " << layer << endl;
    int avg_par_layer = ceil(static_cast<double>(NumberOfParticles) / layer);
    cout << "avg_par_layer " << avg_par_layer << endl;
    double z_gap = (BoxLength - layer * Height) / layer;
    cout << "z_gap " << z_gap << endl;*/

    double PositionX = Rx[0];
    for (int i = 1; i < NumberOfParticles; ++i) {
        std::cout << "Insertion Of Molecule " << i + 1 << " Successful\n";
        double CubeSideLength1 = Sizes[i] * AASigma;
        int numy = floor(BoxLength/Length);
        int doublenumy = 2 * numy;
        int fournumy = 4 * numy;
        int totnumy = floor(NumberOfParticles / numy);
        int numx = floor(totnumy / numy);

        double gapy = (BoxLength - numy * Length) / numy;
        double xleft = BoxLength / numx - (Length-AASigma) * sqrt(3) / 12.0 + 0.25 * AASigma;
        double xright = BoxLength * 4.0 / numx - (Length - AASigma) * sqrt(3) / 3.0 - AASigma - xleft;
        int i_remain = i % totnumy;
        double Alpha = 0.0;
        double Beta = 0.0;
        double Gamma = 0.0;

        if (i_remain % fournumy == 0) {
           PositionX = Rx[0] + (BoxLength * 4.0 / numx) * floor(i_remain / fournumy);
        }

        if (i_remain % fournumy < doublenumy && i_remain % fournumy > 0) {
           if (i_remain % 2 == 1) {
               PositionX = PositionX - xleft;
               Gamma = 180.0;
           } else {
               PositionX = PositionX + xleft;
               Gamma = 0.0;
           }
        } else if (i_remain % fournumy >= doublenumy) {
          if (i_remain % fournumy == doublenumy) {
               PositionX = PositionX + xleft + xright;
               Gamma = 180.0;
          }

          if (i_remain % 2 == 1) {
             PositionX = PositionX + xleft;
             Gamma = 0.0;
          } else {
             PositionX = PositionX - xleft;
             Gamma = 180.0;
          }
        }

        double PositionY = Ry[0] + 0.5 * (Length + gapy) * (i_remain % doublenumy);
        double PositionZ = Rz[0] + Length * aspectratio * floor(i / totnumy);

        bool Repeat = true;
        while (Repeat) {
            Repeat = false;

            ParticleRotation({PositionX, PositionY, PositionZ},
                 VectorX[i], VectorY[i], VectorZ[i],
                 {PositionX, PositionY, PositionZ},
                 {Alpha, Beta, Gamma},
                 RotationCenter_To_CubeCentroid,
                 CubeVectors);

            for (int K = 0; K < i; ++K) {
                //cout <<"K " << K << endl;
                //cout << "Rx " << Rx[K] << " " << Ry[K] << " " << Rz[K] << endl;
                //cout << "Position " << PositionX << " " << PositionY << " " << PositionZ << endl;
                //cout << "index " << index << endl;
                //cout << "selectedPositions[index] " << selectedPositions[index].first << " " << selectedPositions[index].second << endl;
                //cout << "selectedPositions[K] " << selectedPositions[K].first << " " << selectedPositions[K].second << endl;
                double CubeSideLength2 = Sizes[K] * AASigma;
                double DistanceFromAnOldAtomX = PositionX - Rx[K];
                double DistanceFromAnOldAtomY = PositionY - Ry[K];
                double DistanceFromAnOldAtomZ = PositionZ - Rz[K];

                DistanceFromAnOldAtomX -= BoxLength * round(DistanceFromAnOldAtomX / BoxLength);
                DistanceFromAnOldAtomY -= BoxLength * round(DistanceFromAnOldAtomY / BoxLength);
                DistanceFromAnOldAtomZ -= BoxLength * round(DistanceFromAnOldAtomZ / BoxLength);

                Vector3 ParticleiCentroid = {0.0, 0.0, 0.0};
                Vector3 ParticlejCentroid = {DistanceFromAnOldAtomX, DistanceFromAnOldAtomY, DistanceFromAnOldAtomZ};

                double DistanceSquare = pow(DistanceFromAnOldAtomX, 2)+pow(DistanceFromAnOldAtomY, 2)+pow(DistanceFromAnOldAtomZ, 2);
                double maxDistanceSquare = 2.0*pow(CubeSideLength2,2)+pow(CubeSideLength2*aspectratio,2);
                //cout << "here" << endl;
                if (DistanceSquare > maxDistanceSquare) {
                    Repeat = false;
                }
                //else if (CheckOverlapForTwoCubes(0,1,BoxLength,AASigma,CubeSideLength1,CubeSideLength2,CubeSideLength1*aspectratio,CubeSideLength2*aspectratio,ParticleiCentroid,ParticlejCentroid,VectorX[K],VectorY[K],VectorZ[K],CubeVectors[0],CubeVectors[1],CubeVectors[2])) {
                //    Repeat = true;
                //    break;
                //}
            }
            //cout << "Repeat " << Repeat << endl;

            for (int j = 0; j < 3; ++j) {
                VectorX[i][j] = CubeVectors[0][j];
                VectorY[i][j] = CubeVectors[1][j];
                VectorZ[i][j] = CubeVectors[2][j];
            }

            Rx[i] = PositionX;
            Ry[i] = PositionY;
            Rz[i] = PositionZ;
        }
    }
    std::cout << "Start from random configuration\n";
}

// Volume expansion: rescaling of the center-of-mass coordinates of the NPs relative to the center of the simulation box after each volume change.
void VolumeMove(int NumberOfParticles, double BoxLength, double AASigma, double aspectratio, double volrate,
                          unordered_map<double, double>& SizeRatioDictionary, string shape,
                          vector<double>& Sizes, vector<double>& Rx, vector<double>& Ry, vector<double>& Rz,
                          double (*VectorX)[3], double (*VectorY)[3], double (*VectorZ)[3]) {
/*
    Rescales the positions of particles after a volume change, effectively simulating volume expansion or compression in a
    simulation box.

    The function applies a volume expansion (or compression) by rescaling the center-of-mass coordinates of the particles.
    The particles are relocated proportionally to the current expansion or compression rate (`volrate`) relative to the
    box size. Additionally, the function checks and prevents particle overlaps during this adjustment by evaluating
    distances between all pairs of particles.

    Inputs:
      NumberOfParticles       - Total number of particles in the system
      BoxLength               - Length of the simulation box in each direction
      AASigma                 - Atom diameter
      aspectratio             - Aspect ratio of the particles (height vs side length)
      volrate                 - Volume scaling factor, i.e., the expansion or compression rate
      SizeRatioDictionary     - Dictionary of particle sizes and their relative ratios
      shape                   - Shape of the particles (e.g., "cuboid", "triangular prism", etc.)
      Sizes                   - A vector storing the sizes of the particles
      Rx, Ry, Rz              - Vectors storing the x, y, and z positions of the particles
      VectorX, VectorY, VectorZ - Orientation vectors representing the orientations of the particles

    Output:
      None (the function modifies the `Rx`, `Ry`, `Rz`, `VectorX`, `VectorY`, and `VectorZ` vectors in place)*/

    vector<double> keys_list;
    vector<double> values_list;

    for (const auto& entry : SizeRatioDictionary) {
        keys_list.push_back(entry.first);
        values_list.push_back(entry.second);
    }

    double valuesSum = accumulate(values_list.begin(), values_list.end(), 0.0);
    int ParticleID = 0;
    vector<int> CumulativeNumberOfParticlesList;

    for (double size : keys_list) {
        int NumberOfCurrentSize = round(SizeRatioDictionary[size] / valuesSum * NumberOfParticles);
        CumulativeNumberOfParticlesList.push_back(NumberOfCurrentSize);
    }

    for (size_t i = 0; i < CumulativeNumberOfParticlesList.size(); ++i) {
        if (ParticleID < NumberOfParticles) {
            int NumberOfCurrentSize = CumulativeNumberOfParticlesList[i];
            for (int j = 0; j < NumberOfCurrentSize; ++j) {
                Sizes[ParticleID] = keys_list[i];
                ParticleID++;
            }
        } else {
            break; // Exit the loop if ParticleID exceeds NumberOfParticles
        }
    }

    double BoxLengthHalf = BoxLength / 2.0;
    double Length = Sizes[0] * AASigma;
    double Height = Length * aspectratio;

    /*for (int i = 0; i < NumberOfParticles; ++i) {
        cout << "Particle " << i << " " << Rx[i] << " " << Ry[i] << " " << Rz[i] << endl;
    }*/

    double FirstPositionX = Rx[0] * volrate / BoxLength; // scaling their x position based on expansion rate
    double FirstPositionY = Ry[0] * volrate / BoxLength; // scaling their y position based on expansion rate
    double FirstPositionZ = Rz[0] * volrate / BoxLength; // scaling their z position based on expansion rate
    double FirstCubeAlpha = 0.0; // you don't change the orientation of these cubes, that is why you don't perform rotations on them
    double FirstCubeBeta = 0.0;
    double FirstCubeGamma = 0.0;

    Vector3 RotationCenter_To_CubeCentroid;
    Vector3 FirstCubeVectors[3];
    Vector3 CubeVectors[3];

    ParticleRotation({FirstPositionX, FirstPositionY, FirstPositionZ},
             VectorX[0],VectorY[0],VectorZ[0],
             {FirstPositionX, FirstPositionY, FirstPositionZ},
             {FirstCubeAlpha, FirstCubeBeta, FirstCubeGamma},
             RotationCenter_To_CubeCentroid,
             FirstCubeVectors);

    for (int i = 0; i < 3; ++i) {
        VectorX[0][i] = FirstCubeVectors[0][i];
        VectorY[0][i] = FirstCubeVectors[1][i];
        VectorZ[0][i] = FirstCubeVectors[2][i];
    }

    Rx[0] = BoxLength * FirstPositionX;
    Ry[0] = BoxLength * FirstPositionY;
    Rz[0] = BoxLength * FirstPositionZ;

    /*int layer = floor(BoxLength / Height);
    cout << "NumberOfParticles " << NumberOfParticles << endl;
    cout <<"layer " << layer << endl;
    int avg_par_layer = ceil(static_cast<double>(NumberOfParticles) / layer);
    cout << "avg_par_layer " << avg_par_layer << endl;
    double z_gap = (BoxLength - layer * Height) / layer;
    cout << "z_gap " << z_gap << endl;*/


    for (int i = 1; i < NumberOfParticles; ++i) {
        std::cout << "Insertion Of Molecule " << i + 1 << " Successful\n";
        double ParticleSideLength1 = Sizes[i] * AASigma;
        double PositionX = Rx[i] * volrate;
        double PositionY = Ry[i] * volrate;
        double PositionZ = Rz[i] * volrate;
        //cout << "New particle " << i << " " << PositionX << " " << PositionY << " " << PositionZ << endl;
        bool Repeat = true;
        int repeat_num = 0; // you don't change orientation of these cubes
        while (Repeat) {
            //cout << "repeat_num " << repeat_num << endl;
            Repeat = false;
            double Alpha = 0.0;
            double Beta = 0.0;
            double Gamma = 0.0;
            if (repeat_num > 0) { //this is actually unnecessary for volume expansion, as you don't need to change the cube orientations in volume expansion
               // but I included this in case one would like to do volume compression simulations, then this could be useful
               VectorX[i][0] = 1.0; VectorX[i][1] = 0.0; VectorX[i][2] = 0.0;
               VectorY[i][0] = 0.0; VectorY[i][1] = 1.0; VectorY[i][2] = 0.0;
               VectorZ[i][0] = 0.0; VectorZ[i][1] = 0.0; VectorZ[i][2] = 1.0;
               if (NumberOfParticles*pow(ParticleSideLength1,3)*aspectratio/pow(BoxLength,3) < 0.3) {
                  Alpha = 0.0;
                  Beta = 0.0;
                  Gamma = ((double)rand()/(RAND_MAX)) * 90;
               } else if (NumberOfParticles*pow(ParticleSideLength1,3)*aspectratio/pow(BoxLength,3) < 0.45) {
                  Alpha = (2.0 * ((double)rand()/(RAND_MAX)) - 1.0) * 90.0;
                  Beta = 0.0;//(2.0 * ((double)rand()/(RAND_MAX)) - 1.0) * 30.0;
                  Gamma = 0.0;
               } else {
                  Alpha = (2.0 * ((double)rand()/(RAND_MAX)) - 1.0) * 20.0;
                  Beta = 0.0;
                  Gamma = 0.0;
               }
            }

            if (repeat_num > 20) { // this is unnecessary for volume expansion simulations but could be useful for volume compression simulations
               PositionX = PositionX + ParticleSideLength1 * (2.0 * ((double)rand()/(RAND_MAX)) - 1.0);
               PositionX = PositionX - BoxLength * round(PositionX / BoxLength);
               PositionY = PositionY + ParticleSideLength1 * (2.0 * ((double)rand()/(RAND_MAX)) - 1.0);
               PositionY = PositionY - BoxLength * round(PositionY / BoxLength);
               PositionZ = PositionZ + ParticleSideLength1 * aspectratio * (2.0 * ((double)rand()/(RAND_MAX)) - 1.0);
               PositionZ = PositionZ - BoxLength * round(PositionZ / BoxLength);
            }

            ParticleRotation({PositionX, PositionY, PositionZ},
                 VectorX[i], VectorY[i], VectorZ[i],
                 {PositionX, PositionY, PositionZ},
                 {Alpha, Beta, Gamma},
                 RotationCenter_To_CubeCentroid,
                 CubeVectors);

            for (int K = 0; K < i; ++K) {
                //cout << "New particle " << i << " " << PositionX << " " << PositionY << " " << PositionZ << endl;
                //cout <<"K " << K << endl;
                //cout << "Rx " << Rx[K] << " " << Ry[K] << " " << Rz[K] << endl;
                //cout << "Position " << PositionX << " " << PositionY << " " << PositionZ << endl;
                //cout << "index " << index << endl;
                //cout << "selectedPositions[index] " << selectedPositions[index].first << " " << selectedPositions[index].second << endl;
                //cout << "selectedPositions[K] " << selectedPositions[K].first << " " << selectedPositions[K].second << endl;
                double ParticleSideLength2 = Sizes[K] * AASigma;
                double DistanceFromAnOldAtomX = PositionX - Rx[K];
                double DistanceFromAnOldAtomY = PositionY - Ry[K];
                double DistanceFromAnOldAtomZ = PositionZ - Rz[K];

                DistanceFromAnOldAtomX -= BoxLength * round(DistanceFromAnOldAtomX / BoxLength);
                DistanceFromAnOldAtomY -= BoxLength * round(DistanceFromAnOldAtomY / BoxLength);
                DistanceFromAnOldAtomZ -= BoxLength * round(DistanceFromAnOldAtomZ / BoxLength);

                Vector3 ParticleiCentroid = {0.0, 0.0, 0.0};
                Vector3 ParticlejCentroid = {DistanceFromAnOldAtomX, DistanceFromAnOldAtomY, DistanceFromAnOldAtomZ};

                double DistanceSquare = pow(DistanceFromAnOldAtomX, 2)+pow(DistanceFromAnOldAtomY, 2)+pow(DistanceFromAnOldAtomZ, 2);
                double maxDistanceSquare = 2.0*pow(ParticleSideLength2,2)+pow(ParticleSideLength2*aspectratio,2);
                if (DistanceSquare > maxDistanceSquare) {
                    Repeat = false;
                }
                else if (CheckOverlapForTwoParticles(0,1,BoxLength,AASigma,ParticleSideLength1,ParticleSideLength2,ParticleSideLength1*aspectratio,ParticleSideLength2*aspectratio,ParticleiCentroid,ParticlejCentroid,VectorX[K],VectorY[K],VectorZ[K],CubeVectors[0],CubeVectors[1],CubeVectors[2],shape)) {
                    //cout << "DistanceFromOldAtom " << DistanceFromAnOldAtomX << " " << DistanceFromAnOldAtomY << " " << DistanceFromAnOldAtomZ << endl;
                    //cout << "Alpha Beta Gamma " << Alpha << " " << Beta << " " << Gamma << endl;
                    Repeat = true;
                    repeat_num = repeat_num + 1;
                    break;
                }
            }
            //cout << "Repeat " << Repeat << endl;

            for (int j = 0; j < 3; ++j) {
                VectorX[i][j] = CubeVectors[0][j];
                VectorY[i][j] = CubeVectors[1][j];
                VectorZ[i][j] = CubeVectors[2][j];
            }

            Rx[i] = PositionX;
            Ry[i] = PositionY;
            Rz[i] = PositionZ;
        }
    }
    std::cout << "Start from random configuration\n";
}


void findAtomPositionsInCube(double CubeSize, double ParticleHeight, double Sigma, double COMX, double COMY, double COMZ, const Vector3 VectorX, const Vector3 VectorY, const Vector3 VectorZ, double atomList[][3]) {
/*
    Computes the atomic positions of a cuboid-shaped particle for use in energy computations (e.g., AACG energy).

    This function calculates the positions of atoms that make up a cuboid particle. It assumes the particle is made of atoms
    arranged in a grid within the cuboid, where the grid's spacing is determined by the particle's size (CubeSize and ParticleHeight).
    The resulting positions are transformed from the local coordinate system (centered at the particle's center of mass) to
    the global coordinate system based on the particle's orientation vectors.

    Inputs:
      CubeSize               - The length of the cuboid along the x and y axes (side length)
      ParticleHeight         - The height of the cuboid along the z-axis
      Sigma                  - Atom diameter
      COMX, COMY, COMZ       - The center of mass coordinates for the particle in the simulation box
      VectorX, VectorY, VectorZ - The orientation vectors of the particle, defining its rotation in space
      atomList               - A 2D array (numAtoms x 3) to store the computed atomic positions
    Output:
      atomList               - A 2D array containing the positions of all atoms in the particle after translation and rotation.
*/

    // Determine coordinates based on CubeSize
    double numCoordinates = (CubeSize - 1) / 2;
    double numCoordinates2 = (ParticleHeight - 1) / 2;

    // Fill the 2D array with coordinates
    int index = 0;
    for (double x = -numCoordinates; x < CubeSize/2; x += 1) {
        for (double y = -numCoordinates; y < CubeSize/2; y += 1) {
            for (double z = -numCoordinates2; z < ParticleHeight/2; z += 1) {
                atomList[index][0] = x * Sigma;
                atomList[index][1] = y * Sigma;
                atomList[index][2] = z * Sigma;
                ++index;
            }
        }
    }

    double x, y, z;
    // Translate and rotate the coordinates
    for (int i = 0; i < index; ++i) {
        x = atomList[i][0];
        y = atomList[i][1];
        z = atomList[i][2];
        atomList[i][0] = COMX + 1.0 * (x * VectorX[0] + y * VectorY[0] + z * VectorZ[0]);
        atomList[i][1] = COMY + 1.0 * (x * VectorX[1] + y * VectorY[1] + z * VectorZ[1]);
        atomList[i][2] = COMZ + 1.0 * (x * VectorX[2] + y * VectorY[2] + z * VectorZ[2]);
    }
}

void findAtomPositionsInTriangularPrism(double CubeSize, double CubeHeight, double Sigma, double COMX, double COMY, double COMZ, const Vector3 VectorX, const Vector3 VectorY, const Vector3 VectorZ, double atomList[][3]) {
/*
    Computes the atomic positions of a triangular prism-shaped particle for use in energy computations (e.g., AACG energy).

    This function calculates the positions of atoms that make up a triangular prism-shaped particle. The atoms are arranged
    in a grid, with the x and y coordinates forming a triangular lattice, and the z-coordinate representing the height of
    the prism. The resulting atomic positions are then transformed from the local coordinate system (centered at the particle's
    center of mass) to the global coordinate system based on the particle's orientation vectors.

    Inputs:
      CubeSize               - The base length of the triangular prism along the x and y axes
      CubeHeight             - The height of the triangular prism along the z-axis
      Sigma                  - The atom diameter
      COMX, COMY, COMZ       - The center of mass coordinates for the particle in the simulation box
      VectorX, VectorY, VectorZ - The orientation vectors of the particle, defining its rotation in space
      atomList               - A 2D array (numAtoms x 3) to store the computed atomic positions

    Output:
      atomList               - A 2D array containing the positions of all atoms in the particle after translation and rotation.
*/
    // Determine coordinates based on CubeSize
    double numCoordinates = (CubeSize - 1) / 2;
    double numCoordinates2 = (CubeHeight - 1) / 2;

    // Fill the 2D array with coordinates
    int index = 0;
    for (double x = 0.0; x < CubeSize; x += 1) {
        for (double y = 0.0; y <= x; y += 1) {
            for (double z = -numCoordinates2; z < CubeHeight/2; z += 1) {
                atomList[index][0] = x * Sigma * sqrt(3) / 2.0 - numCoordinates * 2.0 * Sigma / sqrt(3);
                atomList[index][1] = y * Sigma - x * Sigma / 2.0;
                atomList[index][2] = z * Sigma;
                ++index;
                //cout << "index " << index << endl;
            }
        }
    }
    //cout << "index " << index << endl;
    double x, y, z;
    // Translate and rotate the coordinates
    for (int i = 0; i < index; ++i) {
        x = atomList[i][0];
        y = atomList[i][1];
        z = atomList[i][2];
        atomList[i][0] = COMX + 1.0 * (x * VectorX[0] + y * VectorY[0] + z * VectorZ[0]);
        atomList[i][1] = COMY + 1.0 * (x * VectorX[1] + y * VectorY[1] + z * VectorZ[1]);
        atomList[i][2] = COMZ + 1.0 * (x * VectorX[2] + y * VectorY[2] + z * VectorZ[2]);
    }
}

double calculate_angle_between_vectors(Vector3 vector_a, Vector3 vector_b){
/*
    Calculates the angle between two 3D vectors.

    This function computes the angle between two vectors in 3D space using the dot product formula. The angle is returned in degrees and lies between 0 and 360 degrees. The calculation involves determining the dot product of the two vectors, followed by normalization to compute the cosine of the angle. The `acos` function is then used to find the angle in radians, which is subsequently converted to degrees.

    Inputs:
      vector_a               - The first 3D vector (array of 3 doubles)
      vector_b               - The second 3D vector (array of 3 doubles)

    Output:
      Returns the angle between the two vectors in degrees (ranging from 0 to 360 degrees).
*/
    double dot_product = vector_a[0] * vector_b[0] + vector_a[1] * vector_b[1] + vector_a[2] * vector_b[2];
    double norm_a = sqrt(vector_a[0] * vector_a[0] + vector_a[1] * vector_a[1] + vector_a[2] * vector_a[2]);
    double norm_b = sqrt(vector_b[0] * vector_b[0] + vector_b[1] * vector_b[1] + vector_b[2] * vector_b[2]);
    double angle = acos(dot_product / (norm_a * norm_b));
    double angleDegrees = fmod(abs(round(angle * (180.0 / M_PI))), 360.0);
    return angleDegrees;
}

double wrapAngle(double angle) {
/*
    Wraps an angle to the range of [-90, 270].

    This function takes an angle (in degrees) and modifies it so that it falls within the range of -90 to 270 degrees. The wrapping is done by adjusting angles outside this range to their equivalent within the desired interval. Specifically, angles are mapped to this range by adding or subtracting multiples of 90 degrees, ensuring the result falls within the target interval.

    Inputs:
      angle                  - The input angle in degrees

    Output:
      Returns the wrapped angle in degrees, constrained to the range [-90, 270].
*/
    if (angle <= 90 && angle >= 0) {
        // angle = angle;  // No change needed in this case
    } else if (angle > 90 && angle <= 180) {
        angle -= 90;
    } else if (angle > 180 && angle <= 270) {
        angle -= 180;
    } else if (angle > 270 && angle <= 360) {
        angle -= 270;
    } else if (angle > -90 && angle < 0) {
        angle = 90 + angle;
    } else if (angle > -180 && angle <= -90) {
        angle = 180 + angle;
    } else if (angle > -270 && angle <= -180) {
        angle = 270 + angle;
    } else if (angle > -360 && angle <= -270) {
        angle = 360 + angle;
    }
    return angle;
}

pair<double, bool> EnergyBetweenTwoParticles(string model, double AASigma, double CGSigma, double CGEpsilon, double cutoffCG, double Hamaker, double AtomDensity, int Size_i, int Size_j, double aspectratio, double BoxLength, const Vector3 ParticleiCentroid, const Vector3 ParticlejCentroid, const Vector3 VectorX_i, const Vector3 VectorY_i, const Vector3 VectorZ_i, const Vector3 VectorX_j, const Vector3 VectorY_j, const Vector3 VectorZ_j, string shape) {
/*
    Computes the energy between two particles based on their model and shape.

    This function calculates the interaction energy between two particles using either the "AACG" model or the "vdW" model. It depends on the shape of the particles (e.g., "cuboid" or "triangularprism").The energy calculation can be based on Lennard-Jones (LJ) potential or van der Waals (vdW) potential.

    Inputs:
      model                   - The model type to use for energy calculation ("AACG" or "vdW").
      AASigma                 - The sigma parameter for the AA interaction.
      CGSigma                 - The sigma parameter for the AACG (AA/coarse-grained) interaction.
      CGEpsilon               - The epsilon parameter for the AACG interaction.
      cutoffCG                - The cutoff distance for the AACG potential.
      Hamaker                 - The Hamaker constant for vdW energy calculation.
      AtomDensity             - The atom density used in vdW potential calculation.
      Size_i                  - The size of the first particle.
      Size_j                  - The size of the second particle.
      aspectratio             - The aspect ratio of the particle.
      BoxLength               - The length of the box used for simulation.
      ParticleiCentroid       - The centroid position of the first particle.
      ParticlejCentroid       - The centroid position of the second particle.
      VectorX_i, VectorY_i, VectorZ_i - The orientation vectors of the first particle.
      VectorX_j, VectorY_j, VectorZ_j - The orientation vectors of the second particle.
      shape                   - The shape of the particles ("cuboid" or "triangularprism").

    Output:
      Returns a pair of:
        CubeCubeEnergy         - The computed energy between the two particles.
        error_flag             - A boolean flag indicating whether an error occurred during the calculation.
*/

    double CubeCubeEnergy = 0.0;
    bool error_flag = false;
    if (model == "AACG") {
        // AA/CG LJ energy calculation
        if (shape == "cuboid") {
          double COMX, COMY, COMZ;
          COMX = ParticleiCentroid[0];
          COMY = ParticleiCentroid[1];
          COMZ = ParticleiCentroid[2];
          int NumberofRows = pow(Size_i * AASigma / CGSigma, 3) * aspectratio;
          int CGParticleSize1 = Size_i * AASigma / CGSigma;
          int CGHeight1 = CGParticleSize1 * aspectratio;
          double atomList1[NumberofRows][3];
          findAtomPositionsInCube(CGParticleSize1, CGHeight1, CGSigma, COMX, COMY, COMZ, VectorX_i, VectorY_i, VectorZ_i, atomList1);

          COMX = ParticlejCentroid[0];
          COMY = ParticlejCentroid[1];
          COMZ = ParticlejCentroid[2];
          NumberofRows = pow(Size_j * AASigma / CGSigma, 3) * aspectratio;
          int CGParticleSize2 = Size_j * AASigma / CGSigma;
          int CGHeight2 = CGParticleSize2 * aspectratio;
          double atomList2[NumberofRows][3];
          findAtomPositionsInCube(CGParticleSize2, CGHeight2, CGSigma, COMX, COMY, COMZ, VectorX_j, VectorY_j, VectorZ_j, atomList2);

          // Energy
          CubeCubeEnergy = LJPotentialBetweenTwoParticles(CGSigma, CGEpsilon, cutoffCG, atomList1, CGParticleSize1, atomList2, CGParticleSize2, aspectratio, BoxLength, shape);
        } else if (shape == "triangularprism") {
          double COMX, COMY, COMZ;
          COMX = ParticleiCentroid[0];
          COMY = ParticleiCentroid[1];
          COMZ = ParticleiCentroid[2];
          int NumberofRows = pow(Size_i * AASigma / CGSigma, 2) * aspectratio * (Size_i * AASigma / CGSigma + 1) / 2.0;
          //cout << "NumberofRows " << NumberofRows << endl;
          int CGCubeSize1 = Size_i * AASigma / CGSigma;
          int CGHeight1 = CGCubeSize1 * aspectratio;
          double atomList1[NumberofRows][3];
          findAtomPositionsInTriangularPrism(CGCubeSize1, CGHeight1, CGSigma, COMX, COMY, COMZ, VectorX_i, VectorY_i, VectorZ_i, atomList1);

          COMX = ParticlejCentroid[0];
          COMY = ParticlejCentroid[1];
          COMZ = ParticlejCentroid[2];
          NumberofRows = pow(Size_j * AASigma / CGSigma, 2) * aspectratio * (Size_j * AASigma / CGSigma + 1) / 2.0;
          int CGCubeSize2 = Size_j * AASigma / CGSigma;
          int CGHeight2 = CGCubeSize2 * aspectratio;
          double atomList2[NumberofRows][3];
          findAtomPositionsInTriangularPrism(CGCubeSize2, CGHeight2, CGSigma, COMX, COMY, COMZ, VectorX_j, VectorY_j, VectorZ_j, atomList2);

          // Energy
          CubeCubeEnergy = LJPotentialBetweenTwoParticles(CGSigma, CGEpsilon, cutoffCG, atomList1, CGCubeSize1, atomList2, CGCubeSize2, aspectratio, BoxLength, shape);
        }
    } else if (model == "vdW") {
        // vdW energy calculation
        int i = 0;
        int j = 1;

        double Particle1Centroid[3] = {ParticleiCentroid[0], ParticleiCentroid[1], ParticleiCentroid[2]};
        double Particle2Centroid[3] = {ParticlejCentroid[0], ParticlejCentroid[1], ParticlejCentroid[2]};
        Vector3 VectorX_1,VectorX_2,VectorY_1,VectorY_2,VectorZ_1,VectorZ_2;

        for (int k = 0; k < 3; ++k) {
           VectorX_1[k] = VectorX_i[k];
           VectorY_1[k] = VectorY_i[k];
           VectorZ_1[k] = VectorZ_i[k];
           VectorX_2[k] = VectorX_j[k];
           VectorY_2[k] = VectorY_j[k];
           VectorZ_2[k] = VectorZ_j[k];
        }
        if (shape == "cuboid") {
          double Particle2VerticesEdge[8][3];
          double y_coord11[4];
          double z_coord11[4];
          double ds;
          double lc = Size_i/2.0;
          double Particle2SideLength = Size_j * AASigma;
          double Particle2Height = Particle2SideLength * aspectratio;
          // Reorientation of two particles such that first particle will align along coordinate frame axes
          Reorientation_cuboid(i,j,Particle2SideLength,Particle2Height,BoxLength,Particle1Centroid,Particle2Centroid,VectorX_1,VectorY_1,VectorZ_1,VectorX_2,VectorY_2,VectorZ_2,Particle2VerticesEdge,y_coord11,z_coord11,ds,lc);

          // Parameters
          double catr1 = -7.75596;
          double atre1 = 3.4339;
          double catr2 = -4.65093;
          double atre2 = 2.7167;
          double crep = 3.85898;
          double repe = 10.66482;

          double Epsilon = (Hamaker * (pow(10, -19)) / (4 * M_PI * M_PI * pow(AtomDensity, 2) * pow(AASigma, 6) * 6.9477 * 2 * (pow(10, -21))));
          // Call function to calculate potential energy between two particles
          CubeCubeEnergy = calc_ver_cuboid(Particle2VerticesEdge, AASigma, catr1, atre1, catr2, atre2, crep, repe, Epsilon, lc, aspectratio, y_coord11, z_coord11, ds, error_flag);
          if (std::isinf(std::abs(CubeCubeEnergy))) { // if energy value is infinite, set this 0
              CubeCubeEnergy = 0.0;
          }
        }
        if (shape == "triangularprism") {
          double Particle2VerticesEdge[6][3];
          vector<double> y_coord11;
          vector<double> z_coord11;
          double ds;
          double lc = Size_i/2.0;
          double Cube2SideLength = Size_j * AASigma;
          double Cube2Height = Cube2SideLength * aspectratio;

          Reorientation_triangularprism(i,j,Cube2SideLength,Cube2Height,AASigma,BoxLength,Particle1Centroid,Particle2Centroid,VectorX_1,VectorY_1,VectorZ_1,VectorX_2,VectorY_2,VectorZ_2,Particle2VerticesEdge,y_coord11,z_coord11,ds,lc);
          // Parameters
          double catr1 = -7.75596;
          double atre1 = 3.4339;
          double catr2 = -4.65093;
          double atre2 = 2.7167;
          double crep = 3.85898;
          double repe = 10.66482;

        /*cout << "Reorientation" << endl;
        cout << "Particle2Centroid[0] " << Particle2Centroid[0] << " Particle2Centroid[1] " << Particle2Centroid[1] << " Particle2Centroid[2] " << Particle2Centroid[2] << endl;
        cout << "VectorX_1[0] " << VectorX_1[0] << " VectorX_1[1] " << VectorX_1[1] << " VectorX_1[2] " << VectorX_1[2] << endl;
        cout << "VectorY_1[0] " << VectorY_1[0] << " VectorY_1[1] " << VectorY_1[1] << " VectorY_1[2] " << VectorY_1[2] << endl;
        cout << "VectorZ_1[0] " << VectorZ_1[0] << " VectorZ_1[1] " << VectorZ_1[1] << " VectorZ_1[2] " << VectorZ_1[2] << endl;
        cout << "VectorX_2[0] " << VectorX_2[0] << " VectorX_2[1] " << VectorX_2[1] << " VectorX_2[2] " << VectorX_2[2] << endl;
        cout << "VectorY_2[0] " << VectorY_2[0] << " VectorY_2[1] " << VectorY_2[1] << " VectorY_2[2] " << VectorY_2[2] << endl;
        cout << "VectorZ_2[0] " << VectorZ_2[0] << " VectorZ_2[1] " << VectorZ_2[1] << " VectorZ_2[2] " << VectorZ_2[2] << endl;

        for (int i = 0; i < 6; ++i) {
           for (int j = 0; j < 3; ++j) {
              std::cout << Particle2VerticesEdge[i][j] << " ";
           }
           std::cout << std::endl;
        }*/



          double Epsilon = (Hamaker * (pow(10, -19)) / (4 * M_PI * M_PI * pow(AtomDensity, 2) * pow(AASigma, 6) * 6.9477 * 2 * (pow(10, -21))));
          // Call function to calculate potential energy between two particles
          CubeCubeEnergy = calc_ver_triangularprism(Particle2VerticesEdge, AASigma, catr1, atre1, catr2, atre2, crep, repe, Epsilon, lc, aspectratio, y_coord11, z_coord11, ds, error_flag);
          if (std::isinf(std::abs(CubeCubeEnergy))) { // if energy value is infinite, set this 0
              CubeCubeEnergy = 0.0;
          }
        }


    }

    return make_pair(CubeCubeEnergy, error_flag);
}

pair<double, bool> OneParticleEnergy(string model, int i, double Size_i, double aspectratio, double cutoffCG, double Hamaker, double AtomDensity, const double Rx_i, const double Ry_i, const double Rz_i, 
                                     vector<double>& Sizes, vector<double>& Rx, vector<double>& Ry, vector<double>& Rz,const double VectorX[][3], const double VectorY[][3], const double VectorZ[][3],
                                     int NumberOfParticles, double BoxLength, double CGEpsilon, double CGSigma, double AASigma, double globalminEnergy,
                                     string shape, vector<int>& clusterParticleList){
/*
    Calculates the total energy of one particle with all interacting particles within a cutoff distance.

    This function computes the potential energy of a single particle (specified by the index `i`) due to interactions with other particles. The interactions are evaluated based on the model selected (`AACG` or `vdW`), and the energy is calculated only for particles that lie within a cutoff distance. The function also checks if the particles overlap and computes the energy accordingly. If the particles do overlap, their potential energy is recalculated using the appropriate interaction model.

    Inputs:
      model                   - The model used for energy calculation ("AACG" or "vdW").
      i                        - The index of the particle for which the energy is calculated.
      Size_i                   - The size of the particle `i`.
      aspectratio              - The aspect ratio of the particle.
      cutoffCG                 - The cutoff distance for the CG potential.
      Hamaker                  - The Hamaker constant for vdW energy calculation.
      AtomDensity              - The atom density used in vdW potential calculation.
      Rx_i, Ry_i, Rz_i         - The position coordinates of the particle `i`.
      Sizes                    - A vector containing the sizes of all particles.
      Rx, Ry, Rz               - Vectors containing the positions of all particles.
      VectorX, VectorY, VectorZ - 3D vectors for the orientations of all particles.
      NumberOfParticles        - The total number of particles in the system.
      BoxLength                - The length of the box used in the simulation.
      CGEpsilon                - The epsilon parameter for the CG potential.
      CGSigma                  - The sigma parameter for the CG potential.
      AASigma                  - The sigma parameter for the AA/CG interaction.
      globalminEnergy          - The global minimum energy threshold for checking energy validity.
      shape                    - The shape of the particles ("cuboid" or "triangularprism").
      clusterParticleList      - A list of particle indices that are part of a cluster

    Output:
      Returns a pair:
        CurrentAtomTotalPotentialEnergy - The total energy of the particle `i` with respect to all interacting particles.
        overLapFlag                 - A boolean flag indicating whether the particle overlaps with another particle.
*/

   // Calculate energy of one cube
   bool overLapFlag = false;
   double CurrentAtomTotalPotentialEnergy = 0;
   double CurrentPairPotentialEnergy = 0;
   Vector3 VectorX_i,VectorY_i,VectorZ_i, VectorX_j, VectorY_j, VectorZ_j;

   for (int j = 0; j < NumberOfParticles; ++j) {
       if (j != i && find(clusterParticleList.begin(), clusterParticleList.end(), j) == clusterParticleList.end()) {
          Size_i = Sizes[i];
          //cout << i << " and " << j << " investigation" << endl;
          double Size_j = Sizes[j];
          for (int k = 0; k < 3; ++k) {
             VectorX_i[k] = VectorX[i][k];
             VectorY_i[k] = VectorY[i][k];
             VectorZ_i[k] = VectorZ[i][k];
          }
          VectorX_j[0] = VectorX[j][0];VectorX_j[1] = VectorX[j][1];VectorX_j[2] = VectorX[j][2];
          VectorY_j[0] = VectorY[j][0];VectorY_j[1] = VectorY[j][1];VectorY_j[2] = VectorY[j][2];
          VectorZ_j[0] = VectorZ[j][0];VectorZ_j[1] = VectorZ[j][1];VectorZ_j[2] = VectorZ[j][2];
          double Rx_j = Rx[j];
          double Rx_ij = Rx_i - Rx_j;
          double Ry_j = Ry[j];
          double Ry_ij = Ry_i - Ry_j;
          double Rz_j = Rz[j];
          double Rz_ij = Rz_i - Rz_j;
          /*cout << "i " << i << " and j " << j << endl;
          cout << "VectorX_i[0] " << VectorX_i[0] << " VectorX_i[1] " << VectorX_i[1] << " VectorX_i[2] " << VectorX_i[2] << endl;
          cout << "VectorY_i[0] " << VectorY_i[0] << " VectorY_i[1] " << VectorY_i[1] << " VectorY_i[2] " << VectorY_i[2] << endl;
          cout << "VectorZ_i[0] " << VectorZ_i[0] << " VectorZ_i[1] " << VectorZ_i[1] << " VectorZ_i[2] " << VectorZ_i[2] << endl;
          cout << "VectorX_j[0] " << VectorX_j[0] << " VectorX_j[1] " << VectorX_j[1] << " VectorX_j[2] " << VectorX_j[2] << endl;
          cout << "VectorY_j[0] " << VectorY_j[0] << " VectorY_j[1] " << VectorY_j[1] << " VectorY_j[2] " << VectorY_j[2] << endl;
          cout <<"VectorZ_j[0] " << VectorZ_j[0] << " VectorZ_j[1] " << VectorZ_j[1] << " VectorZ_j[2] " << VectorZ_j[2] << endl;
          cout << "Rxi " << Rx_i << " Ryi " << Ry_i << " Rzi " << Rz_i << endl;
          cout << "Rxj " << Rx_j << " Ryj " << Ry_j << " Rzj " << Rz_j << endl;*/
          Rx_ij = Rx_ij - BoxLength * round(Rx_ij / BoxLength);
          Ry_ij = Ry_ij - BoxLength * round(Ry_ij / BoxLength);
          Rz_ij = Rz_ij - BoxLength * round(Rz_ij / BoxLength);
          double RijSquare = Rx_ij * Rx_ij + Ry_ij * Ry_ij + Rz_ij * Rz_ij;
          Vector3 ParticleiCentroid = {0.0, 0.0, 0.0};
          Vector3 ParticlejCentroid = {-Rx_ij, -Ry_ij, -Rz_ij};
          double MaxSize = max(Size_i,Size_j);
          double CutOff = 2.5 * MaxSize * CGSigma * max(1.0,aspectratio); // cutoff distance for determining interacting particles
          double CutOffSquare = CutOff * CutOff;
          if (RijSquare < CutOffSquare) {
             double ParticleSideLength1 = Size_i*AASigma;
             double ParticleSideLength2 = Size_j*AASigma;
             double ParticleHeight1 = ParticleSideLength1*aspectratio;
             double ParticleHeight2 = ParticleSideLength2*aspectratio;

             if (CheckOverlapForTwoParticles(i,j,BoxLength,AASigma,ParticleSideLength1,ParticleSideLength2,ParticleHeight1,ParticleHeight2,ParticleiCentroid,ParticlejCentroid,VectorX_i,VectorY_i,VectorZ_i,VectorX_j,VectorY_j,VectorZ_j,shape)) {
                pair<double, bool> result = EnergyBetweenTwoParticles(model,AASigma,CGSigma,CGEpsilon,cutoffCG,Hamaker,AtomDensity,Size_i,Size_j,aspectratio,BoxLength,ParticleiCentroid,ParticlejCentroid,VectorX_i,VectorY_i,VectorZ_i,VectorX_j,VectorY_j,VectorZ_j,shape);
                CurrentPairPotentialEnergy = result.first;
                bool error_flag = result.second;
                overLapFlag = true;
                //cout << i << " and " << j << " overlap" << endl;
             } else {
                pair<double, bool> result = EnergyBetweenTwoParticles(model,AASigma,CGSigma,CGEpsilon,cutoffCG,Hamaker,AtomDensity,Size_i,Size_j,aspectratio,BoxLength,ParticleiCentroid,ParticlejCentroid,VectorX_i,VectorY_i,VectorZ_i,VectorX_j,VectorY_j,VectorZ_j,shape);
                CurrentPairPotentialEnergy = result.first;
                bool error_flag = result.second;
                //cout << "error_flag " << error_flag << endl;
                if (CurrentPairPotentialEnergy > 0 || abs(CurrentPairPotentialEnergy) > globalminEnergy || isnan(CurrentPairPotentialEnergy || error_flag)) { // as a measure in case our overlapping algorithm did not detect the overlap
                   overLapFlag = true;
                }
             }
             CurrentAtomTotalPotentialEnergy += CurrentPairPotentialEnergy;
          }
       }
   }
   // cout << "CurrentAtomTotalPotentialEnergy " << CurrentAtomTotalPotentialEnergy << endl;
   return make_pair(CurrentAtomTotalPotentialEnergy, overLapFlag);
}


double TotalEnergy(vector<double> Sizes, vector<double> Rx, vector<double> Ry, vector<double> Rz, double VectorX[][3],double VectorY[][3],double VectorZ[][3], int NumberOfParticles, double BoxLength, double aspectratio, double CGSigma, double CGEpsilon, double AASigma, string model, double cutoffCG, double Hamaker, double AtomDensity, string shape){
/*
    Calculates the total energy of the entire particle system.

    This function computes the total potential energy of a system of particles by iterating over all pairs of particles, checking their distances, and calculating the energy based on a selected interaction model (e.g., "AACG" or "vdW"). It only considers interactions between particles within a defined cutoff distance. The function also checks for overlaps between particles and calculates the corresponding energy based on the overlap. The total energy is then accumulated for all particle pairs.

    Inputs:
      Sizes                    - A vector containing the sizes of all particles in the system.
      Rx, Ry, Rz               - Vectors containing the positions of all particles.
      VectorX, VectorY, VectorZ - 3D orientation vectors for each particle.
      NumberOfParticles        - The total number of particles in the system.
      BoxLength                - The length of the simulation box.
      aspectratio              - The aspect ratio of the particles.
      CGSigma                  - The sigma parameter for the CG potential.
      CGEpsilon                - The epsilon parameter for the CG potential.
      AASigma                  - The sigma parameter for the AA/CG interaction.
      model                    - The interaction model used for energy calculation ("AACG" or "vdW").
      cutoffCG                 - The cutoff distance for the CG potential.
      Hamaker                  - The Hamaker constant for vdW energy calculation.
      AtomDensity              - The atom density used in vdW potential calculation.
      shape                    - The shape of the particles ("cuboid", "sphere", etc.).

    Output:
      Returns the total potential energy of the system as a double.
*/

   //Calculate total energy of the system
   Vector3 VectorX_i,VectorX_j,VectorY_i,VectorY_j,VectorZ_i,VectorZ_j;
   double CurrentPairPotentialEnergy = 0;
   double CurrentAtomTotalPotentialEnergy = 0;
   bool overLapFlag = false;

   for (int i = 0; i < NumberOfParticles; ++i) {
       for (int j = i+1; j < NumberOfParticles; ++j) {
           double Size_i = Sizes[i];
           double Size_j = Sizes[j];
           //cout << i << " and " << j << " investigated" << endl;
           VectorX_i[0] = VectorX[i][0];VectorX_i[1] = VectorX[i][1];VectorX_i[2] = VectorX[i][2];
           VectorY_i[0] = VectorY[i][0];VectorY_i[1] = VectorY[i][1];VectorY_i[2] = VectorY[i][2];
           VectorZ_i[0] = VectorZ[i][0];VectorZ_i[1] = VectorZ[i][1];VectorZ_i[2] = VectorZ[i][2];

           VectorX_j[0] = VectorX[j][0];VectorX_j[1] = VectorX[j][1];VectorX_j[2] = VectorX[j][2];
           VectorY_j[0] = VectorY[j][0];VectorY_j[1] = VectorY[j][1];VectorY_j[2] = VectorY[j][2];
           VectorZ_j[0] = VectorZ[j][0];VectorZ_j[1] = VectorZ[j][1];VectorZ_j[2] = VectorZ[j][2];

           double Rx_ij = Rx[i] - Rx[j];
           double Ry_ij = Ry[i] - Ry[j];
           double Rz_ij = Rz[i] - Rz[j];
           Rx_ij = Rx_ij - BoxLength * round(Rx_ij / BoxLength);
           Ry_ij = Ry_ij - BoxLength * round(Ry_ij / BoxLength);
           Rz_ij = Rz_ij - BoxLength * round(Rz_ij / BoxLength);
           //cout << "Total Energy" << endl;
           //cout << "i " << i << " and j " << j << endl;
           //cout << "VectorX_i[0] " << VectorX_i[0] << " VectorX_i[1] " << VectorX_i[1] << " VectorX_i[2] " << VectorX_i[2] << endl;
           //cout << "VectorY_i[0] " << VectorY_i[0] << " VectorY_i[1] " << VectorY_i[1] << " VectorY_i[2] " << VectorY_i[2] << endl;
           //cout << "VectorZ_i[0] " << VectorZ_i[0] << " VectorZ_i[1] " << VectorZ_i[1] << " VectorZ_i[2] " << VectorZ_i[2] << endl;
           //cout << "VectorX_j[0] " << VectorX_j[0] << " VectorX_j[1] " << VectorX_j[1] << " VectorX_j[2] " << VectorX_j[2] << endl;
           //cout << "VectorY_j[0] " << VectorY_j[0] << " VectorY_j[1] " << VectorY_j[1] << " VectorY_j[2] " << VectorY_j[2] << endl;
           //cout << "VectorZ_j[0] " << VectorZ_j[0] << " VectorZ_j[1] " << VectorZ_j[1] << " VectorZ_j[2] " << VectorZ_j[2] << endl;
           //cout << "Rxj " << Rx[j] << " Ryj " << Ry[j] << " Rzj " << Rz[j] << endl;
           //cout << "Rxi " << Rx[i] << " Ryi " << Ry[i] << " Rzi " << Rz[i] << endl;
           //cout << "Rxj " << Rx[j] << " Ryj " << Ry[j] << " Rzj " << Rz[j] << endl;
           double RijSquare = Rx_ij * Rx_ij + Ry_ij * Ry_ij + Rz_ij * Rz_ij;
           Vector3 ParticleiCentroid = {0.0, 0.0, 0.0};
           Vector3 ParticlejCentroid = {-Rx_ij, -Ry_ij, -Rz_ij};
           //cout << "Here totatl energy part" << endl;
           //cout << ParticlejCentroid[0] << " " << ParticlejCentroid[1] << " " << ParticlejCentroid[2] << endl;
           double MaxSize = max(Size_i,Size_j);
           double CutOff =  2.5 * MaxSize * CGSigma * max(1.0,aspectratio); // cutoff distance for determining interaction between two particles
           double CutOffSquare = CutOff * CutOff;
           if (RijSquare < CutOffSquare) {
              double ParticleSideLength1 = Size_i*AASigma;
              double ParticleSideLength2 = Size_j*AASigma;
              double ParticleHeight1 = ParticleSideLength1*aspectratio;
              double ParticleHeight2 = ParticleSideLength2*aspectratio;

              if (CheckOverlapForTwoParticles(i,j,BoxLength,AASigma,ParticleSideLength1,ParticleSideLength2,ParticleHeight1,ParticleHeight2,ParticleiCentroid,ParticlejCentroid,VectorX_i,VectorY_i,VectorZ_i,VectorX_j,VectorY_j,VectorZ_j,shape)) {
                 pair<double, bool> result = EnergyBetweenTwoParticles(model,AASigma,CGSigma,CGEpsilon,cutoffCG,Hamaker,AtomDensity,Size_i,Size_j,aspectratio,BoxLength,ParticleiCentroid,ParticlejCentroid,VectorX_i,VectorY_i,VectorZ_i,VectorX_j,VectorY_j,VectorZ_j,shape);
                 CurrentPairPotentialEnergy = result.first;
                 bool error_flag = result.second;
                 overLapFlag = true;
                 //cout << i << " and " << j << " overlap" << endl;
              } else {
                 //cout << ParticlejCentroid[0] << " " << ParticlejCentroid[1] << " " << ParticlejCentroid[2] << endl;
                 pair<double, bool> result = EnergyBetweenTwoParticles(model,AASigma,CGSigma,CGEpsilon,cutoffCG,Hamaker,AtomDensity,Size_i,Size_j,aspectratio,BoxLength,ParticleiCentroid,ParticlejCentroid,VectorX_i,VectorY_i,VectorZ_i,VectorX_j,VectorY_j,VectorZ_j,shape);
                 CurrentPairPotentialEnergy = result.first;
              }
              CurrentAtomTotalPotentialEnergy += CurrentPairPotentialEnergy;
              //cout << "CurrentPairPotentialEnergy " << CurrentPairPotentialEnergy << endl;
           }
       }
   }
   //cout << "CurrentAtomTotalPotentialEnergy" << endl;

   return CurrentAtomTotalPotentialEnergy;
}

void apply_periodic_boundary(Vector3 position, double BoxLength, Vector3 result){
/*
    Applies periodic boundary conditions to a particle's position.

    This function adjusts a particle's position to ensure it remains within the boundaries of a periodic simulation box.

    Inputs:
      position    - The current position of the particle (in 3D space).
      BoxLength   - The length of the simulation box (assuming a cubic box).

    Output:
      result      - The updated position of the particle after applying the periodic boundary conditions.
*/
   for (int i = 0; i < 3; ++i) {
       result[i] = position[i] - BoxLength * round(position[i]/ BoxLength);
   }
}

void calculate_cluster_centroid(vector<array<double, 3>> particles, double BoxLength, Vector3 centroid, double ClusterTemp[][3]){
/*
    Finds the coordinates of the centroid (center) of a cluster of particles.

    This function calculates the centroid of a cluster by first adjusting each particle's position with periodic boundary conditions, then computing the mean position of all particles in the cluster.

    Inputs:
      particles      - A vector of arrays representing the positions of the particles in the cluster.
      BoxLength      - The length of the simulation box (assuming a cubic box).

    Output:
      centroid       - The updated position of the centroid after computing the mean position of all particles, adjusted for periodic boundaries.
      ClusterTemp    - A temporary storage for the adjusted coordinates of particles in the cluster after applying periodic boundary conditions.
*/

   size_t numRows = particles.size();
   Vector3 result;
   Vector3 convertedElement;
   array<double, 3> element;

   for (size_t i = 0; i < numRows; ++i) {
       element = particles[i];
       copy(element.begin(), element.end(), convertedElement);
       apply_periodic_boundary(convertedElement, BoxLength, result);
       ClusterTemp[i][0] = result[0]; // temporary cluster after periodic boundary condition is applied
       ClusterTemp[i][1] = result[1];
       ClusterTemp[i][2] = result[2];
   }

   double X_sum = ClusterTemp[0][0];
   double Y_sum = ClusterTemp[0][1];
   double Z_sum = ClusterTemp[0][2];
   centroid[0] = X_sum;
   centroid[1] = Y_sum;
   centroid[2] = Z_sum;
   double X_dist, Y_dist, Z_dist;
   // find the mean of coordiantes of all particles in this cluster
   for (int i = 1; i < numRows; ++i) {
       X_dist = ClusterTemp[i][0]-centroid[0];
       Y_dist = ClusterTemp[i][1]-centroid[1];
       Z_dist = ClusterTemp[i][2]-centroid[2];
       X_dist = X_dist - BoxLength * round(X_dist / BoxLength);
       Y_dist = Y_dist - BoxLength * round(Y_dist / BoxLength);
       Z_dist = Z_dist - BoxLength * round(Z_dist / BoxLength);
       X_sum += centroid[0]+X_dist;
       Y_sum += centroid[1]+Y_dist;
       Z_sum += centroid[2]+Z_dist;
       ClusterTemp[i][0] = centroid[0]+X_dist;
       ClusterTemp[i][1] = centroid[1]+Y_dist;
       ClusterTemp[i][2] = centroid[2]+Z_dist;
   }

   double Xaverage = X_sum / numRows;
   double Yaverage = Y_sum / numRows;
   double Zaverage = Z_sum / numRows;

   centroid[0] = Xaverage;
   centroid[1] = Yaverage;
   centroid[2] = Zaverage;
   apply_periodic_boundary(centroid, BoxLength, centroid);
}

// Function to compute quaternion from rotation matrix
std::array<double, 4> rotationMatrixToQuaternion(const double R[3][3]) {
/*
    Computes the quaternion representation from a 3x3 rotation matrix.

    This function converts a 3x3 rotation matrix into a unit quaternion that represents the same rotation. The conversion is done by first computing the components of the quaternion based on the elements of the matrix, and then determining the correct sign and order of the quaternion components by examining the maximum value.

    The resulting quaternion can be used to represent the same rotation as the matrix but in a more compact form.

    Inputs:
      R             - A 3x3 rotation matrix representing a rotation.

    Output:
      Returns a quaternion (represented as an array of 4 doubles) that corresponds to the same rotation as the input matrix.
*/


Ask ChatGPT

    std::array<double, 4> q;

    q[0] = sqrt((1.0+R[0][0]+R[1][1]+R[2][2])/4.0);
    q[1] = sqrt((1.0+R[0][0]-R[1][1]-R[2][2])/4.0);
    q[2] = sqrt((1.0-R[0][0]+R[1][1]-R[2][2])/4.0);
    q[3] = sqrt((1.0-R[0][0]-R[1][1]+R[2][2])/4.0);

    double max_value = q[0];
    int max_index = 0;

    for (int i = 1; i < q.size(); ++i) {
        if (q[i] > max_value) {
            max_value = q[i];
            max_index = i;
        }
    }

    if (max_index == 0) {
       q[1] = -(R[2][1]-R[1][2])/(4.0*q[0]);
       q[2] = -(R[0][2]-R[2][0])/(4.0*q[0]);
       q[3] = -(R[1][0]-R[0][1])/(4.0*q[0]);
    } else if (max_index == 1) {
       q[0] = -(R[2][1]-R[1][2])/(4.0*q[1]);
       q[2] = (R[0][1]+R[1][0])/(4.0*q[1]);
       q[3] = (R[0][2]+R[2][0])/(4.0*q[1]);
    } else if (max_index == 2) {
       q[0] = -(R[0][2]-R[2][0])/(4.0*q[2]);
       q[1] = (R[0][1]+R[1][0])/(4.0*q[2]);
       q[3] = (R[1][2]+R[2][1])/(4.0*q[2]);
    } else {
       q[0] = -(R[1][0]-R[0][1])/(4.0*q[3]);
       q[1] = (R[0][2]+R[2][0])/(4.0*q[3]);
       q[2] = (R[1][2]+R[2][1])/(4.0*q[3]);
    }

    return q;
}

// Writes trajectories into LAMMPSTrajectory.lammpstrj for visualization
void writeTrajectory(int TrajectoryInterval, string& Style, int Step, int NumberOfParticles, double AASigma, double CGParticleSize, double aspectratio, double CGSigma, double BoxLengthHalf,
                     string shape, vector<double>& Sizes, vector<double>& Rx, vector<double>& Ry, vector<double>& Rz,
                     double VectorX[][3], double VectorY[][3], double VectorZ[][3]) {
/*
    Writes the particle trajectory data to a LAMMPS-style trajectory file.

    This function writes the simulation data (positions, orientations, and sizes) of particles at a given time step into a LAMMPS trajectory file (`LAMMPSTrajectory.lammpstrj`) for visualization purposes. The format of the output depends on the style of the simulation (e.g., "Virtual", "AA", "Vertex") and the shape of the particles (e.g., "cuboid", "triangularprism"). For each style, different information is written, such as particle positions, quaternion rotations, and atomic radii. The function supports periodic boundary conditions and ensures that particle positions are adjusted accordingly.

    Inputs:
      TrajectoryInterval - The frequency at which to write the trajectory (in terms of simulation steps).
      Style              - The style of the visualization, determining the data format ("Virtual", "AA", "Vertex").
      Step               - The current time step in the simulation.
      NumberOfParticles  - The total number of particles in the system.
      AASigma            - Atom diameter.
      CGParticleSize     - The size of the coarse-grained particles.
      aspectratio        - The aspect ratio of the particle shape (height/width for non-spherical shapes).
      CGSigma            - CG-bead diameter.
      BoxLengthHalf      - Half the length of the simulation box (used for periodic boundary conditions).
      shape              - The shape of the particles (e.g., "cuboid", "triangularprism").
      Sizes              - A vector containing the size of each particle.
      Rx, Ry, Rz         - Vectors containing the x, y, and z positions of each particle.
      VectorX, VectorY, VectorZ - 2D arrays containing the orientation vectors (rotation matrices) for each particle.

    Output:
      Writes the trajectory data to a LAMMPS trajectory file (`LAMMPSTrajectory.lammpstrj`).
*/

    if (Step % TrajectoryInterval == 0) {
        if (Style == "Virtual") {
            std::ofstream LAMMPSTrajectoryFile("LAMMPSTrajectory.lammpstrj", std::ios::app);
            LAMMPSTrajectoryFile << "ITEM: TIMESTEP\n" << Step << "\n";
            LAMMPSTrajectoryFile << "ITEM: NUMBER OF ATOMS\n" << NumberOfParticles << "\n";
            LAMMPSTrajectoryFile << "ITEM: BOX BOUNDS pp pp pp\n";
            LAMMPSTrajectoryFile << -BoxLengthHalf << " " << BoxLengthHalf << "\n";
            LAMMPSTrajectoryFile << -BoxLengthHalf << " " << BoxLengthHalf << "\n";
            LAMMPSTrajectoryFile << -BoxLengthHalf << " " << BoxLengthHalf << "\n";
            LAMMPSTrajectoryFile << "ITEM: ATOMS mol id type x y z quatw quati quatj quatk shapex shapey shapez\n";

            for (int i = 0; i < NumberOfParticles; ++i) {
                double ParticleSideLength = Sizes[i] * AASigma;
                double ParticleHeight = ParticleSideLength * aspectratio;
                // Construct rotation matrix from the vector
                double RotationMatrix[3][3] = {
                    {VectorX[i][0], VectorX[i][1], VectorX[i][2]},
                    {VectorY[i][0], VectorY[i][1], VectorY[i][2]},
                    {VectorZ[i][0], VectorZ[i][1], VectorZ[i][2]}
                };

                // Convert to quaternion
                auto Quaternion = rotationMatrixToQuaternion(RotationMatrix);

                // Write to file (note the order of quaternion elements)
                LAMMPSTrajectoryFile << i << " " << i << " 1 " << Rx[i] << " "
                                     << Ry[i] << " " << Rz[i] << " "
                                     << Quaternion[0] << " "
                                     << Quaternion[1] << " "
                                     << Quaternion[2] << " "
                                     << Quaternion[3] << " "
                                     << ParticleSideLength << " "
                                     << ParticleSideLength << " "
                                     << ParticleHeight << "\n";
            }
        }


        if (Style == "AA" && shape == "cuboid") {
            int numAtoms = CGParticleSize * CGParticleSize * CGParticleSize * aspectratio * NumberOfParticles;

            ofstream LAMMPSTrajectoryFile("LAMMPSTrajectory.lammpstrj", std::ios::app);
            LAMMPSTrajectoryFile << "ITEM: TIMESTEP\n" << Step << "\n";
            LAMMPSTrajectoryFile << "ITEM: NUMBER OF ATOMS\n" << numAtoms << "\n";
            LAMMPSTrajectoryFile << "ITEM: BOX BOUNDS pp pp pp\n" << -BoxLengthHalf << " " << BoxLengthHalf << "\n";
            LAMMPSTrajectoryFile << -BoxLengthHalf << " " << BoxLengthHalf << "\n";
            LAMMPSTrajectoryFile << -BoxLengthHalf << " " << BoxLengthHalf << "\n";
            LAMMPSTrajectoryFile << "ITEM: ATOMS mol id type x y z radius\n";

            double COMX,COMY,COMZ;
            int NumAtomsPerParticle = pow(CGParticleSize, 3)*aspectratio;
            double atomList[NumAtomsPerParticle][3];

            for (int i = 0; i < NumberOfParticles; ++i) {
                COMX = Rx[i];
                COMY = Ry[i];
                COMZ = Rz[i];
                double CGHeight = CGParticleSize*aspectratio;
                findAtomPositionsInCube(CGParticleSize, CGHeight, CGSigma, COMX, COMY, COMZ, VectorX[i], VectorY[i], VectorZ[i], atomList);

                for (int j = 0; j < NumAtomsPerParticle; ++j) {
                    double atomX = atomList[j][0];
                    double atomY = atomList[j][1];
                    double atomZ = atomList[j][2];

                    LAMMPSTrajectoryFile << i << " " << i * NumAtomsPerParticle + j + 1 << " 2 " << atomX << " " << atomY << " " << atomZ << " " << CGSigma / 2 << "\n";
                }
            }
        } else if (Style == "Vertex" && shape == "cuboid") {
            std::ofstream LAMMPSTrajectoryFile("LAMMPSTrajectory.lammpstrj", std::ios::app);
            LAMMPSTrajectoryFile << "ITEM: TIMESTEP\n" << Step << "\n";
            LAMMPSTrajectoryFile << "ITEM: NUMBER OF ATOMS\n" << 8 * NumberOfParticles << "\n";
            LAMMPSTrajectoryFile << "ITEM: BOX BOUNDS pp pp pp\n" << -BoxLengthHalf << " " << BoxLengthHalf << "\n";
            LAMMPSTrajectoryFile << -BoxLengthHalf << " " << BoxLengthHalf << "\n";
            LAMMPSTrajectoryFile << -BoxLengthHalf << " " << BoxLengthHalf << "\n";
            LAMMPSTrajectoryFile << "ITEM: ATOMS mol id type x y z radius\n";
            double Vertex[8][3];
            double COMX,COMY,COMZ;

            for (int i = 0; i < NumberOfParticles; ++i) {
                double ParticleSideLength = Sizes[i] * AASigma;
                double ParticleHeight = ParticleSideLength * aspectratio;
                COMX = Rx[i];
                COMY = Ry[i];
                COMZ = Rz[i];
                findVerticesOfCube(ParticleSideLength, ParticleHeight, COMX, COMY, COMZ, VectorX[i], VectorY[i], VectorZ[i], Vertex);
                for (int j = 0; j < 8; ++j) {
                    double atomX = Vertex[j][0];
                    double atomY = Vertex[j][1];
                    double atomZ = Vertex[j][2];

                    LAMMPSTrajectoryFile << i << " " << i * 8 + j + 1 << " 2 " << atomX << " " << atomY << " " << atomZ << " " << AASigma / 2 << "\n";
                }
            }
        } else if (Style == "AA" && shape == "triangularprism") {
            int numAtoms = CGParticleSize * CGParticleSize * (CGParticleSize + 1) * aspectratio * NumberOfParticles / 2.0;

            ofstream LAMMPSTrajectoryFile("LAMMPSTrajectory.lammpstrj", std::ios::app);
            LAMMPSTrajectoryFile << "ITEM: TIMESTEP\n" << Step << "\n";
            LAMMPSTrajectoryFile << "ITEM: NUMBER OF ATOMS\n" << numAtoms << "\n";
            LAMMPSTrajectoryFile << "ITEM: BOX BOUNDS pp pp pp\n" << -BoxLengthHalf << " " << BoxLengthHalf << "\n";
            LAMMPSTrajectoryFile << -BoxLengthHalf << " " << BoxLengthHalf << "\n";
            LAMMPSTrajectoryFile << -BoxLengthHalf << " " << BoxLengthHalf << "\n";
            LAMMPSTrajectoryFile << "ITEM: ATOMS mol id type x y z radius\n";

            double COMX,COMY,COMZ;
            int NumAtomsPerCube = pow(CGParticleSize, 2)*aspectratio*(CGParticleSize+1)/2.0;
            double atomList[NumAtomsPerCube][3];

            for (int i = 0; i < NumberOfParticles; ++i) {
                COMX = Rx[i];
                COMY = Ry[i];
                COMZ = Rz[i];
                double CGHeight = CGParticleSize*aspectratio;
                findAtomPositionsInTriangularPrism(CGParticleSize, CGHeight, CGSigma, COMX, COMY, COMZ, VectorX[i], VectorY[i], VectorZ[i], atomList);

                for (int j = 0; j < NumAtomsPerCube; ++j) {
                    double atomX = atomList[j][0];
                    double atomY = atomList[j][1];
                    double atomZ = atomList[j][2];

                    LAMMPSTrajectoryFile << i << " " << i * NumAtomsPerCube + j + 1 << " 2 " << atomX << " " << atomY << " " << atomZ << " " << CGSigma / 2 << "\n";
                }
            }
        }
    }
}

double ComputeMinimumEnergy(
    const std::string& model,
    double AASigma,
    double CGSigma,
    double CGEpsilon,
    double cutoffCG,
    double Hamaker,
    double AtomDensity,
    int Size_i,
    int Size_j,
    double aspectratio,
    double BoxLength,
    string shape) {

/*
    Computes the global minimum energy for determining the necessary scaling factor to match the minima of the atomistic model.

    This function calculates the global minimum energy for a pair of particles in a specified configuration (parallel) and shape (e.g., "cuboid", "triangularprism"). It iterates over a range of particle separations (based on the particle sizes and scaling factors) and computes the energy between two particles at each separation. The goal is to find the separation that minimizes the interaction energy between the particles, which is used to determine the scaling factor required for matching the atomistic model's energy profile.

    The energy is computed using the `EnergyBetweenTwoParticles` function, which takes into account various physical parameters such as particle size, epsilon (interaction strength), Hamaker constant, and the atom density, among others.

    Inputs:
      model            - The model of the system (e.g., "AACG" or "vdW").
      AASigma          - The atom diameter.
      CGSigma          - The CG-bead diameter.
      CGEpsilon        - The interaction strength between coarse-grained particles.
      cutoffCG         - The cutoff distance for the coarse-grained model interactions.
      Hamaker          - The Hamaker constant for van der Waals interactions.
      AtomDensity      - The atomic density of the material.
      Size_i, Size_j   - The number of atoms per side.
      aspectratio      - The aspect ratio of the particle shape (height/width).
      BoxLength        - The length of the simulation box (used for periodic boundary conditions).
      shape            - The shape of the particles ("cuboid" or "triangularprism").

    Output:
      Returns the minimum energy value found for the particle configuration, which corresponds to the separation at which the particles experience the least interaction energy.
*/

    cout << "Scale factor for matching with atomistic global minimum energy is being computed" << endl;
    double min_energy = 100000; // very big value to initialize
    if (shape == "cuboid") {
      // Fixed orientation vectors (parallel configuration)
      Vector3 VectorX_i = {1.0, 0.0, 0.0};
      Vector3 VectorY_i = {0.0, 1.0, 0.0};
      Vector3 VectorZ_i = {0.0, 0.0, 1.0};
      Vector3 VectorX_j = {1.0, 0.0, 0.0};
      Vector3 VectorY_j = {0.0, 1.0, 0.0};
      Vector3 VectorZ_j = {0.0, 0.0, 1.0};

      Vector3 ParticleiCentroid = {0.0, 0.0, 0.0}; // particle 1 at origin

      double x_start = Size_i * AASigma; // two particles touch each other, and you vary x from there to x_end to find the minimum energy value
      double x_end = x_start + CGSigma;
      Vector3 ParticlejCentroid = {x_start, 0.0, 0.0};

      for (double x = x_start; x <= x_end; x += 0.01) {
          ParticlejCentroid[0] = x; // changing x coordinate of particle 2

          auto result = EnergyBetweenTwoParticles(
              model,
              AASigma,
              CGSigma,
              CGEpsilon,
              cutoffCG,
              Hamaker,
              AtomDensity,
              Size_i,
              Size_j,
              aspectratio,
              BoxLength,
              ParticleiCentroid,
              ParticlejCentroid,
              VectorX_i,
              VectorY_i,
              VectorZ_i,
              VectorX_j,
              VectorY_j,
              VectorZ_j,
              shape
          );

          double energy = result.first;
          bool error = result.second;

          if (!error && energy < min_energy && !std::isinf(std::abs(energy))) {
              min_energy = energy;
          }
      }
    } else if (shape == "triangularprism") {
      Vector3 VectorX_i = {1.0, 0.0, 0.0};
      Vector3 VectorY_i = {0.0, 1.0, 0.0};
      Vector3 VectorZ_i = {0.0, 0.0, 1.0};
      Vector3 VectorX_j = {-1.0, 0.0, 0.0};
      Vector3 VectorY_j = {0.0, -1.0, 0.0};
      Vector3 VectorZ_j = {0.0, 0.0, 1.0};

      Vector3 ParticleiCentroid = {0.0, 0.0, 0.0}; // particle 1 at origin

      double x_start = (Size_i * AASigma / CGSigma - 1) * CGSigma / sqrt(3) + CGSigma; // two particles touch each other
      double x_end = x_start + CGSigma;
      Vector3 ParticlejCentroid = {x_start, 0.0, 0.0};

      for (double x = x_start; x <= x_end; x += 0.01) {
          ParticlejCentroid[0] = x; // changing x coordinate of particle 2

          auto result = EnergyBetweenTwoParticles(
            model,
            AASigma,
            CGSigma,
            CGEpsilon,
            cutoffCG,
            Hamaker,
            AtomDensity,
            Size_i,
            Size_j,
            aspectratio,
            BoxLength,
            ParticleiCentroid,
            ParticlejCentroid,
            VectorX_i,
            VectorY_i,
            VectorZ_i,
            VectorX_j,
            VectorY_j,
            VectorZ_j,
            shape
          );

          double energy = result.first;
          bool error = result.second;

          if (!error && energy < min_energy && !std::isinf(std::abs(energy))) {
              min_energy = energy;
          }
      }

    }
    return min_energy;
}

// Restart from saved .csv file
void readRestart(const std::string& restartFile, int NumberOfParticles,
                 std::vector<double>& Sizes, std::vector<double>& Rx,
                 std::vector<double>& Ry, std::vector<double>& Rz,
                 double VectorX[][3], double VectorY[][3], double VectorZ[][3],
                 int& restartStep) {
/* Reads the restart data from a CSV file and restores the particle configuration (sizes, positions, and orientation vectors).

    This function reads the saved simulation state from a CSV file and restores the particle properties. The data is parsed from the file, and for each particle, the size, position (Rx, Ry, Rz), and orientation vectors (VectorX, VectorY, VectorZ) are restored into the corresponding vectors.

    Inputs:
      restartFile        - The path to the CSV file containing the saved restart data.
      NumberOfParticles  - The number of particles in the simulation.
      Sizes              - A vector to store the particle sizes.
      Rx, Ry, Rz         - Vectors to store the x, y, z coordinates (positions) of the particles.
      VectorX, VectorY, VectorZ - 2D arrays to store the orientation vectors (X, Y, Z) for each particle.
      restartStep        - A reference to the step number at which the restart data was saved.

    Outputs:
      This function modifies the input vectors (`Sizes`, `Rx`, `Ry`, `Rz`, `VectorX`, `VectorY`, `VectorZ`) with the values from the restart file. The `restartStep` is updated to reflect the last step from the file.
*/
    std::ifstream inputFile(restartFile);

    if (!inputFile.is_open()) {
        std::cerr << "Error: Unable to open file " << restartFile << std::endl;
        return;
    }

    std::string line;
    std::getline(inputFile, line);  // Read and discard the header line
    int num = 0;

    // Iterate through each line in the CSV file
    while (std::getline(inputFile, line)) {

        // Use a string stream to parse the values
        std::stringstream ss(line);
        std::string value;
        std::vector<double> values;

        // Iterate through each value in the line
        while (std::getline(ss, value, ',')) {
            values.push_back(std::stod(value));
        }

        // Now you have the values in the 'values' vector
        // Example usage:
        double size = values[0];
        double rx = values[1];
        double ry = values[2];
        double rz = values[3];
        double vxx = values[4];
        double vxy = values[5];
        double vxz = values[6];
        double vyx = values[7];
        double vyy = values[8];
        double vyz = values[9];
        double vzx = values[10];
        double vzy = values[11];
        double vzz = values[12];
        Sizes[num] = size;
        Rx[num] = rx;
        Ry[num] = ry;
        Rz[num] = rz;
        VectorX[num][0] = vxx;
        VectorX[num][1] = vxy;
        VectorX[num][2] = vxz;
        VectorY[num][0] = vyx;
        VectorY[num][1] = vyy;
        VectorY[num][2] = vyz;
        VectorZ[num][0] = vzx;
        VectorZ[num][1] = vzy;
        VectorZ[num][2] = vzz;

        num += 1;
    }

    // Close the input file
    inputFile.close();
}

// Save the center coordinates and orientation vectors of particle (that can be utilized later on to restart the simulation from this saved point)
int writeRestart(int RestartFileInterval, int lastRestartStep, vector<double>& Sizes,
                 vector<double>& Rx, vector<double>& Ry, vector<double>& Rz,
                 double VectorX[][3], double VectorY[][3], double VectorZ[][3], int Step) {
/*
    Writes the current simulation state (particle sizes, positions, and orientation vectors) to a CSV file for later restart.

    This function saves the current particle configuration (sizes, positions, and orientation vectors) to a CSV file. The data is written every `RestartFileInterval` steps. The CSV file is named based on the current simulation step, and it stores particle data in the format: Sizes, Rx, Ry, Rz, VectorX, VectorY, VectorZ.

    Inputs:
      RestartFileInterval - The interval (in steps) at which to write the restart data.
      lastRestartStep      - The step number of the last restart, used to name the previous restart file.
      Sizes                - A vector containing the sizes of the particles.
      Rx, Ry, Rz           - Vectors containing the x, y, z positions of the particles.
      VectorX, VectorY, VectorZ - 2D arrays containing the orientation vectors of the particles.
      Step                 - The current simulation step.

    Outputs:
      This function writes a CSV file containing the particle configuration at the current step. It returns the updated `lastRestartStep` value to keep track of the last saved restart step.
*/

    // write restart file
    if (Step % RestartFileInterval == 0) {
        string lastRestartFile = "restart_" + std::to_string(lastRestartStep) + ".csv";
        ofstream outputFile("restart_" + std::to_string(Step) + ".csv");

        if (outputFile.is_open()) {
            outputFile << "Sizes,Rx,Ry,Rz,VectorXx,VectorXy,VectorXz,VectorYx,VectorYy,VectorYz,VectorZx,VectorZy,VectorZz\n";

            for (size_t i = 0; i < Sizes.size(); ++i) {
                outputFile << Sizes[i] << "," << Rx[i] << "," << Ry[i] << "," << Rz[i] << ","
                           << VectorX[i][0] << "," << VectorX[i][1] << "," << VectorX[i][2] << ","
                           << VectorY[i][0] << "," << VectorY[i][1] << "," << VectorY[i][2] << ","
                           << VectorZ[i][0] << "," << VectorZ[i][1] << "," << VectorZ[i][2] << "\n";
            }

            outputFile.close();

            //if (std::filesystem::exists(lastRestartFile)) {
            //    std::filesystem::remove(lastRestartFile);
            //}

            lastRestartStep = Step;
        }
    }

    return lastRestartStep;
}

void writeEnergy(int EnergyOutputInterval, int Step, double SystemPotentialEnergy, ofstream& EnergyFile) {
/*
    Writes the system's potential energy to a file at specified intervals.

    This function writes the system's potential energy to a specified output file (`EnergyFile`). It only writes the energy at intervals defined by `EnergyOutputInterval`, ensuring that energy data is recorded periodically without overwhelming the output file.

    Inputs:
      EnergyOutputInterval - The interval (in steps) at which to write the energy data.
      Step                 - The current simulation step.
      SystemPotentialEnergy - The current potential energy of the system.
      EnergyFile          - The output file stream to which the energy data will be written.

    Outputs:
      This function writes the potential energy to the output file at the specified intervals. The energy is recorded alongside the simulation step number for reference.
*/

    // Write out energies every EnergyOutputInterval steps
    if (Step % EnergyOutputInterval == 0) {
       EnergyFile << Step << " " << SystemPotentialEnergy << "\n" << flush;
    }
}

// Main function (you need to change only this section where you set the parameters of simulations)
int main() {

    // Variables to store parameters
    string Style, model, shape;
    int NumberOfParticles, EquilibrationSteps, ProductionSteps, EnergyOutputInterval, TrajectoryInterval;
    int RestartFileInterval, volumeMoveInterval;
    double BoxLength, Temperature, kBT, globalminEnergy, AASigma, AAEpsilon, MaxRotation, MaxMoveSingle, MaxMoveCluster, oldvolfrac;
    double SingleMoveProbability, ParticleSideAASize, aspectratio, CGSigma, Hamaker, AtomDensity;
    bool volumeexpansion;

    double kB = 0.0019872;  // Boltzmann constant, kcal/(mol*K)

    // Open the file
    ifstream inputFile("params.txt");
    string line;

    if (!inputFile.is_open()) {
        cerr << "Unable to open file" << endl;
        return 1;
    }

    // Read the file line by line
    while (getline(inputFile, line)) {
        // Remove leading/trailing whitespace
        line.erase(0, line.find_first_not_of(" \t"));
        line.erase(line.find_last_not_of(" \t") + 1);

        // Skip empty lines or lines that start with '#'
        if (line.empty() || line[0] == '#') continue;

        // Split by '=' and extract variable name and value
        size_t delimiterPos = line.find("=");
        if (delimiterPos != string::npos) {
            string varName = line.substr(0, delimiterPos);
            string value = line.substr(delimiterPos + 1);

            // Remove any extra spaces from the variable and value
            varName.erase(0, varName.find_first_not_of(" \t"));
            varName.erase(varName.find_last_not_of(" \t") + 1);
            value.erase(0, value.find_first_not_of(" \t"));
            value.erase(value.find_last_not_of(" \t") + 1);

            // Assign values based on the variable names
            if (varName == "Style") {
                Style = value;
            } else if (varName == "model") {
                model = value;
            } else if (varName == "NumberOfParticles") {
                NumberOfParticles = stoi(value);
            } else if (varName == "BoxLength") {
                BoxLength = stod(value);
            } else if (varName == "Temperature") {
                Temperature = stod(value);
                kBT = kB * Temperature;
            } else if (varName == "globalminEnergy") {
                globalminEnergy = stod(value) * kBT;
            } else if (varName == "AASigma") {
                AASigma = stod(value);
            } else if (varName == "AAEpsilon") {
                AAEpsilon = stod(value);
            } else if (varName == "MaxRotation") {
                MaxRotation = stod(value) * M_PI / 180;
            } else if (varName == "MaxMoveSingle") {
                MaxMoveSingle = stod(value) * ParticleSideAASize * AASigma;
            } else if (varName == "MaxMoveCluster") {
                MaxMoveCluster = stod(value) * ParticleSideAASize * AASigma;
            } else if (varName == "SingleMoveProbability") {
                SingleMoveProbability = stod(value);
            } else if (varName == "EquilibrationSteps") {
                EquilibrationSteps = stoi(value);
            } else if (varName == "ProductionSteps") {
                ProductionSteps = stoi(value);
            } else if (varName == "EnergyOutputInterval") {
                EnergyOutputInterval = stoi(value);
            } else if (varName == "TrajectoryInterval") {
                TrajectoryInterval = stoi(value);
            } else if (varName == "RestartFileInterval") {
                RestartFileInterval = stoi(value);
            } else if (varName == "volumeMoveInterval") {
                volumeMoveInterval = stoi(value);
            } else if (varName == "ParticleSideAASize") {
                ParticleSideAASize = stod(value);
            } else if (varName == "aspectratio") {
                aspectratio = stod(value);
            } else if (varName == "volumeexpansion") {
                volumeexpansion = (value == "true");
            } else if (varName == "CGSigma") {
                CGSigma = stod(value);
            } else if (varName == "Hamaker") {
                Hamaker = stod(value);
            } else if (varName == "AtomDensity") {
                AtomDensity = stod(value);
            } else if (varName == "shape") {
                shape = value;  // Store the shape type
            }
        }
    }

    // Close the file
    inputFile.close();

    // Debug print: Check if the values are correctly assigned
    cout << "Style: " << Style << "\n";
    cout << "Model: " << model << "\n";
    cout << "NumberOfParticles: " << NumberOfParticles << "\n";
    cout << "BoxLength: " << BoxLength << "\n";
    cout << "Temperature: " << Temperature << " K\n";
    cout << "Global Minimum Energy: " << globalminEnergy << "\n";
    cout << "AASigma: " << AASigma << "\n";
    cout << "MaxRotation: " << MaxRotation << "\n";
    cout << "MaxMoveSingle: " << MaxMoveSingle << "\n";
    cout << "EquilibrationSteps: " << EquilibrationSteps << "\n";
    cout << "AtomDensity: " << AtomDensity << "\n";
    double BoxLengthHalf = BoxLength/2.0;
    int lastRestartStep = 0;
    if (shape == "cuboid") {
      oldvolfrac = NumberOfParticles*pow(ParticleSideAASize*AASigma,3)*aspectratio/pow(BoxLength,3); // volume fraction of the system
    } else if (shape == "triangularprism") {
      oldvolfrac = NumberOfParticles*pow(ParticleSideAASize*AASigma,3)*sqrt(3)*0.25*aspectratio/pow(BoxLength,3);
    }
    unordered_map<double, double> SizeRatioDictionary = {{ParticleSideAASize, 1.0}};
    double CGEpsilon, CGParticleSize, cutoffCG, MaxMove;

    if (model == "AACG") {
       cutoffCG = 4.0 * CGSigma; // Cutoff distance for inter-bead energy computation
       CGParticleSize = ParticleSideAASize * AASigma / CGSigma; // Number of CG beads per side
       double min_energy_AA = ComputeMinimumEnergy(model,AASigma,AASigma,AAEpsilon,4.0*AASigma,0.0,0.0,ParticleSideAASize,ParticleSideAASize,aspectratio,BoxLength,shape);
       double min_energy_CG = ComputeMinimumEnergy(model,AASigma,CGSigma,AAEpsilon,4.0*CGSigma,0.0,0.0,ParticleSideAASize,ParticleSideAASize,aspectratio,BoxLength,shape);
       double ScaleFactor = min_energy_AA / min_energy_CG; //Scale epsilon of AA to match minimum of CG energy with minimum of AA energy
       //cout << "ScaleFactor " << ScaleFactor << endl;
       AAEpsilon = AAEpsilon * (globalminEnergy / abs(min_energy_AA)); // Necessary AA Epsilon value to have global minimum energy
       // Scale epsilon of AA to match minimum of CG energy with minimum of AA energy
       CGEpsilon = AAEpsilon * ScaleFactor; // Epsilon of value of CG model
    } else if (model == "vdW") {
       CGSigma = 1 * AASigma; // vdW is not CG model, this CG size is just for trajectory file outputing AA style purposes
       CGParticleSize = ParticleSideAASize; // vdW is not CG model, this CG size is just for trajectory file outputing AA style purposes
       double min_energy_AA = ComputeMinimumEnergy("AACG",AASigma,AASigma,AAEpsilon,4.0*AASigma,0.0,0.0,ParticleSideAASize,ParticleSideAASize,aspectratio,BoxLength,shape);
       double min_energy_vdW = ComputeMinimumEnergy(model,AASigma,CGSigma,AAEpsilon,4.0*AASigma,Hamaker,AtomDensity,ParticleSideAASize,ParticleSideAASize,aspectratio,BoxLength,shape);
       double ScaleFactor = min_energy_AA / min_energy_vdW; // scale analytical vdW Hamaker constant to match the minimum of AA energy
       //cout << "ScaleFactor " << ScaleFactor << endl;
       ScaleFactor = ScaleFactor * (globalminEnergy / abs(min_energy_AA)); // Necessary scale factor to reach the determined global minimum energy
       Hamaker = Hamaker * ScaleFactor; // Final Hamaker constant yielding a global minimum energy set in the main function
       //cout << "Hamaker " << Hamaker << endl;
    }

    vector<double> Sizes(NumberOfParticles, 0.0);
    vector<double> Rx(NumberOfParticles, 0.0);
    vector<double> Ry(NumberOfParticles, 0.0);
    vector<double> Rz(NumberOfParticles, 0.0);
    double VectorX[NumberOfParticles][3];
    double VectorY[NumberOfParticles][3];
    double VectorZ[NumberOfParticles][3];

    int restartStep = 0;
    int Step = 0;

    ofstream LAMMPSTrajectoryFile("LAMMPSTrajectory.lammpstrj");
    ofstream EnergyFile("Energies.out");
    ofstream TimeFile("Time.out");
    //ofstream RDFFile("rdf.out");
    //ofstream MeanRDFFile("mean_rdf.out");
    //ofstream MeanAggFile("meanagg.out");
    //RDFFile.flush();
    //MeanRDFFile.flush();
    //MeanAggFile.flush();
    // Call the function to create initial configuration of the system at selected number of particles and simulation box size
    if (shape == "cuboid") {
      InitialConfigurationCuboid(NumberOfParticles, BoxLength, AASigma, aspectratio, shape, SizeRatioDictionary,
                        Sizes, Rx, Ry, Rz, VectorX, VectorY, VectorZ,
                        restartStep, LAMMPSTrajectoryFile, EnergyFile, TimeFile);
    } else if (shape == "triangularprism") {
       InitialConfigurationTriangularPrism(NumberOfParticles, BoxLength, AASigma, aspectratio, shape, SizeRatioDictionary,
                        Sizes, Rx, Ry, Rz, VectorX, VectorY, VectorZ,
                        restartStep, LAMMPSTrajectoryFile, EnergyFile, TimeFile);
    }
    //string restartFile = "restart_25000000.csv";
    //readRestart(restartFile,NumberOfParticles,Sizes,Rx,Ry,Rz,VectorX,VectorY,VectorZ,restartStep);

    writeTrajectory(TrajectoryInterval, Style, Step, NumberOfParticles,AASigma,CGParticleSize,aspectratio,CGSigma,BoxLengthHalf,shape,Sizes,Rx,Ry,Rz,VectorX,VectorY,VectorZ);
    lastRestartStep = writeRestart(RestartFileInterval,lastRestartStep,Sizes,Rx,Ry,Rz,VectorX,VectorY,VectorZ,Step);
    ofstream ProbFile("prob.out");
    ProbFile.flush();
    // Initial total energy of the system
    double SystemPotentialEnergy = TotalEnergy(Sizes,Rx,Ry,Rz,VectorX,VectorY,VectorZ,NumberOfParticles,BoxLength,aspectratio,CGSigma,CGEpsilon,AASigma,model,cutoffCG,Hamaker,AtomDensity,shape);
    cout << "Total energy " << SystemPotentialEnergy << endl;
    double Attempt = 0.0;
    double Accept = 0.0;
    int CountBoltzmann = 0;
    double SingleMoveAccept = 0.0;
    double SingleMoveNumber = 0.0;
    double ClusterMoveAccept = 0.0;
    double ClusterMoveNumber = 0.0;

    // For rdf calculation
    //int num_bins = 1000;
    //double normalized_bin_counts[num_bins];
    //double mean_bin_counts[num_bins];
    //for (int k = 0; k < num_bins; ++k) {
      // mean_bin_counts[k] = 0.0;
    //}
    //double bin_width = 0.1;
    //double normalization_constant = NumberOfParticles * (NumberOfParticles / pow(BoxLength,3.0)) * 4.0 * M_PI * bin_width;

    for (int Step = 1 + restartStep; Step <= restartStep + EquilibrationSteps + ProductionSteps; ++Step) {
        //double mean_agg = 0.0;
        //for (int k = 0; k < num_bins; ++k) {
          // normalized_bin_counts[k] = 0.0;
        //}
        bool EarlyTermination = false;
        bool overLapFlag = false;
        //ProbFile << " " << endl;
        //ProbFile << "Step " << Step << endl;
        //cout << "Step " << Step << endl;
        //cout << endl;
        if (Step == 1) {
           cout << "\nStarting Equilibration\n";
        }
        if (Step == EquilibrationSteps + 1) {
           cout << "\nStarting Production\n";
        }

        if (Step % 100 == 0) {
           if (Step <= EquilibrationSteps) {
              cout << "Equilibration Step " << Step << "\n";
           } else {
              cout << "Production Step " << Step << "\n";
           }
        }

        // Choose a random particle i as seed
        int i = static_cast<int>(NumberOfParticles * ((double)rand()/(RAND_MAX)));
        //ProbFile << "Seed particle: " << i << endl;

        // Position and orientation of the seed
        double seedX = Rx[i]; // center of mass in x, Angstrom
        double seedY = Ry[i]; // center of mass in y, Angstrom
        double seedZ = Rz[i]; // center of mass in z, Angstrom
        Vector3 seedCentroid = {seedX, seedY, seedZ};
        Vector3 seedVectorX, seedVectorY, seedVectorZ;

        for (int k = 0; k < 3; ++k) {
             seedVectorX[k] = VectorX[i][k]; // x component of orientation vector
             seedVectorY[k] = VectorY[i][k]; // y component of orientation vector
             seedVectorZ[k] = VectorZ[i][k]; // z component of orientation vector
        }

        // Create a list to recruit particles into a cluster
        vector<int> clusterParticleList;
        // Recruit the seed particle as the first member of the cluster
        clusterParticleList.push_back(i);
        // Create a list for particles haven't been recruited to the cluster
        vector<int> ParticlesNotInClusterList(NumberOfParticles);
        for (int ParticleID = 0; ParticleID < NumberOfParticles; ++ParticleID) {
            ParticlesNotInClusterList[ParticleID] = ParticleID;
        }
        // Particles that for sure will be recruited to cluster but due to multiple of such particles bonding with a linker at the same time
        // They will be revisited later after going through one branch
        // Each of these particles will begin a new branch of recruitment
        vector<int> ToBeRecruitedParticlesList;
        // Temporary List for the j(s) linked with current i
        vector<int> TemperaryjList;

        // Removing the seed particle from ParticlesNotInClusterList
        ParticlesNotInClusterList.erase(remove(ParticlesNotInClusterList.begin(), ParticlesNotInClusterList.end(), i), ParticlesNotInClusterList.end());

        vector<double> reverseMoveProbabilityList, forwardMoveProbabilityList, unacceptedMoveProbabilityList, unacceptedReverseMoveProbabilityList;
        vector<int> linkConnectionList;

        double MaxClusterSize;
        double singleOrCluster = ((double)rand()/(RAND_MAX)); //to determine if move is single or cluster move
        if (oldvolfrac > 0.4) { // at higher volume fractions, only single moves are utilized to reduce computational cost
           singleOrCluster = 0.0;
        }
        if (singleOrCluster < SingleMoveProbability) {
           MaxClusterSize = 1; // single move with 50% probability
           MaxMove = MaxMoveSingle;
           //ProbFile << "Single move selected" << endl;
        } else {
           MaxClusterSize = NumberOfParticles / 6.0; // After this number of particles, it is time consuming to create clusters and then reject
           MaxMove = MaxMoveCluster;
        }

        double tranOrRot = ((double)rand()/(RAND_MAX)); // to determine if the move is translation or rotation
        bool isTranslation = false;
        bool isRotation = false;
        double VirtualMoveX, VirtualMoveY, VirtualMoveZ, VirtualMoveAlpha, VirtualMoveBeta, VirtualMoveGamma;

        if (tranOrRot <= 0.5) { // 50% chance of translation or rotation
           // Translation
           isTranslation = true;
           VirtualMoveX = (2 * ((double)rand()/(RAND_MAX)) - 1) * MaxMove;
           VirtualMoveY = (2 * ((double)rand()/(RAND_MAX)) - 1) * MaxMove;
           VirtualMoveZ = (2 * ((double)rand()/(RAND_MAX)) - 1) * MaxMove;
           VirtualMoveAlpha = 0.0;
           VirtualMoveBeta = 0.0;
           VirtualMoveGamma = 0.0;
        }else {
           // Rotation
           isRotation = true;
           VirtualMoveX = 0;
           VirtualMoveY = 0;
           VirtualMoveZ = 0;
           double randomAxis = ((double)rand()/(RAND_MAX));
           // randomly choose an axis from x,y,z to rotate
           if (randomAxis <= 0.333) { // pick x-axis to rotate, the other two rotation angles are zero
              VirtualMoveAlpha = (2 * ((double)rand()/(RAND_MAX)) - 1) * MaxRotation;
              VirtualMoveBeta = 0.0;
              VirtualMoveGamma = 0.0;
           } else if (randomAxis <= 0.666){
              VirtualMoveAlpha = 0.0;
              VirtualMoveBeta = (2 * ((double)rand()/(RAND_MAX)) - 1) * MaxRotation;
              VirtualMoveGamma = 0.0;
           } else {
              VirtualMoveAlpha = 0.0;
              VirtualMoveBeta = 0.0;
              VirtualMoveGamma = (2 * ((double)rand()/(RAND_MAX)) - 1) * MaxRotation;
           }
        }
        // Rotation matrices based on these rotation angles
        Matrix3x3 VirtualMoveRotationX = {{1, 0, 0}, {0, cos(VirtualMoveAlpha), -sin(VirtualMoveAlpha)}, {0, sin(VirtualMoveAlpha), cos(VirtualMoveAlpha)}};
        Matrix3x3 VirtualMoveRotationY = {{cos(VirtualMoveBeta), 0, sin(VirtualMoveBeta)}, {0, 1, 0}, {-sin(VirtualMoveBeta), 0, cos(VirtualMoveBeta)}};
        Matrix3x3 VirtualMoveRotationZ = {{cos(VirtualMoveGamma), -sin(VirtualMoveGamma), 0}, {sin(VirtualMoveGamma), cos(VirtualMoveGamma), 0}, {0, 0, 1}};

        Matrix3x3 ReverseMoveRotationX = {{1, 0, 0}, {0, cos(-VirtualMoveAlpha), -sin(-VirtualMoveAlpha)}, {0, sin(-VirtualMoveAlpha), cos(-VirtualMoveAlpha)}};
        Matrix3x3 ReverseMoveRotationY = {{cos(-VirtualMoveBeta), 0, sin(-VirtualMoveBeta)}, {0, 1, 0}, {-sin(-VirtualMoveBeta), 0, cos(-VirtualMoveBeta)}};
        Matrix3x3 ReverseMoveRotationZ = {{cos(-VirtualMoveGamma), -sin(-VirtualMoveGamma), 0}, {sin(-VirtualMoveGamma), cos(-VirtualMoveGamma), 0}, {0, 0, 1}};

        //ProbFile << "VirtualMoveX " << VirtualMoveX << endl;
        //ProbFile << "VirtualMoveY " << VirtualMoveY << endl;
        //ProbFile << "VirtualMoveZ " << VirtualMoveZ << endl;
        //ProbFile << "VirtualMoveAlpha " << VirtualMoveAlpha << " VirtualMoveBeta " << VirtualMoveBeta << " VirtualMoveGamma " << VirtualMoveGamma << endl;

        // if maximum cluster size is 1, then there is no need to build cluster
        // because the particle itself is considered a cluster
        bool Loop;
        double overallReverseProbability, overallUnacceptedReverseProbability, overallMoveProbability, overallUnacceptedMoveProbability;
        if (MaxClusterSize > 1) {
           Loop = true; // continue to loop if Loop is True
        } else {
           Loop = false; // continue to loop if Loop is True
           overallReverseProbability = 1.0;
           overallUnacceptedReverseProbability = 1.0;
           overallMoveProbability = 1.0;
           overallUnacceptedMoveProbability = 1.0;
        }
        // to regulate cutoff for cluster size
        double Nc = ((double)rand()/(RAND_MAX));
        //ProbFile << "Nc " << (1/Nc) << endl;

        //ProbFile << "clusterParticleList: " << endl;
        //for (int particleID : clusterParticleList) {
        //   ProbFile << particleID << " ";
        //}
        //ProbFile << endl;
        while (Loop) {
           // check if size of cluster exceeds maximum
           if (clusterParticleList.size() >= MaxClusterSize) {
              Loop = false;
              //ProbFile << "Cluster exceeds maximum cluster size" << endl;
              break;
           }
           // to abort the link formation procedure if the cluster size exceeds Nc
           if (1.0 / clusterParticleList.size() < Nc) {
              //ProbFile << "Early termination" << endl;
              //ProbFile << "clusterParticleList: " << endl;
              //for (int particleID : clusterParticleList) {
              //   ProbFile << particleID << " ";
              //}
              //ProbFile << endl;
              Loop = false; // End recruiting particles to cluster
              EarlyTermination = true;
              break;
           }
           double Size_i = Sizes[i];
           double Rx_i_old = Rx[i];
           double Ry_i_old = Ry[i];
           double Rz_i_old = Rz[i];
           //ProbFile << "Rx_i_old " << Rx_i_old << " Ry_i_old " << Ry_i_old << " Rz_i_old " << Rz_i_old << endl;
           Vector3 ParticleiCentroid_old = {Rx_i_old, Ry_i_old, Rz_i_old};
           Vector3 VectorX_i_old, VectorY_i_old, VectorZ_i_old, VectorX_j_old, VectorY_j_old, VectorZ_j_old;
           for (int k = 0; k < 3; ++k) {
              VectorX_i_old[k] = VectorX[i][k]; // x component of orientation vector
              VectorY_i_old[k] = VectorY[i][k]; // y component of orientation vector
              VectorZ_i_old[k] = VectorZ[i][k]; // z component of orientation vector
           }
           // loop over all j(s) in ParticlesNotInClusterList
           // to find all j(s) that link to the current i
           for (int ParticleNotInCluster : ParticlesNotInClusterList) {
              // Go through each particle in ParticlesNotInClusterList
              // to check i-j interaction to determine if j should be recruited to cluster
              int j = ParticleNotInCluster;
              //ProbFile << i << " and " << j << " investigation" << endl;

              //ProbFile << "Rx_i_old " << Rx_i_old << " Ry_i_old " << Ry_i_old << " Rz_i_old " << Rz_i_old << endl;
              double Size_i = Sizes[i];
              double Size_j = Sizes[j];
              double Rx_j_old = Rx[j];  // center of mass in x, Angstrom
              double Ry_j_old = Ry[j];  // center of mass in y, Angstrom
              double Rz_j_old = Rz[j];  // center of mass in z, Angstrom
              //ProbFile << "Rx_j_old " << Rx_j_old << " Ry_j_old " << Ry_j_old << " Rz_j_old " << Rz_j_old << endl;
              for (int k = 0; k < 3; ++k) {
                 VectorX_j_old[k] = VectorX[j][k]; // x component of orientation vector
                 VectorY_j_old[k] = VectorY[j][k]; // y component of orientation vector
                 VectorZ_j_old[k] = VectorZ[j][k]; // z component of orientation vector
              }

              //ProbFile << "VectorX_i_old[0] " << VectorX_i_old[0] << " VectorX_i_old[1] " << VectorX_i_old[1] << " VectorX_i_old[2] " << VectorX_i_old[2] <<endl;
              //ProbFile << "VectorY_i_old[0] " << VectorY_i_old[0] << " VectorY_i_old[1] " << VectorY_i_old[1] << " VectorY_i_old[2] " << VectorY_i_old[2] <<endl;
              //ProbFile << "VectorZ_i_old[0] " << VectorZ_i_old[0] << " VectorZ_i_old[1] " << VectorZ_i_old[1] << " VectorZ_i_old[2] " << VectorZ_i_old[2] <<endl;
              //ProbFile << "VectorX_j_old[0] " << VectorX_j_old[0] << " VectorX_j_old[1] " << VectorX_j_old[1] << " VectorX_j_old[2] " << VectorX_j_old[2] <<endl;
              //ProbFile << "VectorY_j_old[0] " << VectorY_j_old[0] << " VectorY_j_old[1] " << VectorY_j_old[1] << " VectorY_j_old[2] " << VectorY_j_old[2] <<endl;
              //ProbFile << "VectorZ_j_old[0] " << VectorZ_j_old[0] << " VectorZ_j_old[1] " << VectorZ_j_old[1] << " VectorZ_j_old[2] " << VectorZ_j_old[2] <<endl;
              // relative distance between i and j
              double Rx_ij_old = Rx_i_old - Rx_j_old;
              double Ry_ij_old = Ry_i_old - Ry_j_old;
              double Rz_ij_old = Rz_i_old - Rz_j_old;
              // smallest relative distance between i and j considering periodic boundary conditions
              Rx_ij_old = Rx_ij_old - BoxLength * round(Rx_ij_old/BoxLength);
              Ry_ij_old = Ry_ij_old - BoxLength * round(Ry_ij_old/BoxLength);
              Rz_ij_old = Rz_ij_old - BoxLength * round(Rz_ij_old/BoxLength);
              // absolute distance between i and j
              double ijDistanceSquare_old = Rx_ij_old * Rx_ij_old + Ry_ij_old * Ry_ij_old + Rz_ij_old * Rz_ij_old;
              double ijDistance_old = pow(ijDistanceSquare_old, 0.5);
              // i as reference, move j to the image where ij distance is minimal
              Rx_j_old = Rx_i_old - Rx_ij_old;
              Ry_j_old = Ry_i_old - Ry_ij_old;
              Rz_j_old = Rz_i_old - Rz_ij_old;
              Vector3 ParticlejCentroid_old = {Rx_j_old, Ry_j_old, Rz_j_old};
              //ProbFile << "ParticleiCentroid_old[0] " << ParticleiCentroid_old[0] << " ParticleiCentroid_old[1] " << ParticleiCentroid_old[1] << " ParticleiCentroid_old[2] " << ParticleiCentroid_old[2] << endl;
              //ProbFile << "ParticlejCentroid_old[0] " << ParticlejCentroid_old[0] << " ParticlejCentroid_old[1] " << ParticlejCentroid_old[1] << " ParticlejCentroid_old[2] " << ParticlejCentroid_old[2] << endl;

              // larger size
              double MaxSize = max(Size_i, Size_j);
              double CutOff = 1.5 * MaxSize * AASigma * max(1.0,aspectratio);
              double CutOffSquare = CutOff*CutOff;
              double ijEnergy_old, ijEnergy_new, ParticleSideLength1, ParticleSideLength2, ParticleHeight1, ParticleHeight2, ijEnergy_reverse;
              double virtualMoveProbability,reverseMoveProbability;
              bool error_flag;

              if (ijDistanceSquare_old <= CutOffSquare) {
                 //ProbFile << i << " and " << j << " investigation " << endl;
                 // check energy between i and j before vitual move
                 pair<double, bool> result_old = EnergyBetweenTwoParticles(model,AASigma,CGSigma,CGEpsilon,cutoffCG,Hamaker,AtomDensity,Size_i,Size_j,aspectratio,BoxLength,ParticleiCentroid_old,ParticlejCentroid_old,VectorX_i_old, VectorY_i_old, VectorZ_i_old, VectorX_j_old, VectorY_j_old, VectorZ_j_old,shape);
                 ijEnergy_old = result_old.first;
                 error_flag = result_old.second;
                 //cout << "error_flag_old " << error_flag << endl;
                 if (error_flag) {
                    Loop = false;
                    EarlyTermination = true;
                    break;
                 }
                 ParticleSideLength1 = Size_i * AASigma;
                 ParticleSideLength2 = Size_j * AASigma;
                 ParticleHeight1 = ParticleSideLength1*aspectratio;
                 ParticleHeight2 = ParticleSideLength2*aspectratio;

                 // particle centroid in new coordinate system, with seed particle as the centroid
                 Vector3 ParticleCentroid = {Rx[i]-seedX,Ry[i]-seedY,Rz[i]-seedZ};
                 apply_periodic_boundary(ParticleCentroid,BoxLength,ParticleCentroid);
                 Vector3 ParticleVectorX, ParticleVectorY, ParticleVectorZ;
                 for (int k = 0; k < 3; ++k) {
                    ParticleVectorX[k] = VectorX[i][k]; // x component of orientation vector
                    ParticleVectorY[k] = VectorY[i][k]; // y component of orientation vector
                    ParticleVectorZ[k] = VectorZ[i][k]; // z component of orientation vector
                    //ProbFile << "ParticleVectorX[" << k << "] " << ParticleVectorX[k] << "ParticleVectorY[" << k << "] " << ParticleVectorY[k]<< "ParticleVectorZ[" << k << "] " << ParticleVectorZ[k] << endl;
                 }
                 // rotate each particle centroid and vectors about the centroid of the seed
                 Vector3 tempresult1, tempresult2;
                 Vector3 RotatedParticleCentroid, VectorX_i_new, VectorY_i_new, VectorZ_i_new;
                 MultiplyMatrixVector(VirtualMoveRotationX, ParticleCentroid, tempresult1);
                 MultiplyMatrixVector(VirtualMoveRotationY, tempresult1, tempresult2);
                 MultiplyMatrixVector(VirtualMoveRotationZ, tempresult2, RotatedParticleCentroid);

                 MultiplyMatrixVector(VirtualMoveRotationX, ParticleVectorX, tempresult1);
                 MultiplyMatrixVector(VirtualMoveRotationY, tempresult1, tempresult2);
                 MultiplyMatrixVector(VirtualMoveRotationZ, tempresult2, VectorX_i_new);

                 MultiplyMatrixVector(VirtualMoveRotationX, ParticleVectorY, tempresult1);
                 MultiplyMatrixVector(VirtualMoveRotationY, tempresult1, tempresult2);
                 MultiplyMatrixVector(VirtualMoveRotationZ, tempresult2, VectorY_i_new);

                 MultiplyMatrixVector(VirtualMoveRotationX, ParticleVectorZ, tempresult1);
                 MultiplyMatrixVector(VirtualMoveRotationY, tempresult1, tempresult2);
                 MultiplyMatrixVector(VirtualMoveRotationZ, tempresult2, VectorZ_i_new);

                 // shift the particle coordinates back to the original coordinate system
                 double Rx_i_new = RotatedParticleCentroid[0] + seedX;
                 double Ry_i_new = RotatedParticleCentroid[1] + seedY;
                 double Rz_i_new = RotatedParticleCentroid[2] + seedZ;

                 // new position of i after virtual move
                 Rx_i_new = Rx_i_new + VirtualMoveX;
                 Ry_i_new = Ry_i_new + VirtualMoveY;
                 Rz_i_new = Rz_i_new + VirtualMoveZ;

                 // apply periodic boundary conditions to positions
                 Vector3 R_i_new = {Rx_i_new,Ry_i_new,Rz_i_new};
                 apply_periodic_boundary(R_i_new, BoxLength,R_i_new);
                 Rx_i_new = R_i_new[0];
                 Ry_i_new = R_i_new[1];
                 Rz_i_new = R_i_new[2];

                 // relatvie distance between new i and old j
                 double Rx_ij_new = Rx_i_new - Rx_j_old;
                 double Ry_ij_new = Ry_i_new - Ry_j_old;
                 double Rz_ij_new = Rz_i_new - Rz_j_old;

                 // smallest relative distance between new i and old j considering periodic boundary conditions
                 Rx_ij_new = Rx_ij_new - BoxLength * round(Rx_ij_new/BoxLength);
                 Ry_ij_new = Ry_ij_new - BoxLength * round(Ry_ij_new/BoxLength);
                 Rz_ij_new = Rz_ij_new - BoxLength * round(Rz_ij_new/BoxLength);

                 // new i as reference, move j to the image where ij distance is minimal
                 Rx_j_old = Rx_i_new - Rx_ij_new;
                 Ry_j_old = Ry_i_new - Ry_ij_new;
                 Rz_j_old = Rz_i_new - Rz_ij_new;

                 Vector3 ParticleiCentroid_new = {Rx_i_new, Ry_i_new, Rz_i_new};
                 Vector3 ParticlejCentroid_minimal = {Rx_j_old,Ry_j_old,Rz_j_old};

                 // check energy between i and j after virtual move
                 if (CheckOverlapForTwoParticles(i,j,BoxLength,AASigma,ParticleSideLength1,ParticleSideLength2,ParticleHeight1,ParticleHeight2,ParticleiCentroid_new,ParticlejCentroid_minimal,VectorX_i_new,VectorY_i_new,VectorZ_i_new,VectorX_j_old,VectorY_j_old,VectorZ_j_old,shape)) {
                    pair<double, bool> result_new = EnergyBetweenTwoParticles(model,AASigma,CGSigma,CGEpsilon,cutoffCG,Hamaker,AtomDensity,Size_i,Size_j,aspectratio,BoxLength,ParticleiCentroid_new, ParticlejCentroid_old,VectorX_i_new, VectorY_i_new, VectorZ_i_new, VectorX_j_old, VectorY_j_old, VectorZ_j_old,shape);
                    ijEnergy_new = result_new.first;
                    error_flag = result_new.second;
                    overLapFlag = true;
                 } else {
                    pair<double, bool> result_new = EnergyBetweenTwoParticles(model,AASigma,CGSigma,CGEpsilon,cutoffCG,Hamaker,AtomDensity,Size_i,Size_j,aspectratio,BoxLength,ParticleiCentroid_new, ParticlejCentroid_old,VectorX_i_new, VectorY_i_new, VectorZ_i_new, VectorX_j_old, VectorY_j_old, VectorZ_j_old,shape);
                    ijEnergy_new = result_new.first;
                    error_flag = result_new.second;
                 }
                 //cout << "error_flag_new " << error_flag << endl;
                 if (error_flag) {
                    Loop = false;
                    EarlyTermination = true;
                    break;
                 }
                 try {
                    // Calculate virtualMoveProbability, considering the possibility of OverflowError
                    virtualMoveProbability = max(0.0, 1.0 - exp((ijEnergy_old - ijEnergy_new) / kBT));
                 } catch (const std::overflow_error& e) {
                    virtualMoveProbability = 0.0;
                 }
                 //ProbFile << "ijEnergy_old " << ijEnergy_old << endl;
                 //ProbFile << "ijEnergy_new " << ijEnergy_new << endl;
                 //ProbFile << "virtualMoveProbability " << virtualMoveProbability << endl;
                 // check acceptance probability of this virtual move
                 double rand_num = ((double)rand()/(RAND_MAX));
                 if (virtualMoveProbability >= rand_num) {
                    //ProbFile << "acceptedVirtualMoveProbability " << virtualMoveProbability << endl;
                    // add j to TemperaryjList for later use
                    TemperaryjList.push_back(j);
                    //ProbFile << "TemperaryjList j " << j << endl;
                    // calculate reverse move acceptance probability
                    // particle centroid in new coordinate system, with seed particle as the centroid
                    ParticleCentroid[0] = Rx[i] - seedX;
                    ParticleCentroid[1] = Ry[i] - seedY;
                    ParticleCentroid[2] = Rz[i] - seedZ;
                    apply_periodic_boundary(ParticleCentroid, BoxLength, ParticleCentroid);
                    for (int k = 0; k < 3; ++k) {
                       ParticleVectorX[k] = VectorX[i][k]; // x component of orientation vector
                       ParticleVectorY[k] = VectorY[i][k]; // y component of orientation vector
                       ParticleVectorZ[k] = VectorZ[i][k]; // z component of orientation vector
                    }
                    MultiplyMatrixVector(ReverseMoveRotationX, ParticleCentroid, tempresult1);
                    MultiplyMatrixVector(ReverseMoveRotationY, tempresult1, tempresult2);
                    MultiplyMatrixVector(ReverseMoveRotationZ, tempresult2, RotatedParticleCentroid);

                    // shift the particle coordinates back to the original coordinate system
                    double Rx_i_reverse = RotatedParticleCentroid[0] + seedX;
                    double Ry_i_reverse = RotatedParticleCentroid[1] + seedY;
                    double Rz_i_reverse = RotatedParticleCentroid[2] + seedZ;

                    // new position of i after reverse move
                    Rx_i_reverse = Rx_i_reverse - VirtualMoveX;
                    Ry_i_reverse = Ry_i_reverse - VirtualMoveY;
                    Rz_i_reverse = Rz_i_reverse - VirtualMoveZ;

                    Vector3 VectorX_i_reverse, VectorY_i_reverse, VectorZ_i_reverse;
                    MultiplyMatrixVector(ReverseMoveRotationX, VectorX_i_old, tempresult1);
                    MultiplyMatrixVector(ReverseMoveRotationY, tempresult1, tempresult2);
                    MultiplyMatrixVector(ReverseMoveRotationZ, tempresult2, VectorX_i_reverse);

                    MultiplyMatrixVector(ReverseMoveRotationX, VectorY_i_old, tempresult1);
                    MultiplyMatrixVector(ReverseMoveRotationY, tempresult1, tempresult2);
                    MultiplyMatrixVector(ReverseMoveRotationZ, tempresult2, VectorY_i_reverse);

                    MultiplyMatrixVector(ReverseMoveRotationX, VectorZ_i_old, tempresult1);
                    MultiplyMatrixVector(ReverseMoveRotationY, tempresult1, tempresult2);
                    MultiplyMatrixVector(ReverseMoveRotationZ, tempresult2, VectorZ_i_reverse);

                    // pickup the central image for the new location of i after reverse move
                    Rx_i_reverse = Rx_i_reverse - BoxLength*round(Rx_i_reverse/BoxLength);
                    Ry_i_reverse = Ry_i_reverse - BoxLength*round(Ry_i_reverse/BoxLength);
                    Rz_i_reverse = Rz_i_reverse - BoxLength*round(Rz_i_reverse/BoxLength);

                    // relative distance between reversed i and old j
                    double Rx_ij_reverse = Rx_i_reverse - Rx_j_old;
                    double Ry_ij_reverse = Ry_i_reverse - Ry_j_old;
                    double Rz_ij_reverse = Rz_i_reverse - Rz_j_old;

                    Rx_ij_reverse = Rx_ij_reverse - BoxLength * round(Rx_ij_reverse/BoxLength);
                    Ry_ij_reverse = Ry_ij_reverse - BoxLength * round(Ry_ij_reverse/BoxLength);
                    Rz_ij_reverse = Rz_ij_reverse - BoxLength * round(Rz_ij_reverse/BoxLength);

                    // reversed i as reference, move j to the image where ij distance is minimal
                    Rx_j_old = Rx_i_reverse - Rx_ij_reverse;
                    Ry_j_old = Ry_i_reverse - Ry_ij_reverse;
                    Rz_j_old = Rz_i_reverse - Rz_ij_reverse;

                    ParticleSideLength1 = Sizes[i] * AASigma;
                    ParticleSideLength2 = Sizes[j] * AASigma;
                    ParticleHeight1 = ParticleSideLength1*aspectratio;
                    ParticleHeight2 = ParticleSideLength2*aspectratio;

                    Vector3 ParticleiCentroid_reverse = {Rx_i_reverse,Ry_i_reverse,Rz_i_reverse};
                    // check energy between i and j after reverse move
                    if (CheckOverlapForTwoParticles(i,j,BoxLength,AASigma,ParticleSideLength1,ParticleSideLength2,ParticleHeight1,ParticleHeight2,ParticleiCentroid_reverse,ParticlejCentroid_old,VectorX_i_reverse,VectorY_i_reverse,VectorZ_i_reverse, VectorX_j_old, VectorY_j_old, VectorZ_j_old,shape)) {
                       pair<double, bool> result_reverse = EnergyBetweenTwoParticles(model,AASigma,CGSigma,CGEpsilon,cutoffCG,Hamaker,AtomDensity,Size_i,Size_j,aspectratio,BoxLength,ParticleiCentroid_reverse,ParticlejCentroid_old,VectorX_i_reverse,VectorY_i_reverse,VectorZ_i_reverse, VectorX_j_old, VectorY_j_old, VectorZ_j_old,shape);
                       ijEnergy_reverse = result_reverse.first;
                       error_flag = result_reverse.second;
                       //ProbFile << i << " and " << j << " overlap in reverse move" << endl;
                    } else {
                       pair<double, bool> result_reverse = EnergyBetweenTwoParticles(model,AASigma,CGSigma,CGEpsilon,cutoffCG,Hamaker,AtomDensity,Size_i,Size_j,aspectratio,BoxLength,ParticleiCentroid_reverse,ParticlejCentroid_old,VectorX_i_reverse,VectorY_i_reverse,VectorZ_i_reverse, VectorX_j_old, VectorY_j_old, VectorZ_j_old,shape);
                       ijEnergy_reverse = result_reverse.first;
                       error_flag = result_reverse.second;
                    }
                    // calculate reverse move acceptance probability
                    if (error_flag) {
                       Loop = false;
                       EarlyTermination = true;
                       break;
                    }
                    try {
                       reverseMoveProbability = max(0.0, 1.0 - exp((ijEnergy_old - ijEnergy_reverse) / kBT));
                    } catch (const std::overflow_error& e) {
                       reverseMoveProbability = 0.0;
                    }
                    reverseMoveProbabilityList.push_back(reverseMoveProbability);
                    forwardMoveProbabilityList.push_back(virtualMoveProbability);
                    //ProbFile << "virtualMoveProbability " << virtualMoveProbability << endl;
                    //ProbFile << "reverseMoveProbability " << reverseMoveProbability << endl;
                 } else {
                    double qij = 1-virtualMoveProbability;
                    unacceptedMoveProbabilityList.push_back(qij);
                    //ProbFile << "unacceptedMoveProbability " << qij << endl;
                    // calculate reverse move acceptance probability
                    ParticleCentroid[0] = Rx[i] - seedX;
                    ParticleCentroid[1] = Ry[i] - seedY;
                    ParticleCentroid[2] = Rz[i] - seedZ;
                    apply_periodic_boundary(ParticleCentroid, BoxLength, ParticleCentroid);
                    MultiplyMatrixVector(ReverseMoveRotationX, ParticleCentroid, tempresult1);
                    MultiplyMatrixVector(ReverseMoveRotationY, tempresult1, tempresult2);
                    MultiplyMatrixVector(ReverseMoveRotationZ, tempresult2, RotatedParticleCentroid);

                    // shift the particle coordinates back to the original coordinate system
                    double Rx_i_reverse = RotatedParticleCentroid[0] + seedX;
                    double Ry_i_reverse = RotatedParticleCentroid[1] + seedY;
                    double Rz_i_reverse = RotatedParticleCentroid[2] + seedZ;

                    // new position of i after reverse move
                    Rx_i_reverse = Rx_i_reverse - VirtualMoveX;
                    Ry_i_reverse = Ry_i_reverse - VirtualMoveY;
                    Rz_i_reverse = Rz_i_reverse - VirtualMoveZ;

                    Vector3 VectorX_i_reverse, VectorY_i_reverse, VectorZ_i_reverse;
                    MultiplyMatrixVector(ReverseMoveRotationX, VectorX_i_old, tempresult1);
                    MultiplyMatrixVector(ReverseMoveRotationY, tempresult1, tempresult2);
                    MultiplyMatrixVector(ReverseMoveRotationZ, tempresult2, VectorX_i_reverse);

                    MultiplyMatrixVector(ReverseMoveRotationX, VectorY_i_old, tempresult1);
                    MultiplyMatrixVector(ReverseMoveRotationY, tempresult1, tempresult2);
                    MultiplyMatrixVector(ReverseMoveRotationZ, tempresult2, VectorY_i_reverse);

                    MultiplyMatrixVector(ReverseMoveRotationX, VectorZ_i_old, tempresult1);
                    MultiplyMatrixVector(ReverseMoveRotationY, tempresult1, tempresult2);
                    MultiplyMatrixVector(ReverseMoveRotationZ, tempresult2, VectorZ_i_reverse);

                    // pickup the central image for the new location of i after reverse move
                    Rx_i_reverse = Rx_i_reverse - BoxLength*round(Rx_i_reverse/BoxLength);
                    Ry_i_reverse = Ry_i_reverse - BoxLength*round(Ry_i_reverse/BoxLength);
                    Rz_i_reverse = Rz_i_reverse - BoxLength*round(Rz_i_reverse/BoxLength);

                    // relative distance between reversed i and old j
                    double Rx_ij_reverse = Rx_i_reverse - Rx_j_old;
                    double Ry_ij_reverse = Ry_i_reverse - Ry_j_old;
                    double Rz_ij_reverse = Rz_i_reverse - Rz_j_old;

                    Rx_ij_reverse = Rx_ij_reverse - BoxLength * round(Rx_ij_reverse/BoxLength);
                    Ry_ij_reverse = Ry_ij_reverse - BoxLength * round(Ry_ij_reverse/BoxLength);
                    Rz_ij_reverse = Rz_ij_reverse - BoxLength * round(Rz_ij_reverse/BoxLength);

                    // reversed i as reference, move j to the image where ij distance is minimal
                    Rx_j_old = Rx_i_reverse - Rx_ij_reverse;
                    Ry_j_old = Ry_i_reverse - Ry_ij_reverse;
                    Rz_j_old = Rz_i_reverse - Rz_ij_reverse;

                    ParticleSideLength1 = Sizes[i] * AASigma;
                    ParticleSideLength2 = Sizes[j] * AASigma;
                    ParticleHeight1 = ParticleSideLength1*aspectratio;
                    ParticleHeight2 = ParticleSideLength2*aspectratio;
                    Vector3 ParticleiCentroid_reverse = {Rx_i_reverse,Ry_i_reverse,Rz_i_reverse};
                    // check energy between i and j after reverse move
                    if (CheckOverlapForTwoParticles(i,j,BoxLength,AASigma,ParticleSideLength1,ParticleSideLength2,ParticleHeight1,ParticleHeight2,ParticleiCentroid_reverse,ParticlejCentroid_old,VectorX_i_reverse,VectorY_i_reverse,VectorZ_i_reverse, VectorX_j_old, VectorY_j_old, VectorZ_j_old,shape)) {
                       pair<double, bool> result_reverse = EnergyBetweenTwoParticles(model,AASigma,CGSigma,CGEpsilon,cutoffCG,Hamaker,AtomDensity,Size_i,Size_j,aspectratio,BoxLength,ParticleiCentroid_reverse,ParticlejCentroid_old,VectorX_i_reverse,VectorY_i_reverse,VectorZ_i_reverse, VectorX_j_old, VectorY_j_old, VectorZ_j_old,shape);
                       ijEnergy_reverse = result_reverse.first;
                       error_flag = result_reverse.second;
                       //ProbFile << i << " and " << j << " overlap in reverse move" << endl;
                    } else {
                       pair<double, bool> result_reverse = EnergyBetweenTwoParticles(model,AASigma,CGSigma,CGEpsilon,cutoffCG,Hamaker,AtomDensity,Size_i,Size_j,aspectratio,BoxLength,ParticleiCentroid_reverse,ParticlejCentroid_old,VectorX_i_reverse,VectorY_i_reverse,VectorZ_i_reverse, VectorX_j_old, VectorY_j_old, VectorZ_j_old,shape);
                       ijEnergy_reverse = result_reverse.first;
                       error_flag = result_reverse.second;
                    }
                    if (error_flag) {
                       Loop = false;
                       EarlyTermination = true;
                       break;
                    }
                    // calculate reverse move acceptance probability
                    try {
                       reverseMoveProbability = max(0.0, 1.0 - exp((ijEnergy_old - ijEnergy_reverse) / kBT));
                    } catch (const std::overflow_error& e) {
                       reverseMoveProbability = 0.0;
                    }
                    double qji = 1-reverseMoveProbability;
                    unacceptedReverseMoveProbabilityList.push_back(qji);
                    linkConnectionList.push_back(j);
                    //ProbFile << "unacceptedReverseMoveProbability " << qji << endl;
                 }
              }
           }
           overallReverseProbability = 1.0;
           overallMoveProbability = 1.0;
           overallUnacceptedMoveProbability = 1.0;
           overallUnacceptedReverseProbability = 1.0;

           if (reverseMoveProbabilityList.size() > 0) {
              for (double reverseProbability : reverseMoveProbabilityList) {
                 overallReverseProbability *= reverseProbability;
                 //ProbFile << "overallReverseProbability " << overallReverseProbability << endl;
              }
           }

           if (unacceptedReverseMoveProbabilityList.size() > 0) {
              for (size_t index = 0; index < unacceptedReverseMoveProbabilityList.size(); ++index) {
                 if (find(clusterParticleList.begin(), clusterParticleList.end(), linkConnectionList[index]) != clusterParticleList.end()) {
                    // Unformed link internal to the cluster
                    double unacceptedReverseProbability = unacceptedReverseMoveProbabilityList[index];
                    overallUnacceptedReverseProbability *= unacceptedReverseProbability;
                    //ProbFile << "overallUnacceptedReverseProbability " << overallUnacceptedReverseProbability << endl;
                 }
              }
           }

           if (forwardMoveProbabilityList.size() > 0) {
              for (double forwardMoveProbability : forwardMoveProbabilityList) {
                 overallMoveProbability *= forwardMoveProbability;
                 //ProbFile << "overallMoveProbability " << overallMoveProbability << endl;
              }
           }

           if (unacceptedMoveProbabilityList.size() > 0) {
              for (size_t unacceptedindex = 0; unacceptedindex < unacceptedMoveProbabilityList.size(); ++unacceptedindex) {
                  if (find(clusterParticleList.begin(), clusterParticleList.end(), linkConnectionList[unacceptedindex]) != clusterParticleList.end()) {
                     // Unformed link internal to the cluster
                     double unacceptedMoveProbability = unacceptedMoveProbabilityList[unacceptedindex];
                     overallUnacceptedMoveProbability *= unacceptedMoveProbability;
                     //ProbFile << "overallUnacceptedMoveProbability " << overallUnacceptedMoveProbability << endl;
                  }
              }
           }

           if (TemperaryjList.size() == 0) {
              if (ToBeRecruitedParticlesList.size() == 0) {
                 Loop = false; // End recruiting particles to cluster
              } else {
                 // Pick a particle from ToBeRecruitedParticlesList as the new i
                 int rand_index = ((double)rand()/(RAND_MAX)) * (ToBeRecruitedParticlesList.size() - 1);
                 i = ToBeRecruitedParticlesList[rand_index];

                 // Add the new i to clusterParticleList
                 if (find(clusterParticleList.begin(), clusterParticleList.end(), i) == clusterParticleList.end()) {
                    clusterParticleList.push_back(i);
                    //ProbFile << "clusterParticleList: ";
                    //for (int particleID : clusterParticleList) {
                    //   ProbFile << particleID << " ";
                    //}
                    //ProbFile << " " << endl;
                 }

                 // Delete the new i from ParticlesNotInClusterList
                 ParticlesNotInClusterList.erase(remove(ParticlesNotInClusterList.begin(), ParticlesNotInClusterList.end(), i), ParticlesNotInClusterList.end());
                 // Delete the new i from ToBeRecruitedParticlesList
                 ToBeRecruitedParticlesList.erase(remove(ToBeRecruitedParticlesList.begin(), ToBeRecruitedParticlesList.end(), i), ToBeRecruitedParticlesList.end());
              }
           } else { // There are some j(s) that link to i
              // Pick a particle from TemperaryjList as the new i
              int rand_ind = ((double)rand()/(RAND_MAX)) * (TemperaryjList.size() - 1);
              i = TemperaryjList[rand_ind];

              // Delete the new i from TemperaryjList
              TemperaryjList.erase(remove(TemperaryjList.begin(), TemperaryjList.end(), i), TemperaryjList.end());

              // Add the new i to clusterParticleList
             if (find(clusterParticleList.begin(), clusterParticleList.end(), i) == clusterParticleList.end()) {
                clusterParticleList.push_back(i);
                //ProbFile << "clusterParticleList: ";
                //for (int particleID : clusterParticleList) {
                //   ProbFile << particleID << " ";
                //}
                //ProbFile << " " << endl;
             }

             // Delete the new i from ParticlesNotInClusterList
             ParticlesNotInClusterList.erase(remove(ParticlesNotInClusterList.begin(), ParticlesNotInClusterList.end(), i), ParticlesNotInClusterList.end());

             // Remove the new i from ToBeRecruitedParticlesList if it exists there
             ToBeRecruitedParticlesList.erase(remove(ToBeRecruitedParticlesList.begin(), ToBeRecruitedParticlesList.end(), i), ToBeRecruitedParticlesList.end());

             // Add the rest particles (if any) in TemperaryjList to ToBeRecruitedParticlesList
             if (TemperaryjList.size() > 0) {
                for (int Temperaryj : TemperaryjList) {
                   if (find(ToBeRecruitedParticlesList.begin(), ToBeRecruitedParticlesList.end(), Temperaryj) == ToBeRecruitedParticlesList.end()) {
                      ToBeRecruitedParticlesList.push_back(Temperaryj);
                   }
                }
             }
             // Make TemperaryjList empty
             TemperaryjList.clear();
           }
        }
        //////////////////////////  Building Cluster Ends  /////////////////
        // now the cluster is built
        if (EarlyTermination == false) {
           double clusterEnergy_old = 0.0;
           double clusterEnergy_new = 0.0;
           vector<array<double, 3>> particlesPositionList;
           array<double, 3> Ri;
           vector<int> pairList;
           double currentParticleOldEnergy,currentParticleNewEnergy;;
           bool overLapFlag, overLapFlagNew, overLapFlagOld;
           int clusterSize = clusterParticleList.size();
           for (int particleID : clusterParticleList) {
              i = particleID;
              pairList.push_back(i);
              //ProbFile << "pairList values: ";
              //for (int m : pairList) {
              //   ProbFile << m << " ";
              //}
              //ProbFile << std::endl;
              double Size_i = Sizes[i];

              //cout << "Rx old: ";
              //for (double m : Rx) {
              //   cout << m << endl;
              //}
              //cout << "Ry old: ";
              //for (double m : Ry) {
              //   cout << m << endl;
              //}
              //cout << "Rz old: ";
              //for (double m : Rz) {
              //   cout << m << endl;
              //}

              // old position and orientation of this particle before virtual move
              double Rx_i_old = Rx[i];
              double Ry_i_old = Ry[i];
              double Rz_i_old = Rz[i];
              //ProbFile <<"i " << i << " Rx_i_old " << Rx_i_old << " Ry_i_old " << Ry_i_old << " Rz_i_old " << Rz_i_old << endl;
              Ri[0] = Rx_i_old;
              Ri[1] = Ry_i_old;
              Ri[2] = Rz_i_old;
              Vector3 VectorX_i_old,VectorY_i_old,VectorZ_i_old;
              //for (int k = 0; k < 3; ++k) {
              //   VectorX_i_old[k] = VectorX[i][k]; // x component of orientation vector
              //   VectorY_i_old[k] = VectorY[i][k]; // y component of orientation vector
              //   VectorZ_i_old[k] = VectorZ[i][k]; // z component of orientation vector
              //}
              // energy of cluster before cluster move
              pair<double, bool> result_initial = OneParticleEnergy(model,i,Size_i,aspectratio,cutoffCG,Hamaker,AtomDensity,Rx_i_old,Ry_i_old,Rz_i_old,Sizes,Rx,Ry,Rz,VectorX,VectorY,VectorZ,NumberOfParticles,BoxLength,CGEpsilon,CGSigma,AASigma,globalminEnergy,shape,pairList);
              currentParticleOldEnergy = result_initial.first;
              overLapFlagOld = result_initial.second;
              clusterEnergy_old = clusterEnergy_old + currentParticleOldEnergy;
              //ProbFile << "currentParticleOldEnergy " << currentParticleOldEnergy << endl;
              //cout << "overLapFlagOld " << overLapFlagOld << endl;
              particlesPositionList.push_back(Ri);
           }

           // find centroid of the cluster
           Vector3 centroid;
           double ClusterTemp[clusterSize][3];
           calculate_cluster_centroid(particlesPositionList,BoxLength,centroid,ClusterTemp);
           Vector3 n_vector;
           double D_value = 1.0; // initialization
           // damping (based on Stokes' law, please check Whitelam & Geissler, 2008 paper to understand its implementation)
           if (isTranslation == true && MaxClusterSize > 1) { // damping coefficient calculation for clusters
              n_vector[0] = VirtualMoveX;
              n_vector[1] = VirtualMoveY;
              n_vector[2] = VirtualMoveZ;
              double n_vector_norm = sqrt(n_vector[0]*n_vector[0] + n_vector[1]*n_vector[1] + n_vector[2]*n_vector[2]);
              n_vector[0] /= n_vector_norm;
              n_vector[1] /= n_vector_norm;
              n_vector[2] /= n_vector_norm;
              vector<double> RiMinusRcCrossN_vectorNormSquareList;
              vector<double> RadiiOfCubeInCluster;

              for (int i : clusterParticleList) {
                 RadiiOfCubeInCluster.push_back(Sizes[i]*AASigma);
                 auto it = find(clusterParticleList.begin(), clusterParticleList.end(), i);
                 if (it != clusterParticleList.end()) {
                    size_t index = std::distance(clusterParticleList.begin(), it);
                    double Rx_i_old = ClusterTemp[index][0];
                    double Ry_i_old = ClusterTemp[index][1];
                    double Rz_i_old = ClusterTemp[index][2];
                    Vector3 Ri = {Rx_i_old, Ry_i_old, Rz_i_old};
                    Vector3 Rc = {centroid[0], centroid[1], centroid[2]};
                    Vector3 RiminusRc = {Ri[0] - Rc[0], Ri[1] - Rc[1], Ri[2] - Rc[2]};
                    Vector3 RiMinusRcCrossN_vector = {RiminusRc[1]*n_vector[2]-RiminusRc[2]*n_vector[1],RiminusRc[2]*n_vector[0]-RiminusRc[0]*n_vector[2],RiminusRc[0]*n_vector[1]-RiminusRc[1]*n_vector[0]};
                    double RiMinusRcCrossN_vectorNorm = sqrt(RiMinusRcCrossN_vector[0]*RiMinusRcCrossN_vector[0]+RiMinusRcCrossN_vector[1]*RiMinusRcCrossN_vector[1]+RiMinusRcCrossN_vector[2]*RiMinusRcCrossN_vector[2]);
                    double RiMinusRcCrossN_vectorNormSquare = RiMinusRcCrossN_vectorNorm * RiMinusRcCrossN_vectorNorm;
                    RiMinusRcCrossN_vectorNormSquareList.push_back(RiMinusRcCrossN_vectorNormSquare);
                 }
              }
              double sumRadiiOfCubeInCluster = accumulate(RadiiOfCubeInCluster.begin(), RadiiOfCubeInCluster.end(), 0.0);
              double AverageRadiusOfCubeInCluster = sumRadiiOfCubeInCluster / RadiiOfCubeInCluster.size();
              double sumRiMinusRcCrossN_vectorNormSquare = accumulate(RiMinusRcCrossN_vectorNormSquareList.begin(),RiMinusRcCrossN_vectorNormSquareList.end(),0.0);
              double RiMinusRcCrossN_vectorNormSquareAverage = sumRiMinusRcCrossN_vectorNormSquare / RiMinusRcCrossN_vectorNormSquareList.size();
              double EffectiveHydrodynamicRadius = sqrt(RiMinusRcCrossN_vectorNormSquareAverage) + AverageRadiusOfCubeInCluster;
              double Dt = AverageRadiusOfCubeInCluster / EffectiveHydrodynamicRadius; // translational damping
              D_value = Dt;
           } else if (isRotation == true && MaxClusterSize > 1) { // for rotation
              if (VirtualMoveAlpha != 0) {
                 n_vector[0] = 1.0;
                 n_vector[1] = 0.0;
                 n_vector[2] = 0.0;
              } else if (VirtualMoveBeta != 0) {
                 n_vector[0] = 0.0;
                 n_vector[1] = 1.0;
                 n_vector[2] = 0.0;
              } else if (VirtualMoveGamma != 0) {
                 n_vector[0] = 0.0;
                 n_vector[1] = 0.0;
                 n_vector[2] = 1.0;
              } else { // if the rotation is set to 0
                 n_vector[0] = 0.0;
                 n_vector[1] = 0.0;
                 n_vector[2] = 1.0;
              }
              vector<double> RiMinusRcCrossN_vectorNormSquareList;
              vector<double> RadiiOfCubeInCluster;

              for (int i : clusterParticleList) {
                 RadiiOfCubeInCluster.push_back(Sizes[i]*AASigma);
                 auto it = find(clusterParticleList.begin(), clusterParticleList.end(), i);
                 if (it != clusterParticleList.end()) {
                    size_t index = std::distance(clusterParticleList.begin(), it);
                    double Rx_i_old = ClusterTemp[index][0];
                    double Ry_i_old = ClusterTemp[index][1];
                    double Rz_i_old = ClusterTemp[index][2];
                    Vector3 Ri = {Rx_i_old, Ry_i_old, Rz_i_old};
                    Vector3 Rc = {centroid[0], centroid[1], centroid[2]};
                    Vector3 RiminusRc = {Ri[0] - Rc[0], Ri[1] - Rc[1], Ri[2] - Rc[2]};
                    Vector3 RiMinusRcCrossN_vector = {RiminusRc[1]*n_vector[2]-RiminusRc[2]*n_vector[1],RiminusRc[2]*n_vector[0]-RiminusRc[0]*n_vector[2],RiminusRc[0]*n_vector[1]-RiminusRc[1]*n_vector[0]};
                    double RiMinusRcCrossN_vectorNorm = sqrt(RiMinusRcCrossN_vector[0]*RiMinusRcCrossN_vector[0]+RiMinusRcCrossN_vector[1]*RiMinusRcCrossN_vector[1]+RiMinusRcCrossN_vector[2]*RiMinusRcCrossN_vector[2]);
                    double RiMinusRcCrossN_vectorNormSquare = RiMinusRcCrossN_vectorNorm * RiMinusRcCrossN_vectorNorm;
                    RiMinusRcCrossN_vectorNormSquareList.push_back(RiMinusRcCrossN_vectorNormSquare);
                 }
              }
              double sumRadiiOfCubeInCluster = accumulate(RadiiOfCubeInCluster.begin(), RadiiOfCubeInCluster.end(), 0.0);
              double AverageRadiusOfCubeInCluster = sumRadiiOfCubeInCluster / RadiiOfCubeInCluster.size();
              double sumRiMinusRcCrossN_vectorNormSquare = accumulate(RiMinusRcCrossN_vectorNormSquareList.begin(),RiMinusRcCrossN_vectorNormSquareList.end(),0.0);
              double RiMinusRcCrossN_vectorNormSquareAverage = sumRiMinusRcCrossN_vectorNormSquare / RiMinusRcCrossN_vectorNormSquareList.size();
              double EffectiveHydrodynamicRadius = sqrt(RiMinusRcCrossN_vectorNormSquareAverage) + AverageRadiusOfCubeInCluster;
              double Dr = pow(AverageRadiusOfCubeInCluster / EffectiveHydrodynamicRadius,3); // rotation damping
              D_value = Dr;
           }
           vector<double> Rx_All_temp;
           vector<double> Ry_All_temp;
           vector<double> Rz_All_temp;
           double VectorX_All_temp[NumberOfParticles][3];
           double VectorY_All_temp[NumberOfParticles][3];
           double VectorZ_All_temp[NumberOfParticles][3];
           // Here you temporarily store the positions of particles after the proposed move to be able to compute the energy after proposed move
           for (int k = 0; k < NumberOfParticles; ++k) {
              for (int m = 0; m < 3; ++m) {
                 VectorX_All_temp[k][m] = VectorX[k][m];
                 VectorY_All_temp[k][m] = VectorY[k][m];
                 VectorZ_All_temp[k][m] = VectorZ[k][m];
              }
              Rx_All_temp.push_back(Rx[k]);
              Ry_All_temp.push_back(Ry[k]);
              Rz_All_temp.push_back(Rz[k]);
           }

           // Now, you perform the selected move on the particles within the formed cluster
           Vector3 ParticleCentroid;
           for (int particleID : clusterParticleList) {
              i = particleID;
              ParticleCentroid[0] = Rx[i] - seedX;
              ParticleCentroid[1] = Ry[i] - seedY;
              ParticleCentroid[2] = Rz[i] - seedZ;
              apply_periodic_boundary(ParticleCentroid, BoxLength, ParticleCentroid);

              Vector3 ParticleVectorX, ParticleVectorY, ParticleVectorZ;
              for (int k = 0; k < 3; ++k) {
                 ParticleVectorX[k] = VectorX[i][k]; // x component of orientation vector
                 ParticleVectorY[k] = VectorY[i][k]; // y component of orientation vector
                 ParticleVectorZ[k] = VectorZ[i][k]; // z component of orientation vector
              }

              Matrix3x3 VirtualMoveRotationX = {{1, 0, 0}, {0, cos(VirtualMoveAlpha), -sin(VirtualMoveAlpha)}, {0, sin(VirtualMoveAlpha), cos(VirtualMoveAlpha)}};
              Matrix3x3 VirtualMoveRotationY = {{cos(VirtualMoveBeta), 0, sin(VirtualMoveBeta)}, {0, 1, 0}, {-sin(VirtualMoveBeta), 0, cos(VirtualMoveBeta)}};
              Matrix3x3 VirtualMoveRotationZ = {{cos(VirtualMoveGamma), -sin(VirtualMoveGamma), 0}, {sin(VirtualMoveGamma), cos(VirtualMoveGamma), 0}, {0, 0, 1}};

              Vector3 tempresult1, tempresult2;
              Vector3 RotatedParticleCentroid, VectorX_i_new, VectorY_i_new, VectorZ_i_new;
              MultiplyMatrixVector(VirtualMoveRotationX, ParticleCentroid, tempresult1);
              MultiplyMatrixVector(VirtualMoveRotationY, tempresult1, tempresult2);
              MultiplyMatrixVector(VirtualMoveRotationZ, tempresult2, RotatedParticleCentroid);

              MultiplyMatrixVector(VirtualMoveRotationX, ParticleVectorX, tempresult1);
              MultiplyMatrixVector(VirtualMoveRotationY, tempresult1, tempresult2);
              MultiplyMatrixVector(VirtualMoveRotationZ, tempresult2, VectorX_i_new);

              MultiplyMatrixVector(VirtualMoveRotationX, ParticleVectorY, tempresult1);
              MultiplyMatrixVector(VirtualMoveRotationY, tempresult1, tempresult2);
              MultiplyMatrixVector(VirtualMoveRotationZ, tempresult2, VectorY_i_new);

              MultiplyMatrixVector(VirtualMoveRotationX, ParticleVectorZ, tempresult1);
              MultiplyMatrixVector(VirtualMoveRotationY, tempresult1, tempresult2);
              MultiplyMatrixVector(VirtualMoveRotationZ, tempresult2, VectorZ_i_new);

              double Rx_i_new = RotatedParticleCentroid[0] + seedX;
              double Ry_i_new = RotatedParticleCentroid[1] + seedY;
              double Rz_i_new = RotatedParticleCentroid[2] + seedZ;

              Rx_i_new = Rx_i_new + VirtualMoveX;
              Ry_i_new = Ry_i_new + VirtualMoveY;
              Rz_i_new = Rz_i_new + VirtualMoveZ;
              Vector3 R_i_new = {Rx_i_new,Ry_i_new,Rz_i_new};
              apply_periodic_boundary(R_i_new, BoxLength,R_i_new);
              Rx_i_new = R_i_new[0];
              Ry_i_new = R_i_new[1];
              Rz_i_new = R_i_new[2];

              Rx_All_temp[i] = Rx_i_new;
              Ry_All_temp[i] = Ry_i_new;
              Rz_All_temp[i] = Rz_i_new;
              for (int k = 0; k < 3; ++k) {
                 //cout << "Here" << endl;
                 //cout << VectorX_i_new[0] << " " << VectorX_i_new[1] << " " << VectorX_i_new[2] << endl;
                 VectorX_All_temp[i][k] = VectorX_i_new[k]; // x component of orientation vector
                 VectorY_All_temp[i][k] = VectorY_i_new[k]; // y component of orientation vector
                 VectorZ_All_temp[i][k] = VectorZ_i_new[k]; // z component of orientation vector
              }
           }
           //cout << "Rx temp: ";
           //for (int m : Rx_All_temp) {
           //    cout << m << endl;
           //}
           //cout << "Ry temp: ";
           //for (int m : Ry_All_temp) {
           //   cout << m << endl;
           //}
           //cout << "Rz temp: ";
           //for (int m : Rz_All_temp) {
           //   cout << m << endl;
           //}
           bool clusterOverlapAfterMove = false;
           pairList.clear();
           // energy of cluster after cluster move
           for (int particleID : clusterParticleList) {
              i = particleID;
              pairList.push_back(i);
              double Size_i = Sizes[i];
              double Rx_i_new = Rx_All_temp[i];
              double Ry_i_new = Ry_All_temp[i];
              double Rz_i_new = Rz_All_temp[i];
              Vector3 VectorX_i_new, VectorY_i_new, VectorZ_i_new;
              for (int k = 0; k < 3; ++k) {
                 VectorX_i_new[k] = VectorX_All_temp[i][k]; // x component of orientation vector
                 VectorY_i_new[k] = VectorY_All_temp[i][k]; // y component of orientation vector
                 VectorZ_i_new[k] = VectorZ_All_temp[i][k]; // z component of orientation vector
              }
              pair<double,bool> result_second = OneParticleEnergy(model,i,Size_i,aspectratio,cutoffCG,Hamaker,AtomDensity,Rx_i_new,Ry_i_new,Rz_i_new,Sizes,Rx_All_temp,Ry_All_temp,Rz_All_temp,VectorX_All_temp,VectorY_All_temp,VectorZ_All_temp,NumberOfParticles,BoxLength,CGEpsilon,CGSigma,AASigma,globalminEnergy,shape,pairList);
              currentParticleNewEnergy = result_second.first;
              overLapFlagNew = result_second.second;
              //cout << "overLapFlagNew " << overLapFlagNew << endl;
              if (overLapFlagNew || overLapFlagOld) {
                 clusterOverlapAfterMove = true;
              }
              clusterEnergy_new = clusterEnergy_new + currentParticleNewEnergy;
              //ProbFile << "clusterEnergy_new " << clusterEnergy_new << endl;
           }
           double ClusterEnergyChange = clusterEnergy_new - clusterEnergy_old; // delta_V, difference between new and old total energy
           //ProbFile << "clusterEnergy_new " << clusterEnergy_new << endl;
           //ProbFile << "clusterEnergy_old " << clusterEnergy_old << endl;
           //ProbFile << "ClusterEnergyChange " << ClusterEnergyChange << endl;
           //cout << "ClusterEnergyChange " << ClusterEnergyChange << endl;
           //double PotentialEnergy = TotalEnergy(Sizes,Rx,Ry,Rz,VectorX,VectorY,VectorZ,NumberOfParticles,BoxLength,CGSigma,CGEpsilon,AASigma,model,cutoffCG, Hamaker, AtomDensity);
           //ProbFile << "RealEnergy Old " << PotentialEnergy << endl;
           //cout << "RealEnergy Old" << PotentialEnergy << endl;
           //cout << "OldEnergy " << SystemPotentialEnergy << endl;
           Attempt += 1;
           clusterSize = clusterParticleList.size();
           double RandomNumber = ((double)rand()/(RAND_MAX));
           //ProbFile << "RandomNumber " << RandomNumber << endl;
           double BoltzmannFactor = exp(-ClusterEnergyChange / kBT);
           //ProbFile << "BoltzmannFactor " << BoltzmannFactor << endl;
           if (MaxClusterSize == 1) {
              SingleMoveNumber += 1;
           } else {
              ClusterMoveNumber += 1;
           }
           //ProbFile << "overallReverseProbability " << overallReverseProbability << endl;
           //ProbFile << "overallUnacceptedReverseProbability " << overallUnacceptedReverseProbability << endl;
           //ProbFile << "overallMoveProbability " << overallMoveProbability << endl;
           //ProbFile << "overallUnacceptedMoveProbability " << overallUnacceptedMoveProbability << endl;
           //ProbFile << "D_value " << D_value << endl;

           double overallReverse = overallReverseProbability * overallUnacceptedReverseProbability;
           double overallMove = overallMoveProbability * overallUnacceptedMoveProbability;
           double W_accept = D_value * min(1.0, overallReverse / overallMove); // acceptance criterion of entire cluster's move
           if (MaxClusterSize == 1) { // if single-move is selected, you can calculate this by Metropolis criterion
              W_accept = min(1.0,BoltzmannFactor);
              //ProbFile << "W_accept " << W_accept << endl;
              //ProbFile << "Single move case" << endl;
           }
           //ProbFile << "W_accept " << W_accept << endl;
           //ProbFile << "clusterOverlapAfterMove " << clusterOverlapAfterMove << endl;

           if (RandomNumber <= W_accept && !clusterOverlapAfterMove) {  // accept this move if it satisfies Boltzmann factor criteria, make sure no overlap
              SystemPotentialEnergy += ClusterEnergyChange;
              for (int particleID : clusterParticleList) {
                 i = particleID;
                 Rx[i] = Rx_All_temp[i];
                 Ry[i] = Ry_All_temp[i];
                 Rz[i] = Rz_All_temp[i];
                 for (int k = 0; k < 3; ++k) {
                    VectorX[i][k] = VectorX_All_temp[i][k];
                    VectorY[i][k] = VectorY_All_temp[i][k];
                    VectorZ[i][k] = VectorZ_All_temp[i][k];
                 }
              }
              Accept += 1;
              if (MaxClusterSize == 1) {
                 SingleMoveAccept += 1;
                 //ProbFile << "Single move accepted!!!" << endl;
              } else {
                 ClusterMoveAccept += 1;
                 //ProbFile << "Cluster move accepted!!!" << endl;
              }
              //ProbFile << "Accept Move!!!" << endl;
           } else {
              //ProbFile << "Reject cluster Move~~~" << endl;
           }
        }
        //cout << "Rx new: ";
        //for (double m : Rx) {
        //   cout << m << endl;
        //}
        //cout << "Ry new: ";
        //for (double m : Ry) {
        //   cout << m << endl;
        //}
        //cout << "Rz new: ";
        //for (double m : Rz) {
        //  cout << m << endl;
        //}
        lastRestartStep = writeRestart(RestartFileInterval,lastRestartStep,Sizes,Rx,Ry,Rz,VectorX,VectorY,VectorZ,Step);
        if (Step > EquilibrationSteps) {
          writeEnergy(EnergyOutputInterval, Step, SystemPotentialEnergy, EnergyFile);
          writeTrajectory(TrajectoryInterval,Style,Step,NumberOfParticles,AASigma,CGParticleSize,aspectratio,CGSigma,BoxLengthHalf,shape,Sizes,Rx,Ry,Rz,VectorX,VectorY,VectorZ);
        }
        //ProbFile << "New Energy " << SystemPotentialEnergy << endl;
        // cout <<  "New Energy " << SystemPotentialEnergy << endl;
        //double PotentialEnergy = TotalEnergy(Sizes,Rx,Ry,Rz,VectorX,VectorY,VectorZ,NumberOfParticles,BoxLength,CGSigma,CGEpsilon,AASigma,model,cutoffCG, Hamaker,AtomDensity);
        //ProbFile << "Real Energy New " << PotentialEnergy << endl;
        // cout <<  "Real Energy New" << PotentialEnergy << endl;
        double accept_ratio = Accept/Attempt;
        double cluster_accept_ratio = ClusterMoveAccept / ClusterMoveNumber;
        double single_accept_ratio = SingleMoveAccept / SingleMoveNumber;
        //ProbFile << "Fraction of accepting moves " << accept_ratio << endl;
        if (ClusterMoveNumber > 0) {
           //ProbFile << "Fraction of accepting cluster move " << cluster_accept_ratio << endl;
        }
        if (SingleMoveNumber>0) {
           //ProbFile << "Fraction of accepting single move " << single_accept_ratio << endl;
        }
        cout << "Step-EquilibrationSteps " << Step-EquilibrationSteps << endl;
        if (((Step-EquilibrationSteps) % volumeMoveInterval == 0 || oldvolfrac == 1) && volumeexpansion && Step > EquilibrationSteps) { //first compression step is very early as it starts from 100% volume fraction (no need to run long equilibration for volume fraction=1 case, as there is no space to move)
           double newvolfrac = oldvolfrac - 0.05;
           double volrate = pow(oldvolfrac/newvolfrac,1.0/3.0);
           //cout << "oldvolfrac " << oldvolfrac << endl;
           //cout << "rate " << volrate << endl;
           BoxLength = BoxLength * volrate;
           BoxLengthHalf = BoxLength / 2.0;
           VolumeMove(NumberOfParticles, BoxLength, AASigma, aspectratio, volrate,SizeRatioDictionary,shape,Sizes,Rx,Ry,Rz,VectorX,VectorY,VectorZ);
           oldvolfrac = newvolfrac;
        }

        /*
        vector<bool> inCluster(NumberOfParticles, false);
        for (int k = 0; k < NumberOfParticles; ++k) {
           for (int m = k+1; m < NumberOfParticles; ++m) {
              double pos_x = Rx[k]-Rx[m];
              double pos_y = Ry[k]-Ry[m];
              double pos_z = Rz[k]-Rz[m];
              pos_x = pos_x-BoxLength*round(pos_x/BoxLength);
              pos_y = pos_y-BoxLength*round(pos_y/BoxLength);
              pos_z = pos_z-BoxLength*round(pos_z/BoxLength);
              double dist = sqrt(pos_x*pos_x+pos_y*pos_y+pos_z*pos_z);
              //cout << "dist " << dist << endl;
              int bin_index = static_cast<int>(dist / bin_width);
              if (bin_index < num_bins) {
                 double bin_center = (bin_index+0.5)*bin_width;
                 //cout << "bin_center " << bin_center << endl;
                 //cout << "test" << normalized_bin_counts[bin_index] << endl;
                 normalized_bin_counts[bin_index] += 2.0/(normalization_constant*bin_center*bin_center);
              }
              if (dist < (CubeSideAASize*AASigma*sqrt(3))) {
                 mean_agg += 1;
              }
           }
        }*/

        /*MeanAggFile << mean_agg << endl;
        for (int k = 0; k < num_bins; ++k) {
           mean_bin_counts[k] += normalized_bin_counts[k];
           mean_bin_counts[k] = mean_bin_counts[k]/Step;
           RDFFile << normalized_bin_counts[k] << " " ;
           MeanRDFFile << mean_bin_counts[k] << " ";
        }
        RDFFile << endl;
        MeanRDFFile << endl;*/

    }
    LAMMPSTrajectoryFile.close();
    EnergyFile.close();
    return 0;
}



