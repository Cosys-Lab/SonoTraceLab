
/*
 * Fast detection test if a ray intersects with a set of triangles using CUDA accelerated kernels based on the work of Möller and Trumbore (1997).
 *
 * Arguments: 'vertices' [(3 * 3 * numTriangles),type float32]
                These are interleaved, example:                 
                vertices=[T0_V0_X, T0_V0_Y, T0_V0_Z;
                         T0_V1_X, T0_V1_Y, T0_V1_Z;
                         T0_V2_X, T0_V2_Y, T0_V2_Z;
                         T1_V0_X, T1_V0_Y, T1_V0_Z;
                         T1_V1_X, T1_V1_Y, T1_V1_Z;
                         T1_V2_X, T1_V2_Y, T1_V2_Z;
                         ...]
                 with Ti: triangles, Vj:vertices 1 to 3 of triangle and X/Y/Z the 3 coordinates
              'raysFrom' [(3 * numRays),type float32]
                These are the starting points of each ray in X Y Z
              'directions' [(3 * numRays),type float32]
                These are the directions vectors of each ray in X Y Z
              
 * Returns:   'results' [numTriangles * numRays, type int32]
                These are the binary test results to see if there is an intersection.
                1 = intersection found, 0 = no intersection
                
 * Compile with 'mexcuda -v mollertrumbore.cu' (-v for extra details for debugging)
 * Requires CUDA toolkit and a C compiler.
 * Make sure to correctly set c++ compiler with 'mex -setup c++' and clicking on the link of the version you want if asked.
 * And make sure to set the CUDA enviroment variable correctly with
 * 'setenv('MW_NVCC_PATH','/usr/local/cuda-X/bin')' on Linux
 * 'setenv('MW_NVCC_PATH','C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\vX\bin')' on Windows
 * 
 * Original CUDA implementation of Möller and Trumbore is by Raymond Leung (2022), 
 * 'GPU implementation of a ray-surface intersection algorithm in CUDA',
 * arXiv e-print 2209.02878, 2022.
 * Source code available at: https://github.com/raymondleung8/gpu-ray-surface-intersection-in-cuda
 * and is made under a BSD 3 license.
 *
 * MATLAB implementation by Wouter Jansen, Cosys-Lab, University of Antwerp
 */

#define EPSILON 0.0000000001
#include <string>   

__device__ void subtract(const float *a, const float *b, float *out)
{
    out[0] = a[0] - b[0];
    out[1] = a[1] - b[1];
    out[2] = a[2] - b[2];
}

__device__ void dot(const float *a, const float *b, float &out)
{
    out = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

__device__ void cross(const float *a, const float *b, float *out)
{
    out[0] = a[1]*b[2] - a[2]*b[1];
    out[1] = a[2]*b[0] - a[0]*b[2];
    out[2] = a[0]*b[1] - a[1]*b[0];
}

__device__ float tolerance(const float *d, const float *eAB, const float *eAC)
{
    float scaling = sqrtf((d[0]*d[0] + d[1]*d[1] + d[2]*d[2]) *
                          (eAB[0]*eAB[0] + eAB[1]*eAB[1] + eAB[2]*eAB[2]) *
                          (eAC[0]*eAC[0] + eAC[1]*eAC[1] + eAC[2]*eAC[2]));
    return (scaling > 1)? scaling * EPSILON : EPSILON;
}

// Implement the Moller-Trumbore ray-triangle intersection algorithm
// - Ray model: R(t) = Q0 + t *(dir), where Q0 denote segment start point
// - and dir the direction of the segment
// - Point on triangle: T(u,v) = (1-u-v)*V0 + u*V1 + v*V2
//
__device__ int intersectMoller(
                const float *v0, const float *v1, const float *v2,
                const float *edge1, const float *edge2,
                const float *q0, const float *dir)
{
    float avec[3], bvec[3], tvec[3], t, u, v, det, inv_det;
    cross(dir, edge2, avec);
    dot(avec, edge1, det);
    float epsilon = tolerance(dir, edge1, edge2);
    if (det > epsilon) {
        subtract(q0, v0, tvec);
        dot(avec, tvec, u);
        if (u < 0 || u > det){
            return 0;
        }
        cross(tvec, edge1, bvec);
        dot(bvec, dir, v);
        if (v < 0 || u + v > det){
            return 0;
        }
    }
    else if (det < -epsilon) {
        subtract(q0, v0, tvec);
        dot(avec, tvec, u);
        if (u > 0 || u < det){
            return 0;
        }
        cross(tvec, edge1, bvec);
        dot(bvec, dir, v);
        if (v > 0 || u + v < det){
            return 0;
        }
    }
    else{
        return 0;
    }
    inv_det = 1.0 / det;
    dot(bvec, edge2, t);
    t *= inv_det;
    if (t < 0) {
        return 0;
    }
    else {
        return 1;
    }
}

__device__ void checkRayTriangleIntersection(const float* __restrict__ vertices,
                                             const float* __restrict__ raysFrom,
                                             const float* __restrict__ directions,
                                             int* __restrict__ results,
                                             int triangleIdx, int rayIdx,
                                             int numTriangles, int numRays)
{
    float triangleVerts[9], edge1[3], edge2[3];
    const float *v0 = &triangleVerts[0],
                *v1 = &triangleVerts[3],
                *v2 = &triangleVerts[6];
    for(int j = 0; j < 3; j++) { // loop by triangle v0, v1, v2
        for (int k = 0; k < 3; k++) { // loop by x, y and z
            triangleVerts[3*j + k] = vertices[9*(triangleIdx) + (3*j) + k];
        }
    }
    subtract(v1, v0, edge1);
    subtract(v2, v0, edge2);
   
    // Apply Moller-Trumbore ray-triangle intersection test
    const float *start = &raysFrom[3*rayIdx], *dir = &directions[3*rayIdx];
    if (intersectMoller(v0, v1, v2, edge1, edge2, start, dir)) {
        // printf("HIT! triangleIdx: %i rayIdx: %i\n", triangleIdx, rayIdx);
        int tempRayIdx = rayIdx / 32;
        int bitOffset = rayIdx % 32;
        atomicOr(&results[tempRayIdx * numTriangles + triangleIdx], 1 << bitOffset);
    }
}

__global__ void checkRayTriangleIntersectionKernel(const float* vertices,
                                                   const float* raysFrom,
                                                   const float* directions,
                                                   int* results,
                                                   int numTriangles, int numRays)
{
    int triangleIdx = blockIdx.x * blockDim.x + threadIdx.x;
    int rayIdx = blockIdx.y * blockDim.y + threadIdx.y;
    if(triangleIdx < numTriangles && rayIdx < numRays){
        checkRayTriangleIntersection(vertices, raysFrom, directions, results, triangleIdx, rayIdx, numTriangles, numRays);
    }
}

#include "mex.h"
#include "gpu/mxGPUArray.h"   

/*
 * Host code for CPU
 */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    /* Declare all variables.*/
    mxGPUArray const *deviceVertices;
    mxGPUArray const *deviceRaysFrom;
    mxGPUArray const *deviceDirections;
    mxGPUArray *deviceResults;
    float const *d_vertices;
    float const *d_raysFrom;
    float const *d_directions;
    int *d_results;
    int numTriangles;
    int numRays;

    /* Initialize the MathWorks GPU API. */
    mxInitGPU();

    /* Throw an error if the input are not a CPU arrays. */
    if ((mxIsGPUArray(prhs[0])) || (mxIsGPUArray(prhs[1])) || (mxIsGPUArray(prhs[2]))) {
        mexErrMsgIdAndTxt("parallel:gpu:mexGPUExample:InvalidInput", "The input matrices have to be normal CPU arrays, not GPUArrays.\n");
    }

    /* Throw an error if the input are not the correct datatype. */
    if ( mxGetClassID(prhs[0]) != mxSINGLE_CLASS) {
        mexErrMsgIdAndTxt("parallel:gpu:mexGPUExample:InvalidInput", "The vertices data matrix has to be of datatype 'single'.\n");
    }
    if ( mxGetClassID(prhs[1]) != mxSINGLE_CLASS) {
        mexErrMsgIdAndTxt("parallel:gpu:mexGPUExample:InvalidInput", "The ray-from data matrix has to be of datatype 'single'.\n");
    }
    if ( mxGetClassID(prhs[2]) != mxSINGLE_CLASS) {
        mexErrMsgIdAndTxt("parallel:gpu:mexGPUExample:InvalidInput", "The ray-to data matrix has to be of datatype 'single'.\n");
    }

    deviceVertices = mxGPUCreateFromMxArray(prhs[0]);
    deviceRaysFrom = mxGPUCreateFromMxArray(prhs[1]);
    deviceDirections = mxGPUCreateFromMxArray(prhs[2]);
    numTriangles = (int)((float)mxGPUGetDimensions(deviceVertices)[1]) / 3;
    numRays = mxGPUGetDimensions(deviceRaysFrom)[1];

    /* Extract a pointer to the input data on the device. */
    d_vertices = (float const *)(mxGPUGetDataReadOnly(deviceVertices));
    d_raysFrom = (float const *)(mxGPUGetDataReadOnly(deviceRaysFrom));
    d_directions = (float const *)(mxGPUGetDataReadOnly(deviceDirections));

    // printf("numTriangles:%i numRays:%i \n", numTriangles, numRays);

    /* Create a GPUArray to hold the result and get its underlying pointer. */
    int numRaysLogicalStorage = ceil((float)numRays / (float)32);
    mwSize dims[2] = {numTriangles, numRaysLogicalStorage};
    // printf("Storing results in array of %ix%i\n", numTriangles, numRaysLogicalStorage);
    deviceResults = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(deviceVertices),
                                        dims,
                                        mxINT32_CLASS,
                                        mxREAL,
                                        MX_GPU_INITIALIZE_VALUES);
    d_results = (int *)(mxGPUGetData(deviceResults));

    /* Execute the kernel. */
    dim3 threadsPerBlock(256, 4, 1);
    int gridX = (int)ceil(numTriangles / (threadsPerBlock.x*1.0));
    int gridY = (int)ceil(numRays / (threadsPerBlock.y*1.0));
    dim3 numBlocks(gridX, gridY, 1);
    checkRayTriangleIntersectionKernel<<<numBlocks, threadsPerBlock>>>(d_vertices, d_raysFrom, d_directions, d_results, numTriangles, numRays);

    /* Wrap the result up as a MATLAB gpuArray for return. */
    plhs[0] = mxGPUCreateMxArrayOnCPU(deviceResults);

    /*
     * The mxGPUArray pointers are host-side structures that refer to device
     * data. These must be destroyed before leaving the mex function.
     */
    mxGPUDestroyGPUArray(deviceVertices);
    mxGPUDestroyGPUArray(deviceRaysFrom);
    mxGPUDestroyGPUArray(deviceDirections);
    mxGPUDestroyGPUArray(deviceResults);
}
