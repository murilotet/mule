/* 
 * File:   Hierarchical4DEncoder.h
 * Author: murilo
 *
 * Created on December 28, 2017, 11:41 AM
 */
#include "Block4D.h"
#include "ABACoder.h"
#include "ProbabilityModel.h"

#ifndef HIERARCHICAL4DENCODER_H
#define HIERARCHICAL4DENCODER_H

//#define MAX_DEPH_CONDICIONING 9
#define BITPLANE_BYPASS -1
#define BITPLANE_BYPASS_FLAGS -1
#define SYMBOL_PROBABILITY_MODEL_INDEX 1
#define SEGMENTATION_PROB_MODEL_INDEX 32
#define NUMBER_OF_MODELS 161

class Hierarchical4DEncoder {
public:
    Block4D mSubbandLF;   
    ABACoder mEntropyCoder;
    ProbabilityModel *mPmodel;
    ProbabilityModel *mOptimizationPmodel;
    int mSuperiorBitPlane, mInferiorBitPlane;
    int mSegmentationFlagProbabilityModelIndex;
    int mSymbolProbabilityModelIndex;
    int mPreSegmentation;
    char *mSegmentationTreeCodeBuffer;
    long int mSegmentationTreeCodeBufferSize;
    Hierarchical4DEncoder(void);
    ~Hierarchical4DEncoder(void);
    void StartEncoder(FILE *outputFilePointer);
    void RestartProbabilisticModel(void);
    void EncodeBlock(int position_t, int position_s, int position_v, int position_u, int length_t, int length_s, int length_v, int length_u, int bitplane);
    void EncodeCoefficient(int coefficient, int bitplane);
    void EncodeSegmentationFlag(int flag, int bitplane);
    void EncodePartitionFlag(int flag);
    void EncodeInteger(int integerValue, int precision);
    void EncodeAll(double lambda, int inferiorBitPlane);
    void EncodeSubblock(double lambda);
    double RdOptimizeHexadecaTree(int position_t, int position_s, int position_v, int position_u, int length_t, int length_s, int length_v, int length_u, double lambda, int bitplane, char **codeString, double &signalEnergy);
    void RdEncodeHexadecatree(int position_t, int position_s, int position_v, int position_u, int length_t, int length_s, int length_v, int length_u, int bitplane, int &flagIndex);
    void DoneEncoding(void);
    void SetDimension(int length_t, int length_s, int length_v, int length_u);
    int OptimumBitplane(double lambda);
    void LoadOptimizerState(void);
    void GetOptimizerProbabilisticModelState(ProbabilityModel **state);
    void SetOptimizerProbabilisticModelState(ProbabilityModel *state);
    void DeleteProbabilisticModelState(ProbabilityModel *state);
};
#endif /* HIERARCHICAL4DENCODER_H */

