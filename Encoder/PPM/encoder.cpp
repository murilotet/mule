// No subbands version (encoding is performed on a block-by-block basis)
#include "LightField.h"
#include "Hierarchical4DEncoder.h"
#include "MultiscaleTransform.h"
#include "TransformPartition.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

enum ExtensionType {ZERO_PADDING, REPEAT_LAST, CYCLIC, EXTENDED_DCT, TRANSFORM_DOMAIN_ZERO_PADDING, MEAN_VALUE, NONE};
void ExtendDCT(Matrix &extendedDCT, ExtensionType extensionMethod, int transformLength, int extensionLength);
void ExtendBlock4D(Block4D &extendedblock, ExtensionType extensionMethod, int extensionLength, char direction);
void RGB2YCbCr_BT601(Block4D &Y, Block4D &Cb, Block4D &Cr, Block4D const &R, Block4D const &G, Block4D const &B, int Scale);
void RGB2YCoCg(Block4D &Y, Block4D &Co, Block4D &Cg, Block4D const &R, Block4D const &G, Block4D const &B, int Scale);

void BigEndianUnsignedIntegerWrite(unsigned long int value, int precision, FILE *outputFilePointer);

class EncoderParameters {
public:
    double Lambda;
    int transformLength_t;
    int transformLength_s;
    int transformLength_v;
    int transformLength_u;
    int min_transformLength_t;
    int min_transformLength_s;
    int min_transformLength_v;
    int min_transformLength_u;
    int inputNumberOfVerticalViews;
    int inputNumberOfHorizontalViews;
    int firstHorizontalViewNumber;
    int firstVerticalViewNumber;
    char isLenslet13x13;
    char inputDirectory[1024];
    char outputFileName[1024];
    char configurationFileName[1024];        
    ExtensionType extensionMethod;    
    double transform_scale_t;
    double transform_scale_s;
    double transform_scale_v;
    double transform_scale_u;    
    int colorTransformType;
    int verbosity;
    void ReadConfigurationFile(char *parametersFileName); 
    void DisplayConfiguration(void);
};
void EncoderParameters :: ReadConfigurationFile(char *parametersFileName) {
    char command[128];
    FILE *parametersFilepointer;
    
    if((parametersFilepointer = fopen(parametersFileName, "r")) == NULL) {
        printf("ERROR: unable to open configuration file %s\n", parametersFileName);
        exit(0);
    } 
    
    fscanf(parametersFilepointer, "%s", command);
    while(feof(parametersFilepointer) == 0) {
        if(strcmp(command, "-lambda") == 0) {
            fscanf(parametersFilepointer, "%lf", &Lambda);
        }
        if(strcmp(command, "-l") == 0) {
            fscanf(parametersFilepointer, "%d", &transformLength_u);
            fscanf(parametersFilepointer, "%d", &transformLength_v);
        }
        if(strcmp(command, "-u") == 0) fscanf(parametersFilepointer, "%d", &transformLength_u);
        if(strcmp(command, "-v") == 0) fscanf(parametersFilepointer, "%d", &transformLength_v);
        if(strcmp(command, "-k") == 0) {
            fscanf(parametersFilepointer, "%d", &transformLength_t);
            fscanf(parametersFilepointer, "%d", &transformLength_s);
        }
        if(strcmp(command, "-s") == 0) fscanf(parametersFilepointer, "%d", &transformLength_s);
        if(strcmp(command, "-t") == 0) fscanf(parametersFilepointer, "%d", &transformLength_t);
        if(strcmp(command, "-lf") == 0) fscanf(parametersFilepointer, "%s", &inputDirectory);
        if(strcmp(command, "-o") == 0) fscanf(parametersFilepointer, "%s", &outputFileName);
        if(strcmp(command, "-nv") == 0) fscanf(parametersFilepointer, "%d", &inputNumberOfVerticalViews);
        if(strcmp(command, "-nh") == 0) fscanf(parametersFilepointer, "%d", &inputNumberOfHorizontalViews);
        if(strcmp(command, "-off_h") == 0) fscanf(parametersFilepointer, "%d", &firstHorizontalViewNumber);
        if(strcmp(command, "-off_v") == 0) fscanf(parametersFilepointer, "%d", &firstVerticalViewNumber);
        if(strcmp(command, "-lenslet13x13") == 0) isLenslet13x13 = 1;
        if(strcmp(command, "-extension_repeat") == 0) extensionMethod = REPEAT_LAST;        
        if(strcmp(command, "-extension_dct_ext") == 0) extensionMethod = EXTENDED_DCT;        
        if(strcmp(command, "-extension_dct_zero") == 0) extensionMethod = TRANSFORM_DOMAIN_ZERO_PADDING;        
        if(strcmp(command, "-extension_zero") == 0) extensionMethod = ZERO_PADDING;        
        if(strcmp(command, "-extension_mean") == 0) extensionMethod = MEAN_VALUE;        
        if(strcmp(command, "-extension_none") == 0) extensionMethod = NONE;        
        if(strcmp(command, "-extension_cyclic") == 0) extensionMethod = CYCLIC;        
        if(strcmp(command, "-t_scale") == 0) fscanf(parametersFilepointer, "%lf", &transform_scale_t);        
        if(strcmp(command, "-s_scale") == 0) fscanf(parametersFilepointer, "%lf", &transform_scale_s);        
        if(strcmp(command, "-v_scale") == 0) fscanf(parametersFilepointer, "%lf", &transform_scale_v);        
        if(strcmp(command, "-u_scale") == 0) fscanf(parametersFilepointer, "%lf", &transform_scale_u);        
        if(strcmp(command, "-v0") == 0) verbosity = 0;        
        if(strcmp(command, "-min_u") == 0) fscanf(parametersFilepointer, "%d", &min_transformLength_u);
        if(strcmp(command, "-min_v") == 0) fscanf(parametersFilepointer, "%d", &min_transformLength_v);
        if(strcmp(command, "-min_s") == 0) fscanf(parametersFilepointer, "%d", &min_transformLength_s);
        if(strcmp(command, "-min_t") == 0) fscanf(parametersFilepointer, "%d", &min_transformLength_t);        
        if(strcmp(command, "-bt601") == 0) colorTransformType = 0;        
        if(strcmp(command, "-ycocg") == 0) colorTransformType = 1;        
        fscanf(parametersFilepointer, "%s", command);
    }
}

void EncoderParameters :: DisplayConfiguration(void) {
    printf("Lambda = %f\n", Lambda);
    printf("transformLength_t = %d\n", transformLength_t);
    printf("transformLength_s = %d\n", transformLength_s);
    printf("transformLength_v = %d\n", transformLength_v);
    printf("transformLength_u = %d\n", transformLength_u);
    printf("min_transformLength_t = %d\n", min_transformLength_t);
    printf("min_transformLength_s = %d\n", min_transformLength_s);
    printf("min_transformLength_v = %d\n", min_transformLength_v);
    printf("min_transformLength_u = %d\n", min_transformLength_u);
    printf("inputNumberOfVerticalViews = %d\n", inputNumberOfVerticalViews);
    printf("inputNumberOfHorizontalViews = %d\n", inputNumberOfHorizontalViews);
    printf("firstHorizontalViewNumber = %d\n", firstHorizontalViewNumber);
    printf("firstVerticalViewNumber = %d\n", firstVerticalViewNumber);
    printf("isLenslet13x13 = %d\n", isLenslet13x13);
    printf("inputDirectory = %s\n", inputDirectory);
    printf("outputFileName = %s\n", outputFileName);
    printf("configurationFileName = %s\n", configurationFileName);        
    printf("extensionMethod = %d\n", extensionMethod);    
    printf("transform_scale_t = %f\n", transform_scale_t);
    printf("transform_scale_s = %f\n", transform_scale_s);
    printf("transform_scale_v = %f\n", transform_scale_v);
    printf("transform_scale_u = %f\n", transform_scale_u);    
    printf("colorTransformType = %d\n", colorTransformType);    
    printf("verbosity = %d\n", verbosity);    
}

int main(int argc, char **argv) {
  
    EncoderParameters par;
    
    //Set default parameters
    par.Lambda = 1024;
    par.transformLength_t=13;
    par.transformLength_s=13;
    par.transformLength_v=15;
    par.transformLength_u=15;
    par.min_transformLength_t=4;
    par.min_transformLength_s=4;
    par.min_transformLength_v=4;
    par.min_transformLength_u=4;
    par.inputNumberOfVerticalViews=13;
    par.inputNumberOfHorizontalViews=13;
    par.firstHorizontalViewNumber=0;
    par.firstVerticalViewNumber=0;
    par.isLenslet13x13=0;
        
    strcpy(par.inputDirectory, "./");
    strcpy(par.outputFileName, "out.comp");
    strcpy(par.configurationFileName, "");
    
    par.extensionMethod=REPEAT_LAST;
    
    par.transform_scale_t=1.0;
    par.transform_scale_s=1.0;
    par.transform_scale_v=1.0;
    par.transform_scale_u=1.0;
    
    par.colorTransformType=0;
    
    par.verbosity=1;
    
    for(int n = 0; n < argc; n++) {
        if(strcmp(argv[n], "-cf") == 0) {
            strcpy(par.configurationFileName, argv[n+1]);
            par.ReadConfigurationFile(par.configurationFileName);
        }
    }
    for(int n = 0; n < argc; n++) {
        if(strcmp(argv[n], "-lambda") == 0) {
            par.Lambda = atof(argv[n+1]);
        }
        if(strcmp(argv[n], "-l") == 0) {
            par.transformLength_u = atoi(argv[n+1]);
            par.transformLength_v = atoi(argv[n+1]);
        }
        if(strcmp(argv[n], "-u") == 0) par.transformLength_u = atoi(argv[n+1]);
        if(strcmp(argv[n], "-v") == 0) par.transformLength_v = atoi(argv[n+1]);
        if(strcmp(argv[n], "-k") == 0) {
            par.transformLength_t = atoi(argv[n+1]);
            par.transformLength_s = atoi(argv[n+1]);
        }
        if(strcmp(argv[n], "-s") == 0) par.transformLength_s = atoi(argv[n+1]);
        if(strcmp(argv[n], "-t") == 0) par.transformLength_t = atoi(argv[n+1]);
        if(strcmp(argv[n], "-lf") == 0) strcpy(par.inputDirectory, argv[n+1]);
        if(strcmp(argv[n], "-o") == 0) strcpy(par.outputFileName, argv[n+1]);
        if(strcmp(argv[n], "-nv") == 0) par.inputNumberOfVerticalViews = atoi(argv[n+1]);
        if(strcmp(argv[n], "-nh") == 0) par.inputNumberOfHorizontalViews = atoi(argv[n+1]);
        if(strcmp(argv[n], "-off_h") == 0) par.firstHorizontalViewNumber = atoi(argv[n+1]);
        if(strcmp(argv[n], "-off_v") == 0) par.firstVerticalViewNumber = atoi(argv[n+1]);
        if(strcmp(argv[n], "-lenslet13x13") == 0) par.isLenslet13x13 = 1;
        if(strcmp(argv[n], "-extension_repeat") == 0) par.extensionMethod = REPEAT_LAST;        
        if(strcmp(argv[n], "-extension_dct_ext") == 0) par.extensionMethod = EXTENDED_DCT;        
        if(strcmp(argv[n], "-extension_dct_zero") == 0) par.extensionMethod = TRANSFORM_DOMAIN_ZERO_PADDING;        
        if(strcmp(argv[n], "-extension_zero") == 0) par.extensionMethod = ZERO_PADDING;        
        if(strcmp(argv[n], "-extension_mean") == 0) par.extensionMethod = MEAN_VALUE;        
        if(strcmp(argv[n], "-extension_none") == 0) par.extensionMethod = NONE;        
        if(strcmp(argv[n], "-extension_cyclic") == 0) par.extensionMethod = CYCLIC;        
        if(strcmp(argv[n], "-t_scale") == 0) par.transform_scale_t = atof(argv[n+1]);        
        if(strcmp(argv[n], "-s_scale") == 0) par.transform_scale_s = atof(argv[n+1]);        
        if(strcmp(argv[n], "-v_scale") == 0) par.transform_scale_v = atof(argv[n+1]);        
        if(strcmp(argv[n], "-u_scale") == 0) par.transform_scale_u = atof(argv[n+1]);        
        if(strcmp(argv[n], "-v0") == 0) par.verbosity = 0;        
        if(strcmp(argv[n], "-min_u") == 0) par.min_transformLength_u = atoi(argv[n+1]);
        if(strcmp(argv[n], "-min_v") == 0) par.min_transformLength_v = atoi(argv[n+1]);
        if(strcmp(argv[n], "-min_s") == 0) par.min_transformLength_s = atoi(argv[n+1]);
        if(strcmp(argv[n], "-min_t") == 0) par.min_transformLength_t = atoi(argv[n+1]);
        if(strcmp(argv[n], "-bt601") == 0) par.colorTransformType = 0;        
        if(strcmp(argv[n], "-ycocg") == 0) par.colorTransformType = 1;        
    }
    
    if(par.min_transformLength_t >  par.transformLength_t)
        par.min_transformLength_t = par.transformLength_t;
    if(par.min_transformLength_s > par.transformLength_s)
        par.min_transformLength_s = par.transformLength_s;
    if(par.min_transformLength_v >  par.transformLength_v)
        par.min_transformLength_v = par.transformLength_v;
    if(par.min_transformLength_u > par.transformLength_u)
        par.min_transformLength_u = par.transformLength_u;
    
    par.DisplayConfiguration();

    int numberOfCacheViewLines = par.transformLength_v;
    LightField inputLF(par.inputNumberOfVerticalViews,par.inputNumberOfHorizontalViews,numberOfCacheViewLines);
    
    Block4D lfBlock, rBlock, gBlock, bBlock, yBlock, cbBlock, crBlock;
    
    Hierarchical4DEncoder hdt;
    
    TransformPartition tp;
    tp.mPartitionData.SetDimension(par.transformLength_t,par.transformLength_s,par.transformLength_v,par.transformLength_u);
    tp.mlength_t_min = par.min_transformLength_t;
    tp.mlength_s_min = par.min_transformLength_s;
    tp.mlength_v_min = par.min_transformLength_v;
    tp.mlength_u_min = par.min_transformLength_u;
   
    lfBlock.SetDimension(par.transformLength_t,par.transformLength_s,par.transformLength_v,par.transformLength_u);
    rBlock.SetDimension(par.transformLength_t,par.transformLength_s,par.transformLength_v,par.transformLength_u);
    gBlock.SetDimension(par.transformLength_t,par.transformLength_s,par.transformLength_v,par.transformLength_u);
    bBlock.SetDimension(par.transformLength_t,par.transformLength_s,par.transformLength_v,par.transformLength_u);
    yBlock.SetDimension(par.transformLength_t,par.transformLength_s,par.transformLength_v,par.transformLength_u);
    cbBlock.SetDimension(par.transformLength_t,par.transformLength_s,par.transformLength_v,par.transformLength_u);
    crBlock.SetDimension(par.transformLength_t,par.transformLength_s,par.transformLength_v,par.transformLength_u);
   
    inputLF.mVerticalViewNumberOffset = par.firstVerticalViewNumber;
    inputLF.mHorizontalViewNumberOffset = par.firstHorizontalViewNumber;
    inputLF.OpenLightFieldPPM(par.inputDirectory, ".ppm", par.inputNumberOfVerticalViews, par.inputNumberOfHorizontalViews, 3, 3, 'r');
       
    int extensionLength_t = par.inputNumberOfVerticalViews%par.transformLength_t;
    int extensionLength_s = par.inputNumberOfHorizontalViews%par.transformLength_s;
    int extensionLength_v = inputLF.mNumberOfViewLines%par.transformLength_v;
    int extensionLength_u = inputLF.mNumberOfViewColumns%par.transformLength_u;   

    MultiscaleTransform DCTarray;
    
    DCTarray.SetDimension(par.transformLength_t, par.transformLength_s, par.transformLength_v, par.transformLength_u);
    DCTarray.mTransformGain_t = par.transform_scale_t;
    DCTarray.mTransformGain_s = par.transform_scale_s;
    DCTarray.mTransformGain_v = par.transform_scale_v;
    DCTarray.mTransformGain_u = par.transform_scale_u;
       
    FILE *outputFileNamePointer;
    if((outputFileNamePointer = fopen(par.outputFileName, "wb")) == NULL) {
        printf("Error: input file %s not found\n", par.outputFileName);
        exit(0);
    }
    //writes the superior bit plane value
    BigEndianUnsignedIntegerWrite(hdt.mSuperiorBitPlane, 2, outputFileNamePointer);
    
    //writes the maximum transform sizes
    BigEndianUnsignedIntegerWrite(par.transformLength_t, 2, outputFileNamePointer);
    BigEndianUnsignedIntegerWrite(par.transformLength_s, 2, outputFileNamePointer);
    BigEndianUnsignedIntegerWrite(par.transformLength_v, 2, outputFileNamePointer);
    BigEndianUnsignedIntegerWrite(par.transformLength_u, 2, outputFileNamePointer);
        
    //writes the number of views of the lightfield
    BigEndianUnsignedIntegerWrite(inputLF.mNumberOfVerticalViews, 2, outputFileNamePointer);
    BigEndianUnsignedIntegerWrite(inputLF.mNumberOfHorizontalViews, 2, outputFileNamePointer);
    
    //writes the number of lines and columns of each view
    BigEndianUnsignedIntegerWrite(inputLF.mNumberOfViewLines, 2, outputFileNamePointer);
    BigEndianUnsignedIntegerWrite(inputLF.mNumberOfViewColumns, 2, outputFileNamePointer);
    
    //writes the bit precision of each component of the pixels of the views
    BigEndianUnsignedIntegerWrite(inputLF.mPGMScale, 2, outputFileNamePointer);
    
    hdt.StartEncoder(outputFileNamePointer);
    
    for(int verticalView = 0; verticalView < inputLF.mNumberOfVerticalViews; verticalView += par.transformLength_t) {
        for(int horizontalView = 0; horizontalView < inputLF.mNumberOfHorizontalViews; horizontalView += par.transformLength_s) {
            for(int viewLine = 0; viewLine < inputLF.mNumberOfViewLines; viewLine += par.transformLength_v) {
                for(int viewColumn = 0; viewColumn < inputLF.mNumberOfViewColumns; viewColumn += par.transformLength_u) {
                    if(par.verbosity > 0)
                        printf("transforming the 4D block at position (%d %d %d %d)\n", verticalView, horizontalView, viewLine, viewColumn);
                    rBlock.Zeros();
                    gBlock.Zeros();
                    bBlock.Zeros();
                    inputLF.ReadBlock4DfromLightField(&rBlock, verticalView, horizontalView, viewLine, viewColumn, 0);
                    inputLF.ReadBlock4DfromLightField(&gBlock, verticalView, horizontalView, viewLine, viewColumn, 1);
                    inputLF.ReadBlock4DfromLightField(&bBlock, verticalView, horizontalView, viewLine, viewColumn, 2);
                    if(par.isLenslet13x13 == 1) {
                        if(verticalView == 0) {
                            if(horizontalView == 0) {
                                rBlock.Shift_UVPlane(2, 0, 0);
                                gBlock.Shift_UVPlane(2, 0, 0);
                                bBlock.Shift_UVPlane(2, 0, 0);
                            }
                            if((horizontalView + par.transformLength_s >= inputLF.mNumberOfHorizontalViews)&&(horizontalView <= inputLF.mNumberOfHorizontalViews)) { 
                                rBlock.Shift_UVPlane(2, 0, inputLF.mNumberOfHorizontalViews-horizontalView-1);
                                gBlock.Shift_UVPlane(2, 0, inputLF.mNumberOfHorizontalViews-horizontalView-1);
                                bBlock.Shift_UVPlane(2, 0, inputLF.mNumberOfHorizontalViews-horizontalView-1);
                            }
                        }
                        if((verticalView + par.transformLength_t >= inputLF.mNumberOfVerticalViews)&&(verticalView <= inputLF.mNumberOfVerticalViews)) {
                            if(horizontalView == 0) {
                                rBlock.Shift_UVPlane(2, inputLF.mNumberOfVerticalViews-verticalView-1, 0);
                                gBlock.Shift_UVPlane(2, inputLF.mNumberOfVerticalViews-verticalView-1, 0);
                                bBlock.Shift_UVPlane(2, inputLF.mNumberOfVerticalViews-verticalView-1, 0);
                            }
                            if((horizontalView + par.transformLength_s >= inputLF.mNumberOfHorizontalViews)&&(horizontalView <= inputLF.mNumberOfHorizontalViews)) {
                                rBlock.Shift_UVPlane(2, inputLF.mNumberOfVerticalViews-verticalView-1, inputLF.mNumberOfHorizontalViews-horizontalView-1);
                                gBlock.Shift_UVPlane(2, inputLF.mNumberOfVerticalViews-verticalView-1, inputLF.mNumberOfHorizontalViews-horizontalView-1);
                                bBlock.Shift_UVPlane(2, inputLF.mNumberOfVerticalViews-verticalView-1, inputLF.mNumberOfHorizontalViews-horizontalView-1);
                            }
                        }
                    }
                    if(par.colorTransformType == 0)
                        RGB2YCbCr_BT601(yBlock, cbBlock, crBlock, rBlock, gBlock, bBlock, inputLF.mPGMScale);
                    if(par.colorTransformType == 1)
                        RGB2YCoCg(yBlock, cbBlock, crBlock, rBlock, gBlock, bBlock, inputLF.mPGMScale);
                    for(int spectralComponent = 0; spectralComponent < 3; spectralComponent++) {
                        printf("\nProcessing spectral component %d\n", spectralComponent);
                        if(spectralComponent == 0)
                            lfBlock.CopySubblockFrom(yBlock, 0, 0, 0, 0);
                        if(spectralComponent == 1)
                            lfBlock.CopySubblockFrom(cbBlock, 0, 0, 0, 0);
                        if(spectralComponent == 2)
                            lfBlock.CopySubblockFrom(crBlock, 0, 0, 0, 0);

                        lfBlock = lfBlock - (inputLF.mPGMScale+1)/2;
                                    
                        if(viewColumn + par.transformLength_u > inputLF.mNumberOfViewColumns)
                            ExtendBlock4D(lfBlock, par.extensionMethod, extensionLength_u, 'u');                                   
                        if(viewLine + par.transformLength_v > inputLF.mNumberOfViewLines) 
                            ExtendBlock4D(lfBlock, par.extensionMethod, extensionLength_v, 'v');                                  
                        if(horizontalView + par.transformLength_s > inputLF.mNumberOfHorizontalViews)
                            ExtendBlock4D(lfBlock, par.extensionMethod, extensionLength_s, 's');                                  
                        if(verticalView + par.transformLength_t > inputLF.mNumberOfVerticalViews)
                            ExtendBlock4D(lfBlock, par.extensionMethod, extensionLength_t, 't');                                                                      

                        hdt.RestartProbabilisticModel();
                        tp.RDoptimizeTransform(lfBlock, DCTarray, hdt, par.Lambda);
                        tp.EncodePartition(hdt, par.Lambda);

		    }            
                }
            }
        }
    }
    
    hdt.DoneEncoding();
    
    inputLF.CloseLightField();
    fclose(outputFileNamePointer);
    
    
}

void ExtendDCT(Matrix &extendedDCT, ExtensionType extensionMethod, int transformLength, int extensionLength) {
    
        extendedDCT.SetDimension(transformLength,transformLength);
        if(extensionMethod == EXTENDED_DCT) {
            Matrix C00, C01, C10, C11, INVC11;
            extendedDCT.DCT();
            C00.SetDimension(extensionLength,extensionLength);
            C01.SetDimension(extensionLength,transformLength-extensionLength);
            C10.SetDimension(transformLength-extensionLength,extensionLength);
            C11.SetDimension(transformLength-extensionLength,transformLength-extensionLength);
            INVC11.SetDimension(transformLength-extensionLength,transformLength-extensionLength);
            C00.CopyFrom(extendedDCT);
            C01.CopyFrom(extendedDCT,0,extensionLength);
            C10.CopyFrom(extendedDCT,extensionLength,0);
            C11.CopyFrom(extendedDCT,extensionLength,extensionLength);
            INVC11.Inverse(C11);
            INVC11.Multiply(-1.0);
            C10.PreMultiply(INVC11);    //C10 <- - INVC11 * C10
            C10.PreMultiply(C01);       //C10 <- - C01 * INVC11 * C10
            C00.Add(C10);               //C00 <- C00 - C01 * INVC11 * C10
            extendedDCT.Zeros();
            extendedDCT.CopyFrom(C00);
        }
        if(extensionMethod == REPEAT_LAST) {
            extendedDCT.DCT();
            extendedDCT.AccumulateFromColumn(extensionLength-1);
        }
        if(extensionMethod == TRANSFORM_DOMAIN_ZERO_PADDING) {
            Matrix C00;
            C00.SetDimension(extensionLength,extensionLength);
            C00.DCT();
            C00.Multiply((1.0*transformLength)/extensionLength);
            extendedDCT.Zeros();
            extendedDCT.CopyFrom(C00);
        }
        if(extensionMethod == ZERO_PADDING) {
            Matrix C00;
            C00.SetDimension(transformLength,extensionLength);
            extendedDCT.DCT();
            C00.CopyFrom(extendedDCT);
            extendedDCT.Zeros();
            extendedDCT.CopyFrom(C00);
        }
        if(extensionMethod == NONE) {
            extendedDCT.DCT();
        }
    
}

void ExtendBlock4D(Block4D &extendedBlock, ExtensionType extensionMethod, int extensionLength, char direction) {
    
    if(extensionMethod == REPEAT_LAST) {
        
        if(direction == 't') 
            extendedBlock.Extend_T(extensionLength-1);
        if(direction == 's') 
            extendedBlock.Extend_S(extensionLength-1);
        if(direction == 'v') 
            extendedBlock.Extend_V(extensionLength-1);
        if(direction == 'u') 
            extendedBlock.Extend_U(extensionLength-1);
        
    }
    if(extensionMethod == CYCLIC) {
        
        if(direction == 't') 
            extendedBlock.CopySubblockFrom(extendedBlock, 0, 0, 0, 0, extensionLength, 0, 0, 0);
        if(direction == 's') 
            extendedBlock.CopySubblockFrom(extendedBlock, 0, 0, 0, 0, 0, extensionLength, 0, 0);
        if(direction == 'v') 
            extendedBlock.CopySubblockFrom(extendedBlock, 0, 0, 0, 0, 0, 0, extensionLength, 0);
        if(direction == 'u') 
            extendedBlock.CopySubblockFrom(extendedBlock, 0, 0, 0, 0, 0, 0, 0, extensionLength);
        
    }
    if(extensionMethod == NONE) {
            
    }
}


void RGB2YCbCr_BT601(Block4D &Y, Block4D &Cb, Block4D &Cr, Block4D const &R, Block4D const &G, Block4D const &B, int Scale) {
   
    for(int n = 0; n < R.mlength_t*R.mlength_s*R.mlength_v*R.mlength_u; n++) {
        double pixel =  0.299 * R.mPixelData[n] + 0.587 * G.mPixelData[n] + 0.114 * B.mPixelData[n];
        Y.mPixelData[n] = round(pixel);
        pixel = -0.16875 * R.mPixelData[n] -0.33126 * G.mPixelData[n] + 0.5 * B.mPixelData[n];
        Cb.mPixelData[n] = round(pixel) + (Scale + 1)/2;
        pixel = 0.5 * R.mPixelData[n] -0.41869 * G.mPixelData[n] -0.08131  * B.mPixelData[n];
        Cr.mPixelData[n] = round(pixel) + (Scale + 1)/2;
    }
       
}

void RGB2YCoCg(Block4D &Y, Block4D &Co, Block4D &Cg, Block4D const &R, Block4D const &G, Block4D const &B, int Scale) {
    
    for(int n = 0; n < R.mlength_t*R.mlength_s*R.mlength_v*R.mlength_u; n++) {
        int t;
        Co.mPixelData[n] = R.mPixelData[n] - B.mPixelData[n];
        t = B.mPixelData[n] + (Co.mPixelData[n]>>1);
        Cg.mPixelData[n] = G.mPixelData[n] - t;
        Y.mPixelData[n] = t + (Cg.mPixelData[n]>>1);
        Co.mPixelData[n] += (Scale + 1)/2;
        Cg.mPixelData[n] += (Scale + 1)/2;
    }
        
}

void BigEndianUnsignedIntegerWrite(unsigned long int value, int precision, FILE *outputFilePointer) {

    if(precision > sizeof(long int)) {
        printf("Unsuported %d bytes integer precision\n", precision);
        exit(0);
    }
    
    unsigned char output_byte;
    for(int byte_index = precision-1; byte_index >= 0; byte_index--) {
        output_byte = (value >> (8*byte_index))&255;
        fwrite(&output_byte, 1, 1, outputFilePointer);
    }
    
}
