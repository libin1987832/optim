fprintf('Compiling minFunc files...\n');
mex -compatibleArrayDims minFunc/lbfgsC.c
fprintf('Compiling projection files...\n');
mex -compatibleArrayDims project/projectRandom2C.c
mex -compatibleArrayDims -Iproject project/projectBlockL1.c project/oneProjectorCore.c project/heap.c
mex -compatibleArrayDims -Iproject project/projectBlockL2.c
fprintf('Compiling CRF files (requires a good compilier)...\n\n');
mex -compatibleArrayDims crfChain/crfChain_lossL2C.c

% UGM
fprintf('Compiling UGM files...\n');
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_makeEdgeVEC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Decode_ExactC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Infer_ExactC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Infer_ChainC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_makeClampedPotentialsC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Decode_ICMC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Decode_GraphCutC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Sample_GibbsC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Infer_MFC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Infer_LBPC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Decode_LBPC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Infer_TRBPC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Decode_TRBPC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_CRF_makePotentialsC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_CRF_PseudoNLLC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_LogConfigurationPotentialC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Decode_AlphaExpansionC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Decode_AlphaExpansionBetaShrinkC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_CRF_NLLC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Decode_ChainC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_makeCRFmapsC.c
