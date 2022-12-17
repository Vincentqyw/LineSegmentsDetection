## [DELETE SUBMODULES]
# git rm --cached SOLD2
# git rm --cached DeepLSD
# git rm --cached TP-LSD
# git rm --cached Unified-Line-Segment-Detection
# git rm --cached LETR
# git rm --cached Deep-Hough-Transform-Line-Priors
# git rm --cached F-Clip
# git rm --cached mlsd
# git rm --cached afm_cvpr2019
# git rm --cached lcnn
# git rm --cached hawp

# rm -rf SOLD2
# rm -rf DeepLSD
# rm -rf TP-LSD
# rm -rf Unified-Line-Segment-Detection
# rm -rf LETR
# rm -rf Deep-Hough-Transform-Line-Priors
# rm -rf F-Clip
# rm -rf mlsd
# rm -rf afm_cvpr2019
# rm -rf lcnn
# rm -rf hawp

# rm .gitmodules

## [ADD SUBMODULES]
# mkdir thirdparty
# git submodule add https://github.com/cvg/DeepLSD.git thirdparty/DeepLSD
# git submodule add https://github.com/navervision/mlsd.git thirdparty/MLSD
# git submodule add https://github.com/Delay-Xili/F-Clip.git thirdparty/F-Clip
# git submodule add https://github.com/cvg/SOLD2.git thirdparty/SOLD2
# git submodule add https://github.com/mlpc-ucsd/LETR.git thirdparty/LETR
# git submodule add https://github.com/cherubicXN/hawp.git thirdparty/HAWP
# git submodule add https://github.com/Siyuada7/TP-LSD.git thirdparty/TP-LSD
# git submodule add https://github.com/lh9171338/Unified-Line-Segment-Detection.git thirdparty/ULSD
# git submodule add https://github.com/yanconglin/Deep-Hough-Transform-Line-Priors.git thirdparty/DHTLP
# git submodule add https://github.com/cherubicXN/afm_cvpr2019.git thirdparty/AFM-LSD
# git submodule add https://github.com/zhou13/lcnn.git thirdparty/LCNN
# git submodule add https://github.com/mtamburrano/LBD_Descriptor.git thirdparty/LBD