线段检测算法汇总，请各位自行调试。

## [LSD](http://www.ipol.im/pub/art/2012/gjmr-lsd/)

- 论文标题："LSD: a Line Segment Detector"
- 项目主页：http://www.ipol.im/pub/art/2012/gjmr-lsd/
- 论文地址：http://www.ipol.im/pub/art/2012/gjmr-lsd/article.pdf
- 代码地址：http://www.ipol.im/pub/art/2012/gjmr-lsd/lsd_1.5.zip

该方法是目前性价比（速度精度）最好的算法，现已经集成到`opencv`中[`LSDDetector`](https://docs.opencv.org/master/d1/dbd/classcv_1_1line__descriptor_1_1LSDDetector.html)。LSD能够在线性时间内检测到亚像素精度的线段。无需调整参数，适用于各种场景。因为每张图有误检，LSD能够控制误检率。PS：论文此处不介绍了，可以参考[这里](https://blog.csdn.net/chishuideyu/article/details/78081643?locationNum=9&fps=1)。

## TP-LSD

- 论文标题："TP-LSD: Tri-Points Based Line Segment Detector", ECCV 2020.
- 论文地址：https://arxiv.org/abs/2009.05505
- 代码地址：https://github.com/Siyuada7/TP-LSD


## ULSD-ISPRS

- 论文标题："ULSD: Unified Line Segment Detection across Pinhole, Fisheye, and Spherical Cameras",  arXiv 2021.
- 论文地址：https://arxiv.org/abs/2011.03174
- 代码地址：https://github.com/lh9171338/Unified-Line-Segment-Detection

## LETR

- 论文标题："Line Segment Detection Using Transformers without Edges",  CVPR 2021.
- 论文地址：[https://arxiv.org/abs/2101.01909](https://arxiv.org/abs/2101.01909)
- 代码地址：https://github.com/mlpc-ucsd/LETR


## Deep-Hough-Transform-Line-Priors

- 论文标题："Deep Hough-Transform Line Priors",  ECCV 2020.
- 论文地址：[https://arxiv.org/pdf/2007.09493.pdf](https://arxiv.org/pdf/2007.09493.pdf)
- 代码地址：https://github.com/yanconglin/Deep-Hough-Transform-Line-Priors


## F-Clip

- 论文标题："Fully Convolutional Line Parsing",  arXiv 2021.
- 论文地址：[https://arxiv.org/pdf/2104.11207v2.pdf](https://arxiv.org/pdf/2104.11207v2.pdf)
- 代码地址：https://github.com/Delay-Xili/F-Clip


## M-LSD

- 论文标题："M-LSD: Towards Light-weight and Real-time Line Segment Detection",  arXiv 2021.
- 论文地址：[https://arxiv.org/pdf/2106.00186.pdf](https://arxiv.org/pdf/2106.00186.pdf)
- 代码地址：https://github.com/navervision/mlsd


## AFM-LSD

- 论文标题："Learning Attraction Field Representation for Robust Line Segment Detection",  CVPR 2019.
- 论文地址：[https://arxiv.org/pdf/1812.02122.pdf](https://arxiv.org/pdf/1812.02122.pdf)
- 代码地址：https://github.com/cherubicXN/afm_cvpr2019


## L-CNN

- 论文标题："End-to-End Wireframe Parsing",  ICCV 2019.
- 论文地址：[https://arxiv.org/pdf/1905.03246.pdf](https://arxiv.org/pdf/1905.03246.pdf)
- 代码地址：https://github.com/zhou13/lcnn


## LSWMS

- 论文标题："Line segment detection using weighted mean shift procedures on a 2D slice sampling strategy"
- 论文地址：[https://www.researchgate.net/LSWMS.pdf](https://www.researchgate.net/profile/Marcos_Nieto3/publication/220654859_Line_segment_detection_using_weighted_mean_shift_procedures_on_a_2D_slice_sampling_strategy/links/56a5d56a08aef91c8c16b1ac.pdf?inViewer=0&origin=publication_detail&pdfJsDownload=0)
- 代码地址：https://sourceforge.net/projects/lswms/

## EDline（ED: Edge Drawing）
- 论文标题："Edge Drawing: A Combined Real-Time Edge and Segment Detector"
- 论文地址：https://www.sciencedirect.com/science/article/abs/pii/S1047320312000831
- 代码地址：https://github.com/mtamburrano/LBD_Descriptor

## [CannyLines](http://cvrs.whu.edu.cn/projects/cannyLines/)

- 论文标题："CannyLines: A Parameter-Free Line Segment Detector"
- 项目主页：http://cvrs.whu.edu.cn/projects/cannyLines/
- 论文地址：http://cvrs.whu.edu.cn/projects/cannyLines/papers/CannyLines-ICIP2015.pdf
- 代码地址：http://cvrs.whu.edu.cn/projects/cannyLines/codes/CannyLines-v3.rar

本文提出了一种鲁棒的线段检测算法来有效地检测来自输入图像的线段。首先，文章提出了一种无参数的Canny算子，称为Canny（PANYPF），通过自适应地设置传统Canny算子的低阈值和高阈值来鲁棒地从输入图像中提取边缘。然后，提出了有效的像素连接和分割技术，直接从边缘图中收集共线点群，用于基于最小二乘拟合方法拟合初始线段。第三，通过有效的扩展和合并，产生更长、更完整的线段。最后，利用亥姆霍兹原理（Helmholtz Principle）对所有的检测线段检测，主要考虑梯度方向和幅度信息。该算法能够在人工场景中能够获得比LSD以及EDline精度更高以及平均长度更高的线段。

## MCMLSD

- 论文标题："MCMLSD: A Dynamic Programming Approach to Line Segment Detection"
- 论文地址：http://openaccess.thecvf.com/content_cvpr_2017/papers/Almazan_MCMLSD_A_Dynamic_CVPR_2017_paper.pdf
- 代码地址：http://www.elderlab.yorku.ca/?smd_process_download=1&download_id=8423
