<p align="center">
  <h1 align="center"><br><ins>Line Segment Detection</ins><br>A collection of line segment detection algorithms</h1>
 
</p>

## 

This repository hosts the papers with code for line segment detection, enjoy yourself. 
(线段检测算法汇总，请各位自行调试)

<!-- ![](assets/images/line-segment-detection.jpg) -->

![SOLD2](assets/images/demo_moving_camera.gif)

<!-- ![demo_deeplsd](assets/images/demo_deeplsd.gif) -->

## Usage

``` bash
git clone --recursive https://github.com/Vincentqyw/LineSegmentsDetection.git
```

## DeepLSD

- Title: "DeepLSD: Line Segment Detection and Refinement with Deep Image Gradients", Arxiv 2022.
- Paper: https://arxiv.org/abs/2212.07766
- Code: https://github.com/cvg/DeepLSD

## M-LSD

- Title: "M-LSD: Towards Light-weight and Real-time Line Segment Detection", AAAI 2022.
- Paper: [https://arxiv.org/abs/2106.00186](https://arxiv.org/abs/2106.00186)
- Code: https://github.com/navervision/mlsd

## F-Clip

- Title: "Fully Convolutional Line Parsing", Neurocomputing 2022.
- Paper: [https://arxiv.org/abs/2104.11207](https://arxiv.org/abs/2104.11207)
- Code: https://github.com/Delay-Xili/F-Clip

## SOLD2

- Title: "SOLD2: Self-supervised Occlusion-aware Line Description and Detection", CVPR 2021.
- Paper: https://arxiv.org/abs/2104.03362
- Code: https://github.com/cvg/SOLD2, [[Kornia Tutorial]](https://kornia-tutorials.readthedocs.io/en/latest/line_detection_and_matching_sold2.html)

## LETR

- Title: "Line Segment Detection Using Transformers without Edges", CVPR 2021.
- Paper: [https://arxiv.org/abs/2101.01909](https://arxiv.org/abs/2101.01909)
- Code: https://github.com/mlpc-ucsd/LETR

## HAWP

- Title: "Holistically-Attracted Wireframe Parsing: From Supervised Learning to Self-Supervised Learning", CVPR 2020.
- Paper: https://arxiv.org/abs/2210.12971
- Code: https://github.com/cherubicXN/hawp

## TP-LSD

- Title: "TP-LSD: Tri-Points Based Line Segment Detector", ECCV 2020.
- Paper: https://arxiv.org/abs/2009.05505
- Code: https://github.com/Siyuada7/TP-LSD

## ULSD-ISPRS

- Title: "ULSD: Unified Line Segment Detection across Pinhole, Fisheye, and Spherical Cameras", ISPRS 2020.
- Paper: https://arxiv.org/abs/2011.03174
- Code: https://github.com/lh9171338/Unified-Line-Segment-Detection

## Deep-Hough-Transform-Line-Priors

- Title: "Deep Hough-Transform Line Priors", ECCV 2020.
- Paper: [https://arxiv.org/abs/2007.09493](https://arxiv.org/abs/2007.09493)
- Code: https://github.com/yanconglin/Deep-Hough-Transform-Line-Priors

## AFM-LSD

- Title: "Learning Attraction Field Representation for Robust Line Segment Detection", CVPR 2019.
- Paper: [https://arxiv.org/abs/1812.02122](https://arxiv.org/abs/1812.02122)
- Code: https://github.com/cherubicXN/afm_cvpr2019

## L-CNN

- Title: "End-to-End Wireframe Parsing", ICCV 2019.
- Paper: [https://arxiv.org/abs/1905.03246](https://arxiv.org/abs/1905.03246)
- Code: https://github.com/zhou13/lcnn

## MCMLSD

- Title: "MCMLSD: A Dynamic Programming Approach to Line Segment Detection", CVPR 2017.
- Paper: http://openaccess.thecvf.com/content_cvpr_2017/papers/Almazan_MCMLSD_A_Dynamic_CVPR_2017_paper.pdf
- Code: http://www.elderlab.yorku.ca/?smd_process_download=1&download_id=8423

## CannyLines

- Title: "CannyLines: A Parameter-Free Line Segment Detector", ICIP 2015.
- Project Page: https://cvrs.whu.edu.cn/cannylines
- Paper: http://cvrs.whu.edu.cn/projects/cannyLines/papers/CannyLines-ICIP2015.pdf
- Code: http://cvrs.whu.edu.cn/projects/cannyLines/codes/CannyLines-v3.rar

> 本文提出了一种鲁棒的线段检测算法来有效地检测来自输入图像的线段。首先，文章提出了一种无参数的Canny算子，称为Canny（PANYPF），通过自适应地设置传统Canny算子的低阈值和高阈值来鲁棒地从输入图像中提取边缘。然后，提出了有效的像素连接和分割技术，直接从边缘图中收集共线点群，用于基于最小二乘拟合方法拟合初始线段。第三，通过有效的扩展和合并，产生更长、更完整的线段。最后，利用亥姆霍兹原理（Helmholtz Principle）对所有的检测线段检测，主要考虑梯度方向和幅度信息。该算法能够在人工场景中能够获得比LSD以及EDline精度更高以及平均长度更高的线段。

## EDline（ED: Edge Drawing）

- Title: "Edge Drawing: A Combined Real-Time Edge and Segment Detector", JVCIR 2012.
- Paper: https://www.sciencedirect.com/science/article/abs/pii/S1047320312000831
- Code: https://github.com/mtamburrano/LBD_Descriptor

## [LSD](http://www.ipol.im/pub/art/2012/gjmr-lsd/)

- Title: "LSD: a Line Segment Detector", Image Processing On Line 2012
- Project Page: http://www.ipol.im/pub/art/2012/gjmr-lsd/
- Paper: http://www.ipol.im/pub/art/2012/gjmr-lsd/article.pdf
- Code: http://www.ipol.im/pub/art/2012/gjmr-lsd/lsd_1.5.zip

> 该方法是目前性价比（速度精度）最好的算法，现已经集成到`opencv`中[`LSDDetector`](https://docs.opencv.org/master/d1/dbd/classcv_1_1line__descriptor_1_1LSDDetector.html)。LSD能够在线性时间内检测到亚像素精度的线段。无需调整参数，适用于各种场景。因为每张图有误检，LSD能够控制误检率。PS: 论文此处不介绍了，可以参考[这里](https://blog.csdn.net/chishuideyu/article/details/78081643?locationNum=9&fps=1)。

## LSWMS

- Title: "Line segment detection using weighted mean shift procedures on a 2D slice sampling strategy", Pattern Analysis and Applications 2011.
- Paper: [https://www.researchgate.net/LSWMS.pdf](https://www.researchgate.net/profile/Marcos_Nieto3/publication/220654859_Line_segment_detection_using_weighted_mean_shift_procedures_on_a_2D_slice_sampling_strategy/links/56a5d56a08aef91c8c16b1ac.pdf?inViewer=0&origin=publication_detail&pdfJsDownload=0)
- Code: https://sourceforge.net/projects/lswms/