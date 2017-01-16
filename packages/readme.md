# Packages applied in 3DCS project
The packages used in the three-dimensioanal compressive sensing (3DCS) project.
## TVAL3
Total Variation minimization by Augmented Lagrangian and ALternating direction ALgorithms ([TVAL3](http://www.caam.rice.edu/~optimization/L1/TVAL3/ "TVAL3")) is adapted from the open-access TV-minimization solver by C. Li, W. Yin and Y. Zhang, Department of Computational and Applied Mathematics, Rice University, Houston, Texas 77005.

The version applied here is the beta2.4 version released on Nov. 14, 2010.

Two major references are:
[1]   Y. Wang, J. Yang, W. Yin, and Y. Zhang, “A New Alternating Minimization Algorithm for Total Variation Image Reconstruction,” *SIAM Journal on Imaging Sciences*, vol. 1, no. 3, pp. 248–272, 2008.

[2]   C. Li, W. Yin, H. Jiang, and Y. Zhang, “An efficient augmented Lagrangian method with applications to total variation minimization,”
*Computational Optimization and Applications*, vol. 56, no. 3, pp. 507–
530, 2013.

## wavelet
The wavelet directory contains the MATLAB source code of constructing general wavelet bases, such as Haar, Daubechies, Symmlet. This directory contains three MATLAB scripts, the first one `MakeONFilter.m` Copy Right by [Waevelab](http://statweb.stanford.edu/~wavelab/ "Wavelab"), the other two `formhg.m` Copy Right by [Xuejun Liao](http://people.ee.duke.edu/~xjliao/ "Xuejun Liao, Duke University") and the rests Copy Right by [Xin Yuan](https://sites.google.com/site/eiexyuan/ "Xin Yuan, Bell lab").
Four major references are:
[1]   S. Mallat, A wavelet tour of signal processing. Academic press, 1999.
[2]   X. Liao, H. Li, and L. Carin, “Generalized alternating projection for
weighted-ℓ 2,1 minimization with applications to model-based compressive sensing,” SIAM Journal on Imaging Sciences, vol. 7, no. 2, pp.797–823, 2014.
[3]   X. Yuan, P. Llull, X. Liao, J. Yang, D. J. Brady, G. Sapiro, and L. Carin, “Low-cost compressive sensing for color video and depth,” 2014 IEEE
Conference on Computer Vision and Pattern Recognition, pp. 3318–3325, 2014.
[4]    X. Yuan, “Generalized alternating projection based total variation minimization for compressive sensing,” in 2016 IEEE International Conference on Image Processing (ICIP), Conference Proceedings, pp. 2539–2543.

## YUV2Image
The MATLAB source codes of [converting YUV CIF 4:2:0 video file to image files](http://www.mathworks.com/matlabcentral/fileexchange/6318-convert-yuv-cif-4-2-0-video-file-to-image-files "File Exchange"). This package is downloaded from MathWorks(R) File Exchange, Copyright(C) by Da Yu.