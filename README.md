# (H)GIR-MRF
Structured graph based image regression for unsupervised multimodal change detection

## Introduction
MATLAB Code: (H)GIR-MRF - 2021
This is a test program for the graph based image regression and MRF segmentation method (GIR-MRF) for multimodal change detection problem.

GIR-MRF is an unsupervised image regression method based on the inherent structure consistency between heterogeneous images, which learns a structured graph and computes the regression image by graph projection. Firstly, the proposed method uses the self-expression property to preserve the global structure of image and uses the adaptive neighbor approach to capture the local structure of image in the graph learning process. Then, with the learned graph, two types of structure constraints are introduced into the regression model: one corresponds to the global selfexpression constraint and the other corresponds to the local similarity constraint, which can be further implemented by using graph or hypergraph Laplacian based regularization. Finally, a Markov segmentation model is designed to calculate the binary change map, which combines the change information and spatial information to improve the detection accuracy.

Please refer to the paper for details. You are more than welcome to use the code! 

===================================================

## Available datasets

#2-Texas is download from Professor Michele Volpi's webpage at https://sites.google.com/site/michelevolpiresearch/home.

#6-California is download from Dr. Luigi Tommaso Luppino's webpage (https://sites.google.com/view/luppino/data) and it was downsampled to 875*500 as shown in our paper.

===================================================

## Citation

If you use this code for your research, please cite our paper. Thank you!

@article{SUN202216,  
title = {Structured graph based image regression for unsupervised multimodal change detection},  
journal = {ISPRS Journal of Photogrammetry and Remote Sensing},  
volume = {185},  
pages = {16-31},  
year = {2022},  
issn = {0924-2716},  
doi = {https://doi.org/10.1016/j.isprsjprs.2022.01.004},  
url = {https://www.sciencedirect.com/science/article/pii/S0924271622000089},  
author = {Yuli Sun and Lin Lei and Xiang Tan and Dongdong Guan and Junzheng Wu and Gangyao Kuang}} 

## Future work

Our future work is to improve the computation effciency and design an effective fusion strategy to fuse the forward and backward transformations, thus improving the CD performance.

## Q & A

If you have any queries, please do not hesitate to contact me (sunyuli@mail.ustc.edu.cn ).
