# spiral_aliasing_reduction

This is the repository for code used in the paper "*Ye Tian, Yongwan Lim, Ziwei Zhao, Dani Byrd, Shrikanth Narayanan, Krishna S. Nayak. Aliasing Artifact Reduction in Spiral Real-Time MRI. Magn Reson Med. https://doi.org/10.1002/mrm.28746.*"

- `demo_recon.m` reproduces reconstructions without correction, with LF correction, and with ES correction. One example dataset is provided under /rawdata to reproduce **Figure 3** in the paper.

- `demo_simu.m` reproduces simulation shown in **Figure 2** in the paper.

## Requirements
- **Orchestra Software Development Kit** is required. This can be found at https://collaborate.mr.gehealthcare.com
- **Variable-Density Spiral Design Functions** is required. This can be found at: https://mrsrl.stanford.edu/~brian/vdspiral/

## Citation
Please cite this paper if you use the code in this repository:

Tian, Y, Lim, Y, Zhao, Z, Byrd, D, Narayanan, S, Nayak, KS. Aliasing artifact reduction in spiral real‐time MRI. Magn Reson Med. https://doi.org/10.1002/mrm.28746

## Data
One example dataset is provided in this repository. All of the datasets used in this work are part of the USC 75-Speaker Speech MRI Database: https://sail.usc.edu/span/75speakers/. More detailed description can be found in this paper: https://arxiv.org/abs/2102.07896.
