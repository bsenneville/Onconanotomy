# Onconanotomy
Deciphering Tumour Tissue Organization by 3D Electron Microscopy

This code is designed to identify bioarchitectural parameters that shape the internal and spatial organization of tumours.

This code contains 3 tutorials:

- "test_cells_alignement.m": In this tutorial, each cell is characterized using its main axis. To this end, a Principal Component Analysis (PCA) is applied on the voxel coordinates within the segmented mask of a structure of interest. The main axes of cells is used to calculate the best 2D “alignment plane” (in the least-squares sense: the sum of squared differences is minimized between the observed main axes of cells, and the fitted value provides by a 2D plane equation). Angles between the main axis of each cell and the alignment plane are then calculated.

- "test_cells_polarisation.m": In this tutorial, each cell is characterized using its main axis. To this end, a Principal Component Analysis (PCA) is applied on the voxel coordinates within the segmented mask of a structure of interest. 3D virtual rays are emitted by each cell along its main axis in both directions (the ray radius is a user-defined input parameter for the algorithm). The accumulated beam intensity is then calculated on a voxel-by-voxel basis.

- 'test_volume_versus_vessel.m': In this tutorial, the size of different cell components is compared to the distance to a blood capillary.

This code has been written by Baudouin Denis de Senneville.

This code was developed under the commercial software Matlab ©1994-2021 The MathWorks, Inc.

Please cite the following paper if you are using this code:

Baudouin Denis de Senneville, Fatma Zohra Khoubai, Marc Bevilacqua, Alexandre Labedade, Kathleen Flosseau, Christophe Chardot, Sophie Branchereau, Jean Ripoche, Stefano Cairo, Etienne Gontier, Christophe F. Grosset,  Deciphering Tumour Tissue Organization by 3D Electron Microscopy and machine learning, bioRxiv 2021.06.15.446847; doi: https://doi.org/10.1101/2021.06.15.446847. [Download](https://www.biorxiv.org/content/10.1101/2021.06.15.446847v1.abstract)
