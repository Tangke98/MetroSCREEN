.. highlight:: shell

.. role:: bash(code)
   :language: bash

Mouse olfactory bulb (Stereo-seq)
---------------------------------

Here we use a mouse olfactory bulb dataset profiled by Stereo-seq (`Chen et al., Cell, 2022 <https://www.sciencedirect.com/science/article/pii/S0092867422003993>`_) to demonstrate the usage of Cellist. The original study provides spatial expression and a corresponding ssDNA staining image. Users can download the processed data from `here <https://github.com/wanglabtongji/Cellist/tree/main/test/Stereoseq_Mouse_OB>`_.

Step 1 Registration between staining and spatial expression profile
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

While the images are initially aligned with RNA coordinates, minor misalignments may still occur. The :bash:`align` function could be used to refine the alignment, which integrates the :bash:`refine_alignment` function with the rigid mode in `Spateo <https://spateo-release.readthedocs.io/en/latest/technicals/cell_segmentation.html#alignment-of-stain-and-rna-coordinates>`_. 
::

   cellist align --gem Data/DP8400013846TR_F5.bin1.olfactorybulb_cropped.gem \
   --tif Data/ssDNA_cropped_3625_9545_950_5630.tiff \
   --nworkers 8 \
   --outdir Result/Alignment \
   --outprefix DP8400013846TR_F5

.. image:: ../_static/img/DP8400013846TR_F5_alignment_compare.png
   :width: 100%
   :align: center

Step 2 Watershed segmentation of nucleus
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

The initial nucleus segmentation is required for refined cell segmentation by Cellist. In Cellist, we utilize the watershed algorithm to segment nuclei in the ssDNA image, which is implemented by the function of :bash:`watershed`. 

::

   cellist watershed --gem Data/DP8400013846TR_F5.bin1.olfactorybulb_cropped.gem \
   --tif Result/Alignment/DP8400013846TR_F5_regist_transposed_aligned_by_Spateo.tiff \
   --min-distance 6 \
   --outdir Result/Watershed \
   --outprefix DP8400013846TR_F5

.. image:: ../_static/img/DP8400013846TR_F5_cell_boundary.png
   :width: 100%
   :align: center

Step 3 Cell segmentation by Cellist
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

With nucleus segmentation completed, the next step is to expand the nucleus labels to include the cytoplasm, namely, cell segmentation. In cellist, we take both expression similarity and spatial proximity into consideration when assigning non-nucleus spots to labelled nuclei. 

::

   cellist seg --platform barcoding \
   --resolution 0.5 \
   --gem Data/DP8400013846TR_F5.bin1.olfactorybulb_cropped.gem \
   --spot-count-h5 Result/Watershed/DP8400013846TR_F5_bin1.h5 \
   --nuclei-prop Result/Watershed/DP8400013846TR_F5_watershed_nucleus_property.txt \
   --nuclei-count-h5 Result/Watershed/DP8400013846TR_F5_waterhsed_segmentation_cell_count.h5 \
   --watershed-seg Result/Watershed/DP8400013846TR_F5_watershed_nucleus_coord.txt \
   --nworkers 16 \
   --cell-radius 15 \
   --spot-imputation-distance 2.5 \
   --prob-cutoff 0.6 \
   --outdir Result/Cellist \
   --outprefix DP8400013846TR_F5

The results of :bash:`seg` will be stored in the :bash:`Result/Cellist` floder, and the detailed descritions are shown as below.

+-----------------------------------------------+-------------------------------------------------------------------------------+
| File                                          | Description                                                                   |
+===============================================+===============================================================================+
| Data_HVG/                                     | The directory stores small patches cropped from the slide.                    |
+-----------------------------------------------+-------------------------------------------------------------------------------+
| {outprefix}_segmentation.txt                  | The spot-level cell segmentation result where each row represents a spot.     |
+-----------------------------------------------+-------------------------------------------------------------------------------+
| {outprefix}_segmentation_cell_count.h5        | The aggrefated cell-level expression matrix, stored in the format of h5,      |
|                                               | where each row represents a gene and each column represents a cell.           |
+-----------------------------------------------+-------------------------------------------------------------------------------+
| {outprefix}_segmentation_cell_coord.txt       | The spatial coordinates of the segmented cells, which correspond to the cells |
|                                               | in the above expression file.                                                 |
+-----------------------------------------------+-------------------------------------------------------------------------------+
| {outprefix}_segmentation_plot.pdf             | Visualization of the cell segmentation results.                               |
+-----------------------------------------------+-------------------------------------------------------------------------------+
| {outprefix}_cellist_corr_nucl_cyto_df.txt     | The correlation of expression between nucleus and cytoplasm within each cell. |
+-----------------------------------------------+-------------------------------------------------------------------------------+
| parameters.json                               | Parameters to run :bash:`cellist` and statistics of the segmentation results. |
+-----------------------------------------------+-------------------------------------------------------------------------------+

Step 4 Spatially-aware expression imputation at the cell level (Optional)
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

In certain cases, the gene coverage within each cell may still be insufficient for in-depth single-cell analyses. To mitigate this issue, Cellist offers an optional imputation function that recovers missing gene expression at the cell level, leveraging information from neighboring cells close in both physical space and low dimensional embedding learned from expression data.

::

   cellist impute --expr Result/Cellist/alpha_0.8_sigma_1.0_beta_10_gene_HVG_dist_15_iter_False_prob_0.6_neigh_2.5/DP8400013846TR_F5_segmentation_cell_count.h5 \
   --spatial Result/Cellist/alpha_0.8_sigma_1.0_beta_10_gene_HVG_dist_15_iter_False_prob_0.6_neigh_2.5/DP8400013846TR_F5_segmentation_cell_coord.txt \
   --nworkers 8 \
   --outdir Result/Imputation_louvain \
   --outprefix DP8400013846TR_F5


