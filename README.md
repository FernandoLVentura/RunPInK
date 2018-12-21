# RunPInK
Scripts for running PInK that:

1. Take a FITS catalogue and the directory of associated images.
2. Find suitable sources based on peak signal strength
3. Find the associated FITS images if they exist. If they do not an alternative is searched for and the source omitted if not found.
4. Create cutouts of some given size centred on the source and then compile them into a binary for PINK and new catalogue for us.
5. Use PINK to create a SOM based on the images.
6. Uses SDSS VAC data to create heatmaps of the extracted data.
7. There is also a script to get the indices of the worst fitting sources. These indices may be used to view them from the collection of cutouts created in step 3.

Made for the FIRST (Faint Images of the Radio Sky at Twenty cm) catalogue.

Instructions:
-------------
Ensure that all required variables are set correctly and point to a FIRST catalogue and a directory of images.
Run Philatelist.py (the stamp a collector). This should place individual cutouts in FITS format in a directory called stamps and in binary format, along with a catalogue of indices in CSV format, in the output directory.

In the Output directory, run PINK on the stamps file using the command:

Pink --cuda-off True --seed 42 --dist-func gaussian 0.8 0.05 --inter-store overwrite --neuron-dimension 64 --numrot 360 --num-iter 10 --progress 0.01 --som-width 11 --som-height 11 --layout hexagonal --train stamps.bin result.bin

This creates the hexagonal SOM, 11x11, stored in binary format in the file "result.bin." It uses 360 rotations and 10 iterations. Increase iterations for a better SOM, but increased processing time.

The images may now be mapped onto the SOM by PINK using the command:

Pink --seed 42 --inter-store overwrite --dist-func gaussian 0.8 0.05 --neuron-dimension 64 --numrot 360 --num-iter 10 --progress 0.01 --som-width 11 --som-height 11 --layout hexagonal --map stamps.bin stampmap150k.bin result.bin 

This maps the images in stamps.bin onto the SOM in result.bin, storing the mapping in stampmap150k.bin.
The mapped images may now be processed.

Ensure that all variables are set correctly and run mapDataToSOM.py. This will generate the heatmaps and store the means and standard deviations in CSV format.

The SOM and heatmaps may now be visualised using the scripts in visualiseResults. These script to visualise the SOM came with PINK and has not been modified. The script to visualise the heatmaps consists of the same script to view the SOM, marginally modified to handle RGBa data instead of monochrome images. displayMaps.py contains functions that allow for easy use of these other scripts.

