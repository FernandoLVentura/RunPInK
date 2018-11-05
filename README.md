# RunPInK
Scripts for running PInK that:

-Take a FITS catalogue and the directory of associated images.
-Find suitable sources based on peak signal strength
-Find the associated FITS images if they exist. If they do not an alternative is searched for and the source omitted if not found.
-Create cutouts of some given size centred on the source and then compile them into a binary for PInK and new catalogue for us.

Made for the FIRST (Faint Images of the Radio Sky at Twenty cm) catalogue.
Outputs a single .bin for PInK and a CSV catalogue of the included sources.
