
 Updated 1 hour ago
Surface tomography, version 3.0. This code calculates three dimensional point-clouds from tomographic tilt series of micrographs, using methods from differential geometry, as opposed to conventional back-projection. In some common and useful forms of electron microscopy, particularly phase contrast imaging, there is rarely an operative projection law and so the Radon transform is not reliable for reconstructing three dimensional (3D) shape. Sparse features exhibited in rotating micrographs can nonetheless be tracked as a smooth function of specimen orientation or rocking of the incident beam tilt, the continuous variation of which can be harnessed for differential triangulation to recover shape information.

StomoVer3 builds upon two previous versions (published in Computer Physics Communications), which find topological features such as edges, ridges and/or valleys, to track the motion of common sets of points between micrographs in a tilt series, through variants of Canny edge detection. The interior and exterior morphology of scattering surfaces are then reconstructed by fitting polynomials to the angular variation of these features, to recover the intersection points of local ray bundles with these surfaces and also estimate the surface normal. The persistence of such features across several micrographs need only be nascent for robust tracking, permitting shape reconstruction across a broad class of specimen types for which strong scattering may invalidate the requisite projection law of computer automated tomography.

By forgoing reconstruction of the full 3D scattering density, the resulting sparse and unstructured point clouds are orders of magnitude smaller in size than conventional tomograms and are largely undistorted by systematic errors such as the missing wedge problem. The "surface tomography" of StomoVer3 is not tomography per-se, as back-projection is not used and the algorithm is based upon smooth differentiation of localised features, rather than integration of a distributed scattering density. As such, for specimens and imagining modalities amenable to reliable back-projection (such as high angle annular dark field), StomoVer3 will likely produce noisier reconstructions yet is expected to provide a complementary and light-weight quantification of specimen morphology.
