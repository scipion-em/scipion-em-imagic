=============
IMAGIC plugin
=============

This plugin includes two protocols to provide wrappers around Multivariate Statistical Analysis (MSA) module of `IMAGIC <https://www.imagescience.de/imagic.html>`_ software suite. IMAGIC is a licensed software, not distributed with Scipion and has to be installed by user.

Installation
------------

You will need to use `3.0 <https://github.com/I2PC/scipion/releases/tag/V3.0.0>`_ version of Scipion to be able to run these protocols. To install the plugin, you have two options:

a) Stable version

.. code-block::

    scipion installp -p scipion-em-imagic

b) Developer's version

    * download repository

    .. code-block::

        git clone https://github.com/scipion-em/scipion-em-imagic.git

    * install

    .. code-block::

        scipion installp -p path_to_scipion-em-imagic --devel

Also, you need a working IMAGIC installation. Default installation path assumed is ``software/em/imagic-180311``, if you want to change it, set *IMAGIC_HOME* in ``scipion.conf`` file to the folder where the IMAGIC is installed (it is the same as *IMAGIC_ROOT* variable in your shell environment). If you want to use MPI-based parallel job execution, make sure you have `openmpi` directory inside IMAGIC installation folder.
To check the installation, simply run the following Scipion test:

``scipion tests tests.em.workflows.test_workflow_imagicMSA.TestImagicWorkflow``

Supported versions
------------------

Since almost with every version of IMAGIC software user interaction with IMAGIC programs is changed, we came up with a way to provide multiple version support. In `imagic/scripts` directory you have a folder for each corresponding version, that contains batch scripts similar to those used by IMAGIC. This way one can create a similar script specific to a certain version. At this moment, versions 110308 (Mar 2011), 160418 (Apr 2016) and 180311 (Mar 2018) are supported. If you experience any problems or need help with adapting scripts for your IMAGIC version, do not hesitate to `create an issue on Github <https://github.com/scipion-em/scipion-em-imagic/issues/new>`_. Besides editing scripts directory you would need to add version number to `_supportedVersions` list in file ``imagic/__init__.py`` and add IMAGIC_HOME variable to ``scipion.conf`` if necessary.

Protocols
---------

    * imagic - msa

        Multivariate Statistical Analysis (MSA) is a powerful technique that allows to identify largest variations in a big data set. It was originally introduced to discriminate between various classes of molecular projections prior to averaging. In the MSA approach, aligned molecular images are submitted to correspondence analysis (CA), that determines the main (orthogonal) directions of inter-image variance and calculates the image coordinates in a system spanned by these newly determined axes. Since this new coordinate system is adapted to the general behavior of the image data, a large reduction in the total amount of data can be obtained: for example, instead of 64x64=4096 density values (pixels) per image, each image is now characterized by the first eight factorial-axis coordinates at the most! With this large data reduction, the classification of the images becomes much simpler.

        To launch MSA protocol, you have to provide an aligned (at least, centered) `SetOfParticles`, number of factors (eigenvectors), maximum number of iterations for algorithm to converge and a mask if you want to analyze variance withing specific area of your particles (fig. 1). Usually 20-25 factors and similar number of iterations are enough even for large data sets.

        .. figure:: https://user-images.githubusercontent.com/6952870/50742308-79dc5800-1209-11e9-843e-f3a6afbacc26.png
           :align: left
           :alt: GUI input form of the imagic - msa protocol

        If you want to play with advanced parameters, select *Advanced* expert level and look at the Help message for any particular option.

        .. figure:: https://user-images.githubusercontent.com/6952870/50742309-7c3eb200-1209-11e9-8121-c358e6893a71.png
           :align: left
           :alt: Advanced protocol parameters

        This protocol does not generate any results except eigenimages. Eigenimages represent eigenvectors in the image space and account for major density variations in the data set (fig. 3). The very first eigenimage is a total sum of all particles. The following eigenimages show data set variance in a decreasing order. Last eigenimages are usually very noisy and can be discarded from further analysis.

        .. figure:: https://user-images.githubusercontent.com/6952870/50742310-7d6fdf00-1209-11e9-8b19-2d888bdcce48.png
           :align: left
           :alt: Displaying results of MSA protocol

    * imagic - msa-classify

        After MSA analysis you can use a subset of eigenimages for clustering original images (that will be reconstructed from a linear combination of selected eigenvectors) into groups. IMAGIC MSA module implements hierarchical ascendant classification (HAC) that tries to merge images into clusters by minimizing intra-class variance and maximizing inter-class variance between different clusters.

        The msa-classify protocol requires the `SetOfParticles` from the previous run of msa, a number of factors to use for analysis and a number of classes. At this moment only *first N eigenimages* can be chosen for MSA-based classification. In the future versions of the protocol it will be possible to select eigenimages independently and also assign weighting coefficients for more advanced image analysis.

        .. figure:: https://user-images.githubusercontent.com/6952870/50742311-7ea10c00-1209-11e9-86ad-80a8aac6bc1a.png
           :align: left
           :alt: GUI input form of the imagic - msa classify protocol

        As always, if you want to play with advanced parameters, select *Advanced* expert level and look at the Help message for any particular option.

        .. figure:: https://user-images.githubusercontent.com/6952870/50742315-819bfc80-1209-11e9-83fb-21336230eeee.png
           :align: left
           :alt: Advanced protocol parameters

        In the end you will obtain 2D classes that will most likely display what kind of heterogeneity you have in your data set.

        .. figure:: https://user-images.githubusercontent.com/6952870/50742313-806acf80-1209-11e9-8694-bbab48b8f296.png
           :align: left
           :alt: Output 2D classes

References
----------

1. M van Heel and W Keegstra (1981). IMAGIC: A fast, flexible and friendly image analysis software system. Ultramicroscopy 7: 113-130.
2. M van Heel, G Harauz, EV Orlova, R Schmidt and M Schatz (1996). A new generation of the IMAGIC image processing system. J. Struct. Biol. 116: 17-24.
3. M van Heel, R Portugal, A Rohou, C Linnemayr, C Bebeacua, R Schmidt, T Grant and M Schatz (2012). Four-Dimensional Cryo Electron Microscopy at Quasi Atomic Resolution: "IMAGIC 4D”. International Tables for Crystallography, vol. F, ch. 19.9: 624-628.
4. M van Heel (1984). Multivariate statistical classification of noisy images (randomly oriented biological macromolecules). Ultramicroscopy 13(1-2): 165-183.
5. Lisa Borland and Marin van Heel (1990). Classification of image data in conjugate representation spaces. Journal of the Optical Society of America A 7(4): 601-610.
