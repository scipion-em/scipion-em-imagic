# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *              J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os

from pyworkflow.protocol.params import (IntParam, PointerParam,
                                        EnumParam, FloatParam)
from pyworkflow.constants import PROD
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pwem.emlib.image import ImageHandler
import pyworkflow.utils as pwutils

from ..constants import MODULATION
from .protocol_base import ImagicProtocol


class ImagicProtMSA(ImagicProtocol):
    """
    Performs multivariate statistical analysis of aligned cryo-EM particle
    images in order to identify the principal sources of structural
    variability within a dataset.

    AI Generated:

    Multivariate Statistical Analysis (ImagicProtMSA) — User Manual
        Overview

        The Multivariate Statistical Analysis protocol applies IMAGIC-based
        statistical dimensionality reduction methods to a set of aligned
        particle images. Its primary purpose is to extract the dominant
        patterns of variability present in the dataset and represent them
        through a reduced number of eigenimages and associated factors.

        In cryo-EM workflows, this protocol is commonly used as a preparatory
        stage before particle classification. By transforming high-dimensional
        particle images into a smaller statistical representation, it becomes
        possible to identify meaningful structural differences while reducing
        the influence of noise and redundant information. More info:
        https://imagic4d.readthedocs.io/en/latest/

        Biological Context and Purpose

        Biological macromolecules often display conformational flexibility,
        compositional heterogeneity, preferred orientations, and varying image
        quality. Multivariate statistical analysis helps reveal these sources
        of variability by identifying the dominant trends present across the
        particle population.

        From a biological perspective, the extracted eigenimages may capture
        conformational transitions, orientation-dependent differences, or
        systematic experimental variations. This information provides a more
        interpretable representation of the dataset and improves the ability
        to separate structurally distinct particle populations during later
        classification stages.

        Inputs and General Workflow

        The protocol requires a set of aligned particle images. Proper alignment
        before analysis is essential because the statistical decomposition
        assumes that corresponding structural regions occupy similar spatial
        positions across all particles.

        During execution, the protocol computes a reduced representation of
        the dataset in terms of eigenimages and eigenvalues. Each particle is
        projected into this reduced feature space, allowing the dominant modes
        of structural variation to be analyzed more efficiently than in the
        original image space.

        Choice of Statistical Metric

        The protocol provides several statistical distance metrics that define
        how similarities between particles are evaluated. Different metrics may
        emphasize different aspects of the data and can influence the biological
        interpretation of variability.

        The Euclidian metric corresponds to classical principal component
        analysis and is appropriate for many general-purpose analyses. The
        Chi-square metric is related to correspondence analysis and may provide
        advantages when relative intensity distributions are important. The
        Modulation metric is commonly recommended for cryo-EM applications
        because it often provides improved robustness for noisy experimental
        data.

        In practical workflows, users frequently begin with the Modulation
        metric and later compare results using alternative metrics if the
        dataset contains unusual heterogeneity or challenging signal-to-noise
        conditions.

        Number of Factors and Dimensionality Reduction

        One of the most important biological parameters is the number of
        eigenimages retained during analysis. These factors represent the main
        systematic variations present in the dataset.

        Retaining too few factors may discard meaningful structural variability,
        while retaining too many may introduce noise into subsequent analyses.
        In many cryo-EM datasets, the first factors capture major conformational
        differences or dominant orientation changes, whereas higher-order
        factors progressively contain weaker information.

        Biological users should consider the structural complexity of the
        specimen when selecting the dimensionality. Flexible complexes or
        heterogeneous assemblies often require more factors than highly rigid
        particles.

        Iterative Convergence and Stability

        The eigenimage estimation process is iterative and gradually converges
        toward a stable statistical representation of the dataset. The number
        of iterations controls how extensively the protocol refines this
        solution.

        An overcorrection factor controls the convergence behavior of the
        iterative optimization. Moderate values generally improve convergence
        speed, whereas excessively large values may destabilize the process and
        produce oscillatory behavior. Careful parameter selection therefore
        contributes to stable and biologically meaningful results.

        Masking and Biological Focus

        Masking is a particularly important aspect of multivariate analysis
        because it determines which image regions contribute to the statistical
        decomposition.

        Circular masks are convenient for globular particles and exploratory
        analyses. They provide a simple way to exclude large solvent regions
        that would otherwise dominate the variance calculations. Custom masks
        are more appropriate for irregular or flexible biological assemblies,
        especially when the user wishes to focus the analysis on a stable core
        region while excluding highly mobile domains or surrounding noise.

        From a biological standpoint, a well-designed mask should include the
        structurally informative regions of the particle while minimizing
        irrelevant background variation. Poor masking may reduce the quality
        and interpretability of the resulting eigenimages.

        Outputs and Interpretation

        The protocol produces eigenimages, eigenvalues, and statistical
        coordinate representations describing the principal modes of variation
        within the dataset. These outputs are commonly used for downstream
        classification and heterogeneity analysis.

        Eigenimages should not necessarily be interpreted as direct physical
        structures. Instead, they represent statistical patterns that describe
        correlated variability among particles. Nevertheless, biologically
        meaningful conformational transitions often become visible through the
        dominant factors.

        The reduced statistical coordinates generated by the protocol provide
        the foundation for clustering and classification methods that separate
        particles into structurally coherent groups.

        Practical Recommendations

        In routine cryo-EM workflows, it is often advisable to begin with a
        moderate number of factors and inspect the resulting classifications
        obtained downstream. Increasing the number of factors may reveal subtle
        variability, but excessive dimensionality can amplify noise and reduce
        classification stability.

        Careful masking is frequently one of the most important determinants of
        success. Flexible regions, detergent micelles, or large solvent areas
        may dominate the variance if not properly excluded. For difficult
        datasets, testing several masks and comparing the biological coherence
        of the resulting classifications is often beneficial.

        Users should also monitor convergence behavior when selecting the
        overcorrection factor and iteration count. Stable convergence generally
        produces more reliable statistical representations for subsequent
        biological interpretation.

        Final Perspective

        For most cryo-EM analyses, multivariate statistical analysis is not
        merely a mathematical compression technique but a biologically relevant
        strategy for uncovering structural variability hidden within noisy
        particle datasets. Appropriate selection of statistical metrics,
        dimensionality, convergence settings, and masking approaches strongly
        influences the quality of downstream classification and the reliability
        of structural interpretation.
    """
    _label = 'msa'
    MSA_DIR = 'MSA'
    _devStatus = PROD

    def __init__(self, **kwargs):
        ImagicProtocol.__init__(self, **kwargs)

        self._params = {'eigen_img': os.path.join(self.MSA_DIR, 'eigen_img'),
                        'msa_pixvec_coord': os.path.join(self.MSA_DIR, 'msa_pixvec_coord'),
                        'msa_eigen_pixel': os.path.join(self.MSA_DIR, 'msa_eigen_pixel')
                        }

# --------------------------- DEFINE param functions --------------------------

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputParticles', PointerParam,
                      label="Input particles", important=True,
                      pointerClass='SetOfParticles',
                      help='Select the input particles to perform MSA.')
        form.addParam('distanceType', EnumParam,
                      default=MODULATION,
                      choices=['EUCLIDIAN', 'CHISQUARE', 'MODULATION'],
                      label='MSA distance type',
                      expertLevel=LEVEL_ADVANCED,
                      help='Select general metric for square distance '
                           'between images:\n\n'
                           'a. Euclidian metric (Principal Components '
                           'Analysis: PCA)\n'
                           'b. Chi-square metric (Correspondence Analysis: CA)\n'
                           'c. Modulation metric (Modulation Analysis: MA, '
                           'recommended)')
        form.addParam('numberOfFactors', IntParam, default=25,
                      label='Number of factors (eigenimages)',
                      help='A 64x64 image can be expressed as a vector of '
                           '4096 dimensions. In this step, we will reduce '
                           'this number of dimensions to the number of factors '
                           'specified here. These factors will represent the '
                           'largest systematic variations in the data.\n\n'
                           'The number of eigenimages that should be '
                           'used depends on the '
                           'complexity of the input data.')
        form.addParam('numberOfIterations', IntParam, default=25,
                      label='Number of iterations',
                      help='The calculation of the eigenimages and the '
                           'related eigenvalues are determined in an iterative '
                           'process.\nPlease give the number of iterations '
                           'wanted.\n\nNOTE: The iterations will stop '
                           'automatically if the eigenimage (eigenvector '
                           'eigenvalue) calculations are converging')
        form.addParam('overcorrectionFactor', FloatParam,
                      default=0.8, expertLevel=LEVEL_ADVANCED,
                      label='Overcorrection factor [0 - 0.9]',
                      help='The overcorrection factor is a very important '
                           'parameter in the MSA program. It determines '
                           'the convergence speed of the Eigenvector '
                           'Eigenvalue algorithm. However, if a too large '
                           'overcorrection is chosen, the algorithm may '
                           'start oscillating. Oscillations of the algorithm '
                           'may be observed in the plot of the sum of the '
                           'eigenvalues versus iteration number which is '
                           'part of the output of this program. Divergence '
                           'may thus only be detected a posteriori.\n\n'
                           'The accepted values for OVER_CORRECTION lie '
                           'between 0 and 0.9.')
        form.addParam('maskType', EnumParam,
                      choices=['circular', 'object'], default=0,
                      display=EnumParam.DISPLAY_HLIST,
                      label='Mask type',
                      help='Select which type of mask do you want to apply. '
                           'Only the pixels beneath this mask will be analyzed. '
                           'In the simplest case, a circular mask can be used. '
                           'Alternatively, a custom mask can be used '
                           'which follows the contour of the particle '
                           '(but not too tightly).')
        form.addParam('radius', IntParam, default=-1,
                      label='Mask radius (px)', condition='maskType==0',
                      help='If -1, the entire image (in pixels) will be considered.')
        form.addParam('maskImage', PointerParam,
                      label="Mask image", condition='maskType==1',
                      pointerClass='Mask',
                      help="Select a mask file")

        form.addParallelSection(threads=0, mpi=1)

    # --------------------------- INSERT steps functions ----------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep', needsGPU=False)

        if self.maskType > 0:
            self._insertFunctionStep('convertMaskStep',
                                     self.maskImage.get().getObjId(),
                                     needsGPU=False)
        else:
            self._insertFunctionStep('createMaskStep', needsGPU=False)

        self._insertFunctionStep('msaStep', needsGPU=False)

    # --------------------------- STEPS functions -----------------------------

    def convertInputStep(self):
        # we need to put all images into a single stack
        # TODO: skip writeStack, convert directly via e2proc2d.py
        inputParticles = self.inputParticles.get()
        tmpStack = self._getTmpPath('input_particles.stk')
        inputParticles.writeStack(tmpStack, applyTransform=True)
        ImageHandler().convert(tmpStack, self.getParticlesStack())

    def convertMaskStep(self, maskType):
        """ Convert the input mask to Imagic. """
        if maskType > 0:  # mask from file
            maskFn = self._getTmpPath('mask.img')
            ImageHandler().convert(self.maskImage.get(), maskFn)

    def createMaskStep(self):
        """ Create a circular mask in Imagic format. """
        inputParticles = self.inputParticles.get()
        radius = self.radius.get()

        if self.maskType.get() == 0:
            if radius < 0:  # usually -1
                radiusMask = inputParticles.getDim()[0] / 2
                # use half of input dim
            else:
                radiusMask = radius
            outMask = self._getTmpPath('mask.img')
            ih = ImageHandler()
            ih.createCircularMask(radiusMask, inputParticles.getFirstItem(),
                                  outMask)

    def msaStep(self):
        """ Run MSA on input particles. """

        distances = ['EUCLIDIAN', 'CHISQUARE', 'MODULATION']
        distance_name = distances[self.distanceType.get()]

        self._params.update({'msa_dir': self.MSA_DIR,
                             'msa_distance': distance_name,
                             'num_factors': self.numberOfFactors.get(),
                             'num_iter': self.numberOfIterations.get(),
                             'overcorrectionFactor': self.overcorrectionFactor.get(),
                             'mpi_procs': self.numberOfMpi.get()
                             })

        msaDir = self._getPath(self.MSA_DIR)
        pwutils.cleanPath(msaDir)
        pwutils.makePath(msaDir)

        self.runTemplate('msa/msa-run.b', self._params)

    # --------------------------- INFO functions ------------------------------

    def _validate(self):
        errors = []

        if self.maskImage.get():
            # check pixel size and image size
            pixel_inp = self.inputParticles.get().getSamplingRate()
            pixel_mask = self.maskImage.get().getSamplingRate()
            if pixel_inp != pixel_mask:
                errors.append('Pixel sizes of input images and mask '
                              'should be the same!')

            if self.maskImage.get().getDim()[0] != self.inputParticles.get().getDim()[0]:
                errors.append('Image size of input images and mask is not the same!')

        # check overcorrection factor value
        value = round(self.overcorrectionFactor.get(), 1)
        if not 0.0 <= value <= 0.9:
            errors.append('Overcorrection factor value should be '
                          'in [0 - 0.9] range!')

        # for compatibility with old IMAGIC versions
        if not (self.numberOfFactors.get() <= 64 and
                self.numberOfIterations.get() <= 64):
            errors.append('For compatibility with old IMAGIC versions, '
                          'number of eigenimages and iterations should be <= 64')

        # check radius size vs input particles
        radiusmax = self.inputParticles.get().getDim()[0] / 2
        if not self.radius.get() <= radiusmax:
            errors.append('Radius cannot be bigger than half-size '
                          'of input images!')

        return errors

    def _citations(self):
        return ['Borland1990']

    def _summary(self):
        summary = list()
        summary.append('This protocol generates only eigenimages '
                       '(factors), that will be used later for'
                       ' MSA-based classification.')

        if self.distanceType == 0:
            summary.append('Distance type: *Euclidian*')
        if self.distanceType == 1:
            summary.append('Distance type: *ChiSquare*')
        if self.distanceType == 2:
            summary.append('Distance type: *Modulation*')

        summary.append('Number of factors: *%s*' % self.numberOfFactors)

        if self.maskType == 0:  # circular mask
            if self.radius == -1:
                summary.append('Mask: *Circular, of radius 1/2 image dimension*')
            else:
                summary.append('Mask: *Circular, of radius: %s*' % self.radius)
        else:  # custom mask
            summary.append('Mask: *Custom file*')

        return summary

    def _methods(self):

        msg = "\nInput particles %s were subjected to MSA using " %\
              self.getObjectTag('inputParticles')

        if self.distanceType == 0:
            msg += "Euclidian metric, "
        if self.distanceType == 1:
            msg += "ChiSquare metric, "
        if self.distanceType == 2:
            msg += "Modulation metric, "

        msg += "computing %s factors, and using a " % self.numberOfFactors

        if self.maskType == 0:  # circular mask
            if self.radius == -1:
                msg += "circular mask of radius half the image dimension."
            else:
                msg += "circular mask of radius %s pixels." % self.radius
        else:  # custom mask
            msg += "custom mask %s." % self.getObjectTag('maskImage')

        return [msg]

    # --------------------------- UTILS functions -----------------------------

    def getParticlesStack(self):
        return self._getPath('input_particles.img')

    def _getOutputPath(self, fn):
        """ Return the output file from the run directory and the MSA dir. """
        return self._getPath(self.MSA_DIR, fn)

    def getOutputEigenImages(self):
        return self._getOutputPath('eigen_img.img')

    def getOutputLis(self):
        return self._getOutputPath('msa.lis')

    def getOutputPlt(self):
        return self._getOutputPath('msa.plt')
