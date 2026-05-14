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

from enum import Enum

from pwem.protocols import ProtClassify2D
from pwem.objects import SetOfClasses2D
from pyworkflow.object import Float
from pyworkflow.protocol.params import PointerParam, IntParam, BooleanParam
from pyworkflow.protocol.constants import LEVEL_ADVANCED
import pyworkflow.utils as pwutils
from pyworkflow.constants import PROD

from ..utils import ImagicPltFile, ImagicLisFile
from .protocol_base import ImagicProtocol


class outputs(Enum):
    outputClasses = SetOfClasses2D


class ImagicProtMSAClassify(ProtClassify2D, ImagicProtocol):
    """
    Performs two-dimensional particle classification using the IMAGIC
    MSA-CLASSIFY workflow based on multivariate statistical analysis and
    hierarchical ascendant classification methods.

    AI Generated:

    MSA Classification (ImagicProtMSAClassify) — User Manual
        Overview

        The MSA Classification protocol organizes cryo-EM particle images into
        groups of structurally similar views using statistical classification
        methods implemented in IMAGIC. Its primary objective is to separate
        heterogeneous particle populations into more homogeneous classes that
        can later be interpreted biologically, refined independently, or used
        for quality assessment and downstream reconstruction workflows.

        In single-particle cryo-EM, classification is one of the most important
        stages for identifying structural variability, removing poor-quality
        particles, and improving the interpretability of experimental data.
        This protocol is particularly useful after multivariate statistical
        analysis, where particles have already been projected into a reduced
        dimensional space defined by eigenimages.

        Biological Context and Purpose

        Biological macromolecules often exist in multiple conformational or
        compositional states. Even carefully prepared datasets may contain
        damaged particles, contaminants, aggregation artifacts, or preferred
        orientations. Classification helps distinguish these populations and
        provides a clearer representation of the underlying structural states.

        For biological interpretation, the resulting class averages can reveal
        dominant conformations, identify rare structural states, or highlight
        flexible regions. Well-defined classes frequently improve confidence in
        subsequent three-dimensional analysis and may help determine whether
        the sample contains biologically meaningful heterogeneity.

        Inputs and General Workflow

        The protocol requires particles that have already undergone multivariate
        statistical analysis. These particles are represented within a reduced
        feature space where the main structural variability is captured by
        eigenimages. Classification is then performed using hierarchical
        clustering methods that group particles according to their statistical
        similarity.

        The workflow generates a defined number of particle classes together
        with representative class averages. Each particle becomes associated
        with one class, allowing users to inspect the internal consistency of
        the dataset and identify structurally coherent subsets.

        Choice of Eigenimages

        One of the most biologically important parameters is the number of
        eigenimages used during classification. Eigenimages represent the main
        directions of variability present in the dataset. Selecting too few may
        ignore meaningful structural information, while selecting too many may
        introduce noise and reduce classification stability.

        In practical cryo-EM workflows, the first eigenimages often capture
        major conformational differences or dominant view orientations. Higher
        order eigenimages may progressively contain more noise. Biological users
        generally benefit from testing several values and visually inspecting
        whether class averages become more coherent or more fragmented.

        Number of Classes

        The number of requested classes strongly influences the granularity of
        the final analysis. A small number of classes produces broader,
        more inclusive particle groups that may combine related conformations.
        A larger number separates subtle differences more effectively but may
        also generate sparse or noisy classes.

        From a biological perspective, exploratory analyses often begin with
        moderate class counts to identify major structural trends. Subsequent
        refinement may increase the number of classes to resolve finer
        conformational variability or isolate rare states.

        Ignoring Poor-Quality Particles

        The protocol allows exclusion of a fraction of particles during the
        hierarchical classification process. This option is particularly useful
        for noisy experimental datasets where some particles contribute little
        meaningful structural information.

        Biologically, moderate exclusion values may improve the quality of class
        averages by reducing the influence of damaged particles, contaminants,
        or poorly aligned images. However, excessive exclusion may inadvertently
        remove rare but biologically relevant conformations. Careful visual
        inspection of the resulting classes is therefore recommended.

        Downweighting Small Classes

        Small classes can arise from noise, rare orientations, or genuine
        structural heterogeneity. The protocol optionally reduces the influence
        of these classes during classification to stabilize the overall
        partitioning process.

        In many biological datasets this improves robustness by preventing
        isolated outliers from dominating the classification structure.
        Nevertheless, users studying rare conformational states should evaluate
        the results carefully, since biologically important minor populations
        may also appear as small classes.

        Refinement of Class Quality

        An additional polishing stage can remove the weakest members from each
        class according to their contribution to intra-class variability. This
        process often produces cleaner and more interpretable class averages.

        From a biological perspective, this option is especially valuable when
        preparing class averages for visualization, particle cleaning, or
        subsequent structural refinement. Removing inconsistent particles may
        enhance structural detail and reduce blurring caused by heterogeneity.

        Outputs and Interpretation

        The protocol generates a set of two-dimensional classes together with
        representative averages for each class. These averages summarize the
        common structural features shared by the particles assigned to that
        group.

        Additional quality measurements describing variance and class coherence
        help evaluate the reliability of the classification. Classes with low
        internal variability and clear structural features are generally more
        suitable for downstream refinement and interpretation.

        Biological users should visually inspect class averages to determine
        whether the classes represent meaningful structural states, preferred
        orientations, contaminants, or alignment artifacts.

        Practical Recommendations

        In routine cryo-EM analysis, it is often advisable to begin with a
        moderate number of eigenimages and classes, followed by careful visual
        inspection of the resulting averages. If classes appear overly broad,
        increasing the number of classes may reveal hidden heterogeneity. If
        classes become noisy or fragmented, reducing the dimensionality or class
        count may improve stability.

        Datasets with substantial noise frequently benefit from excluding a
        modest fraction of poor-quality particles or weak class members.
        However, aggressive filtering should be avoided unless the biological
        objective specifically prioritizes highly homogeneous subsets.

        Final Perspective

        For most cryo-EM workflows, classification is not simply a statistical
        grouping operation but a biologically meaningful step that reveals the
        structural organization of the dataset. Appropriate selection of the
        number of eigenimages, class granularity, and quality filtering options
        can substantially influence the interpretation of conformational
        variability and the success of downstream structural analysis.
    """
    _label = 'msa-classify'
    CLASS_DIR = 'MSA-cls'
    _devStatus = PROD
    _possibleOutputs = outputs

    def __init__(self, **kwargs):
        ImagicProtocol.__init__(self, **kwargs)

        self._params = {'cls_dir': self.CLASS_DIR,
                        'msa_cls_img': 'classes'}

# --------------------------- DEFINE param functions --------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputMSA', PointerParam,
                      label="Input particles", important=True,
                      pointerClass='ImagicProtMSA',
                      help='Input images after MSA')
        form.addParam('numberOfFactors', IntParam, default=15,
                      label='Number of eigenimages to use',
                      help='Select the first N eigenimages to use for '
                           'classification.\nTypically all but the first '
                           'few are noisy.')
        form.addParam('numberOfClasses', IntParam, default=10,
                      label='Number of classes',
                      help='Desired number of classes.')
        form.addParam('percentIgnore', IntParam, default=0,
                      expertLevel=LEVEL_ADVANCED,
                      label='Percent of images to ignore',
                      help='This option allows for a percentage of the '
                           'original images to be ignored. The last individual '
                           'images to be merged into a class are set inactive '
                           'in the HAC algorithm. For noisy raw data a value '
                           'of 15% could be tried, for example (this is '
                           'STATISTICS, remember?).')
        form.addParam('doDownweight', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Downweight small classes?',
                      help='A consequence of downweighting small classes '
                           'is that classes with only '
                           'one member will contain only zeroes')
        form.addParam('percentIgnoreBad', IntParam, default=0,
                      expertLevel=LEVEL_ADVANCED,
                      label='Percent of worst class members to ignore',
                      help='Here you get a final chance to polish your '
                           'classes. Since in the CLS file '
                           'the sequence of images in a class is sorted by '
                           'their contribution to the internal variance of '
                           'that class, then  we can enhance the class '
                           'qualities by ignoring the last images of each '
                           'class. The fraction of the images to be ignored '
                           'is what you are supposed to specify here.')

    # --------------------------- INSERT steps functions ----------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('classifyStep', needsGPU=False)
        self._insertFunctionStep('createOutputStep', needsGPU=False)

    # --------------------------- STEPS functions -----------------------------

    def classifyStep(self):
        """ Run MSA-CL and MSA-SUM from IMAGIC. """
        inputFile = self.inputMSA.get().getParticlesStack()
        inputFileBase = pwutils.removeExt(inputFile)
        inputFileImg = inputFileBase + '.img'
        inputFileHed = inputFileBase + '.hed'

        pwutils.createLink(inputFileImg, self._getTmpPath("particles.img"))
        pwutils.createLink(inputFileHed, self._getTmpPath("particles.hed"))
        inputFn = "tmp/particles"

        if self.doDownweight.get():
            downweight = 'YES'
        else:
            downweight = 'NO'

        self._params.update({'particles': inputFn,
                             'eigs_num': self.numberOfFactors.get(),
                             'cls_num': self.numberOfClasses.get(),
                             'perc_ign': self.percentIgnore.get(),
                             'downweight': downweight,
                             'perc_ign_bad': self.percentIgnoreBad.get()
                             })

        classDir = self._getPath(self.CLASS_DIR)
        pwutils.cleanPath(classDir)
        pwutils.makePath(classDir)

        self.runTemplate('msa/msa-cls.b', self._params)

    def createOutputStep(self):
        """ Create the SetOfClass from the cls file with the images-class
        assignment and the averages for each class.
        """
        particles = self.inputMSA.get().inputParticles.get()
        classes2D = self._createSetOfClasses2D(particles)
        # Load the class assignment file from results
        plt = ImagicPltFile(self._getPath(self.CLASS_DIR, 'class_assignment.plt'))
        self._loadClassInfo(self.numberOfClasses.get())

        # Here we are assuming that the order of the class assignment rows
        # is the same for the input particles and the generated img stack
        classes2D.classifyItems(updateItemCallback=self._updateParticle,
                                updateClassCallback=self._updateClass,
                                itemDataIterator=plt.iterRows())

        self._defineOutputs(**{outputs.outputClasses.name: classes2D})
        self._defineSourceRelation(particles, classes2D)

    # --------------------------- INFO functions ------------------------------

    def _validate(self):
        errors = []
        return errors

    def _citations(self):
        return ['vanHeel1984', 'vanHeel1989', 'Borland1990']

    def _summary(self):
        summary = list()
        summary.append('Number of classes: *%s*' % self.numberOfClasses.get())
        summary.append('Number of eigenimages: *%s*' % self.numberOfFactors.get())
        return summary

    def _methods(self):
        msg = "\nInput particles after MSA run were divided into "
        msg += "%s classes by hierarchical ascendant classification (HAC) " %\
               self.numberOfClasses.get()
        msg += "using first %s eigenimages." % self.numberOfFactors.get()
        return [msg]

    # --------------------------- UTILS functions -----------------------------

    def _getOutputPath(self, fn):
        """ Return the output file from the run directory and the CLASS dir. """
        return self._getPath(self.CLASS_DIR, fn)

    def getOutputLis(self):
        return self._getOutputPath('classes.lis')

    def _updateParticle(self, item, row):
        _, classNum = row
        item.setClassId(classNum)

    def _updateClass(self, item):
        classId = item.getObjId()
        avgFile = self._getPath(self.CLASS_DIR,
                                self._params['msa_cls_img'] + '_avg.img')
        rep = item.getRepresentative()
        rep.setSamplingRate(item.getSamplingRate())
        rep.setLocation(classId, avgFile)

        item._intraClassVariance = Float(self.varianceDict[classId])
        item._representationQuality = Float(self.quality1Dict[classId])
        item._overallQuality = Float(self.quality2Dict[classId])

    def _loadClassInfo(self, cls):
        fn = self.getOutputLis()
        self.varianceDict, self.quality1Dict, self.quality2Dict = ImagicLisFile(fn, cls).getParams()
