 
import os

from WtTestSuite.case import TestCase
from WtTestSuite.imageTools import SequenceInfo, FrameRange
from WtTestSuite import mayaTools

class PlayblastTest( TestCase ):
    """
    Playblast example
    
    Demonstrates:
    * WtTestSuite.mayaTools.playblastFile
    * WtTestSuite.case.TestCase.assertSequencesSimilar
    """
    def testcase01( self ):
        """Runs a playblast. Compare the result against a reference"""

        # Prepare input parameters        
        playblastOutputFolder = os.path.join( self.outputDir(), 'playblast' )
        playblastFileBase     = 'case1Playblast'
        playblastBase         = os.path.join( playblastOutputFolder, '%s.#.png'%playblastFileBase )
        frameRange            = FrameRange( 1, 24 )
        
        # Run the playblast
        mayaTools.playblastFile(
            film          = 'weta',
            scene         = 'regression',
            shot          = 'figaro',
            mayaSceneFile = os.path.abspath( './testing/testCases/case1.mb' ),
            camera        = 'playblastCam',
            outFolder     = playblastOutputFolder,
            outBaseName   = playblastFileBase,
            startTime     = frameRange.start(),
            endTime       = frameRange.end(),
            local         = True,
            shell         = self.shell
        )
        
        # Compare the results
        playblastSequence = SequenceInfo( base=playblastBase,             frameRange=frameRange, label='Playblast output' )
        refSequence       = SequenceInfo( base='./testing/testPlayblasts/ref/testCase01/case1Playblast.#.png', frameRange=frameRange, label='Reference playblast' )
        
        self.assertSequencesSimilar( playblastSequence, refSequence, threshold=1e-3 )           
    
    def testcase01a( self ):
        """Runs a playblast. Compare the result against a reference"""

        # Prepare input parameters        
        playblastOutputFolder = os.path.join( self.outputDir(), 'playblast' )
        playblastFileBase     = 'case1aPlayblast'
        playblastBase         = os.path.join( playblastOutputFolder, '%s.#.png'%playblastFileBase )
        frameRange            = FrameRange( 1, 24 )
        
        # Run the playblast
        mayaTools.playblastFile(
            film          = 'weta',
            scene         = 'regression',
            shot          = 'figaro',
            mayaSceneFile = os.path.abspath( './testing/testCases/case1a.mb' ),
            camera        = 'playblastCam',
            outFolder     = playblastOutputFolder,
            outBaseName   = playblastFileBase,
            startTime     = frameRange.start(),
            endTime       = frameRange.end(),
            local         = True,
            shell         = self.shell
        )
        
        # Compare the results
        playblastSequence = SequenceInfo( base=playblastBase,             frameRange=frameRange, label='Playblast output' )
        refSequence       = SequenceInfo( base='./testing/testPlayblasts/ref/testCase01a/case1aPlayblast.#.png', frameRange=frameRange, label='Reference playblast' )
        
        self.assertSequencesSimilar( playblastSequence, refSequence, threshold=1e-3 )           
  
    def testcase01b( self ):
        """Runs a playblast. Compare the result against a reference"""

        # Prepare input parameters        
        playblastOutputFolder = os.path.join( self.outputDir(), 'playblast' )
        playblastFileBase     = 'case1bPlayblast'
        playblastBase         = os.path.join( playblastOutputFolder, '%s.#.png'%playblastFileBase )
        frameRange            = FrameRange( 1, 24 )
        
        # Run the playblast
        mayaTools.playblastFile(
            film          = 'weta',
            scene         = 'regression',
            shot          = 'figaro',
            mayaSceneFile = os.path.abspath( './testing/testCases/case1b.mb' ),
            camera        = 'playblastCam',
            outFolder     = playblastOutputFolder,
            outBaseName   = playblastFileBase,
            startTime     = frameRange.start(),
            endTime       = frameRange.end(),
            local         = True,
            shell         = self.shell
        )
        
        # Compare the results
        playblastSequence = SequenceInfo( base=playblastBase,             frameRange=frameRange, label='Playblast output' )
        refSequence       = SequenceInfo( base='./testing/testPlayblasts/ref/testCase01b/case1bPlayblast.#.png', frameRange=frameRange, label='Reference playblast' )
        
        self.assertSequencesSimilar( playblastSequence, refSequence, threshold=1e-3 )           
  
    
    def testcase02( self ):
        """Runs a playblast. Compare the result against a reference"""

        # Prepare input parameters        
        playblastOutputFolder = os.path.join( self.outputDir(), 'playblast' )
        playblastFileBase     = 'case2Playblast'
        playblastBase         = os.path.join( playblastOutputFolder, '%s.#.png'%playblastFileBase )
        frameRange            = FrameRange( 1, 24 )
        
        # Run the playblast
        mayaTools.playblastFile(
            film          = 'weta',
            scene         = 'regression',
            shot          = 'figaro',
            mayaSceneFile = os.path.abspath( './testing/testCases/case2.mb' ),
            camera        = 'playblastCam',
            outFolder     = playblastOutputFolder,
            outBaseName   = playblastFileBase,
            startTime     = frameRange.start(),
            endTime       = frameRange.end(),
            local         = True,
            shell         = self.shell
        )
        
        # Compare the results
        playblastSequence = SequenceInfo( base=playblastBase,             frameRange=frameRange, label='Playblast output' )
        refSequence       = SequenceInfo( base='./testing/testPlayblasts/ref/testCase02/case2Playblast.#.png', frameRange=frameRange, label='Reference playblast' )
        
        self.assertSequencesSimilar( playblastSequence, refSequence, threshold=1e-3 )       
	
    def testcase03( self ):
        """Runs a playblast. Compare the result against a reference"""

        # Prepare input parameters        
        playblastOutputFolder = os.path.join( self.outputDir(), 'playblast' )
        playblastFileBase     = 'case3Playblast'
        playblastBase         = os.path.join( playblastOutputFolder, '%s.#.png'%playblastFileBase )
        frameRange            = FrameRange( 1, 24 )
        
        # Run the playblast
        mayaTools.playblastFile(
            film          = 'weta',
            scene         = 'regression',
            shot          = 'figaro',
            mayaSceneFile = os.path.abspath( './testing/testCases/case3.mb' ),
            camera        = 'playblastCam',
            outFolder     = playblastOutputFolder,
            outBaseName   = playblastFileBase,
            startTime     = frameRange.start(),
            endTime       = frameRange.end(),
            local         = True,
            shell         = self.shell
        )
        
        # Compare the results
        playblastSequence = SequenceInfo( base=playblastBase,             frameRange=frameRange, label='Playblast output' )
        refSequence       = SequenceInfo( base='./testing/testPlayblasts/ref/testCase03/case3Playblast.#.png', frameRange=frameRange, label='Reference playblast' )
        
        self.assertSequencesSimilar( playblastSequence, refSequence, threshold=1e-3 )

 
    def testcase04( self ):
        """Runs a playblast. Compare the result against a reference"""

        # Prepare input parameters        
        playblastOutputFolder = os.path.join( self.outputDir(), 'playblast' )
        playblastFileBase     = 'case4Playblast'
        playblastBase         = os.path.join( playblastOutputFolder, '%s.#.png'%playblastFileBase )
        frameRange            = FrameRange( 1, 24 )
        
        # Run the playblast
        mayaTools.playblastFile(
            film          = 'weta',
            scene         = 'regression',
            shot          = 'figaro',
            mayaSceneFile = os.path.abspath( './testing/testCases/case4.mb' ),
            camera        = 'playblastCam',
            outFolder     = playblastOutputFolder,
            outBaseName   = playblastFileBase,
            startTime     = frameRange.start(),
            endTime       = frameRange.end(),
            local         = True,
            shell         = self.shell
        )
        
        # Compare the results
        playblastSequence = SequenceInfo( base=playblastBase,             frameRange=frameRange, label='Playblast output' )
        refSequence       = SequenceInfo( base='./testing/testPlayblasts/ref/testCase04/case4Playblast.#.png', frameRange=frameRange, label='Reference playblast' )
        
        self.assertSequencesSimilar( playblastSequence, refSequence, threshold=1e-3 )

    def testcase05( self ):
        """Runs a playblast. Compare the result against a reference"""

        # Prepare input parameters        
        playblastOutputFolder = os.path.join( self.outputDir(), 'playblast' )
        playblastFileBase     = 'case5Playblast'
        playblastBase         = os.path.join( playblastOutputFolder, '%s.#.png'%playblastFileBase )
        frameRange            = FrameRange( 1, 24 )
        
        # Run the playblast
        mayaTools.playblastFile(
            film          = 'weta',
            scene         = 'regression',
            shot          = 'figaro',
            mayaSceneFile = os.path.abspath( './testing/testCases/case5.mb' ),
            camera        = 'playblastCam',
            outFolder     = playblastOutputFolder,
            outBaseName   = playblastFileBase,
            startTime     = frameRange.start(),
            endTime       = frameRange.end(),
            local         = True,
            shell         = self.shell
        )
        
        # Compare the results
        playblastSequence = SequenceInfo( base=playblastBase,             frameRange=frameRange, label='Playblast output' )
        refSequence       = SequenceInfo( base='./testing/testPlayblasts/ref/testCase05/case5Playblast.#.png', frameRange=frameRange, label='Reference playblast' )
        
        self.assertSequencesSimilar( playblastSequence, refSequence, threshold=1e-3 )

    def testcase05a( self ):
        """Runs a playblast. Compare the result against a reference"""

        # Prepare input parameters        
        playblastOutputFolder = os.path.join( self.outputDir(), 'playblast' )
        playblastFileBase     = 'case5aPlayblast'
        playblastBase         = os.path.join( playblastOutputFolder, '%s.#.png'%playblastFileBase )
        frameRange            = FrameRange( 1, 24 )
        
        # Run the playblast
        mayaTools.playblastFile(
            film          = 'weta',
            scene         = 'regression',
            shot          = 'figaro',
            mayaSceneFile = os.path.abspath( './testing/testCases/case5a.mb' ),
            camera        = 'playblastCam',
            outFolder     = playblastOutputFolder,
            outBaseName   = playblastFileBase,
            startTime     = frameRange.start(),
            endTime       = frameRange.end(),
            local         = True,
            shell         = self.shell
        )
        
        # Compare the results
        playblastSequence = SequenceInfo( base=playblastBase,             frameRange=frameRange, label='Playblast output' )
        refSequence       = SequenceInfo( base='./testing/testPlayblasts/ref/testCase05a/case5aPlayblast.#.png', frameRange=frameRange, label='Reference playblast' )
        
        self.assertSequencesSimilar( playblastSequence, refSequence, threshold=1e-3 )

    def testcase06( self ):
        """Runs a playblast. Compare the result against a reference"""

        # Prepare input parameters        
        playblastOutputFolder = os.path.join( self.outputDir(), 'playblast' )
        playblastFileBase     = 'case6Playblast'
        playblastBase         = os.path.join( playblastOutputFolder, '%s.#.png'%playblastFileBase )
        frameRange            = FrameRange( 1, 24 )
        
        # Run the playblast
        mayaTools.playblastFile(
            film          = 'weta',
            scene         = 'regression',
            shot          = 'figaro',
            mayaSceneFile = os.path.abspath( './testing/testCases/case6.mb' ),
            camera        = 'playblastCam',
            outFolder     = playblastOutputFolder,
            outBaseName   = playblastFileBase,
            startTime     = frameRange.start(),
            endTime       = frameRange.end(),
            local         = True,
            shell         = self.shell
        )
        
        # Compare the results
        playblastSequence = SequenceInfo( base=playblastBase,             frameRange=frameRange, label='Playblast output' )
        refSequence       = SequenceInfo( base='./testing/testPlayblasts/ref/testCase06/case6Playblast.#.png', frameRange=frameRange, label='Reference playblast' )
        
        self.assertSequencesSimilar( playblastSequence, refSequence, threshold=1e-3 )

    def testcase07( self ):
        """Runs a playblast. Compare the result against a reference"""

        # Prepare input parameters        
        playblastOutputFolder = os.path.join( self.outputDir(), 'playblast' )
        playblastFileBase     = 'case7Playblast'
        playblastBase         = os.path.join( playblastOutputFolder, '%s.#.png'%playblastFileBase )
        frameRange            = FrameRange( 1, 24 )
        
        # Run the playblast
        mayaTools.playblastFile(
            film          = 'weta',
            scene         = 'regression',
            shot          = 'figaro',
            mayaSceneFile = os.path.abspath( './testing/testCases/case7.mb' ),
            camera        = 'playblastCam',
            outFolder     = playblastOutputFolder,
            outBaseName   = playblastFileBase,
            startTime     = frameRange.start(),
            endTime       = frameRange.end(),
            local         = True,
            shell         = self.shell
        )
        
        # Compare the results
        playblastSequence = SequenceInfo( base=playblastBase,             frameRange=frameRange, label='Playblast output' )
        refSequence       = SequenceInfo( base='./testing/testPlayblasts/ref/testCase07/case7Playblast.#.png', frameRange=frameRange, label='Reference playblast' )
        
        self.assertSequencesSimilar( playblastSequence, refSequence, threshold=1e-3 )

    def testcase08( self ):
        """Runs a playblast. Compare the result against a reference"""

        # Prepare input parameters        
        playblastOutputFolder = os.path.join( self.outputDir(), 'playblast' )
        playblastFileBase     = 'case8Playblast'
        playblastBase         = os.path.join( playblastOutputFolder, '%s.#.png'%playblastFileBase )
        frameRange            = FrameRange( 1, 24 )
        
        # Run the playblast
        mayaTools.playblastFile(
            film          = 'weta',
            scene         = 'regression',
            shot          = 'figaro',
            mayaSceneFile = os.path.abspath( './testing/testCases/case8.mb' ),
            camera        = 'playblastCam',
            outFolder     = playblastOutputFolder,
            outBaseName   = playblastFileBase,
            startTime     = frameRange.start(),
            endTime       = frameRange.end(),
            local         = True,
            shell         = self.shell
        )
        
        # Compare the results
        playblastSequence = SequenceInfo( base=playblastBase,             frameRange=frameRange, label='Playblast output' )
        refSequence       = SequenceInfo( base='./testing/testPlayblasts/ref/testCase08/case8Playblast.#.png', frameRange=frameRange, label='Reference playblast' )
        
        self.assertSequencesSimilar( playblastSequence, refSequence, threshold=1e-3 )

    def testcase09( self ):
        """Runs a playblast. Compare the result against a reference"""

        # Prepare input parameters        
        playblastOutputFolder = os.path.join( self.outputDir(), 'playblast' )
        playblastFileBase     = 'case9Playblast'
        playblastBase         = os.path.join( playblastOutputFolder, '%s.#.png'%playblastFileBase )
        frameRange            = FrameRange( 1, 24 )
        
        # Run the playblast
        mayaTools.playblastFile(
            film          = 'weta',
            scene         = 'regression',
            shot          = 'figaro',
            mayaSceneFile = os.path.abspath( './testing/testCases/case9.mb' ),
            camera        = 'playblastCam',
            outFolder     = playblastOutputFolder,
            outBaseName   = playblastFileBase,
            startTime     = frameRange.start(),
            endTime       = frameRange.end(),
            local         = True,
            shell         = self.shell
        )
        
        # Compare the results
        playblastSequence = SequenceInfo( base=playblastBase,             frameRange=frameRange, label='Playblast output' )
        refSequence       = SequenceInfo( base='./testing/testPlayblasts/ref/testCase09/case9Playblast.#.png', frameRange=frameRange, label='Reference playblast' )
        
        self.assertSequencesSimilar( playblastSequence, refSequence, threshold=1e-3 )

    def testcase09a( self ):
        """Runs a playblast. Compare the result against a reference"""

        # Prepare input parameters        
        playblastOutputFolder = os.path.join( self.outputDir(), 'playblast' )
        playblastFileBase     = 'case9aPlayblast'
        playblastBase         = os.path.join( playblastOutputFolder, '%s.#.png'%playblastFileBase )
        frameRange            = FrameRange( 1, 24 )
        
        # Run the playblast
        mayaTools.playblastFile(
            film          = 'weta',
            scene         = 'regression',
            shot          = 'figaro',
            mayaSceneFile = os.path.abspath( './testing/testCases/case9a.mb' ),
            camera        = 'playblastCam',
            outFolder     = playblastOutputFolder,
            outBaseName   = playblastFileBase,
            startTime     = frameRange.start(),
            endTime       = frameRange.end(),
            local         = True,
            shell         = self.shell
        )
        
        # Compare the results
        playblastSequence = SequenceInfo( base=playblastBase,             frameRange=frameRange, label='Playblast output' )
        refSequence       = SequenceInfo( base='./testing/testPlayblasts/ref/testCase09a/case9aPlayblast.#.png', frameRange=frameRange, label='Reference playblast' )
        
        self.assertSequencesSimilar( playblastSequence, refSequence, threshold=1e-3 )

    def testcase09b( self ):
        """Runs a playblast. Compare the result against a reference"""

        # Prepare input parameters        
        playblastOutputFolder = os.path.join( self.outputDir(), 'playblast' )
        playblastFileBase     = 'case9bPlayblast'
        playblastBase         = os.path.join( playblastOutputFolder, '%s.#.png'%playblastFileBase )
        frameRange            = FrameRange( 1, 24 )
        
        # Run the playblast
        mayaTools.playblastFile(
            film          = 'weta',
            scene         = 'regression',
            shot          = 'figaro',
            mayaSceneFile = os.path.abspath( './testing/testCases/case9b.mb' ),
            camera        = 'playblastCam',
            outFolder     = playblastOutputFolder,
            outBaseName   = playblastFileBase,
            startTime     = frameRange.start(),
            endTime       = frameRange.end(),
            local         = True,
            shell         = self.shell
        )
        
        # Compare the results
        playblastSequence = SequenceInfo( base=playblastBase,             frameRange=frameRange, label='Playblast output' )
        #refSequence       = SequenceInfo( base='./testing/testPlayblasts/ref/testCase09b/case9bPlayblast.#.png', frameRange=frameRange, label='Reference playblast' )
        
        #self.assertSequencesSimilar( playblastSequence, #refSequence, threshold=1e-3 )

    def testcase09c( self ):
        """Runs a playblast. Compare the result against a reference"""

        # Prepare input parameters        
        playblastOutputFolder = os.path.join( self.outputDir(), 'playblast' )
        playblastFileBase     = 'case9cPlayblast'
        playblastBase         = os.path.join( playblastOutputFolder, '%s.#.png'%playblastFileBase )
        frameRange            = FrameRange( 1, 24 )
        
        # Run the playblast
        mayaTools.playblastFile(
            film          = 'weta',
            scene         = 'regression',
            shot          = 'figaro',
            mayaSceneFile = os.path.abspath( './testing/testCases/case9c.mb' ),
            camera        = 'playblastCam',
            outFolder     = playblastOutputFolder,
            outBaseName   = playblastFileBase,
            startTime     = frameRange.start(),
            endTime       = frameRange.end(),
            local         = True,
            shell         = self.shell
        )
        
        # Compare the results
        playblastSequence = SequenceInfo( base=playblastBase,             frameRange=frameRange, label='Playblast output' )
        #refSequence       = SequenceInfo( base='./testing/testPlayblasts/ref/testCase09c/case9cPlayblast.#.png', frameRange=frameRange, label='Reference playblast' )
        
        #self.assertSequencesSimilar( playblastSequence, #refSequence, threshold=1e-3 )

    def testcase10( self ):
        """Runs a playblast. Compare the result against a reference"""

        # Prepare input parameters        
        playblastOutputFolder = os.path.join( self.outputDir(), 'playblast' )
        playblastFileBase     = 'case10Playblast'
        playblastBase         = os.path.join( playblastOutputFolder, '%s.#.png'%playblastFileBase )
        frameRange            = FrameRange( 1, 24 )
        
        # Run the playblast
        mayaTools.playblastFile(
            film          = 'weta',
            scene         = 'regression',
            shot          = 'figaro',
            mayaSceneFile = os.path.abspath( './testing/testCases/case10.mb' ),
            camera        = 'playblastCam',
            outFolder     = playblastOutputFolder,
            outBaseName   = playblastFileBase,
            startTime     = frameRange.start(),
            endTime       = frameRange.end(),
            local         = True,
            shell         = self.shell
        )
        
        # Compare the results
        playblastSequence = SequenceInfo( base=playblastBase,             frameRange=frameRange, label='Playblast output' )
        #refSequence       = SequenceInfo( base='./testing/testPlayblasts/ref/testCase10/case10Playblast.#.png', frameRange=frameRange, label='Reference playblast' )
        
        #self.assertSequencesSimilar( playblastSequence, #refSequence, threshold=1e-3 )

    def testcase10a( self ):
        """Runs a playblast. Compare the result against a reference"""

        # Prepare input parameters        
        playblastOutputFolder = os.path.join( self.outputDir(), 'playblast' )
        playblastFileBase     = 'case10aPlayblast'
        playblastBase         = os.path.join( playblastOutputFolder, '%s.#.png'%playblastFileBase )
        frameRange            = FrameRange( 1, 24 )
        
        # Run the playblast
        mayaTools.playblastFile(
            film          = 'weta',
            scene         = 'regression',
            shot          = 'figaro',
            mayaSceneFile = os.path.abspath( './testing/testCases/case10a.mb' ),
            camera        = 'playblastCam',
            outFolder     = playblastOutputFolder,
            outBaseName   = playblastFileBase,
            startTime     = frameRange.start(),
            endTime       = frameRange.end(),
            local         = True,
            shell         = self.shell
        )
        
        # Compare the results
        playblastSequence = SequenceInfo( base=playblastBase,             frameRange=frameRange, label='Playblast output' )
        #refSequence       = SequenceInfo( base='./testing/testPlayblasts/ref/testCase10a/case10aPlayblast.#.png', frameRange=frameRange, label='Reference playblast' )
        
        #self.assertSequencesSimilar( playblastSequence, #refSequence, threshold=1e-3 )
	
    def testcase10b( self ):
        """Runs a playblast. Compare the result against a reference"""

        # Prepare input parameters        
        playblastOutputFolder = os.path.join( self.outputDir(), 'playblast' )
        playblastFileBase     = 'case10bPlayblast'
        playblastBase         = os.path.join( playblastOutputFolder, '%s.#.png'%playblastFileBase )
        frameRange            = FrameRange( 1, 24 )
        
        # Run the playblast
        mayaTools.playblastFile(
            film          = 'weta',
            scene         = 'regression',
            shot          = 'figaro',
            mayaSceneFile = os.path.abspath( './testing/testCases/case10b.mb' ),
            camera        = 'playblastCam',
            outFolder     = playblastOutputFolder,
            outBaseName   = playblastFileBase,
            startTime     = frameRange.start(),
            endTime       = frameRange.end(),
            local         = True,
            shell         = self.shell
        )
        
        # Compare the results
        playblastSequence = SequenceInfo( base=playblastBase,             frameRange=frameRange, label='Playblast output' )
        #refSequence       = SequenceInfo( base='./testing/testPlayblasts/ref/testCase10b/case10bPlayblast.#.png', frameRange=frameRange, label='Reference playblast' )
        
        #self.assertSequencesSimilar( playblastSequence, #refSequence, threshold=1e-3 )	
	
	
    def testcase11( self ):
        """Runs a playblast. Compare the result against a reference"""

        # Prepare input parameters        
        playblastOutputFolder = os.path.join( self.outputDir(), 'playblast' )
        playblastFileBase     = 'case11Playblast'
        playblastBase         = os.path.join( playblastOutputFolder, '%s.#.png'%playblastFileBase )
        frameRange            = FrameRange( 1, 24 )
        
        # Run the playblast
        mayaTools.playblastFile(
            film          = 'weta',
            scene         = 'regression',
            shot          = 'figaro',
            mayaSceneFile = os.path.abspath( './testing/testCases/case11.mb' ),
            camera        = 'playblastCam',
            outFolder     = playblastOutputFolder,
            outBaseName   = playblastFileBase,
            startTime     = frameRange.start(),
            endTime       = frameRange.end(),
            local         = True,
            shell         = self.shell
        )
        
        # Compare the results
        playblastSequence = SequenceInfo( base=playblastBase,             frameRange=frameRange, label='Playblast output' )
        #refSequence       = SequenceInfo( base='./testing/testPlayblasts/ref/testCase11/case11Playblast.#.png', frameRange=frameRange, label='Reference playblast' )
        
        #self.assertSequencesSimilar( playblastSequence, #refSequence, threshold=1e-3 )		
	