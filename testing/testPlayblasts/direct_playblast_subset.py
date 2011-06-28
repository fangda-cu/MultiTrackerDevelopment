 
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

    def testcase02( self ):
        """Runs a playblast. Compare the result against a reference"""

        # Prepare input parameters        
        playblastOutputFolder =  self.outputDir()
        playblastFileBase     = 'case2Playblast'
        playblastBase         = os.path.join( playblastOutputFolder, '%s.#.png'%playblastFileBase )
        frameRange            = FrameRange( 1, 24 )
        
	self.shell.call( 'ssn -f weta regression figaro' )
        self.shell.call( 'setenv DISPLAY `hostname`:0' )

	self.shell.call( 'bob_check_opt' )
        self.shell.call( 'need maya-2010_64' )
	self.shell.call( 'need bob-nover' )
	self.shell.call( 'need codetools-nover' )
	self.shell.call( 'need motionDev-nover' )
	self.shell.call( 'need prodeng-official' )
	self.shell.call( 'need wtAlfJob-1.0')
	
        # Run the playblast
        with self.shell.mayaSession( batch=False ) as maya:
            maya.pbTask = mayaTools.getPlayblastTask( 
            	film          = 'weta',
            	scene         = 'regression',
            	shot          = 'figaro',
            	camera        = 'playblastCam',
            	outFolder     = playblastOutputFolder,
            	outBaseName   = playblastFileBase,
            	mayaSceneFile = os.path.abspath( './testing/testCases/case2.mb' ),
            	startTime     = frameRange.start(),
            	endTime       = frameRange.end(),
		widthHeight   = (1024, 768),
	        percent       = 100,
            	wall          = False
            )
            maya.pbTask.execute()

        # Compare the results
        playblastSequence = SequenceInfo( base=playblastBase,             frameRange=frameRange, label='Playblast output' )
        refSequence       = SequenceInfo( base='./testing/testPlayblasts/ref/testCase02/case2Playblast.#.png',frameRange=frameRange, label='Reference playblast' )
        
        self.assertSequencesSimilar( playblastSequence, refSequence, threshold=1e-3 )       
	
    def testcase03( self ):
        """Runs a playblast. Compare the result against a reference"""

        # Prepare input parameters        
        playblastOutputFolder =  self.outputDir()
        playblastFileBase     = 'case3Playblast'
        playblastBase         = os.path.join( playblastOutputFolder, '%s.#.png'%playblastFileBase )
        frameRange            = FrameRange( 1, 24 )
        
	self.shell.call( 'ssn -f weta regression figaro' )
        self.shell.call( 'setenv DISPLAY `hostname`:0' )

	self.shell.call( 'bob_check_opt' )
        self.shell.call( 'need maya-2010_64' )
	self.shell.call( 'need bob-nover' )
	self.shell.call( 'need codetools-nover' )
	self.shell.call( 'need motionDev-nover' )
	self.shell.call( 'need prodeng-official' )
	self.shell.call( 'need wtAlfJob-1.0')
	
        # Run the playblast
        with self.shell.mayaSession( batch=False ) as maya:
            maya.pbTask = mayaTools.getPlayblastTask( 
            	film          = 'weta',
            	scene         = 'regression',
            	shot          = 'figaro',
            	camera        = 'playblastCam',
            	outFolder     = playblastOutputFolder,
            	outBaseName   = playblastFileBase,
            	mayaSceneFile = os.path.abspath( './testing/testCases/case3.mb' ),
            	startTime     = frameRange.start(),
            	endTime       = frameRange.end(),
		widthHeight   = (1024, 768),
	        percent       = 100,
            	wall          = False
            )
            maya.pbTask.execute()

        # Compare the results
        playblastSequence = SequenceInfo( base=playblastBase,             frameRange=frameRange, label='Playblast output' )
        refSequence       = SequenceInfo( base='./testing/testPlayblasts/ref/testCase03/case3Playblast.#.png',frameRange=frameRange, label='Reference playblast' )
        
        self.assertSequencesSimilar( playblastSequence, refSequence, threshold=1e-3 )        	
	
    def testcase04( self ):
        """Runs a playblast. Compare the result against a reference"""

        # Prepare input parameters        
        playblastOutputFolder =  self.outputDir()
        playblastFileBase     = 'case4Playblast'
        playblastBase         = os.path.join( playblastOutputFolder, '%s.#.png'%playblastFileBase )
        frameRange            = FrameRange( 1, 24 )
        
	self.shell.call( 'ssn -f weta regression figaro' )
        self.shell.call( 'setenv DISPLAY `hostname`:0' )

	self.shell.call( 'bob_check_opt' )
        self.shell.call( 'need maya-2010_64' )
	self.shell.call( 'need bob-nover' )
	self.shell.call( 'need codetools-nover' )
	self.shell.call( 'need motionDev-nover' )
	self.shell.call( 'need prodeng-official' )
	self.shell.call( 'need wtAlfJob-1.0')
	
        # Run the playblast
        with self.shell.mayaSession( batch=False ) as maya:
            maya.pbTask = mayaTools.getPlayblastTask( 
            	film          = 'weta',
            	scene         = 'regression',
            	shot          = 'figaro',
            	camera        = 'playblastCam',
            	outFolder     = playblastOutputFolder,
            	outBaseName   = playblastFileBase,
            	mayaSceneFile = os.path.abspath( './testing/testCases/case4.mb' ),
            	startTime     = frameRange.start(),
            	endTime       = frameRange.end(),
		widthHeight   = (1024, 768),
	        percent       = 100,
            	wall          = False
            )
            maya.pbTask.execute()

        # Compare the results
        playblastSequence = SequenceInfo( base=playblastBase,             frameRange=frameRange, label='Playblast output' )
        refSequence       = SequenceInfo( base='./testing/testPlayblasts/ref/testCase04/case4Playblast.#.png',frameRange=frameRange, label='Reference playblast' )
        
        self.assertSequencesSimilar( playblastSequence, refSequence, threshold=1e-3 )       
		
    def testcase05( self ):
        """Runs a playblast. Compare the result against a reference"""

        # Prepare input parameters        
        playblastOutputFolder =  self.outputDir()
        playblastFileBase     = 'case5Playblast'
        playblastBase         = os.path.join( playblastOutputFolder, '%s.#.png'%playblastFileBase )
        frameRange            = FrameRange( 1, 24 )
        
	self.shell.call( 'ssn -f weta regression figaro' )
        self.shell.call( 'setenv DISPLAY `hostname`:0' )

	self.shell.call( 'bob_check_opt' )
        self.shell.call( 'need maya-2010_64' )
	self.shell.call( 'need bob-nover' )
	self.shell.call( 'need codetools-nover' )
	self.shell.call( 'need motionDev-nover' )
	self.shell.call( 'need prodeng-official' )
	self.shell.call( 'need wtAlfJob-1.0')
	
        # Run the playblast
        with self.shell.mayaSession( batch=False ) as maya:
            maya.pbTask = mayaTools.getPlayblastTask( 
            	film          = 'weta',
            	scene         = 'regression',
            	shot          = 'figaro',
            	camera        = 'playblastCam',
            	outFolder     = playblastOutputFolder,
            	outBaseName   = playblastFileBase,
            	mayaSceneFile = os.path.abspath( './testing/testCases/case5.mb' ),
            	startTime     = frameRange.start(),
            	endTime       = frameRange.end(),
		widthHeight   = (1024, 768),
	        percent       = 100,
            	wall          = False
            )
            maya.pbTask.execute()

        # Compare the results
        playblastSequence = SequenceInfo( base=playblastBase,             frameRange=frameRange, label='Playblast output' )
        refSequence       = SequenceInfo( base='./testing/testPlayblasts/ref/testCase05/case5Playblast.#.png',frameRange=frameRange, label='Reference playblast' )
        
        self.assertSequencesSimilar( playblastSequence, refSequence, threshold=1e-3 )       
	
    def testcase05a( self ):
        """Runs a playblast. Compare the result against a reference"""

        # Prepare input parameters        
        playblastOutputFolder =  self.outputDir()
        playblastFileBase     = 'case5aPlayblast'
        playblastBase         = os.path.join( playblastOutputFolder, '%s.#.png'%playblastFileBase )
        frameRange            = FrameRange( 1, 24 )
        
	self.shell.call( 'ssn -f weta regression figaro' )
        self.shell.call( 'setenv DISPLAY `hostname`:0' )

	self.shell.call( 'bob_check_opt' )
        self.shell.call( 'need maya-2010_64' )
	self.shell.call( 'need bob-nover' )
	self.shell.call( 'need codetools-nover' )
	self.shell.call( 'need motionDev-nover' )
	self.shell.call( 'need prodeng-official' )
	self.shell.call( 'need wtAlfJob-1.0')
	
        # Run the playblast
        with self.shell.mayaSession( batch=False ) as maya:
            maya.pbTask = mayaTools.getPlayblastTask( 
            	film          = 'weta',
            	scene         = 'regression',
            	shot          = 'figaro',
            	camera        = 'playblastCam',
            	outFolder     = playblastOutputFolder,
            	outBaseName   = playblastFileBase,
            	mayaSceneFile = os.path.abspath( './testing/testCases/case5a.mb' ),
            	startTime     = frameRange.start(),
            	endTime       = frameRange.end(),
		widthHeight   = (1024, 768),
	        percent       = 100,
            	wall          = False
            )
            maya.pbTask.execute()

        # Compare the results
        playblastSequence = SequenceInfo( base=playblastBase,             frameRange=frameRange, label='Playblast output' )
        refSequence       = SequenceInfo( base='./testing/testPlayblasts/ref/testCase05a/case5aPlayblast.#.png',frameRange=frameRange, label='Reference playblast' )
        
        self.assertSequencesSimilar( playblastSequence, refSequence, threshold=1e-3 )       
	
    def testcase06( self ):
        """Runs a playblast. Compare the result against a reference"""

        # Prepare input parameters        
        playblastOutputFolder =  self.outputDir()
        playblastFileBase     = 'case6Playblast'
        playblastBase         = os.path.join( playblastOutputFolder, '%s.#.png'%playblastFileBase )
        frameRange            = FrameRange( 1, 24 )
        
	self.shell.call( 'ssn -f weta regression figaro' )
        self.shell.call( 'setenv DISPLAY `hostname`:0' )

	self.shell.call( 'bob_check_opt' )
        self.shell.call( 'need maya-2010_64' )
	self.shell.call( 'need bob-nover' )
	self.shell.call( 'need codetools-nover' )
	self.shell.call( 'need motionDev-nover' )
	self.shell.call( 'need prodeng-official' )
	self.shell.call( 'need wtAlfJob-1.0')
	
        # Run the playblast
        with self.shell.mayaSession( batch=False ) as maya:
            maya.pbTask = mayaTools.getPlayblastTask( 
            	film          = 'weta',
            	scene         = 'regression',
            	shot          = 'figaro',
            	camera        = 'playblastCam',
            	outFolder     = playblastOutputFolder,
            	outBaseName   = playblastFileBase,
            	mayaSceneFile = os.path.abspath( './testing/testCases/case6.mb' ),
            	startTime     = frameRange.start(),
            	endTime       = frameRange.end(),
		widthHeight   = (1024, 768),
	        percent       = 100,
            	wall          = False
            )
            maya.pbTask.execute()

        # Compare the results
        playblastSequence = SequenceInfo( base=playblastBase,             frameRange=frameRange, label='Playblast output' )
        refSequence       = SequenceInfo( base='./testing/testPlayblasts/ref/testCase06/case6Playblast.#.png',frameRange=frameRange, label='Reference playblast' )
        
        self.assertSequencesSimilar( playblastSequence, refSequence, threshold=1e-3 )       
	
    def testcase07( self ):
        """Runs a playblast. Compare the result against a reference"""

        # Prepare input parameters        
        playblastOutputFolder =  self.outputDir()
        playblastFileBase     = 'case7Playblast'
        playblastBase         = os.path.join( playblastOutputFolder, '%s.#.png'%playblastFileBase )
        frameRange            = FrameRange( 1, 24 )
        
	self.shell.call( 'ssn -f weta regression figaro' )
        self.shell.call( 'setenv DISPLAY `hostname`:0' )

	self.shell.call( 'bob_check_opt' )
        self.shell.call( 'need maya-2010_64' )
	self.shell.call( 'need bob-nover' )
	self.shell.call( 'need codetools-nover' )
	self.shell.call( 'need motionDev-nover' )
	self.shell.call( 'need prodeng-official' )
	self.shell.call( 'need wtAlfJob-1.0')
	
        # Run the playblast
        with self.shell.mayaSession( batch=False ) as maya:
            maya.pbTask = mayaTools.getPlayblastTask( 
            	film          = 'weta',
            	scene         = 'regression',
            	shot          = 'figaro',
            	camera        = 'playblastCam',
            	outFolder     = playblastOutputFolder,
            	outBaseName   = playblastFileBase,
            	mayaSceneFile = os.path.abspath( './testing/testCases/case7.mb' ),
            	startTime     = frameRange.start(),
            	endTime       = frameRange.end(),
		widthHeight   = (1024, 768),
	        percent       = 100,
            	wall          = False
            )
            maya.pbTask.execute()

        # Compare the results
        playblastSequence = SequenceInfo( base=playblastBase,             frameRange=frameRange, label='Playblast output' )
        refSequence       = SequenceInfo( base='./testing/testPlayblasts/ref/testCase07/case7Playblast.#.png',frameRange=frameRange, label='Reference playblast' )
        
        self.assertSequencesSimilar( playblastSequence, refSequence, threshold=1e-3 )       
	
    def testcase08( self ):
        """Runs a playblast. Compare the result against a reference"""

        # Prepare input parameters        
        playblastOutputFolder =  self.outputDir()
        playblastFileBase     = 'case8Playblast'
        playblastBase         = os.path.join( playblastOutputFolder, '%s.#.png'%playblastFileBase )
        frameRange            = FrameRange( 1, 24 )
        
	self.shell.call( 'ssn -f weta regression figaro' )
        self.shell.call( 'setenv DISPLAY `hostname`:0' )

	self.shell.call( 'bob_check_opt' )
        self.shell.call( 'need maya-2010_64' )
	self.shell.call( 'need bob-nover' )
	self.shell.call( 'need codetools-nover' )
	self.shell.call( 'need motionDev-nover' )
	self.shell.call( 'need prodeng-official' )
	self.shell.call( 'need wtAlfJob-1.0')
	
        # Run the playblast
        with self.shell.mayaSession( batch=False ) as maya:
            maya.pbTask = mayaTools.getPlayblastTask( 
            	film          = 'weta',
            	scene         = 'regression',
            	shot          = 'figaro',
            	camera        = 'playblastCam',
            	outFolder     = playblastOutputFolder,
            	outBaseName   = playblastFileBase,
            	mayaSceneFile = os.path.abspath( './testing/testCases/case8.mb' ),
            	startTime     = frameRange.start(),
            	endTime       = frameRange.end(),
		widthHeight   = (1024, 768),
	        percent       = 100,
            	wall          = False
            )
            maya.pbTask.execute()

        # Compare the results
        playblastSequence = SequenceInfo( base=playblastBase,             frameRange=frameRange, label='Playblast output' )
        refSequence       = SequenceInfo( base='./testing/testPlayblasts/ref/testCase08/case8Playblast.#.png',frameRange=frameRange, label='Reference playblast' )
        
        self.assertSequencesSimilar( playblastSequence, refSequence, threshold=1e-3 )       
