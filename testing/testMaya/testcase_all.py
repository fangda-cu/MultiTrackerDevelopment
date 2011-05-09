from WtTestSuite.case import TestCase

class MayaTest( TestCase ):
	
	def testcase1a( self ):
               
		self.shell.need( 'maya' , '2010_64' )
		self.shell.call('bob_check_opt' )
		print 'set needs'
		
		with self.shell.mayaSession() as maya:
			maya.cmds.file( '/weta/dev/user/showard/shots/maya/scenes/fur/TEST_CASES/case1a.mb', open=True, force=True )
			for i in range(1,23):
				maya.cmds.currentTime( i, edit=True )
				print( 'time %d'%i )
				print maya.cmds.getAttr( 'wmFigaroNode1.syncAttrs' )
			maya.cmds.currentTime( 1, edit=True )
			print maya.cmds.getAttr( 'wmFigaroNode1.syncAttrs' )
			print "Test Case 1a: worked."


	def testcase1b( self ):
               
		self.shell.need( 'maya' , '2010_64' )
		self.shell.call('bob_check_opt' )
		print 'set needs'
		
		with self.shell.mayaSession() as maya:
			maya.cmds.file( '/weta/dev/user/showard/shots/maya/scenes/fur/TEST_CASES/case1b.mb', open=True, force=True )
			for i in range(1,23):
				maya.cmds.currentTime( i, edit=True )
				print( 'time %d'%i )
				print maya.cmds.getAttr( 'wmFigaroNode1.syncAttrs' )
			maya.cmds.currentTime( 1, edit=True )
			print maya.cmds.getAttr( 'wmFigaroNode1.syncAttrs' )
			print "Test Case 1b: worked."


	def testcase2( self ):
               
		self.shell.need( 'maya' , '2010_64' )
		self.shell.call('bob_check_opt' )
		print 'set needs'
		
		with self.shell.mayaSession() as maya:
			maya.cmds.file( '/weta/dev/user/showard/shots/maya/scenes/fur/TEST_CASES/case2.mb', open=True, force=True )
			for i in range(1,23):
				maya.cmds.currentTime( i, edit=True )
				print( 'time %d'%i )
				print maya.cmds.getAttr( 'wmFigaroNode1.syncAttrs' )
			maya.cmds.currentTime( 1, edit=True )
			print maya.cmds.getAttr( 'wmFigaroNode1.syncAttrs' )
			print "Test Case 2: worked."
	
	
	def testcase3( self ):
               
		self.shell.need( 'maya' , '2010_64' )
		self.shell.call('bob_check_opt' )
		print 'set needs'
		
		with self.shell.mayaSession() as maya:
			maya.cmds.file( '/weta/dev/user/showard/shots/maya/scenes/fur/TEST_CASES/case3.mb', open=True, force=True )
			for i in range(1,23):
				maya.cmds.currentTime( i, edit=True )
				print( 'time %d'%i )
				print maya.cmds.getAttr( 'wmFigaroNode1.syncAttrs' )
			maya.cmds.currentTime( 1, edit=True )
			print maya.cmds.getAttr( 'wmFigaroNode1.syncAttrs' )
			print "Test Case 3: worked."
	
	
	def testcase4( self ):
               
		self.shell.need( 'maya' , '2010_64' )
		self.shell.call('bob_check_opt' )
		print 'set needs'
		
		with self.shell.mayaSession() as maya:
			maya.cmds.file( '/weta/dev/user/showard/shots/maya/scenes/fur/TEST_CASES/case4.mb', open=True, force=True )
			for i in range(1,23):
				maya.cmds.currentTime( i, edit=True )
				print( 'time %d'%i )
				print maya.cmds.getAttr( 'wmFigaroNode1.syncAttrs' )
			maya.cmds.currentTime( 1, edit=True )
			print maya.cmds.getAttr( 'wmFigaroNode1.syncAttrs' )
			print "Test Case 4: worked."


	def testcase5( self ):
               
		self.shell.need( 'maya' , '2010_64' )
		self.shell.call('bob_check_opt' )
		print 'set needs'
		
		with self.shell.mayaSession() as maya:
			maya.cmds.file( '/weta/dev/user/showard/shots/maya/scenes/fur/TEST_CASES/case5.mb', open=True, force=True )
			for i in range(1,23):
				maya.cmds.currentTime( i, edit=True )
				print( 'time %d'%i )
				print maya.cmds.getAttr( 'wmFigaroNode1.syncAttrs' )
			maya.cmds.currentTime( 1, edit=True )
			print maya.cmds.getAttr( 'wmFigaroNode1.syncAttrs' )
			print "Test Case 5: worked."


	def testcase5a( self ):
               
		self.shell.need( 'maya' , '2010_64' )
		self.shell.call('bob_check_opt' )
		print 'set needs'
		
		with self.shell.mayaSession() as maya:
			maya.cmds.file( '/weta/dev/user/showard/shots/maya/scenes/fur/TEST_CASES/case5a.mb', open=True, force=True )
			for i in range(1,23):
				maya.cmds.currentTime( i, edit=True )
				print( 'time %d'%i )
				print maya.cmds.getAttr( 'wmFigaroNode1.syncAttrs' )
			maya.cmds.currentTime( 1, edit=True )
			print maya.cmds.getAttr( 'wmFigaroNode1.syncAttrs' )
			print "Test Case 5a: worked."


	def testcase6( self ):
               
		self.shell.need( 'maya' , '2010_64' )
		self.shell.call('bob_check_opt' )
		print 'set needs'
		
		with self.shell.mayaSession() as maya:
			maya.cmds.file( '/weta/dev/user/showard/shots/maya/scenes/fur/TEST_CASES/case6.mb', open=True, force=True )
			for i in range(1,23):
				maya.cmds.currentTime( i, edit=True )
				print( 'time %d'%i )
				print maya.cmds.getAttr( 'wmFigaroNode1.syncAttrs' )
			maya.cmds.currentTime( 1, edit=True )
			print maya.cmds.getAttr( 'wmFigaroNode1.syncAttrs' )
			print "Test Case 6: worked."


	def testcase7( self ):
               
		self.shell.need( 'maya' , '2010_64' )
		self.shell.call('bob_check_opt' )
		print 'set needs'
		
		with self.shell.mayaSession() as maya:
			maya.cmds.file( '/weta/dev/user/showard/shots/maya/scenes/fur/TEST_CASES/case7.mb', open=True, force=True )
			for i in range(1,23):
				maya.cmds.currentTime( i, edit=True )
				print( 'time %d'%i )
				print maya.cmds.getAttr( 'wmFigaroNode1.syncAttrs' )
			maya.cmds.currentTime( 1, edit=True )
			print maya.cmds.getAttr( 'wmFigaroNode1.syncAttrs' )
			print "Test Case 7: worked."
			

	def testcase8( self ):
               
		self.shell.need( 'maya' , '2010_64' )
		self.shell.call('bob_check_opt' )
		print 'set needs'
		
		with self.shell.mayaSession() as maya:
			maya.cmds.file( '/weta/dev/user/showard/shots/maya/scenes/fur/TEST_CASES/case8.mb', open=True, force=True )
			for i in range(1,23):
				maya.cmds.currentTime( i, edit=True )
				print( 'time %d'%i )
				print maya.cmds.getAttr( 'wmFigaroNode1.syncAttrs' )
			maya.cmds.currentTime( 1, edit=True )
			print maya.cmds.getAttr( 'wmFigaroNode1.syncAttrs' )
			print "Test Case 8: worked."	
	
	def testcase9( self ):
               
		self.shell.need( 'maya' , '2010_64' )
		self.shell.call('bob_check_opt' )
		print 'set needs'
		
		with self.shell.mayaSession() as maya:
			maya.cmds.file( '/weta/dev/user/showard/shots/maya/scenes/fur/TEST_CASES/case9.mb', open=True, force=True )
			for i in range(1,23):
				maya.cmds.currentTime( i, edit=True )
				print( 'time %d'%i )
				print maya.cmds.getAttr( 'wmFigaroNode1.syncAttrs' )
			maya.cmds.currentTime( 1, edit=True )
			print maya.cmds.getAttr( 'wmFigaroNode1.syncAttrs' )
			print "Test Case 9: worked."


	def testcase9a( self ):
               
		self.shell.need( 'maya' , '2010_64' )
		self.shell.call('bob_check_opt' )
		print 'set needs'
		
		with self.shell.mayaSession() as maya:
			maya.cmds.file( '/weta/dev/user/showard/shots/maya/scenes/fur/TEST_CASES/case9a.mb', open=True, force=True )
			for i in range(1,23):
				maya.cmds.currentTime( i, edit=True )
				print( 'time %d'%i )
				print maya.cmds.getAttr( 'wmFigaroNode1.syncAttrs' )
			maya.cmds.currentTime( 1, edit=True )
			print maya.cmds.getAttr( 'wmFigaroNode1.syncAttrs' )
			print "Test Case 9a: worked."
			

	def testcase9b( self ):
               
		self.shell.need( 'maya' , '2010_64' )
		self.shell.call('bob_check_opt' )
		print 'set needs'
		
		with self.shell.mayaSession() as maya:
			maya.cmds.file( '/weta/dev/user/showard/shots/maya/scenes/fur/TEST_CASES/case9b.mb', open=True, force=True )
			for i in range(1,23):
				maya.cmds.currentTime( i, edit=True )
				print( 'time %d'%i )
				print maya.cmds.getAttr( 'wmFigaroNode1.syncAttrs' )
			maya.cmds.currentTime( 1, edit=True )
			print maya.cmds.getAttr( 'wmFigaroNode1.syncAttrs' )
			print "Test Case 9b: worked."
			
			
	def testcase10( self ):
               
		self.shell.need( 'maya' , '2010_64' )
		self.shell.call('bob_check_opt' )
		print 'set needs'
		
		with self.shell.mayaSession() as maya:
			maya.cmds.file( '/weta/dev/user/showard/shots/maya/scenes/fur/TEST_CASES/case10.mb', open=True, force=True )
			for i in range(1,23):
				maya.cmds.currentTime( i, edit=True )
				print( 'time %d'%i )
				print maya.cmds.getAttr( 'wmFigaroNode1.syncAttrs' )
			maya.cmds.currentTime( 1, edit=True )
			print maya.cmds.getAttr( 'wmFigaroNode1.syncAttrs' )
			print "Test Case 10: worked."
			
	def testcase10a( self ):
               
		self.shell.need( 'maya' , '2010_64' )
		self.shell.call('bob_check_opt' )
		print 'set needs'
		
		with self.shell.mayaSession() as maya:
			maya.cmds.file( '/weta/dev/user/showard/shots/maya/scenes/fur/TEST_CASES/case10a.mb', open=True, force=True )
			for i in range(1,23):
				maya.cmds.currentTime( i, edit=True )
				print( 'time %d'%i )
				print maya.cmds.getAttr( 'wmFigaroNode1.syncAttrs' )
			maya.cmds.currentTime( 1, edit=True )
			print maya.cmds.getAttr( 'wmFigaroNode1.syncAttrs' )
			print "Test Case 10a: worked."
			
	def testcase11( self ):
               
		self.shell.need( 'maya' , '2010_64' )
		self.shell.call('bob_check_opt' )
		print 'set needs'
		
		with self.shell.mayaSession() as maya:
			maya.cmds.file( '/weta/dev/user/showard/shots/maya/scenes/fur/TEST_CASES/case11.mb', open=True, force=True )
			for i in range(1,23):
				maya.cmds.currentTime( i, edit=True )
				print( 'time %d'%i )
				print maya.cmds.getAttr( 'wmFigaroNode1.syncAttrs' )
			maya.cmds.currentTime( 1, edit=True )
			print maya.cmds.getAttr( 'wmFigaroNode1.syncAttrs' )
			print "Test Case 11: worked."
			
