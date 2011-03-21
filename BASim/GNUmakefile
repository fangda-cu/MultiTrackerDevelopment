
BUILD := "$(shell perl -e '$$pwd=$$ENV{PWD}; $$time=localtime(); print "\\\"$$pwd $$time\\\"";')"

##############################################################################
include /vol/bob/make/head.mk
##############################################################################

IPATHS 	+= -I$(MKLROOT)/include 
# This is so that the library can access a header in BASimulator app. Nasty, but that's how the code it came
IPATHS 	+= -I../../../

LIB1       := Wfigaro
#LIB1_CF += -DBUILD=${BUILD} -fPIC -DWETA
#LIB1_CF    += -DBUILD=${BUILD} -pthread -fPIC -DHAVE_LAPACK -DWETA
LIB1_CF    += -DBUILD=${BUILD} -openmp -pthread -DHAVE_MKL -fPIC -DHAVE_LAPACK -DWETA
#LIB1_CF    += -DBUILD=${BUILD} -openmp -pthread -fPIC -DHAVE_LAPACK -DWETA

PHYSBAM_LOCATION := /vol/weta/src/amcadams/PhysBAM-svn
PHYSBAM_INCLUDE  := -I$(PHYSBAM_LOCATION)/Public_Library
LIB1_CF += $(PHYSBAM_INCLUDE) -I../../../

LIB1_HSD   := Wfigaro
LIB1_HDIR := src
LIB1_H     := \
Eigen/Array \
Eigen/Cholesky \
Eigen/Core \
Eigen/Dense \
Eigen/Eigen \
Eigen/Geometry \
Eigen/LU \
Eigen/LeastSquares \
Eigen/NewStdVector \
Eigen/QR \
Eigen/QtAlignedMalloc \
Eigen/SVD \
Eigen/Sparse \
Eigen/StdVector \
Eigen/src/Array/BooleanRedux.h \
Eigen/src/Array/CwiseOperators.h \
Eigen/src/Array/Functors.h \
Eigen/src/Array/Norms.h \
Eigen/src/Array/PartialRedux.h \
Eigen/src/Array/Random.h \
Eigen/src/Array/Select.h \
Eigen/src/Cholesky/LDLT.h \
Eigen/src/Cholesky/LLT.h \
Eigen/src/Core/Assign.h \
Eigen/src/Core/Block.h \
Eigen/src/Core/CacheFriendlyProduct.h \
Eigen/src/Core/Coeffs.h \
Eigen/src/Core/CommaInitializer.h \
Eigen/src/Core/Cwise.h \
Eigen/src/Core/CwiseBinaryOp.h \
Eigen/src/Core/CwiseNullaryOp.h \
Eigen/src/Core/CwiseUnaryOp.h \
Eigen/src/Core/DiagonalCoeffs.h \
Eigen/src/Core/DiagonalMatrix.h \
Eigen/src/Core/DiagonalProduct.h \
Eigen/src/Core/Dot.h \
Eigen/src/Core/Flagged.h \
Eigen/src/Core/Functors.h \
Eigen/src/Core/Fuzzy.h \
Eigen/src/Core/GenericPacketMath.h \
Eigen/src/Core/IO.h \
Eigen/src/Core/Map.h \
Eigen/src/Core/MapBase.h \
Eigen/src/Core/MathFunctions.h \
Eigen/src/Core/Matrix.h \
Eigen/src/Core/MatrixBase.h \
Eigen/src/Core/MatrixStorage.h \
Eigen/src/Core/Minor.h \
Eigen/src/Core/NestByValue.h \
Eigen/src/Core/NumTraits.h \
Eigen/src/Core/Part.h \
Eigen/src/Core/Product.h \
Eigen/src/Core/Redux.h \
Eigen/src/Core/SolveTriangular.h \
Eigen/src/Core/Sum.h \
Eigen/src/Core/Swap.h \
Eigen/src/Core/Transpose.h \
Eigen/src/Core/Visitor.h \
Eigen/src/Core/arch/AltiVec/PacketMath.h \
Eigen/src/Core/arch/SSE/PacketMath.h \
Eigen/src/Core/util/Constants.h \
Eigen/src/Core/util/DisableMSVCWarnings.h \
Eigen/src/Core/util/EnableMSVCWarnings.h \
Eigen/src/Core/util/ForwardDeclarations.h \
Eigen/src/Core/util/Macros.h \
Eigen/src/Core/util/Memory.h \
Eigen/src/Core/util/Meta.h \
Eigen/src/Core/util/StaticAssert.h \
Eigen/src/Core/util/XprHelper.h \
Eigen/src/Geometry/AlignedBox.h \
Eigen/src/Geometry/AngleAxis.h \
Eigen/src/Geometry/EulerAngles.h \
Eigen/src/Geometry/Hyperplane.h \
Eigen/src/Geometry/OrthoMethods.h \
Eigen/src/Geometry/ParametrizedLine.h \
Eigen/src/Geometry/Quaternion.h \
Eigen/src/Geometry/Rotation2D.h \
Eigen/src/Geometry/RotationBase.h \
Eigen/src/Geometry/Scaling.h \
Eigen/src/Geometry/Transform.h \
Eigen/src/Geometry/Translation.h \
Eigen/src/LU/Determinant.h \
Eigen/src/LU/Inverse.h \
Eigen/src/LU/LU.h \
Eigen/src/LeastSquares/LeastSquares.h \
Eigen/src/QR/EigenSolver.h \
Eigen/src/QR/HessenbergDecomposition.h \
Eigen/src/QR/QR.h \
Eigen/src/QR/SelfAdjointEigenSolver.h \
Eigen/src/QR/Tridiagonalization.h \
Eigen/src/SVD/SVD.h \
Eigen/src/Sparse/AmbiVector.h \
Eigen/src/Sparse/CholmodSupport.h \
Eigen/src/Sparse/CompressedStorage.h \
Eigen/src/Sparse/CoreIterators.h \
Eigen/src/Sparse/DynamicSparseMatrix.h \
Eigen/src/Sparse/MappedSparseMatrix.h \
Eigen/src/Sparse/RandomSetter.h \
Eigen/src/Sparse/SparseAssign.h \
Eigen/src/Sparse/SparseBlock.h \
Eigen/src/Sparse/SparseCwise.h \
Eigen/src/Sparse/SparseCwiseBinaryOp.h \
Eigen/src/Sparse/SparseCwiseUnaryOp.h \
Eigen/src/Sparse/SparseDiagonalProduct.h \
Eigen/src/Sparse/SparseDot.h \
Eigen/src/Sparse/SparseFlagged.h \
Eigen/src/Sparse/SparseFuzzy.h \
Eigen/src/Sparse/SparseLDLT.h \
Eigen/src/Sparse/SparseLLT.h \
Eigen/src/Sparse/SparseLU.h \
Eigen/src/Sparse/SparseMatrix.h \
Eigen/src/Sparse/SparseMatrixBase.h \
Eigen/src/Sparse/SparseProduct.h \
Eigen/src/Sparse/SparseRedux.h \
Eigen/src/Sparse/SparseTranspose.h \
Eigen/src/Sparse/SparseUtil.h \
Eigen/src/Sparse/SparseVector.h \
Eigen/src/Sparse/SuperLUSupport.h \
Eigen/src/Sparse/TaucsSupport.h \
Eigen/src/Sparse/TriangularSolver.h \
Eigen/src/Sparse/UmfPackSupport.h \
Collisions/BVHAABB.hh \
Collisions/CollisionUtils.hh \
Collisions/BoundingBox.hh \
Collisions/BVH.hh \
Collisions/BVH_inline.hh \
Collisions/BVHNode.hh \
Collisions/CollisionDetector.hh \
Collisions/Geometry.hh \
Core/Definitions.hh \
Core/EigenIncludes.hh \
Core/Handle.hh \
Core/Macros.hh \
Core/ObjectBase.hh \
Core/ObjectControllerBase.hh \
Core/ObjectHandles.hh \
Core/Property.hh \
Core/STLIncludes.hh \
Core/Timer.hh \
Core/TopologicalObject/TopObjHandles.hh \
Core/TopologicalObject/TopObjIterators.hh \
Core/TopologicalObject/TopologicalObject.hh \
Core/TopologicalObject/Topology.hh \
Core/TriangleMesh.hh \
Core/Option.hh \
Core/ScriptingController.hh \
Core/Util.hh \
IO/SolverIO.hh \
IO/Tokenizer.hh \
IO/ObjParser.hh \
IO/RodTextFileParser.hh \
IO/XMLSceneOutputter.hh \
IO/XMLSceneParser.hh \
IO/rapidxml.hpp \
IO/rapidxml_print.hpp \
IO/ObjectSerializer.hh \
IO/SerializationUtils.hh \
Math/BandMatrix.hh \
Math/ConjugateGradient.hh \
Math/DiagonalPreconditioner.hh \
Math/FlensMatrix.hh \
Math/ImplicitEuler.hh \
Math/LinearSolverBase.hh \
Math/MKL/MKLLinearSolver.hh \
Math/Math.hh \
Math/MatrixBase.hh \
Math/MatrixUtils.hh \
Math/Petsc/PetscLinearSolver.hh \
Math/Petsc/PetscMatrix.hh \
Math/Petsc/PetscUtils.hh \
Math/Preconditioner.hh \
Math/SolverUtils.hh \
Math/StaticsSolver.hh \
Math/SymplecticEuler.hh \
Math/TimeSteppingBase.hh \
Math/Pardiso/PardisoLinearSolver.hh \
Physics/DegreeOfFreedom.hh \
Physics/ElasticRods/BridsonStepper.hh \
Physics/ElasticRods/ParallelStepper.hh \
Physics/ElasticRods/AnisotropicRod.hh \
Physics/ElasticRods/AnisotropicRodIncludes.hh \
Physics/ElasticRods/EdgeStencil.hh \
Physics/ElasticRods/ElasticRod.hh \
Physics/ElasticRods/Force.hh \
Physics/ElasticRods/MultipleRodTimeStepper.hh \
Physics/ElasticRods/RodAnisoBending.hh \
Physics/ElasticRods/RodAnisoForce.hh \
Physics/ElasticRods/RodBendingForce.hh \
Physics/ElasticRods/RodBendingForceSym.hh \
Physics/ElasticRods/RodBendingForceIncludes.hh \
Physics/ElasticRods/RodBoundaryCondition.hh \
Physics/ElasticRods/RodExternalForce.hh \
Physics/ElasticRods/RodForce.hh \
Physics/ElasticRods/RodGravity.hh \
Physics/ElasticRods/RodMassDamping.hh \
Physics/ElasticRods/RodPenaltyForce.hh \
Physics/ElasticRods/RodStretchingForce.hh \
Physics/ElasticRods/RodTimeStepper.hh \
Physics/ElasticRods/RodTube.hh \
Physics/ElasticRods/RodTwistingForce.hh \
Physics/ElasticRods/RodTwistingForceSym.hh \
Physics/ElasticRods/RodUtils.hh \
Physics/ElasticRods/Stencil.hh \
Physics/ElasticRods/VertexStencil.hh \
Physics/ElasticRods/ViscousRod.hh \
Physics/ElasticRods/TopologicalObjectSerializer.hh \
Physics/PhysObject.hh \
Physics/ElasticRods/RodRodExternalForce.hh  \
Physics/ElasticRods/RodStateBackup.hh \
Physics/World.hh \
Render/Camera.hh \
Render/Color.hh \
Render/OpenGLDecl.hh \
Render/OpenGLHeaders.hh \
Render/RenderBase.hh \
Render/RenderUtils.hh \
Render/RodRenderer.hh \
Render/Scene.hh \
Render/TrackBall.hh \
Render/Translator.hh \
Render/ViewController.hh \
Render/YImage/YImage.hh \
Render/Zoomer.hh \
Render/TriangleMeshRenderer.hh \
Util/KeyedValue.hh \
Util/StringToolbox.hh \
Util/Utils.hh \
Math/Petsc/PetscMatrix.inl \
Physics/ElasticRods/ElasticRod.inl \
Physics/ElasticRods/Tests/bendingEnergy.inl \
Physics/ElasticRods/Tests/bendingForce.inl \
Physics/ElasticRods/Tests/bendingJacobian.inl \
Physics/ElasticRods/Tests/stretchingEnergy.inl \
Physics/ElasticRods/Tests/stretchingForce.inl \
Physics/ElasticRods/Tests/stretchingJacobian.inl \
Physics/ElasticRods/Tests/twistingEnergy.inl \
Physics/ElasticRods/Tests/twistingForce.inl \
Physics/ElasticRods/Tests/twistingJacobian.inl \
Physics/ElasticRods/MinimalTriangleMeshBackup.hh \
Physics/PhysObject.inl \
Core/StatTracker.hh \
Math/SymmetricImplicitEuler.hh \
Physics/ElasticRods/AdaptiveBinaryStepper.hh \
Physics/ElasticRods/MinimalRodStateBackup.hh \
Physics/ElasticRods/RodReverseSolver.hh \
Core/TopologicalObject/TopologicalObject.inl \
Core/ObjectControllerBase.inl
#Collisions/UniformGrid.hh \
#Collisions/CollisionMeshData.hh \
#Collisions/LevelSet.hh \
#Collisions/RodneyVector.hh \
#Collisions/BridsonArray.hh \
#Collisions/BridsonArray1.hh \
#Collisions/BridsonArray3.hh \
#Collisions/BridsonUtil.hh \
#Collisions/BridsonVec.hh \
#Collisions/CandidateCollision.hh \
#Collisions/Collision.hh \
#Collisions/CollisionObject.hh \
#SceneXML/tinyxml/tinyxml.h \
#SceneXML/tinyxml/tinystr.h \
#SceneXML/SceneXML.hh \


LIB1_S     := \
src/Collisions/BVHAABB.cc \
src/Collisions/Geometry.cc \
src/Collisions/CollisionUtils.cc \
src/Collisions/CollisionDetector.cc \
src/Core/Util.cc \
src/Core/Timer.cc \
src/Core/TriangleMesh.hh \
src/Core/StatTracker.cc \
src/Core/ScriptingController.cc \
src/IO/ObjParser.cc \
src/IO/ObjectSerializer.cc \
src/IO/SerializationUtils.cc \
src/Math/MKL/MKLLinearSolver.cc \
src/Math/SolverUtils.cc \
src/Math/Pardiso/PardisoLinearSolver.cc \
src/Physics/ElasticRods/AnisotropicRod.cc \
src/Physics/ElasticRods/BridsonStepper.cc \
src/Physics/ElasticRods/ElasticRod.cc \
src/Physics/ElasticRods/ParallelStepper.cc \
src/Physics/ElasticRods/RodAnisoBending.cc \
src/Physics/ElasticRods/RodAnisoForce.cc \
src/Physics/ElasticRods/RodBendingForce.cc \
src/Physics/ElasticRods/RodBendingForceSym.cc \
src/Physics/ElasticRods/RodBoundaryCondition.cc \
src/Physics/ElasticRods/RodForce.cc \
src/Physics/ElasticRods/RodPenaltyForce.cc \
src/Physics/ElasticRods/RodStretchingForce.cc \
src/Physics/ElasticRods/RodTube.cc \
src/Physics/ElasticRods/RodTwistingForce.cc \
src/Physics/ElasticRods/RodTwistingForceSym.cc \
src/Physics/ElasticRods/ViscousRod.cc \
src/Physics/ElasticRods/MinimalTriangleMeshBackup.cc \
src/Physics/ElasticRods/TopologicalObjectSerializer.cc \
src/Physics/ElasticRods/RodStateBackup.cc \
src/Render/Camera.cc \
src/Render/RenderBase.cc \
src/Render/RodRenderer.cc \
src/Render/Scene.cc \
src/Render/TrackBall.cc \
src/Render/Translator.cc \
src/Render/TriangleMeshRenderer.cc \
src/Render/ViewController.cc \
src/Render/YImage/YImage.cc \
src/Render/Zoomer.cc \
src/Physics/ElasticRods/MinimalRodStateBackup.cc \
src/Physics/ElasticRods/RodReverseSolver.cc \
src/Physics/World.cc
#src/Collisions/CollisionMeshData.cc \
#src/Collisions/UniformGrid.cc \
#src/Collisions/LevelSet.cc \
#src/Collisions/CandidateCollision.cc \
#src/Collisions/Collision.cc \
#src/SceneXML/tinyxml/tinyxml.cpp \
#src/SceneXML/tinyxml/tinyxmlerror.cpp \
#src/SceneXML/tinyxml/tinyxmlparser.cpp \
#src/SceneXML/tinyxml/tinystr.cpp \
#src/SceneXML/SceneXML.cc \



##############################################################################
include /vol/bob/make/tail.mk
##############################################################################

# Stuff for Emacs:
# Local Variables:
# mode:makefile
# indent-tabs-mode:t
# tab-width:8
# End:

# Stuff for Vim [by default, this needs to be within the last 5 lines]:
# ex:set filetype=make tabstop=8 softtabstop=8 shiftwidth=8 noexpandtab:

