#ifndef WMFIGAROCMD_HH_
#define WMFIGAROCMD_HH_

#include "WmBunsenNode.hh"
#include "WmBunsenRodNode.hh"

#include <map>
#include <set>
#include <string>
#include <vector>

#include <maya/MArgDatabase.h>
#include <maya/MDagModifier.h>
#include <maya/MObject.h>
#include <maya/MPxCommand.h>
#include <maya/MSyntax.h>
#include <maya/MSelectionList.h>
#include <maya/MStringArray.h>
#include <maya/MPointArray.h>
#include <maya/MItCurveCV.h>
#include <maya/MFnTransform.h>

struct WmBunsenHelp
{
    std::string m_longName;
    std::string m_help;
    std::vector<MSyntax::MArgType> m_argTypes;
};

class WmFigaroCmd : public MPxCommand
{
public:    
    /// Constructor
    WmFigaroCmd();
    /// Destructor
    virtual ~WmFigaroCmd();

    /// Static Maya calls to create an instance of the command
    static void *creator();
    /// Static Maya calls to register command syntax
    static MSyntax syntaxCreator();

    /// Entry point the first time the command is invoked
    MStatus doIt( const MArgList &i_mArgList );
    /// Entry point if the command is redone
    MStatus redoIt();
    /// Entry point if the command is undone
    MStatus undoIt();
    /// Maya queries the command to see if it's undoable or not
    bool isUndoable() const { return m_undoable; }
    /// Maya queries the command to see if it has a registered Maya syntax or not
    bool hasSyntax() const { return true; }

    void getNodes( MSelectionList opt_nodes );

    static MString typeName;

protected: 
    static void p_AddFlag( MSyntax &i_mSyntax,
                           const char *const i_shortName,
                           const char *const i_longName,
                           const char *const i_help,
                           const MSyntax::MArgType i_argType1 = MSyntax::kNoArg,
                           const MSyntax::MArgType i_argType2 = MSyntax::kNoArg,
                           const MSyntax::MArgType i_argType3 = MSyntax::kNoArg,
                           const MSyntax::MArgType i_argType4 = MSyntax::kNoArg,
                           const MSyntax::MArgType i_argType5 = MSyntax::kNoArg,
                           const MSyntax::MArgType i_argType6 = MSyntax::kNoArg);

    MStatus createDagNode( const char *transformName, const char *nodeType, MObject &parentObj, 
                           MObject *transformObjP, MObject *shapeObjP, MDagModifier *iDagModifier,
                           MString& o_shapeName );

    static void printHelp();

    static void getSelectedCurves( const MSelectionList &i_selectionList,
                                   MSelectionList &o_meshList );

    void p_PerformConnect();

    void createWmBunsenRodNode( bool useNURBSInput = true, bool i_previewOnly = false, MObject* o_rodNode = NULL );
    void createWmBunsenNode( MObject &o_wmBunsenNodeObj );
    void addCollisionMeshes();
    void attatchEdgeToObject();

public:     // Data
protected:  // Data
  
    void createPreviewNodes();
    static void appendToResultString( MString& i_resultString );
    void quaternionFromMatrix( MMatrix& a, MQuaternion& Q );
    void setColorOfRod();
    
    /// True if the command is undoable, false otherwise
    bool m_undoable;

    enum Operation {
        Error = 0,
        Help,
        CreateRods,
    };

    // The operation we're supposed to be performing
    Operation m_op;
    
    // Context variables that may (or may not) have been provided on the command line
    MString m_node_name;  // What we want a furnode to be called

    int m_cvsPerRod;

    /// Arguments passed the first time the command is called, used on redo
    MArgDatabase *m_mArgDatabase;

    // DG modified for Undo/Redo
    MDGModifier m_dgmodifier;

    /// Dag Modifier For Undo/Redo
    MDagModifier *m_mDagModifier;

    /// Undo For SetAttr
  //  Wmaya::undo_t *m_undo;

    /// Undo Selection List
    MSelectionList m_undoSelectionList;

    MSelectionList m_nurbsCurveList;    // list of nurbs curves
    MSelectionList m_meshList;      // list of meshes
    MSelectionList m_fozzieNodeList;
    MSelectionList m_figRodNodeList;
    MSelectionList m_allOtherTransformNodesList;
    MObject  m_selectedwmBunsenNode;
    MString m_cacheFile;
    int m_rodNumber;
    int m_edgeNumber;
    MColor m_color;
    
    static MStringArray m_results;
        
    static std::map<std::string, WmBunsenHelp> m_help;
    
};

#endif
