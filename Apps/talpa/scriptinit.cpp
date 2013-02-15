// ---------------------------------------------------------
//
//  ScriptInit.cpp
//  Tyson Brochu 2011
//
//  Parse a script text file to initialize the simulation and mesh objects.
//
// ---------------------------------------------------------


#include "scriptinit.h"

#include "drivers/enrightdriver.h"
#include "drivers/faceoff_multi.h"
#include <fstream>
#include "geometryinit.h"
#include <iomesh.h>
#include "drivers/meancurvature.h"
#include "drivers/normaldriver.h"
#include "drivers/sisccurlnoisedriver.h"
#include <subdivisionscheme.h>

using namespace ElTopo;

#ifdef _MSC_VER
#define snprintf _snprintf
#endif


// ---------------------------------------------------------

void ScriptInit::parse_surftrack_parameters( const ParseTree& surftrack_branch )
{
    
    int use_fraction, perform_improvement, topology_changes, collision_safety, non_manifold;
    std::string subdivision_scheme;
    
    surftrack_branch.get_int( "use_fraction", use_fraction );
    surftrack_branch.get_number( "min_edge_length", surf_track_params.m_min_edge_length );
    surftrack_branch.get_number( "max_edge_length", surf_track_params.m_max_edge_length  );
    surftrack_branch.get_number( "max_volume_change", surf_track_params.m_max_volume_change );
    surftrack_branch.get_number( "min_triangle_angle", surf_track_params.m_min_triangle_angle );   
    surftrack_branch.get_number( "max_triangle_angle", surf_track_params.m_max_triangle_angle );      
    surftrack_branch.get_number( "min_triangle_area", surf_track_params.m_min_triangle_area );
    
    int use_curvature_when_splitting;
    if ( surftrack_branch.get_int( "use_curvature_when_splitting", use_curvature_when_splitting ) )
    {
        surf_track_params.m_use_curvature_when_splitting = ( use_curvature_when_splitting != 0 );
    }
    
    int use_curvature_when_collapsing;
    if ( surftrack_branch.get_int( "use_curvature_when_collapsing", use_curvature_when_collapsing ) )
    {
        surf_track_params.m_use_curvature_when_collapsing = ( use_curvature_when_collapsing != 0 );
    }
    
    surftrack_branch.get_number( "min_curvature_multiplier", surf_track_params.m_min_curvature_multiplier );
    surftrack_branch.get_number( "max_curvature_multiplier", surf_track_params.m_max_curvature_multiplier );
    surftrack_branch.get_number( "merge_proximity", surf_track_params.m_merge_proximity_epsilon );
    surftrack_branch.get_number( "repulsion_proximity", surf_track_params.m_proximity_epsilon );
    surftrack_branch.get_number( "friction_coefficient", surf_track_params.m_friction_coefficient );
    surftrack_branch.get_int( "perform_improvement", perform_improvement );
    surftrack_branch.get_int( "allow_topology_changes", topology_changes );   
    surftrack_branch.get_int( "collision_safety", collision_safety );   
    surftrack_branch.get_string( "subdivision_scheme", subdivision_scheme );
    
    surf_track_params.m_use_fraction = ( use_fraction != 0 );   
    surf_track_params.m_perform_improvement = (perform_improvement != 0);
    surf_track_params.m_allow_topology_changes = (topology_changes != 0);
    surf_track_params.m_collision_safety = (collision_safety != 0);   
    
    if ( surftrack_branch.get_int( "allow_non_manifold", non_manifold ) )
    {
      surf_track_params.m_allow_non_manifold = (non_manifold != 0);
    }
    
    if ( strcmp( subdivision_scheme.c_str(), "butterfly" ) == 0 )
    {
        surf_track_params.m_subdivision_scheme = new ButterflyScheme();
    }
    else if( strcmp(subdivision_scheme.c_str(), "modified" ) == 0 || strcmp(subdivision_scheme.c_str(), "modified_butterfly" ) == 0) {
       surf_track_params.m_subdivision_scheme = new ModifiedButterflyScheme();
    }
    else
    {
        surf_track_params.m_subdivision_scheme = new MidpointScheme();
    }
    
    int allow_vertex_movement;
    if ( surftrack_branch.get_int( "allow_vertex_movement", allow_vertex_movement ) )
    {
        surf_track_params.m_allow_vertex_movement_during_collapse = ( allow_vertex_movement != 0 );
        surf_track_params.m_perform_smoothing = (allow_vertex_movement != 0);
    }
    
}

// ---------------------------------------------------------

void ScriptInit::parse_faceoff( const ParseTree& faceoff_sim_branch )
{
    double speed;

    std::vector<std::vector<double>> speed_matrix(region_count, std::vector<double>(region_count, 0));
    if(faceoff_sim_branch.get_number( "speed", speed ) ) { //default two-material case.
       speed_matrix[0][0] = 0;
       speed_matrix[1][1] = 0;

       speed_matrix[0][1] = speed;
       speed_matrix[1][0] = -speed;
    }
    else {
      
      const Array1d* data = faceoff_sim_branch.get_vector( "speedvector");
      if(data->size() != region_count * region_count)
         std::cout << "Error: speed matrix must have length equal to region_count^2\n";
      for(int i = 0; i < region_count; ++i) {
         for(int j = 0; j < region_count; ++j) {
            speed_matrix[i][j] = data->at(i + j*region_count);
         }
      }
      
      for(int i = 0; i < region_count; ++i) {
         assert(speed_matrix[i][i] == 0);
      }
      //could also check that the matrix is skew-symmetric (opposing values across the diagonal are negative to each other)

    }
    
    double reverse_time = 0; //ignored if it is not set
    bool do_reverse = faceoff_sim_branch.get_number("reverse_time", reverse_time);
    if(do_reverse)
       std::cout << "Motion reversal time set to: " << reverse_time << std::endl;


    //check if non-manifold curve should be stationary, for curling up spheres example
    int nmf_stationary = 0;
    bool either_specified = faceoff_sim_branch.get_int("nmf_stationary", nmf_stationary);
    if(either_specified)
       std::cout << "Non-manifold curve is set to stationary, as in two-sphere cyclical invasion test.\n";

    int smooth_w_all = 1; //default to null-space smoothing based on all branches at non-manifold
    faceoff_sim_branch.get_int("smooth_using_all", smooth_w_all);
    bool smooth_using_all = (smooth_w_all == 1);
    if(!smooth_using_all)
       std::cout << "Performing null-space smoothing only on data from the expanding surface at non-manifold vertices.\n";

    //check which is the expanding surface for offsetting in the case where the intersection curve is moving
    int expanding_surf = -1;
    either_specified = either_specified || faceoff_sim_branch.get_int("expanding_surface", expanding_surf);
    if(!either_specified) {
       std::cout << "Must specify either stationary curve or one expanding surface.\n";
       exit(0);
    }
    if(expanding_surf != -1)
       std::cout << "Region " << expanding_surf << " is set to be expanded into.";

    FaceOffMultiDriver* d = new FaceOffMultiDriver( speed_matrix, expanding_surf, nmf_stationary == 1, smooth_using_all );
    if(do_reverse) d->set_reversing(reverse_time);

    //set solution data for the original normal flow tests.
    d->set_solution_data(Vec3d( -0.25, 0.0, 0.0 ), Vec3d( 0.25, 0.0, 0.0 ), 0.4, 0.2);

    driver = d;
    
}

// ---------------------------------------------------------

void ScriptInit::parse_normal( const ParseTree& normal_sim_branch )
{
    double speed;
    normal_sim_branch.get_number( "speed", speed );
    driver = new NormalDriver( speed, Vec3d( -0.25, 0.0, 0.0 ), Vec3d( 0.25, 0.0, 0.0 ), 0.4, 0.2 );
}

// ---------------------------------------------------------

void ScriptInit::parse_mean_curvature( const ParseTree& mean_curvature_sim_branch )
{
    double speed;
    mean_curvature_sim_branch.get_number( "speed", speed );
    
    std::string ground_truth_file;
    mean_curvature_sim_branch.get_string( "ground_truth_file", ground_truth_file );
    
    Array3d sethian_final;
    read_signed_distance( ground_truth_file.c_str(), sethian_final );       
    
    Vec3d phi_domain_low;   
    mean_curvature_sim_branch.get_vec3d( "phi_domain_low", phi_domain_low );
    
    double phi_domain_dx;
    mean_curvature_sim_branch.get_number( "phi_domain_dx", phi_domain_dx );
    
    driver = new MeanCurvatureDriver( speed, sethian_final, phi_domain_low, phi_domain_dx );
}

// ---------------------------------------------------------

void ScriptInit::parse_sisc_curl_noise( const ParseTree&  )
{
    driver = new SISCCurlNoiseDriver( );
}

// ---------------------------------------------------------
void ScriptInit::parse_enright( const ParseTree&  )
{
    driver = new EnrightDriver( );
}

// ---------------------------------------------------------

void ScriptInit::parse_camera( const ParseTree& camera_branch )
{
    camera_branch.get_vec3d( "target", camera_target );
    camera_branch.get_number( "distance", camera_distance );
    camera_branch.get_number( "heading", camera_heading );   
    camera_branch.get_number( "pitch", camera_pitch );
}

// ---------------------------------------------------------

void ScriptInit::parse_sheet( const ParseTree& sheet_branch )
{
    Vec3d lower_corner;
    sheet_branch.get_vec3d( "corner", lower_corner );
    
    double sheet_dx;
    sheet_branch.get_number( "dx", sheet_dx );
    
    int sheet_ni, sheet_nj;
    sheet_branch.get_int( "ni", sheet_ni );
    sheet_branch.get_int( "nj", sheet_nj );
    
    std::vector<Vec3d> sheet_verts;
    std::vector<Vec3st> sheet_tris;
    
    create_sheet( sheet_verts, sheet_tris, lower_corner, sheet_dx, sheet_ni, sheet_nj );      
    
    
    Vec3d rotate_axis;
    double rotate_radians;
    if ( sheet_branch.get_vec3d( "rotate_axis", rotate_axis ) )
    {
        if ( sheet_branch.get_number( "rotate_radians", rotate_radians ) )
        {
            Vec3d centre(0,0,0);
            for ( unsigned int i = 0; i < sheet_verts.size(); ++i )
            {
                centre += sheet_verts[i];
            }
            centre /= static_cast<double>( sheet_verts.size() );
            
            for ( unsigned int i = 0; i < sheet_verts.size(); ++i )
            {
                sheet_verts[i] -= centre;
            }
            
            for ( unsigned int i = 0; i < sheet_verts.size(); ++i )
            {
                sheet_verts[i] = rotate( sheet_verts[i], rotate_radians, rotate_axis );
                sheet_verts[i] += centre;
            }
            
        }
    }
    
    
    int is_solid = 0;
    sheet_branch.get_int( "is_solid", is_solid );  
    
    std::vector<double> sheet_masses;
    
    if ( is_solid )
    {
        sheet_masses.resize( sheet_verts.size(), std::numeric_limits<double>::infinity() );      
    }
    else
    {
        sheet_masses.resize( sheet_verts.size(), 1.0 );
    }
    
    append_mesh( triangles, vertices, labels, masses, sheet_tris, sheet_verts, sheet_masses );
    
}


// ---------------------------------------------------------

void ScriptInit::parse_curved_sheet( const ParseTree& sheet_branch )
{
    Vec3d lower_corner;
    sheet_branch.get_vec3d( "corner", lower_corner );
    
    double sheet_dx;
    sheet_branch.get_number( "dx", sheet_dx );
    
    int sheet_ni, sheet_nj;
    sheet_branch.get_int( "ni", sheet_ni );
    sheet_branch.get_int( "nj", sheet_nj );
    
    std::vector<Vec3d> sheet_verts;
    std::vector<Vec3st> sheet_tris;
    
    create_curved_sheet( sheet_verts, sheet_tris, lower_corner, sheet_dx, sheet_ni, sheet_nj );      
    
    Vec3d rotate_axis;
    double rotate_radians;
    if ( sheet_branch.get_vec3d( "rotate_axis", rotate_axis ) )
    {
        if ( sheet_branch.get_number( "rotate_radians", rotate_radians ) )
        {
            Vec3d centre(0,0,0);
            for ( unsigned int i = 0; i < sheet_verts.size(); ++i )
            {
                centre += sheet_verts[i];
            }
            centre /= static_cast<double>( sheet_verts.size() );
            
            for ( unsigned int i = 0; i < sheet_verts.size(); ++i )
            {
                sheet_verts[i] -= centre;
            }
            
            for ( unsigned int i = 0; i < sheet_verts.size(); ++i )
            {
                sheet_verts[i] = rotate( sheet_verts[i], rotate_radians, rotate_axis );
                sheet_verts[i] += centre;
            }
            
        }
    }
    
    
    std::vector<double> sheet_masses( sheet_verts.size(), 1.0 );
    
    append_mesh( triangles, vertices, labels, masses, sheet_tris, sheet_verts, sheet_masses );
    
}

// ---------------------------------------------------------

void ScriptInit::parse_objfile( const ParseTree& obj_branch) {
   std::string meshpath;
   obj_branch.get_string("filepath", meshpath);
   printf("Got path: %s\n", meshpath.c_str());

   NonDestructiveTriMesh trimesh;

   printf("Reading file\n");

   std::vector<Vec3d> obj_vertices;

   std::vector<Vec3st> obj_triangles;
   read_objfile( obj_triangles, obj_vertices, meshpath.c_str() );

   Vec3d translate;
   if ( obj_branch.get_vec3d("translate", translate) )
   {
      for ( size_t i = 0; i < obj_vertices.size(); ++i )
      {
         obj_vertices[i] += translate;
      }
   }

   std::vector<double> obj_masses(0);
   int is_solid = 0;
   obj_branch.get_int( "is_solid", is_solid );

   if ( is_solid )
   {
      obj_masses.resize( obj_vertices.size(), std::numeric_limits<double>::infinity() );
   }
   else
   {
      obj_masses.resize( obj_vertices.size(), 1.0 );
   }

   int in_label = 0, out_label = 1; //default to the usual thing.
   obj_branch.get_int( "in_label", in_label );
   obj_branch.get_int( "out_label", out_label );
   std::vector<Vec2i> obj_labels(obj_triangles.size(), Vec2i(in_label, out_label));
   std::cout << "read mesh with " << triangles.size() << " triangles and " << vertices.size() << " vertices.\n";
   append_mesh( triangles, vertices, labels, masses, obj_triangles, obj_vertices, obj_labels, obj_masses );
}

// ---------------------------------------------------------

void ScriptInit::parse_sphere( const ParseTree& sphere_branch )
{
    Vec3d sphere_center;
    sphere_branch.get_vec3d( "sphere_center", sphere_center );
    double sphere_radius;
    sphere_branch.get_number( "sphere_radius", sphere_radius );
    
    double dx;
    sphere_branch.get_number( "sphere_dx", dx );
    
    int is_solid = 0;
    sphere_branch.get_int( "is_solid", is_solid );
    
    int in_label = 0, out_label = 1; //default to the usual thing.
    sphere_branch.get_int( "in_label", in_label );
    sphere_branch.get_int( "out_label", out_label );

    std::vector<Vec3d> sphere_vertices;
    std::vector<Vec3st> sphere_triangles;
    std::vector<Vec2i> sphere_labels;

    create_sphere( sphere_center, sphere_radius, dx, sphere_vertices, sphere_triangles );
    sphere_labels.resize( sphere_triangles.size(), Vec2i(in_label, out_label) );

    std::vector<Vec3d> sphere_velocities( sphere_vertices.size(), Vec3d(0) );
    
    std::vector<double> sphere_masses;
    if ( is_solid == 0 )
    {
        sphere_masses.resize( sphere_vertices.size(), 1.0 );
    }
    else
    {
        sphere_masses.resize( sphere_vertices.size(), std::numeric_limits<double>::infinity() );         
    }
    
    append_mesh( triangles, vertices, labels, masses, sphere_triangles, sphere_vertices, sphere_labels, sphere_masses );
}

// ---------------------------------------------------------

void ScriptInit::parse_doublebubble( const ParseTree& bubble_branch )
{
   
   //create attached two-bubble scenario, ported from Fang's BASim code.

   std::vector<Vec3d> bubble_vertices;
   std::vector<Vec3st> bubble_triangles;
   std::vector<Vec2i> bubble_labels;
   
   int N = 20;
   double r = 2.0;

   Vec3d c1 = Vec3d(0.5, 0.5, 0.5 - r * 0.707);
   
   
   bubble_vertices.push_back(Vec3d(c1 + Vec3d(0, 0, r)));

   for (int j = 0; j < N - 1; j++)
   {
      for (int i = 0; i < N * 2; i++)
      {
         
         double theta = (double)i * 2 * M_PI / (N * 2);
         double alpha = (double)(j + 1) * M_PI / N - M_PI / 2;
         bubble_vertices.push_back(c1 + r * Vec3d(cos(alpha) * cos(theta), cos(alpha) * sin(theta), -sin(alpha)));
      }
   }

   bubble_vertices.push_back(Vec3d(c1 - Vec3d(0, 0, r)));
   bubble_vertices.push_back(Vec3d(0.5, 0.5, 0.5));

   Vec3d c2 = Vec3d(0.5, 0.5, 0.5 + r * 0.707);
   
   bubble_vertices.push_back(Vec3d(c2 - Vec3d(0, 0, r)));
   for (int j = 0; j < N - 1; j++)
   {
      for (int i = 0; i < N * 2; i++)
      {
         double theta = (double)i * 2 * M_PI / (N * 2);
         double alpha = (double)(j + 1) * M_PI / N - M_PI / 2;
         bubble_vertices.push_back(c2 + r * Vec3d(cos(alpha) * cos(theta), cos(alpha) * sin(theta), sin(alpha)));
      }
   }
   bubble_vertices.push_back( Vec3d(c2 + Vec3d(0, 0, r)));

   
   
   for (int j = N / 4; j < N; j++)
   {
      for (int i = 0; i < N * 2; i++)
      {
         int v0, v1, v2;
         v0 = (j == 0 ? 0 : 2 * N * (j - 1) + i + 1);
         v1 = (j == 0 ? 0 : 2 * N * (j - 1) + (i + 1) % (N * 2) + 1);
         v2 = (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + (i + 1) % (N * 2) + 1);
         if (!(v0 == v1 || v0 == v2 || v1 == v2))
         {
            bubble_triangles.push_back(Vec3st(v0,v1,v2));
            bubble_labels.push_back(Vec2i(0, 1));
         }

         v0 = (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + (i + 1) % (N * 2) + 1);
         v1 = (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + i + 1);
         v2 = (j == 0 ? 0 : 2 * N * (j - 1) + i + 1);
         if (!(v0 == v1 || v0 == v2 || v1 == v2))
         {
            bubble_triangles.push_back(Vec3st(v0,v1,v2));
            bubble_labels.push_back(Vec2i(0, 1));
         }
      }
   }

   int offset = 2 * N * (N - 1) + 3;
   for (int j = N / 4 + 1; j < N; j++)
   {
      for (int i = 0; i < N * 2; i++)
      {
         int v0, v1, v2;
         v0 = offset + (j == 0 ? 0 : 2 * N * (j - 1) + i + 1);
         v1 = offset + (j == 0 ? 0 : 2 * N * (j - 1) + (i + 1) % (N * 2) + 1);
         v2 = offset + (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + (i + 1) % (N * 2) + 1);
         if (!(v0 == v1 || v0 == v2 || v1 == v2))
         {
            bubble_triangles.push_back(Vec3st(v0,v1,v2));
            bubble_labels.push_back(Vec2i(2, 0));
         }

         v0 = offset + (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + (i + 1) % (N * 2) + 1);
         v1 = offset + (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + i + 1);
         v2 = offset + (j == 0 ? 0 : 2 * N * (j - 1) + i + 1);
         if (!(v0 == v1 || v0 == v2 || v1 == v2))
         {
            bubble_triangles.push_back(Vec3st(v0,v1,v2));
            bubble_labels.push_back(Vec2i(2, 0));
         }
      }
   }

   for (int i = 0; i < N * 2; i++)
   {
      int v0, v1, v2;
      v0 = 2 * N * (N / 4 - 1) + i + 1;
      v1 = 2 * N * (N / 4 - 1) + (i + 1) % (N * 2) + 1;
      v2 = offset + 2 * N * (N / 4) + (i + 1) % (N * 2) + 1;
      if (!(v0 == v1 || v0 == v2 || v1 == v2))
      {
         bubble_triangles.push_back(Vec3st(v0,v1,v2));
         bubble_labels.push_back(Vec2i(2, 0));
      }

      v0 = offset + 2 * N * (N / 4) + (i + 1) % (N * 2) + 1;
      v1 = offset + 2 * N * (N / 4) + i + 1;
      v2 = 2 * N * (N / 4 - 1) + i + 1;
      if (!(v0 == v1 || v0 == v2 || v1 == v2))
      {
         bubble_triangles.push_back(Vec3st(v0,v1,v2));
         bubble_labels.push_back(Vec2i(2, 0));
      }
   }

   for (int i = 0; i < N * 2; i++)
   {
      int v0, v1, v2;
      v0 = 2 * N * (N / 4 - 1) + i + 1;
      v1 = 2 * (N - 1) * N + 2;
      v2 = 2 * N * (N / 4 - 1) + (i + 1) % (N * 2) + 1;
      if (!(v0 == v1 || v0 == v2 || v1 == v2))
      {
         bubble_triangles.push_back(Vec3st(v0,v1,v2));
         bubble_labels.push_back(Vec2i(2, 1));
      }
   }

   /*for (VertexIterator vit = shellObj->vertices_begin(); vit != shellObj->vertices_end(); ++vit)
      if (shellObj->vertexIncidentEdges(*vit) == 0)
         shellObj->deleteVertex(*vit);*/

   std::vector<double> bubble_masses;
   bubble_masses.resize( bubble_vertices.size(), 1.0 );
   
   append_mesh( triangles, vertices, labels, masses, bubble_triangles, bubble_vertices, bubble_labels, bubble_masses );
}

// ---------------------------------------------------------



void ScriptInit::parse_dumbbell( const ParseTree& dumbbell_branch )
{
    double domain_dx;
    dumbbell_branch.get_number( "domain_dx", domain_dx );
    
    Vec3d centre_a;
    dumbbell_branch.get_vec3d( "centre_a", centre_a );
    
    Vec3d centre_b;
    dumbbell_branch.get_vec3d( "centre_b", centre_b );
    
    double sphere_radius;
    dumbbell_branch.get_number( "sphere_radius", sphere_radius );
    
    double handle_width;
    dumbbell_branch.get_number( "handle_width", handle_width );   
    
    Vec3d domain_low = min_union( centre_a, centre_b ) - 2.0 * Vec3d( sphere_radius );
    Vec3d domain_high = max_union( centre_a, centre_b ) + 2.0 * Vec3d( sphere_radius );
    
    Array3d phi;   
    create_dumbbell_signed_distance( centre_a, centre_b, sphere_radius, handle_width, domain_dx, domain_low, domain_high, phi );      
    
    std::vector<Vec3st> new_tris;
    std::vector<Vec3d> new_verts;
    contour_phi( domain_low, domain_dx, phi, new_tris, new_verts );
    project_to_exact_dumbbell( new_verts, centre_a, centre_b, sphere_radius, handle_width );
    
    std::vector<double> new_masses( new_verts.size(), 1.0 );
    std::vector<Vec3d> new_velocities( new_verts.size(), Vec3d(0,0,0) );
    append_mesh( triangles, vertices, labels, masses, new_tris, new_verts, new_masses );
}


// ---------------------------------------------------------

void ScriptInit::parse_script( const char* filename )
{
    
    std::ifstream filestream( filename );
    if ( !filestream.good() )
    {
        std::cerr << "Could not open script file" << std::endl;
        exit(1);
    }
    
    std::cout << "script file: " << filename << std::endl;
    
    ParseTree tree;
    parse_stream( filestream, tree );
    
    
    //
    // Frame stepper
    //
    
    bool ok = tree.get_number( "frame_dt", frame_dt );
    assert( ok );
    
    int num_substeps;
    bool substeps_specified = tree.get_int( "num_substeps", num_substeps );
    if ( substeps_specified )
    {
        sim_dt = frame_dt / (double) num_substeps;
    }
    
    double read_sim_dt;
    if ( tree.get_number( "sim_dt", read_sim_dt ) )
    {
        if ( substeps_specified )
        {
            std::cerr << "WARNING: Both sim_dt and num_substeps specified in config script.  Going with sim_dt." << std::endl;
        }
        
        sim_dt = read_sim_dt;
    }
    
    tree.get_number( "end_sim_t", end_sim_t );
    
    curr_t_specified = tree.get_number( "curr_t", curr_t );
    
    tree.get_int( "region_count", region_count);
    //
    // File output
    //
    
    if ( tree.get_string( "output_path", output_path ) )
    {
        output_path_is_relative = false;
    }
    else if ( tree.get_string( "relative_output_path", output_path ) )
    {
        output_path_is_relative = true;
    }         
    else
    {
        // no path specified
        output_path_is_relative = false;
        output_path = std::string( "./" );
    }
    
    //
    // OpenGL camera
    //
    
    const ParseTree* camera_branch = tree.get_branch( "camera" );
    parse_camera( *camera_branch );
    
    
    //
    // Surface geometry
    //
    
    {
        unsigned int curved_sheet_n = 0;
        char curved_sheet_name[256];
        snprintf( curved_sheet_name, 256, "curved_sheet%d", curved_sheet_n );
        const ParseTree* curved_sheet_branch = tree.get_branch( curved_sheet_name );
        
        while ( curved_sheet_branch != NULL )
        {
            parse_curved_sheet( *curved_sheet_branch );
            curved_sheet_branch = NULL;
            ++curved_sheet_n;
            snprintf( curved_sheet_name, 256, "curved_sheet%d", curved_sheet_n );
            curved_sheet_branch = tree.get_branch( curved_sheet_name );
        }
    }
    
    {
        unsigned int sheet_n = 0;   
        char sheet_name[256];
        snprintf( sheet_name, 256, "sheet%d", sheet_n );
        const ParseTree* sheet_branch = tree.get_branch( sheet_name );
        
        while ( sheet_branch != NULL )
        {
            parse_sheet( *sheet_branch );      
            sheet_branch = NULL;      
            ++sheet_n;
            snprintf( sheet_name, 256, "sheet%d", sheet_n );
            sheet_branch = tree.get_branch( sheet_name );
        }
    }
    
    {
        unsigned int sphere_n = 0;   
        char sphere_name[256];
        snprintf( sphere_name, 256, "sphere%d", sphere_n );
        const ParseTree* sphere_branch = tree.get_branch( sphere_name );
        
        while ( sphere_branch != NULL )
        {
            parse_sphere( *sphere_branch );
            sphere_branch = NULL;
            ++sphere_n;
            snprintf( sphere_name, 256, "sphere%d", sphere_n );
            sphere_branch = tree.get_branch( sphere_name );
        }
    }

    {
      const ParseTree * bubble_branch = tree.get_branch("doublebubble");
      if(bubble_branch != NULL)
         parse_doublebubble( *bubble_branch );
      
    }
    
    
    const ParseTree* sphere_branch = tree.get_branch( "sphere" );
    if ( sphere_branch != NULL )
    {
        parse_sphere( *sphere_branch );
    }      
    
    const ParseTree* dumbbell_branch = tree.get_branch( "dumbbell" );
    if ( dumbbell_branch != NULL )
    {
        parse_dumbbell( *dumbbell_branch );
    }      
    
    const ParseTree* trimesh_branch = tree.get_branch( "trimesh" );
    if ( trimesh_branch != NULL )
    {
        printf("Found trimesh branch\n");
        
        std::string meshpath;
        trimesh_branch->get_string("filepath", meshpath);
        printf("Got path: %s\n", meshpath.c_str());
        
        NonDestructiveTriMesh trimesh;
        
        printf("Reading file\n");
        
        std::vector<Vec3d> input_vertices;
        std::vector<double> in_masses;
        
        read_binary_file( trimesh, input_vertices, in_masses, curr_t, meshpath.c_str() );
        curr_t_specified = true;
        
        int is_solid = 0;
        if ( trimesh_branch->get_int( "is_solid", is_solid ) )
        {
            in_masses.clear();
            if ( is_solid )
            {
                in_masses.resize( input_vertices.size(), std::numeric_limits<double>::infinity() );
            }
            else
            {
                in_masses.resize( input_vertices.size(), 1.0 );
            }
        }
        
        Vec3d translate;
        if ( trimesh_branch->get_vec3d("translate", translate) )
        {
            for ( size_t i = 0; i < input_vertices.size(); ++i )
            {
                input_vertices[i] += translate;
            }
        }
                
        append_mesh( triangles, vertices, labels, masses, trimesh.get_triangles(), input_vertices, in_masses );        
        
        printf("loaded file %s", meshpath.c_str());
    }
    
    {
       unsigned int obj_n = 0;   
       char obj_name[256];
       snprintf( obj_name, 256, "objfile%d", obj_n );
       const ParseTree* obj_branch = tree.get_branch( obj_name );

       while ( obj_branch != NULL )
       {
        printf("Found obj branch\n");
        parse_objfile( *obj_branch );
          obj_branch = NULL;
          ++obj_n;
          snprintf( obj_name, 256, "objfile%d", obj_n );
          obj_branch = tree.get_branch( obj_name );
       }
    
    }
    //
    // SurfTrack parameters
    //
    
    const ParseTree* surftrack_branch = tree.get_branch( "surftrack_parameters" );
    parse_surftrack_parameters( *surftrack_branch );
    
    
    //
    // Mesh drivers
    //
    
    const ParseTree* faceoff_sim_branch = tree.get_branch( "faceoff_simulation" );
    if ( faceoff_sim_branch != NULL )
    {
        parse_faceoff( *faceoff_sim_branch );
    }
    
    const ParseTree* normal_sim_branch = tree.get_branch( "normal_simulation" );
    if ( normal_sim_branch != NULL )
    {
        parse_normal( *normal_sim_branch );
    }
    
    const ParseTree* mean_curvature_sim_branch = tree.get_branch( "mean_curvature_simulation" );
    if ( mean_curvature_sim_branch != NULL )
    {
        parse_mean_curvature( *mean_curvature_sim_branch );
    }
    
    const ParseTree* sisc_curl_noise_sim_branch = tree.get_branch( "sisc_curl_noise_simulation" );
    if ( sisc_curl_noise_sim_branch != NULL )
    {
        parse_sisc_curl_noise( *sisc_curl_noise_sim_branch );
    }
    
    const ParseTree* enright_branch = tree.get_branch( "enright_simulation" );
    if ( enright_branch != NULL )
    {
        parse_enright( *enright_branch );
    }
        
}

