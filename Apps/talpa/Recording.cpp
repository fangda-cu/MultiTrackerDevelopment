//
//  Recording.cpp
//  BASim
//
//  Created by Fang Da on 4/18/13.
//
//

#include "Recording.h"

Recording g_recording;

void Recording::writeSurfTrack(std::ostream & os, const ElTopo::SurfTrack & st)
{
    size_t n;
    n = st.m_mesh.nv();
    os.write((char *)&n, sizeof (size_t));
    for (size_t i = 0; i < n; i++)
    {
        ElTopo::Vec3d x = st.get_position(i);
        os.write((char *)&(x[0]), sizeof (x[0]));
        os.write((char *)&(x[1]), sizeof (x[1]));
        os.write((char *)&(x[2]), sizeof (x[2]));
        
        //    x = st.get_newposition(i);
        //    os.write((char *)&(x[0]), sizeof (x[0]));
        //    os.write((char *)&(x[1]), sizeof (x[1]));
        //    os.write((char *)&(x[2]), sizeof (x[2]));
        //
        //    x = st.get_remesh_velocity(i);
        //    os.write((char *)&(x[0]), sizeof (x[0]));
        //    os.write((char *)&(x[1]), sizeof (x[1]));
        //    os.write((char *)&(x[2]), sizeof (x[2]));
        //
        //    double mass = st.m_masses[i];
        //    os.write((char *)&mass, sizeof (mass));
        //
        //    bool cl = st.m_mesh.m_vertex_constraint_labels[i];
        //    os.write((char *)&cl, sizeof (cl));
    }
    
    n = st.m_mesh.nt();
    os.write((char *)&n, sizeof (size_t));
    for (size_t i = 0; i < n; i++)
    {
        ElTopo::Vec3st t = st.m_mesh.get_triangle(i);
        os.write((char *)&(t[0]), sizeof (t[0]));
        os.write((char *)&(t[1]), sizeof (t[1]));
        os.write((char *)&(t[2]), sizeof (t[2]));
        
        ElTopo::Vec2i l = st.m_mesh.get_triangle_label(i);
        os.write((char *)&(l[0]), sizeof (l[0]));
        os.write((char *)&(l[1]), sizeof (l[1]));
    }
    
}

void Recording::readSurfTrack(std::istream & is, ElTopo::SurfTrack & st)
{
    // clear the mesh
    for (size_t i = 0; i < st.m_mesh.nt(); i++)
    {
        if (st.m_mesh.get_triangle(i)[0] == st.m_mesh.get_triangle(i)[1])
            continue;
        st.remove_triangle(i);
    }
    
    for (size_t i = 0; i < st.m_mesh.nv(); i++)
        st.remove_vertex(i);
    
    size_t n;
    n = st.m_mesh.nv();
    is.read((char *)&n, sizeof (size_t));
    //  std::cout << " nv = " << n << std::endl;
    st.m_mesh.set_num_vertices(n);
    std::vector<ElTopo::Vec3d> pos(n);
    for (size_t i = 0; i < n; i++)
    {
        ElTopo::Vec3d x;
        is.read((char *)&(x[0]), sizeof (x[0]));
        is.read((char *)&(x[1]), sizeof (x[1]));
        is.read((char *)&(x[2]), sizeof (x[2]));
        pos[i] = x;
        
        //    is.read((char *)&(x[0]), sizeof (x[0]));
        //    is.read((char *)&(x[1]), sizeof (x[1]));
        //    is.read((char *)&(x[2]), sizeof (x[2]));
        //    st.set_newposition(i, x);
        //
        //    is.read((char *)&(x[0]), sizeof (x[0]));
        //    is.read((char *)&(x[1]), sizeof (x[1]));
        //    is.read((char *)&(x[2]), sizeof (x[2]));
        //    st.set_remesh_velocity(i, x);
        //
        //    double mass;
        //    is.read((char *)&mass, sizeof (mass));
        //    st.m_masses[i] = mass;
        //
        //    bool cl;
        //    is.read((char *)&cl, sizeof (cl));
        //    st.m_mesh.m_vertex_constraint_labels[i] = cl;
    }
    
    st.m_masses.resize(n);
    for (size_t i = 0; i < n; i++)
        st.m_masses[i] = 1;
    
    st.set_all_positions(pos);
    st.set_all_newpositions(pos);
    st.set_all_remesh_velocities(std::vector<ElTopo::Vec3d>(n, ElTopo::Vec3d(0)));
    
    //  std::cout << "nv: " << pos.size() << " " << st.pm_velocities.size() << " " << st.m_mesh.m_vertex_constraint_labels.size() << " " << st.m_mesh.m_vertex_to_triangle_map.size() << std::endl;
    //  std::cout << "nv: " << st.m_mesh.nv() << std::endl;
    
    n = st.m_mesh.nt();
    is.read((char *)&n, sizeof (size_t));
    //  std::cout << "nt = " << n << std::endl;
    std::vector<ElTopo::Vec3st> tris;
    std::vector<ElTopo::Vec2i> labels;
    for (size_t i = 0; i < n; i++)
    {
        ElTopo::Vec3st t;
        is.read((char *)&(t[0]), sizeof (t[0]));
        is.read((char *)&(t[1]), sizeof (t[1]));
        is.read((char *)&(t[2]), sizeof (t[2]));
        tris.push_back(t);
        
        ElTopo::Vec2i l;
        is.read((char *)&(l[0]), sizeof (l[0]));
        is.read((char *)&(l[1]), sizeof (l[1]));
        labels.push_back(l);
    }
    
    st.m_mesh.replace_all_triangles(tris, labels);
    //  for (size_t i = 0; i < n; i++)
    //    st.m_mesh.set_triangle_label(i, labels[i]);
    
    size_t nv = st.m_mesh.m_vertex_to_triangle_map.size();
    st.pm_positions.resize(nv);
    st.pm_newpositions.resize(nv);
    st.pm_velocities.resize(nv);
    st.m_velocities.resize(nv);
    
    //  std::cout << "nv: " << pos.size() << " " << st.pm_velocities.size() << " " << st.m_mesh.m_vertex_constraint_labels.size() << " " << st.m_mesh.m_vertex_to_triangle_map.size() << std::endl;
    //  std::cout << "nv: " << st.m_mesh.nv() << std::endl;
    
}

void Recording::recordSurfTrack(const ElTopo::SurfTrack & st)
{
    if (!isRecording())
        return;
    
    if (!m_of.is_open())
    {
        std::stringstream filename;
        filename << m_recording_name << "_" << m_current_frame << ".rec";
        m_of.open(filename.str().c_str());
        
        if (!m_of.is_open())
        {
            std::cout << "Cannot open recording frame file " << filename.str() << std::endl;
            return;
        }
    }
    
    m_of.write((char *)&m_current_step, sizeof (m_current_step));
    
    size_t loglen = m_log.str().size();
    m_of.write((char *)&loglen, sizeof (size_t));
    m_of.write((char *)m_log.str().c_str(), loglen);
    m_log.str("");
    
    writeSurfTrack(m_of, st);
    
    m_current_step++;
}

void Recording::loadRecording(ElTopo::SurfTrack & st, int next)
{
    if (!isPlaybackOn())
        return;
    
    if (!m_if.is_open())
    {
        std::stringstream filename;
        filename << m_recording_name << "_" << m_current_frame << ".rec";
        m_if.open(filename.str().c_str());
        
        if (!m_if.is_open())
        {
            std::cout << "Requested recording frame not found!" << std::endl;
            return;
        }
        
        m_step_pos.clear();
        m_step_pos.push_back(m_if.tellg());
        while (!m_if.eof())
        {
            int step = 0;
            m_if.read((char *)&step, sizeof (step));
            /////////
            //      readSurfTrack(m_if, st);
            /////////
            size_t n;
            m_if.read((char *)&n, sizeof (n));  // skip log
            m_if.seekg(n, std::ios_base::cur);
            m_if.read((char *)&n, sizeof (n));  // skip vertices
            m_if.seekg(sizeof (double) * 3 * n, std::ios_base::cur);
            m_if.read((char *)&n, sizeof (n));  // skip faces
            m_if.seekg((sizeof (size_t) * 3 + sizeof (int) * 2) * n, std::ios_base::cur);
            /////////
            m_if.peek();
            if (!m_if.eof())
                m_step_pos.push_back(m_if.tellg());
        }
        std::cout << "Recording file " << filename.str() << " contains " << m_step_pos.size() << " steps." << std::endl;
        
        m_if.seekg(m_step_pos.front());
        m_if.clear();
    }
    
    assert(m_current_step < m_step_pos.size());
    //  std::cout << "Loading step " << (m_current_step + next) % m_step_pos.size() << std::endl;
    m_if.seekg(m_step_pos[(m_current_step + m_step_pos.size() + next) % m_step_pos.size()]);
    m_if.read((char *)&m_current_step, sizeof (m_current_step));
    
    size_t loglen;
    m_if.read((char *)&loglen, sizeof (size_t));
    char * log = new char[loglen + 1];
    m_if.read((char *)log, loglen);
    *(log + loglen) = NULL;
    
    readSurfTrack(m_if, st);
    
    m_if.peek();
    if (m_if.eof())
        m_if.close();
    
    std::cout << "Loaded recording: step " << m_current_step << "/" << m_step_pos.size() << " of frame " << m_current_frame << std::endl;
    std::cout << log << std::endl;
}

