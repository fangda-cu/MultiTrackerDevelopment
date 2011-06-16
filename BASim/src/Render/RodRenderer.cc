
/**
 * \file RodRenderer.cc
 *
 * \author miklos@cs.columbia.edu
 * \date 08/30/2009
 */

#include "RodRenderer.hh"
#include "RenderUtils.hh"

namespace BASim
{

RodRenderer::RodRenderer(ElasticRod& rod) :
    m_rod(rod), m_tube(m_rod), m_mode(SMOOTH), m_drawMaterial(false), m_drawReference(false), m_scaleToRadius(true),
            m_drawArrows(false), m_drawVelocity(false), m_drawResponse(false)
{
    // material colors
    m_palette.push_back(Color(255, 125, 75));
    m_palette.push_back(Color(55, 125, 255));

    // reference frame colors
    m_palette.push_back(Color(200, 200, 200));
    m_palette.push_back(Color(50, 50, 50));

    // curvature binormal color
    m_palette.push_back(Color(102, 204, 51));

    // velocity color
    m_palette.push_back(Color(255,0,0));

    //root and tip color for simple mode
    m_simpleRod.push_back(Color(0, 0, 0));
    m_simpleRod.push_back(Color(155, 200, 100));
}

void RodRenderer::render()
{
    if (m_mode == SMOOTH)
        drawSmoothRod();
    else if (m_mode == SIMPLE)
        drawSimpleRod();

    if (m_drawMaterial)
        drawMaterialFrame();
    if (m_drawReference)
        drawReferenceFrame();

    if (m_drawVelocity)
        drawVelocityVector();
    if (m_drawResponse)
        drawResponseVector();

}

void RodRenderer::drawSimpleRod()
{
    if (m_rod.draw_cl == 0)
    {
        glLineWidth(2);
        glBegin( GL_LINES);

        //const Color& edgeColor = m_palette[1];
        //const Color& fixedColor = m_palette[0];

        ElasticRod::edge_iter eit;
        for (eit = m_rod.edges_begin(); eit != m_rod.edges_end(); ++eit)
        {
            //if (m_rod.edgeFixed(*eit)) OpenGL::color(fixedColor);
            //else OpenGL::color(edgeColor);
            //OpenGL::color(edgeColor);
            ElasticRod::EdgeVertexIter evit = m_rod.ev_iter(*eit);
            for (evit = m_rod.ev_iter(*eit); evit; ++evit)
            {
                Vec3d x = m_rod.getVertex(*evit);
                OpenGL::vertex(x);
            }
        }

        glEnd();
    }
    else
    {
/*
        std::vector<Color> m_paletteq;
        m_paletteq.push_back(Color(200, 200, 100));
        m_paletteq.push_back(Color(155, 200, 100));
        m_paletteq.push_back(Color(255, 000, 000));
        m_paletteq.push_back(Color(105, 100, 100));
        m_paletteq.push_back(Color(255, 200, 100));
        m_paletteq.push_back(Color(200, 000, 000));
        m_paletteq.push_back(Color(255, 100, 000));
        m_paletteq.push_back(Color(200, 100, 000));
        m_paletteq.push_back(Color(155, 000, 000));
        m_paletteq.push_back(Color(105, 000, 000));
        m_paletteq.push_back(Color(155, 100, 000));
        m_paletteq.push_back(Color(155, 000, 100));
        m_paletteq.push_back(Color(105, 000, 100));
        m_paletteq.push_back(Color(105, 100, 000));
        m_paletteq.push_back(Color(255, 200, 000));
        m_paletteq.push_back(Color(200, 200, 000));
        m_paletteq.push_back(Color(155, 200, 000));
        m_paletteq.push_back(Color(105, 200, 000));
        m_paletteq.push_back(Color(255, 000, 100));
        m_paletteq.push_back(Color(200, 000, 100));
        m_paletteq.push_back(Color(255, 100, 100));
        m_paletteq.push_back(Color(200, 100, 100));
        m_paletteq.push_back(Color(155, 100, 100));
        m_paletteq.push_back(Color(105, 200, 100));
*/
        glEnable( GL_COLOR_MATERIAL);
        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

        //OpenGL::color(m_paletteq[(m_rod.draw_cl % (int)m_paletteq.size())]);
        const float lineWidth = m_rod.getRadiusScale();
        glLineWidth(lineWidth);

        glBegin( GL_LINE_STRIP);
        //const Color& edgeColor = m_palette[1];
        //const Color& fixedColor = m_palette[0];

        int numEdges = m_rod.ne();
        int i = 0;

        ElasticRod::edge_iter eit;
        for (eit = m_rod.edges_begin(); eit != m_rod.edges_end(); ++eit)
        {
            ElasticRod::EdgeVertexIter evit = m_rod.ev_iter(*eit);

            // It turns out the old iteration code had redudent display for
            // each rod vertex. It might slow down the display when
            // thousands of rods are involved.
            do
            {
                Vec3d x = m_rod.getVertex(*evit);
                float t = float(i) / float(numEdges - 1);

                Color blended((1 - t) * m_simpleRod[0].data()[0] + t * m_simpleRod[1].data()[0],
                        (1 - t) * m_simpleRod[0].data()[1] + t * m_simpleRod[1].data()[1],
                        (1 - t) * m_simpleRod[0].data()[2] + t * m_simpleRod[1].data()[2]);

                OpenGL::color(blended);
                OpenGL::vertex(x);
                i++;
                ++evit;

            } while (i == numEdges && evit); //only the last edge need to iterate the scond vertex
        }

        glEnd();

        /*
         glPointSize(5);
         glBegin(GL_POINTS);

         ElasticRod::vertex_iter vit;
         for (vit = m_rod.vertices_begin(); vit != m_rod.vertices_end(); ++vit) {
         const Vec3d& x = m_rod.getVertex(*vit);
         OpenGL::vertex(x);
         }

         glEnd();
         */
        glDisable(GL_COLOR_MATERIAL);
        //    glDisable(GL_LIGHTING);


    }

    /* glPointSize(5);
     glBegin(GL_POINTS);

     ElasticRod::vertex_iter vit;
     for (vit = m_rod.vertices_begin(); vit != m_rod.vertices_end(); ++vit) {
     //if (m_rod.vertFixed(*vit)) OpenGL::color(fixedColor);
     //else OpenGL::color(edgeColor);
     OpenGL::color(edgeColor);
     const Vec3d& x = m_rod.getVertex(*vit);
     OpenGL::vertex(x);
     }

     glEnd();*/
}

void RodRenderer::drawSmoothRod()
{
    m_tube.buildTube();

#ifndef WETA
    glEnable( GL_LIGHTING);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glEnable( GL_COLOR_MATERIAL);
#endif

    const Color& color1 = m_palette[0];
    const Color& color2 = m_palette[1];

    int slices = m_tube.getSlices();
    for (int j = 0; j < m_rod.ne(); ++j)
    {
        const Vec3d& v0 = m_rod.getVertex(j);
        const Vec3d& v1 = m_rod.getVertex((j + 1) % m_rod.nv());
        glBegin( GL_QUAD_STRIP);

        for (int k = 0; k <= slices; ++k)
        {
            if ((0 <= k && k < (slices + 3) / 4) || (slices / 2 <= k && k < (3 * slices + 3) / 4))
            {
                OpenGL::color(color1);
            }
            else
            {
                OpenGL::color(color2);
            }

            const Vec3d& x0 = m_tube.getPointAtVert(j, k);
            //Vec3d n = (v0 - x0).normalized();
            Vec3d n = (x0 - v0).normalized();
#ifdef WETA
            n *= -1.0;
#endif
            OpenGL::normal(n);
            OpenGL::vertex(x0);

            const Vec3d& x1 = m_tube.getPointAtVert(j + 1, k);
            //n = (v1 - x1).normalized();
            n = (x1 - v1).normalized();
#ifdef WETA
            n *= -1.0;
#endif
            OpenGL::normal(n);
            OpenGL::vertex(x1);
        }

        glEnd();
    }

#ifndef WETA
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_LIGHTING);
#endif

}

void RodRenderer::drawSmoothPartialRod(const int i_startVertex, const int i_endVertex, const Vec3d i_color)
{
    m_tube.buildTube();

    /*glEnable(GL_LIGHTING);
     glEnable(GL_COLOR_MATERIAL);
     glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);*/

    const Color& color1 = m_palette[0];
    const Color& color2 = m_palette[1];

    int slices = m_tube.getSlices();
    for (int j = i_startVertex; j < i_endVertex; ++j)
    {
        const Vec3d& v0 = m_rod.getVertex(j);
        const Vec3d& v1 = m_rod.getVertex((j + 1) % m_rod.nv());
        glBegin( GL_QUAD_STRIP);

        for (int k = 0; k <= slices; ++k)
        {
            OpenGL::color(i_color[0], i_color[1], i_color[2]);

            const Vec3d& x0 = m_tube.getPointAtVert(j, k);
            //Vec3d n = (v0 - x0).normalized();
            Vec3d n = (x0 - v0).normalized();
#ifdef WETA
            n *= -1.0;
#endif
            OpenGL::normal(n);
            OpenGL::vertex(x0);

            const Vec3d& x1 = m_tube.getPointAtVert(j + 1, k);
            //n = (v1 - x1).normalized();
            n = (x1 - v1).normalized();
#ifdef WETA
            n *= -1.0;
#endif
            OpenGL::normal(n);
            OpenGL::vertex(x1);
        }

        glEnd();
    }

    /*glDisable(GL_COLOR_MATERIAL);
     glDisable(GL_LIGHTING);*/
}

void RodRenderer::drawMaterialFrame()
{
    glLineWidth(2);
    glBegin( GL_LINES);

    const Color& color1 = m_palette[0];
    const Color& color2 = m_palette[1];
    Scalar r = 1.0;

    ElasticRod::edge_iter eit, end = m_rod.edges_end();
    for (eit = m_rod.edges_begin(); eit != end; ++eit)
    {
        ElasticRod::edge_handle& eh = *eit;
        ElasticRod::vertex_handle vh0 = m_rod.fromVertex(eh);
        ElasticRod::vertex_handle vh1 = m_rod.toVertex(eh);

        Vec3d x = (m_rod.getVertex(vh0) + m_rod.getVertex(vh1)) / 2.0;
        //    if (m_scaleToRadius) r = m_rod.radiusA(eh);

        OpenGL::color(color1);
        if (m_drawArrows)
            drawArrow(x, m_rod.getMaterial1(eh), r);
        else
        {
            Vec3d y = x + r * m_rod.getMaterial1(eh);
            OpenGL::vertex(x);
            OpenGL::vertex(y);
        }

        //    if (m_scaleToRadius) r = m_rod.radiusB(eh);

        OpenGL::color(color2);
        if (m_drawArrows)
            drawArrow(x, m_rod.getMaterial2(eh), r);
        else
        {
            Vec3d y = x + r * m_rod.getMaterial2(eh);
            OpenGL::vertex(x);
            OpenGL::vertex(y);
        }
    }

    glEnd();
}

void RodRenderer::drawReferenceFrame()
{
    glLineWidth(2);
    glBegin( GL_LINES);

    const Color& color1 = m_palette[2];
    const Color& color2 = m_palette[3];

    ElasticRod::edge_iter eit, end = m_rod.edges_end();
    for (eit = m_rod.edges_begin(); eit != end; ++eit)
    {
        ElasticRod::edge_handle& eh = *eit;
        ElasticRod::vertex_handle vh0 = m_rod.fromVertex(eh);
        ElasticRod::vertex_handle vh1 = m_rod.toVertex(eh);

        Vec3d x = (m_rod.getVertex(vh0) + m_rod.getVertex(vh1)) / 2.0;

        OpenGL::color(color1);
        if (m_drawArrows)
            drawArrow(x, m_rod.getReferenceDirector1(eh));
        else
        {
            Vec3d y = x + m_rod.getReferenceDirector1(eh);
            OpenGL::vertex(x);
            OpenGL::vertex(y);
        }

        OpenGL::color(color2);
        if (m_drawArrows)
            drawArrow(x, m_rod.getReferenceDirector2(eh));
        else
        {
            Vec3d y = x + m_rod.getReferenceDirector2(eh);
            OpenGL::vertex(x);
            OpenGL::vertex(y);
        }
    }

    glEnd();
    glEnable( GL_LIGHTING);
}

void RodRenderer::drawVelocityVector()
{
    glLineWidth(2);
    glBegin( GL_LINES);

    const Color& color = m_palette[5];

    ElasticRod::edge_iter eit, end = m_rod.edges_end();
    for (eit = m_rod.edges_begin(); eit != end; ++eit)
    {
        ElasticRod::edge_handle& eh = *eit;
        ElasticRod::vertex_handle vh0 = m_rod.fromVertex(eh);

        drawArrow(m_rod.getVertex(vh0), m_rod.getVelocity(vh0), 0.15, false);
    }

    glEnd();
    //glEnable( GL_LIGHTING);
}

void RodRenderer::drawResponseVector()
{
    glLineWidth(2);
    glBegin( GL_LINES);

    const Color& color = m_palette[5];

    ElasticRod::edge_iter eit, end = m_rod.edges_end();
    for (eit = m_rod.edges_begin(); eit != end; ++eit)
    {
        ElasticRod::edge_handle& eh = *eit;
        ElasticRod::vertex_handle vh0 = m_rod.fromVertex(eh);

        drawArrow(m_rod.getVertex(vh0), m_rod.getResponse(vh0), 0.15, false);
    }

    glEnd();
    glEnable( GL_LIGHTING);
}

Vec3d RodRenderer::calculateObjectCenter()
{
    Vec3d center = Vec3d::Zero();

    ElasticRod::vertex_iter vit, end = m_rod.vertices_end();
    for (vit = m_rod.vertices_begin(); vit != end; ++vit)
    {
        center += m_rod.getVertex(*vit);
    }

    center /= m_rod.nv();

    return center;
}

Scalar RodRenderer::calculateObjectBoundingRadius(const Vec3d& center)
{
    Scalar radius = 0.0;

    ElasticRod::vertex_iter vit, end = m_rod.vertices_end();
    for (vit = m_rod.vertices_begin(); vit != end; ++vit)
    {
        radius = std::max(radius, (m_rod.getVertex(*vit) - center).norm());
    }

    return radius;
}

} // namespace BASim
