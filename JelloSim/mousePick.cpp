#include "mousePick.h"
#include "jelloMesh.h"
#include <GL/glut.h>
#include <algorithm>

void mouseRay(double mouseX, double mouseY)
{
	double matModelView[16], matProjection[16]; 
	int viewport[4]; 

	glGetDoublev( GL_MODELVIEW_MATRIX, matModelView ); 
	glGetDoublev( GL_PROJECTION_MATRIX, matProjection ); 
	glGetIntegerv( GL_VIEWPORT, viewport ); 

	double winX = (double)mouseX; 
	double winY = viewport[3] - (double)mouseY; 

	vec3 m_start = vec3(0,0,0);
	vec3 m_end = vec3(0,0,0);

	gluUnProject(winX, winY, 0.0, matModelView, matProjection, viewport, &m_start.n[0], &m_start.n[1], &m_start.n[2]); 
	gluUnProject(winX, winY, 1.0, matModelView, matProjection, viewport, &m_end.n[0], &m_end.n[1], &m_end.n[2]);
}

void testHit(vec3 m_start,vec3 m_end)
{
	//isInterior(int i, int j, int k)
	
	for (int i = 0; i < m_rows+1; i++)
    {
        for (int j = 0; j < m_cols+1; j++)
        {
            for (int k = 0; k < m_stacks+1; k++)
            {
                Particle& p = GetParticle(grid, i,j,k);
				p_vertipoint = p.
                p.force = m_externalForces * p.mass;
            }
        }
    }
}