#include "JelloMesh.h"
#include <GL/glut.h>
#include <algorithm>

// Initialize parameters
double JelloMesh::g_structuralKs = 0.0; 
double JelloMesh::g_structuralKd = 0.0; 
double JelloMesh::g_shearKs = 0.0;
double JelloMesh::g_shearKd = 0.0;
double JelloMesh::g_bendKs = 0.0;
double JelloMesh::g_bendKd = 0.0;
double JelloMesh::g_penaltyKs = 0.0;
double JelloMesh::g_penaltyKd = 0.0;

JelloMesh::JelloMesh() :     
    m_integrationType(JelloMesh::RK4), m_drawflags(MESH | STRUCTURAL),
    m_cols(0), m_rows(0), m_stacks(0), m_width(0.0), m_height(0.0), m_depth(0.0)
{
    SetSize(1.0, 1.0, 1.0);
    SetGridSize(6, 6, 6);
}

JelloMesh::~JelloMesh()
{
}

void JelloMesh::Reset()
{
    InitJelloMesh();
}

JelloMesh::Particle& JelloMesh::GetParticle(JelloMesh::ParticleGrid& grid, int i, int j, int k)
{
    return grid[i][j][k];
}

JelloMesh::Particle& JelloMesh::GetParticle(JelloMesh::ParticleGrid& grid, int idx)
{
    int i,j,k;
    GetCell(idx, i, j, k);
    return GetParticle(grid, i,j,k);
}

const JelloMesh::Particle& JelloMesh::GetParticle(const JelloMesh::ParticleGrid& grid, int i, int j, int k) const
{
    return grid[i][j][k];
}

const JelloMesh::Particle& JelloMesh::GetParticle(const JelloMesh::ParticleGrid& grid, int idx) const
{
    int i,j,k;
    GetCell(idx, i, j, k);
    return GetParticle(grid, i,j,k);
}

bool JelloMesh::isInterior(const JelloMesh::Spring& s) const
{
    int i1,j1,k1,i2,j2,k2;
    GetCell(s.m_p1, i1, j1, k1);
    GetCell(s.m_p2, i2, j2, k2);
    return isInterior(i1,j1,k1) || isInterior(i2,j2,k2);
}


bool JelloMesh::isInterior(int idx) const
{
    int i,j,k;
    GetCell(idx, i, j, k);
    return isInterior(i,j,k);
}

bool JelloMesh::isInterior(int i, int j, int k) const
{
    return (i*j*k*(m_rows-i)*(m_cols-j)*(m_stacks-k) != 0);
}

void JelloMesh::SetGridSize(int cols, int rows, int stacks)
{
    m_cols = cols;
    m_rows = rows;
    m_stacks = stacks;

    if (m_cols > 0 && m_rows > 0 && m_stacks > 0)
    {
        m_vparticles.resize(m_rows+1);
        for (int i = 0; i < m_rows+1; i++)
        {
            m_vparticles[i].resize(m_cols+1);
            for (int j = 0; j < m_cols+1; j++)
            {
                m_vparticles[i][j].resize(m_stacks+1);
            }
        }
    }
    InitJelloMesh();
}

int JelloMesh::GetGridCols() const
{
    return m_cols;
}

int JelloMesh::GetGridRows() const
{
    return m_rows;
}

int JelloMesh::GetGridStacks() const
{
    return m_stacks;
}

void JelloMesh::SetSize(float width, float height, float depth)
{
    m_width = width;
    m_height = height;
    m_depth = depth;
    InitJelloMesh();
}

float JelloMesh::GetWidth() const
{
    return m_width;
}

float JelloMesh::GetHeight() const
{
    return m_height;
}

float JelloMesh::GetDepth() const
{
    return m_depth;
}

int JelloMesh::GetIndex(int i, int j, int k) const
{
    int cols = j;
    int rows = i*(m_cols+1);
    int stacks = k*(m_cols+1)*(m_rows+1);
    return cols + rows + stacks;
}

#define ROUND(x) (floor(x + 0.5))
#define FLOOR(x) (floor(x))
#define FRACT(x) (x - FLOOR(x))
void JelloMesh::GetCell(int idx, int& i, int &j, int& k) const
{
    float rows = m_rows+1;
    float cols = m_cols+1;
    float stacks = m_stacks+1;

    // derived from idx = cols*(rows*k + i) + j
    float tmp = FLOOR(idx/cols);
    j = (int) ROUND(cols*(FRACT(idx/cols)));
    i = (int) ROUND(rows*(FRACT(tmp/rows)));
    k = (int) FLOOR(tmp/rows);
}

void JelloMesh::InitJelloMesh()
{
    m_vsprings.clear();

    if (m_width < 0.01 || m_height < 0.01 || m_depth < 0.01) return;
    if (m_cols < 1 || m_rows < 1 || m_stacks < 1) return;

    // Init particles
    float wcellsize = m_width / m_cols;
    float hcellsize = m_height / m_rows;
    float dcellsize = m_depth / m_stacks;
    
    for (int i = 0; i < m_rows+1; i++)
    {
        for (int j = 0; j < m_cols+1; j++)
        {
            for (int k = 0; k < m_stacks+1; k++)
            {
                float x = -m_width*0.5f + wcellsize*i;
                float y = 0.5 + hcellsize*j; 
                float z = -m_depth*0.5f + dcellsize*k;
                m_vparticles[i][j][k] = Particle(GetIndex(i,j,k), vec3(x, y, z));
            }
        }
    }

    // Setup structural springs
    ParticleGrid& g = m_vparticles;
    for (int i = 0; i < m_rows+1; i++)
    {
        for (int j = 0; j < m_cols+1; j++)
        {
            for (int k = 0; k < m_stacks+1; k++)
            {
                if (j < m_cols) AddStructuralSpring(GetParticle(g,i,j,k), GetParticle(g,i,j+1,k));
                if (i < m_rows) AddStructuralSpring(GetParticle(g,i,j,k), GetParticle(g,i+1,j,k));
                if (k < m_stacks) AddStructuralSpring(GetParticle(g,i,j,k), GetParticle(g,i,j,k+1));
            }
        }
    }

	//TODO
	for (int i = 0; i < m_rows+1; i++)
    {
        for (int j = 0; j < m_cols+1; j++)
        {
            for (int k = 0; k < m_stacks+1; k++)
            {
				//(j < m_cols)   (i < m_rows) (k < m_stacks) 
				if ((j < m_cols)&&(k < m_stacks))
				{
					AddShearSpring(GetParticle(g,i,j,k), GetParticle(g,i,j+1,k+1));
					AddShearSpring(GetParticle(g,i,j,k+1), GetParticle(g,i,j+1,k));
				}
				if((i < m_rows)&&(j < m_cols))
				{
					AddShearSpring(GetParticle(g,i+1,j,k), GetParticle(g,i,j+1,k));
					AddShearSpring(GetParticle(g,i+1,j+1,k), GetParticle(g,i,j,k));
				}
				if((i < m_rows)&&(k < m_stacks))
				{
					AddShearSpring(GetParticle(g,i+1,j,k), GetParticle(g,i,j,k+1));
					AddShearSpring(GetParticle(g,i,j,k), GetParticle(g,i+1,j,k+1));
				}
            }
        }
    }
	
	
	for (int i = 0; i < m_rows+1; i++)
    {
        for (int j = 0; j < m_cols+1; j++)
        {
            for (int k = 0; k < m_stacks+1; k++)
            {
                if (j < m_cols-1) AddBendSpring(GetParticle(g,i,j,k), GetParticle(g,i,j+2,k));
                if (i < m_rows-1) AddBendSpring(GetParticle(g,i,j,k), GetParticle(g,i+2,j,k));
                if (k < m_stacks-1) AddBendSpring(GetParticle(g,i,j,k), GetParticle(g,i,j,k+2));
            }
        }
    }
	//TODO
    // Init mesh geometry
    m_mesh.clear();
    m_mesh.push_back(FaceMesh(*this,XLEFT));
    m_mesh.push_back(FaceMesh(*this,XRIGHT));
    m_mesh.push_back(FaceMesh(*this,YTOP));
    m_mesh.push_back(FaceMesh(*this,YBOTTOM));
    m_mesh.push_back(FaceMesh(*this,ZFRONT));
    m_mesh.push_back(FaceMesh(*this,ZBACK));
}

void JelloMesh::AddStructuralSpring(Particle& p1, Particle& p2)
{
    double restLen = (p1.position - p2.position).Length();
    m_vsprings.push_back(Spring(STRUCTURAL, p1.index, p2.index, g_structuralKs, g_structuralKd, restLen));
}

void JelloMesh::AddBendSpring(JelloMesh::Particle& p1, JelloMesh::Particle& p2)
{
    double restLen = (p1.position - p2.position).Length();
    m_vsprings.push_back(Spring(BEND, p1.index, p2.index, g_bendKs, g_bendKd, restLen));
}

void JelloMesh::AddShearSpring(JelloMesh::Particle& p1, JelloMesh::Particle& p2)
{
    double restLen = (p1.position - p2.position).Length();
    m_vsprings.push_back(Spring(SHEAR, p1.index, p2.index, g_shearKs, g_shearKd, restLen));
}

void JelloMesh::SetIntegrationType(JelloMesh::IntegrationType type)
{
    m_integrationType = type;
}

JelloMesh::IntegrationType JelloMesh::GetIntegrationType() const
{
    return m_integrationType;
}

void JelloMesh::SetDrawFlags(unsigned int flags)
{
    m_drawflags = flags;
}

unsigned int JelloMesh::GetDrawFlags() const
{
    return m_drawflags;
}

void JelloMesh::DrawMesh(const vec3& eyePos)
{
    const ParticleGrid& g = m_vparticles;
    float red[4] = {1.0,0.4,0.4,0.8};
    float white[4] = {1.0,1.0,1.0,1.0};
    float pink[4] = {0.5,0.0,0.0,1.0};
    float black[4] = {0.0,0.0,0.0,1.0};
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, red);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, black);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, black);
    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, pink);

    for (unsigned int i = 0; i < m_mesh.size(); i++)
    {        
       m_mesh[i].CalcDistToEye(*this, eyePos);
    }
    std::sort(m_mesh.begin(), m_mesh.end(), FaceMesh::compare);
    for (unsigned int i = 0; i < m_mesh.size(); i++)
    {        
       m_mesh[i].Draw(*this);
    }

    //glDisable(GL_LIGHTING);
    //for (unsigned int i = 0; i < m_mesh.size(); i++)
    //{        
    //   m_mesh[i].DrawNormals(*this);
    //}

}

void JelloMesh::DrawSprings(double a)
{
    const ParticleGrid& g = m_vparticles;
    glBegin(GL_LINES);
    for (unsigned int i = 0; i < m_vsprings.size(); i++)
    {
        if (!(m_vsprings[i].m_type & m_drawflags)) continue;
        if (isInterior(m_vsprings[i])) continue;

        switch (m_vsprings[i].m_type)
        {
        case BEND:       glColor4f(1.0, 1.0, 0.0, a); break;
        case STRUCTURAL: glColor4f(1.0, 1.0, 0.0, a); break;
        case SHEAR:      glColor4f(0.0, 1.0, 1.0, a); break;
        };

        vec3 p1 = GetParticle(g, m_vsprings[i].m_p1).position;
        vec3 p2 = GetParticle(g, m_vsprings[i].m_p2).position;
        glVertex3f(p1[0], p1[1], p1[2]);
        glVertex3f(p2[0], p2[1], p2[2]);
    }
    glEnd();
}

void JelloMesh::DrawCollisionNormals()
{
    const ParticleGrid& g = m_vparticles;
    glBegin(GL_LINES);
    glColor3f(0.0, 1.0, 0.0);
    for(unsigned int i = 0; i < m_vcollisions.size(); i++)
    {
       Intersection intersection = m_vcollisions[i];
       if (isInterior(intersection.m_p)) continue;

       const Particle& pt = GetParticle(g, intersection.m_p);
       vec3 normal = intersection.m_normal;
       vec3 end = pt.position + 0.2 * normal;
       glVertex3f(pt.position[0], pt.position[1], pt.position[2]);
       glVertex3f(end[0], end[1], end[2]);
    }     
    glEnd();
}

void JelloMesh::DrawForces()
{
    glBegin(GL_LINES);
    glColor3f(1.0, 0.0, 0.0);
    for (int i = 0; i < m_rows+1; i++)
    {
        for (int j = 0; j < m_cols+1; j++)
        {
            for (int k = 0; k < m_stacks+1; k++)
            {
                Particle p = m_vparticles[i][j][k];
                if (isInterior(i,j,k)) continue;

                vec3 normal = p.force.Normalize();
                vec3 end = p.position + 0.1 * normal;
                glVertex3f(p.position[0], p.position[1], p.position[2]);
                glVertex3f(end[0], end[1], end[2]);
            }
        }
    }     
    glEnd();
}

void JelloMesh::Draw(const vec3& eyePos)
{
    if (m_drawflags & MESH) DrawMesh(eyePos);

    if (m_drawflags & (STRUCTURAL|BEND|SHEAR))
    {
        glDisable(GL_LIGHTING);
        glDisable(GL_DEPTH_TEST);
        glLineWidth(1.0);
        DrawSprings(0.2);
        glLineWidth(1.5);
        glEnable(GL_DEPTH_TEST);
        DrawSprings(0.4);
    }

    if (m_drawflags & NORMALS) DrawCollisionNormals();
    if (m_drawflags & FORCES) DrawForces();

    glEnable(GL_LIGHTING);
}

void JelloMesh::Update(double dt, const World& world, const vec3& externalForces)
{
    m_externalForces = externalForces;

	CheckForCollisions(m_vparticles, world);
	ComputeForces(m_vparticles);
	ResolveContacts(m_vparticles);
	ResolveCollisions(m_vparticles);

    switch (m_integrationType)
    {
    case EULER: EulerIntegrate(dt); break;
    case MIDPOINT: MidPointIntegrate(dt); break;
    case RK4: RK4Integrate(dt); break;
    }
}

void JelloMesh::CheckForCollisions(ParticleGrid& grid, const World& world)
{
    m_vcontacts.clear();
    m_vcollisions.clear();

    for (int i = 0; i < m_rows+1; i++)
    {
        for (int j = 0; j < m_cols+1; j++)
        {
            for (int k = 0; k < m_stacks+1; k++)
            {
                Particle& p = GetParticle(grid, i,j,k);

                // 1. Check collisions with world objects 
                for (unsigned int m = 0; m < world.m_shapes.size(); m++)
                {
                    Intersection intersection;

                    if (world.m_shapes[m]->GetType() == World::CYLINDER && 
                        CylinderIntersection(p, (World::Cylinder*) world.m_shapes[m], intersection))
                    {
                        //m_vcontacts.push_back(intersection);
						if(intersection.m_type==CONTACT)
						{
							m_vcontacts.push_back(intersection);
						}
						else if(intersection.m_type==COLLISION)
						{
							m_vcollisions.push_back(intersection);
						}
                    }
                    else if (world.m_shapes[m]->GetType() == World::GROUND && 
                        FloorIntersection(p, intersection))
                    {
                        //m_vcontacts.push_back(intersection);
						if(intersection.m_type==CONTACT)
						{
							m_vcontacts.push_back(intersection);
						}
						else if(intersection.m_type==COLLISION)
						{
							m_vcollisions.push_back(intersection);
						}
                    }
					else if (world.m_shapes[m]->GetType() == World::SPHERE && 
                        SphereIntersection(p, (World::Sphere*) world.m_shapes[m], intersection))
                    {
                        //m_vcontacts.push_back(intersection);
						if(intersection.m_type==CONTACT)
						{
							m_vcontacts.push_back(intersection);
						}
						else if(intersection.m_type==COLLISION)
						{
							m_vcollisions.push_back(intersection);
						}
                    }
					else if (world.m_shapes[m]->GetType() == World::CUBE && 
                        CubeIntersection(p, (World::Cube*) world.m_shapes[m], intersection))
                    {
                        //m_vcontacts.push_back(intersection);
						if(intersection.m_type==CONTACT)
						{
							m_vcontacts.push_back(intersection);
						}
						else if(intersection.m_type==COLLISION)
						{
							m_vcollisions.push_back(intersection);
						}
                    }
                }
            }
        }
    }
}

void JelloMesh::ComputeForces(ParticleGrid& grid)
{
    // Add external froces to all points
    for (int i = 0; i < m_rows+1; i++)
    {
        for (int j = 0; j < m_cols+1; j++)
        {
            for (int k = 0; k < m_stacks+1; k++)
            {
                Particle& p = GetParticle(grid, i,j,k);
                p.force = m_externalForces * p.mass;
            }
        }
    }

    // Update springs
    for(unsigned int i = 0; i < m_vsprings.size(); i++)
    {
        Spring& spring = m_vsprings[i];
        Particle& a = GetParticle(grid, spring.m_p1);
        Particle& b = GetParticle(grid, spring.m_p2);

        // TODO
		float orginalLength = JelloMesh::m_width/JelloMesh::m_cols;
		float olshear = 1.414*orginalLength;
		float olbend = 2*orginalLength;

		if(spring.m_type==0x1)
		{
		
		spring.m_restLen = (a.position-b.position).Length();


		spring.m_Ks = g_structuralKs;
		spring.m_Kd = g_structuralKd;
		a.force = a.force+spring.m_Ks*(spring.m_restLen-orginalLength)/(spring.m_restLen)*(b.position-a.position)+
			spring.m_Kd*(b.velocity-a.velocity)*(a.position-b.position)/spring.m_restLen*(a.position-b.position)/spring.m_restLen;
		b.force = b.force-spring.m_Ks*(spring.m_restLen-orginalLength)/(spring.m_restLen)*(b.position-a.position)-
			spring.m_Kd*(b.velocity-a.velocity)*(a.position-b.position)/spring.m_restLen*(a.position-b.position)/spring.m_restLen;
		}

		if(spring.m_type==0x2)
		{

		spring.m_restLen = (a.position-b.position).Length();

		spring.m_Ks = g_bendKs;
		spring.m_Kd = g_bendKd;
		a.force = a.force+spring.m_Ks*(spring.m_restLen-olbend)/(spring.m_restLen)*(b.position-a.position)+
			spring.m_Kd*(b.velocity-a.velocity)*(a.position-b.position)/spring.m_restLen*(a.position-b.position)/spring.m_restLen;
		b.force = b.force-spring.m_Ks*(spring.m_restLen-olbend)/(spring.m_restLen)*(b.position-a.position)-
			spring.m_Kd*(b.velocity-a.velocity)*(a.position-b.position)/spring.m_restLen*(a.position-b.position)/spring.m_restLen;
		}
		

		if(spring.m_type==0x3)
		{

		spring.m_restLen = (a.position-b.position).Length();

		spring.m_Ks = g_shearKs;
		spring.m_Kd = g_shearKd;
		a.force = a.force+spring.m_Ks*(spring.m_restLen-olshear)/(spring.m_restLen)*(b.position-a.position)+
			spring.m_Kd*(b.velocity-a.velocity)*(a.position-b.position)/spring.m_restLen*(a.position-b.position)/spring.m_restLen;
		b.force = b.force-spring.m_Ks*(spring.m_restLen-olshear)/(spring.m_restLen)*(b.position-a.position)-
			spring.m_Kd*(b.velocity-a.velocity)*(a.position-b.position)/spring.m_restLen*(a.position-b.position)/spring.m_restLen;
		}
		

    }
}

void JelloMesh::ResolveContacts(ParticleGrid& grid)
{
    for (unsigned int i = 0; i < m_vcontacts.size(); i++)
    {
       const Intersection& contact = m_vcontacts[i];
       Particle& p = GetParticle(grid, contact.m_p);
       vec3 normal = contact.m_normal; 

        // TODO
	   double distance = m_vcontacts[i].m_distance;
	   p.position = p.position + normal*distance;

	   float vc = 1;
	   if((p.velocity)*(-normal)>0)
	   {
		vec3 velonormal = ((p.velocity)*(-normal))*(-normal);
		vec3 velotangent = p.velocity-velonormal;
		p.velocity = velotangent-vc*velonormal;
	   }
    }
}

void JelloMesh::ResolveCollisions(ParticleGrid& grid)
{
    for(unsigned int i = 0; i < m_vcollisions.size(); i++)
    {
        Intersection result = m_vcollisions[i];
        Particle& pt = GetParticle(grid, result.m_p);
        vec3 normal = result.m_normal;
        float dist = result.m_distance;

        // TODO
		//float kc = 0.3;
		pt.force = pt.force + dist*g_penaltyKs*result.m_normal + g_penaltyKd*(-pt.velocity);
	}
}

bool JelloMesh::FloorIntersection(Particle& p, Intersection& intersection)
{
    // TODO
	double margin = 0.05;
	vec3 nor = vec3(0,1,0);
	World::Ground ground;
	if (ground.pos.n[1]>=p.position.n[1])
	{
		intersection.m_normal.set(0,1,0);
		intersection.m_p = p.index;
		intersection.m_distance = ground.pos.n[1]-p.position.n[1];
		intersection.m_type = CONTACT;
		return true;
	}
	else if (p.position.n[1]<=(ground.pos.n[1] + margin))
	{
		intersection.m_normal.set(0,1,0);
		intersection.m_p = p.index;
		intersection.m_distance = ground.pos.n[1]-p.position.n[1]+margin;
		intersection.m_type = COLLISION;
		return true;
	}
    return false;
}


bool JelloMesh::CylinderIntersection(Particle& p, World::Cylinder* cylinder, 
                                 JelloMesh::Intersection& result)
{
    vec3 cylinderStart = cylinder->start;
    vec3 cylinderEnd = cylinder->end;
    vec3 cylinderAxis = cylinderEnd - cylinderStart;
    double cylinderRadius = cylinder->r; 

    // TODO
	float shadow = (cylinderAxis.Normalize())*(p.position-cylinderStart);
	vec3 vertipoint = (cylinderAxis.Normalize())*shadow+cylinderStart;
	float vertiline = Distance(vertipoint, p.position);
	double margin = 0.1;

	if((shadow>0)&&(shadow<cylinderAxis.Length())&&(vertiline<cylinderRadius))
	{
		result.m_type = CONTACT;
		result.m_p = p.index;
		result.m_normal = (p.position - vertipoint).Normalize();
		result.m_distance = cylinderRadius + margin - vertiline;
		return true;
	}
	
	if((shadow<=0)&&(shadow>=-margin)&&(vertiline<=cylinderRadius))
	{
		result.m_type = COLLISION;
		result.m_p = p.index;
		result.m_normal = -cylinderAxis.Normalize();
		result.m_distance = (p.position - cylinderStart)*result.m_normal;

		return true;
	}
	
	if((shadow>=cylinderAxis.Length())&&(shadow<=cylinderAxis.Length()+margin)&&(vertiline<=cylinderRadius))
	{
		result.m_type = COLLISION;
		result.m_p = p.index;
		result.m_normal = cylinderAxis.Normalize();
		result.m_distance = (p.position - cylinderEnd)*result.m_normal;
		return true;
	}

	if((shadow<cylinderAxis.Length())&&(shadow>=0)&&(vertiline<=cylinderRadius+margin)&&(vertiline>=cylinderRadius))
	{
		result.m_type = COLLISION;
		result.m_p = p.index;
		result.m_normal = cylinderAxis.Normalize();
		result.m_distance = cylinderRadius + margin -vertiline;
		return true;
	}
    return false;
}
bool JelloMesh::SphereIntersection(Particle& p, World::Sphere* sphere, JelloMesh::Intersection& result)
{
    double radius = sphere->r;
    vec3 location = sphere->pos;
	double margin = 0.1;

	if ((p.position-location).Length()<radius)
	{
		result.m_type = CONTACT;
		result.m_p = p.index;
		result.m_normal = (p.position - location).Normalize();
		result.m_distance = radius + margin - (p.position-location).Length();
		return true;
	}

	else if((p.position-location).Length()<radius+margin)
	{
		result.m_type = COLLISION;
		result.m_p = p.index;
		result.m_normal = (p.position - location).Normalize();
		result.m_distance = radius+ margin - (p.position,location).Length();
		return true;
	}
	return false;

}

bool JelloMesh::CubeIntersection(Particle& p, World::Cube* cube, JelloMesh::Intersection& result)
{
    vec3 location = cube->pos;
	double xmin = location.n[0]-cube->hx/2;
	double xmax = location.n[0]+cube->hx/2;
	double ymin = location.n[1]-cube->hy/2;
	double ymax = location.n[1]+cube->hy/2;
	double zmin = location.n[2]-cube->hz/2;
	double zmax = location.n[2]+cube->hz/2;

	vec3 distance;
	distance[0] = (cube->hx/2-(p.position-location)*vec3(1,0,0));
	distance[1] = (cube->hy/2-(p.position-location)*vec3(0,1,0));
	distance[2] = (cube->hz/2-(p.position-location)*vec3(0,0,1));

	double margin = 0.1;
	double min_distance = distance[0];
	int min_index = 1;

	for(int i = 1;i<=3;i++)
	{
		if (distance[i-1]<min_distance)
		{
			min_distance = distance[i-1]+margin;
			min_index = i;
		}
	}

	if ((p.position.n[0]>xmin)&&(p.position.n[0]<xmax)&&
		(p.position.n[1]>ymin)&&(p.position.n[1]<ymax)&&
		(p.position.n[2]>zmin)&&(p.position.n[2]<zmax))
	{
		//cout<<"detect contact"<<endl;
		result.m_distance = min_distance;
		result.m_type = CONTACT;
		result.m_p = p.index;
		if(min_index == 1)
		{
			result.m_normal = vec3(1,0,0);
		}
		else if (min_index == 2)
		{
			result.m_normal = vec3(0,1,0);
		}
		else if (min_index == 3)
		{
			result.m_normal = vec3(0,0,1);
		}
		return true;
	}

	else if((p.position.n[0]>xmin-margin)&&(p.position.n[0]<xmax+margin)&&
			(p.position.n[1]>ymin-margin)&&(p.position.n[1]<ymax+margin)&&
			(p.position.n[2]>zmin-margin)&&(p.position.n[2]<zmax+margin))
	{
		//cout<<"detect collision"<<endl;
		result.m_type = COLLISION;
		result.m_p = p.index;
		result.m_distance = min_distance;

		if(min_index == 1)
		{
			result.m_normal = vec3(1,0,0);
		}
		else if (min_index == 2)
		{
			result.m_normal = vec3(0,1,0);
		}
		else if (min_index == 3)
		{
			result.m_normal = vec3(0,0,1);
		}
		cout<<result.m_normal<<endl;
		return true;
	}
	return false;

}

void JelloMesh::EulerIntegrate(double dt)
{
    // TODO
	ParticleGrid target = m_vparticles; 
	ParticleGrid& source = m_vparticles; 

	ParticleGrid accum1 = m_vparticles;
	for (int i = 0; i < m_rows+1; i++)
	{
		for (int j = 0; j < m_cols+1; j++)
		{
			for (int k = 0; k < m_stacks+1; k++)
			{
				Particle& s = GetParticle(source, i,j,k);

				Particle& k1 = GetParticle(accum1, i,j,k);
				k1.force = s.force;
				k1.velocity = s.velocity;

				Particle& t = GetParticle(target, i,j,k);
				t.velocity = s.velocity + dt * k1.force * 1/k1.mass;
				t.position = s.position + dt * k1.velocity;
			}
		}
	}
	ComputeForces(target);

	for (int i = 0; i < m_rows+1; i++)
	{
		for (int j = 0; j < m_cols+1; j++)
		{
			for (int k = 0; k < m_stacks+1; k++)
			{
				Particle& p = GetParticle(m_vparticles, i,j,k);
				Particle& k1 = GetParticle(accum1, i,j,k);

				p.velocity = p.velocity + dt*k1.force*1/p.mass;

				p.position = p.position + dt*k1.velocity;
			}
		}
	}
}

void JelloMesh::MidPointIntegrate(double dt)
{
    // TODO
	double halfdt = 0.5 * dt;
	ParticleGrid target = m_vparticles; // target is a copy!
	ParticleGrid& source = m_vparticles; // source is a ptr!
	// Step 1
	ParticleGrid accum1 = m_vparticles;
	for (int i = 0; i < m_rows+1; i++)
	{
		for (int j = 0; j < m_cols+1; j++)
		{
			for (int k = 0; k < m_stacks+1; k++)
			{
				Particle& s = GetParticle(source, i,j,k);

				Particle& k1 = GetParticle(accum1, i,j,k);
				k1.force = s.force;
				k1.velocity = s.velocity;

				Particle& t = GetParticle(target, i,j,k);
				t.velocity = s.velocity + halfdt * k1.force * 1/k1.mass;
				t.position = s.position + halfdt * k1.velocity;
			}
		}
	}

	ComputeForces(target);

	ParticleGrid accum2 = m_vparticles;
	for (int i = 0; i < m_rows+1; i++)
	{
		for (int j = 0; j < m_cols+1; j++)
		{
			for (int k = 0; k < m_stacks+1; k++)
			{
				Particle& t = GetParticle(target, i,j,k);
				Particle& k2 = GetParticle(accum2, i,j,k);

				k2.force = t.force;
				k2.velocity = t.velocity;

				Particle& s = GetParticle(source, i,j,k);
				t.velocity = s.velocity + halfdt * k2.force * 1/k2.mass;
				t.position = s.position + halfdt * k2.velocity;
			}
		}
	}
	ComputeForces(target);

	for (int i = 0; i < m_rows+1; i++)
	{
		for (int j = 0; j < m_cols+1; j++)
		{
			for (int k = 0; k < m_stacks+1; k++)
			{
				Particle& s = GetParticle(source, i,j,k);

				Particle& p = GetParticle(m_vparticles, i,j,k);
				Particle& k1 = GetParticle(accum1, i,j,k);
				Particle& k2 = GetParticle(accum2, i,j,k);
				
				p.velocity = s.velocity + dt*k2.force*1/p.mass;

				p.position = s.position + dt*k2.velocity;
			}
		}
	}

}

void JelloMesh::RK4Integrate(double dt)
{
	double halfdt = 0.5 * dt;
	ParticleGrid target = m_vparticles; // target is a copy!
	ParticleGrid& source = m_vparticles; // source is a ptr!

	// Step 1
	ParticleGrid accum1 = m_vparticles;
	for (int i = 0; i < m_rows+1; i++)
	{
		for (int j = 0; j < m_cols+1; j++)
		{
			for (int k = 0; k < m_stacks+1; k++)
			{
				Particle& s = GetParticle(source, i,j,k);

				Particle& k1 = GetParticle(accum1, i,j,k);
				k1.force = s.force;
				k1.velocity = s.velocity;

				Particle& t = GetParticle(target, i,j,k);
				t.velocity = s.velocity + halfdt * k1.force * 1/k1.mass;
				t.position = s.position + halfdt * k1.velocity;
			}
		}
	}

	ComputeForces(target);

	// Step 2
	ParticleGrid accum2 = m_vparticles;
	for (int i = 0; i < m_rows+1; i++)
	{
		for (int j = 0; j < m_cols+1; j++)
		{
			for (int k = 0; k < m_stacks+1; k++)
			{
				Particle& t = GetParticle(target, i,j,k);
				Particle& k2 = GetParticle(accum2, i,j,k);

				k2.force = t.force;
				k2.velocity = t.velocity;

				Particle& s = GetParticle(source, i,j,k);
				t.velocity = s.velocity + halfdt * k2.force * 1/k2.mass;
				t.position = s.position + halfdt * k2.velocity;
			}
		}
	}

	ComputeForces(target);

	// Step 3
	ParticleGrid accum3 = m_vparticles;
	for (int i = 0; i < m_rows+1; i++)
	{
		for (int j = 0; j < m_cols+1; j++)
		{
			for (int k = 0; k < m_stacks+1; k++)
			{
				Particle& t = GetParticle(target, i,j,k);
				Particle& k3 = GetParticle(accum3, i,j,k);

				k3.force = t.force;
				k3.velocity = t.velocity;

				Particle& s = GetParticle(source, i,j,k);
				t.velocity = s.velocity + dt * k3.force * 1/k3.mass;
				t.position = s.position + dt * k3.velocity;
			}
		}
	}
	ComputeForces(target);

	// Step 4
	ParticleGrid accum4 = m_vparticles;
	for (int i = 0; i < m_rows+1; i++)
	{
		for (int j = 0; j < m_cols+1; j++)
		{
			for (int k = 0; k < m_stacks+1; k++)
			{
				Particle& t = GetParticle(target, i,j,k);
				Particle& k4 = GetParticle(accum4, i,j,k);

				k4.force = t.force;
				k4.velocity = t.velocity;
			}
		}
	}

	// Put it all together
	double asixth = 1/6.0;
	double athird = 1/3.0;
	for (int i = 0; i < m_rows+1; i++)
	{
		for (int j = 0; j < m_cols+1; j++)
		{
			for (int k = 0; k < m_stacks+1; k++)
			{
				Particle& p = GetParticle(m_vparticles, i,j,k);
				Particle& k1 = GetParticle(accum1, i,j,k);
				Particle& k2 = GetParticle(accum2, i,j,k);
				Particle& k3 = GetParticle(accum3, i,j,k);
				Particle& k4 = GetParticle(accum4, i,j,k);

				p.velocity = p.velocity + dt*(asixth * k1.force +
					athird * k2.force + athird * k3.force + asixth * k4.force)*1/p.mass;

				p.position = p.position + dt*(asixth * k1.velocity +
				athird * k2.velocity + athird * k3.velocity + asixth * k4.velocity);
			}
		}
	}
}


//---------------------------------------------------------------------
// Spring
//---------------------------------------------------------------------
JelloMesh::Spring::Spring() :
    m_type(JelloMesh::STRUCTURAL), 
    m_p1(-1), 
    m_p2(-1), 
    m_Ks(1.0), m_Kd(1.0), m_restLen(1.0)
{
}

JelloMesh::Spring::Spring(const JelloMesh::Spring& p) :
    m_type(p.m_type), m_p1(p.m_p1), m_p2(p.m_p2),
    m_Ks(p.m_Ks), m_Kd(p.m_Kd), m_restLen(p.m_restLen)
{
}

JelloMesh::Spring& JelloMesh::Spring::operator=(const JelloMesh::Spring& p)
{
    if (&p == this) return *this;

    m_type = p.m_type;
    m_p1 = p.m_p1;
    m_p2 = p.m_p2;
    m_Ks = p.m_Ks;
    m_Kd = p.m_Kd;
    m_restLen = p.m_restLen;
    return *this;
}

JelloMesh::Spring::Spring(JelloMesh::SpringType t, 
    int p1, int p2, double Ks, double Kd, double restLen) :
    m_type(t), m_Ks(Ks), m_Kd(Kd), m_p1(p1), m_p2(p2), m_restLen(restLen)
{
}

//---------------------------------------------------------------------
// Particle
//---------------------------------------------------------------------

JelloMesh::Particle JelloMesh::Particle::EMPTY;

JelloMesh::Particle::Particle(int idx, const vec3& p, const vec3& v, double m)
{
    index = idx;
    position = p;
    velocity = v;
    force = vec3(0,0,0);
    mass = m;
}

JelloMesh::Particle::Particle() : index(-1), position(0,0,0), velocity(0,0,0), force(0,0,0), mass(1.0)
{
}

JelloMesh::Particle::Particle(const JelloMesh::Particle& p) : 
    index(p.index), position(p.position), velocity(p.velocity), force(p.force), mass(p.mass)
{
}

JelloMesh::Particle& JelloMesh::Particle::operator=(const JelloMesh::Particle& p)
{
    if (&p == this) return *this;

    index = p.index;
    position = p.position;
    velocity = p.velocity;
    force = p.force;
    mass = p.mass;
    return *this;
}

//---------------------------------------------------------------------
// Intersection
//---------------------------------------------------------------------

JelloMesh::Intersection::Intersection() : 
    m_p(-1), m_normal(0,0,0), m_distance(0) , m_type(CONTACT)
{
}

JelloMesh::Intersection::Intersection(const JelloMesh::Intersection& p) :
    m_p(p.m_p), m_normal(p.m_normal), m_distance(p.m_distance), m_type(p.m_type)
{
}

JelloMesh::Intersection& JelloMesh::Intersection::operator=(const JelloMesh::Intersection& p)
{
    if (&p == this) return *this;
    m_p = p.m_p;
    m_normal = p.m_normal;
    m_distance = p.m_distance;
    m_type = p.m_type;
    return *this;
}

JelloMesh::Intersection::Intersection(IntersectionType type, int p, const vec3& normal, double d) :
    m_p(p), m_normal(normal), m_distance(d), m_type(type)
{
}


//---------------------------------------------------------------------
// Drawing
//---------------------------------------------------------------------

void JelloMesh::FaceMesh::Draw(const JelloMesh& m)
{
    const ParticleGrid& g = m.m_vparticles;
    for (unsigned int strip = 0; strip < m_strips.size(); strip++)
    {
        const std::vector<int>& points = m_strips[strip];

        glBegin(GL_TRIANGLE_STRIP);
        for (unsigned int pi = 0; pi < points.size(); pi++)
        {
            int idx = points[pi];
            vec3 p = m.GetParticle(g, idx).position;

            vec3 n(0,0,0);
            const std::vector<int>& neighbors = m_neighbors[idx];
            if (neighbors.size() > 0)
            {
                vec3 pup = m.GetParticle(g, neighbors[0]).position;
                vec3 pdown = m.GetParticle(g, neighbors[1]).position;
                vec3 pleft = m.GetParticle(g, neighbors[2]).position;
                vec3 pright = m.GetParticle(g, neighbors[3]).position;

                vec3 n1 = -((pright - p) ^ (pup - p));
                vec3 n2 = -((pdown - p) ^ (pright - p));
                vec3 n3 = -((pleft - p) ^ (pdown - p));
                vec3 n4 = -((pup - p) ^ (pleft - p));

                n = n1 + n2 + n3 + n4;
                n = n.Normalize();
            }

            glNormal3f(n[0], n[1], n[2]);
            glVertex3f(p[0], p[1], p[2]);
        }
        glEnd();
    }
}

void JelloMesh::FaceMesh::DrawNormals(const JelloMesh& m)
{
    glDisable(GL_LIGHTING);

    glBegin(GL_LINES);
    glColor3f(0.0, 1.0, 0.0);

    const ParticleGrid& g = m.m_vparticles;
    for (unsigned int strip = 0; strip < m_strips.size(); strip++)
    {
        const std::vector<int>& points = m_strips[strip];
        for (unsigned int pi = 0; pi < points.size(); pi++)
        {
            int idx = points[pi];
            vec3 p = m.GetParticle(g, idx).position;

            const std::vector<int>& neighbors = m_neighbors[idx];
            if (neighbors.size() == 0) continue;

            vec3 pup = m.GetParticle(g, neighbors[0]).position;
            vec3 pdown = m.GetParticle(g, neighbors[1]).position;
            vec3 pleft = m.GetParticle(g, neighbors[2]).position;
            vec3 pright = m.GetParticle(g, neighbors[3]).position;

            vec3 n1 = -((pright - p) ^ (pup - p));
            vec3 n2 = -((pdown - p) ^ (pright - p));
            vec3 n3 = -((pleft - p) ^ (pdown - p));
            vec3 n4 = -((pup - p) ^ (pleft - p));

            vec3 n = n1 + n2 + n3 + n4;
            n = n.Normalize();

            vec3 end = p + 0.2 * n;
            glVertex3f(p[0], p[1], p[2]);
            glVertex3f(end[0], end[1], end[2]);
        }
    }

    glEnd();
    glEnable(GL_LIGHTING);
}

#define R(i) max(0, min(i, m.m_rows)) // CLAMP row index
#define C(j) max(0, min(j, m.m_cols)) // CLAMP col index
#define D(j) max(0, min(j, m.m_stacks)) // CLAMP stack index
JelloMesh::FaceMesh::FaceMesh(const JelloMesh& m, JelloMesh::Face f)
{
    const ParticleGrid& g = m.m_vparticles;
    switch(f)
    {
    case ZFRONT:
        m_strips.resize(m.m_rows);
        for (int i = 0; i < m.m_rows+1; i++)
            for (int j = 0; j < m.m_cols+1; j++)
            {
                if (i < m.m_rows)
                {
                    m_strips[i].push_back(m.GetIndex(i+1,j,0));
                    m_strips[i].push_back(m.GetIndex(i,j,0));
                }

                std::vector<int> neighbors;
                neighbors.push_back(m.GetIndex(R(i), C(j+1), D(0)));
                neighbors.push_back(m.GetIndex(R(i), C(j-1), D(0)));
                neighbors.push_back(m.GetIndex(R(i-1), C(j), D(0)));
                neighbors.push_back(m.GetIndex(R(i+1), C(j), D(0)));
                m_neighbors[m.GetIndex(i,j,0)] = neighbors;
            }
        break;
    case ZBACK:
        m_strips.resize(m.m_rows);
        for (int i = 0; i < m.m_rows+1; i++)
            for (int j = 0; j < m.m_cols+1; j++)
            {
                if (i < m.m_rows)
                {
                    m_strips[i].push_back(m.GetIndex(i+1,j,m.m_stacks));
                    m_strips[i].push_back(m.GetIndex(i,j,m.m_stacks));
                }

                std::vector<int> neighbors;
                neighbors.push_back(m.GetIndex(R(i+1), C(j), D(m.m_stacks)));
                neighbors.push_back(m.GetIndex(R(i-1), C(j), D(m.m_stacks)));
                neighbors.push_back(m.GetIndex(R(i), C(j-1), D(m.m_stacks)));
                neighbors.push_back(m.GetIndex(R(i), C(j+1), D(m.m_stacks)));
                m_neighbors[m.GetIndex(i,j,m.m_stacks)] = neighbors;
            }
        break;
    case XLEFT:
        m_strips.resize(m.m_cols);
        for (int j = 0; j < m.m_cols+1; j++)
            for (int k = 0; k < m.m_stacks+1; k++)
            {
                if (j < m.m_cols)
                {
                    m_strips[j].push_back(m.GetIndex(0,j+1,k));
                    m_strips[j].push_back(m.GetIndex(0,j,k));
                }

                std::vector<int> neighbors;
                neighbors.push_back(m.GetIndex(R(0), C(j), D(k+1)));
                neighbors.push_back(m.GetIndex(R(0), C(j), D(k-1)));
                neighbors.push_back(m.GetIndex(R(0), C(j-1), D(k)));
                neighbors.push_back(m.GetIndex(R(0), C(j+1), D(k)));
                m_neighbors[m.GetIndex(0,j,k)] = neighbors;
            }
        break;
    case XRIGHT:
        m_strips.resize(m.m_cols);
        for (int j = 0; j < m.m_cols+1; j++)
            for (int k = 0; k < m.m_stacks+1; k++)
            {
                if (j < m.m_cols)
                {
                    m_strips[j].push_back(m.GetIndex(m.m_rows,j+1,k));
                    m_strips[j].push_back(m.GetIndex(m.m_rows,j,k));
                }

                std::vector<int> neighbors;
                neighbors.push_back(m.GetIndex(R(m.m_rows), C(j+1), D(k)));
                neighbors.push_back(m.GetIndex(R(m.m_rows), C(j-1), D(k)));
                neighbors.push_back(m.GetIndex(R(m.m_rows), C(j), D(k-1)));
                neighbors.push_back(m.GetIndex(R(m.m_rows), C(j), D(k+1)));
                m_neighbors[m.GetIndex(m.m_rows,j,k)] = neighbors;
            }
        break;
    case YBOTTOM:
        m_strips.resize(m.m_rows);
        for (int i = 0; i < m.m_rows+1; i++)
            for (int k = 0; k < m.m_stacks+1; k++)
            {
                if (i < m.m_rows)
                {
                    m_strips[i].push_back(m.GetIndex(i+1,0,k));
                    m_strips[i].push_back(m.GetIndex(i,0,k));
                }

                std::vector<int> neighbors;
                neighbors.push_back(m.GetIndex(R(i+1), C(0), D(k)));
                neighbors.push_back(m.GetIndex(R(i-1), C(0), D(k)));
                neighbors.push_back(m.GetIndex(R(i), C(0), D(k-1)));
                neighbors.push_back(m.GetIndex(R(i), C(0), D(k+1)));
                m_neighbors[m.GetIndex(i,0,k)] = neighbors;
            }
        break;
    case YTOP:
        m_strips.resize(m.m_rows);
        for (int i = 0; i < m.m_rows+1; i++)
            for (int k = 0; k< m.m_stacks+1; k++)
            {
                if (i < m.m_rows)
                {
                    m_strips[i].push_back(m.GetIndex(i+1,m.m_cols,k));
                    m_strips[i].push_back(m.GetIndex(i,m.m_cols,k));
                }

                std::vector<int> neighbors;
                neighbors.push_back(m.GetIndex(R(i), C(m.m_cols), D(k+1)));
                neighbors.push_back(m.GetIndex(R(i), C(m.m_cols), D(k-1)));
                neighbors.push_back(m.GetIndex(R(i-1), C(m.m_cols), D(k)));
                neighbors.push_back(m.GetIndex(R(i+1), C(m.m_cols), D(k)));
                m_neighbors[m.GetIndex(i,m.m_cols,k)] = neighbors;
            }
        break;
    }
}

void JelloMesh::FaceMesh::CalcDistToEye(const JelloMesh& m, const vec3& eyePos)
{
    std::vector<int> points = m_strips[(int) (m_strips.size()*0.5)];
    int idx = points[(int) (points.size()*0.5)];
    vec3 pos = m.GetParticle(m.m_vparticles, idx).position;
    distToEye = (pos - eyePos).Length();
}

bool JelloMesh::FaceMesh::compare(const FaceMesh& one, const FaceMesh& other)
{
    return one.distToEye > other.distToEye;
}

/*
void JelloMesh::rayHit(vec3 m_start,vec3 m_end)
{
	ParticleGrid& g = m_vparticles;
	vec3 normalRay = (m_end-m_start).Normalize();
	double min_length = 100000.0;

	int hit_index = 0;

	Particle& particle = GetParticle(ParticleGrid& grid, i,j,k);

	for (int i = 0; i < rows+1; i++)
    {
		for (int j = 0; j < cols+1; j++)
        {
            for (int k = 0; k < stacks+1; k++)
            {
				
				if(i*j*k*(rows-i)*(cols-j)*(stacks-k)==0)
				{
					
					Particle& p = theJello.GetParticle(g, i,j,k);
					vec3 p_vertipoint = (p.position*normalRay)*normalRay;
					double p_vertilength = (p.position-p_vertipoint).Length();
					double p_longilength = (m_start-p_vertipoint).Length();
					if(p_vertilength<=0.01)
					{
						if(min_length>p_longilength)
						{
							min_length=p_longilength;
							hit_index = p.
						}
					}

				}
				
            }
        }
    }
	cout<<particle.position<<endl;

}
*/
