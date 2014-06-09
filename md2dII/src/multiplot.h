/*
 * multiplot.h
 * Copyright 2002 by Andre Krause.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 * USA.
 *
 * Please report all bugs and problems to "post@andre-krause.net".
 * New versions and bug-fixes under http://www.andre-krause.net
 *	
 */
#ifndef __MULTIPLOT__
#define __MULTIPLOT__


#include <FL/Fl.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/gl.h>
#include <float.h>
#include <time.h>
#include <math.h>
#include <vector>
#include <string>
#include <sstream>


#ifndef min
	#define min(a,b)( (a) < (b) ) ? (a) : (b) )
#endif



using namespace std;

enum MULTIPLOT_GRIDSTYLE
{
	MP_NO_GRID,
	MP_LINEAR_GRID,
	MP_LOG_GRID
};

/**
 * this class defines a plot point, wich essentially consists of its coordinate
 * and its color in rgb-color space.
 */
class PLOT_POINT
{
	public:
		float x;
		float y;
		float r;
		float g;
		float b;
		PLOT_POINT();
		PLOT_POINT(float xx, float yy, float rr=1, float gg=1, float bb=1);
};

/**
 *	this class creates a window to wich you can add an arbitrary number of 
 *	autoscaling traces.
 */
class MULTIPLOT : public Fl_Gl_Window
{
	public:
		/**
		 *	this constructor tells multiplot where to put the window on the 
		 *	desktop in pixel-coordinates(x,y) and with wich width and height (w,h)
		 */
		MULTIPLOT(int x, int y, int w, int h, const char *title="MULTIPLOT - updates on www.andre-krause.net");
		/**
		 *	with this function you cann add a plot-points to a trace. traces are numbered
		 *	from zero to N. memory for the traces is automatically allocated. 
		 */
		void add(unsigned int trace, const PLOT_POINT &p);
		/**
		 *	set the maximum number of points to be plotted. this is useful to
		 *	avoid slow drawing of your trace. if you have 1000 plot-points and 
		 *	set the number of max_points to 100, then only every tenth point gets plotted.
		 */
		void set_max_points(int mx);
		/**
		 *	if you call set_scrolling with a positive number of points to be plotted,
		 *	your graph will scroll left out of the plot-window as you add new plot-points.
		 *	Zero or a negative number disables scrolling.
		 */
		void set_scrolling(int max_points_to_plot);
		
		/**
		 *	call this function if you wish a grid to be plotted in your graph.
		 *	by default, no grids are plotted. call this function with the first 
		 *	two arguments set to either MP_NO_GRID,	MP_LINEAR_GRID or MP_LOG_GRID.
		 *	the next two arguments gridx_step and gridy_step specify the grid spacing.
		 *	Zero or a negative value like -1 enables auto - spacing.
		 *	The last parameter w sets the grid-linewidth. the default is 1 pixel.
		 */
		void set_grid(int gridx, int gridy, float gridx_step=-1.0, float gridy_step=-1.0, float w=1.0);
		
		/**
		 *	this function sets the size of the plot-points. the default value is zero, so 
		 *	no points are drawn at all. if you wish to create a scatter-plot, set the pointsize
		 *  to a value bigger than zero and the linesize to zero.
		 */
		void set_pointsize(unsigned int trace, float s);
		
		/**
		 *	call set_linewidth to change the thickness of the traces. the default
		 *	value is 1 pixel, if you set the linewidth to zero, no lines are drawn. this
		 *	is usefull to create scatter-plots.
		 */
		void set_linewidth(unsigned int trace, float w);
		
		/**
		*	sets the background color
		*/
		void set_bg_color(float r, float g, float b);

		/**
		*	sets the grid color
		*/
		void set_grid_color(float r, float g, float b);

		/**
		 *	this function call simply clears all traces
		 */
		void clear();
		
		/**
		*	this function call clears trace number t
		*/
		void clear(unsigned int t);
	protected:
		string title;		// stores the user-title, so we can add ranges
		vector< vector<PLOT_POINT> > traces;
		PLOT_POINT minimum, maximum;
		PLOT_POINT scale;
		PLOT_POINT offset;
		PLOT_POINT bg_color;
		PLOT_POINT grid_color;
		unsigned int max_points;
		bool scroll;

		// this struct stores information
		// about the traces like line_width, the current pos 
		// in the ring-buffer and the point-size.
		class TRACE
		{
			public:
				unsigned int pos;
				float point_size;
				float line_width;
				TRACE(){ pos=0; point_size=0.0; line_width=1; }
		};
		vector< TRACE > tsettings;

		// grid - vars
		int gridx;
		int gridy;
		float gridx_step, gridy_step, grid_linewidth;		

		void initgl();
		PLOT_POINT draw_grid();
		void draw();
};


// implementation:
inline PLOT_POINT::PLOT_POINT() { x=y=0.0;r=g=b=1; }
inline PLOT_POINT::PLOT_POINT(float xx, float yy, float rr, float gg, float bb) { x=xx;y=yy;r=rr;g=gg;b=bb; }

inline MULTIPLOT::MULTIPLOT(int x, int y, int w, int h, const char *ttitle) : Fl_Gl_Window(x,y,w,h,ttitle) 
{	
	title=string(ttitle);
	//max.x=max.y=-FLT_MAX;
	maximum.x=maximum.y=-FLT_MAX;
	minimum.x=minimum.y=FLT_MAX;
	max_points=max(w,h);
	traces.push_back( vector<PLOT_POINT>() ); // create one trace
	tsettings.push_back(TRACE());

	scroll=false;
	
	gridx=MP_NO_GRID;
	gridy=MP_NO_GRID;
	gridx_step=-1;
	gridy_step=-1;

	bg_color.r=0.0f;
	bg_color.g=0.0f;
	bg_color.b=0.0f;
	
	grid_color.r=0.8f;
	grid_color.g=0.8f;
	grid_color.b=0.8f;
	grid_linewidth=1.0f;

	show();
}

inline void MULTIPLOT::add(unsigned int trace, const PLOT_POINT &p)
{ 
	if(traces.size()<=trace)
	{
		int n=trace-traces.size()+1;
		for(int a=0;a<n;a++)
		{
			traces.push_back(vector<PLOT_POINT>());
			tsettings.push_back(TRACE());			
		}
	}

	if(scroll)
	{
		// this realises a sort of ringbuffer wich is needed for scrolling
		if(tsettings[trace].pos<traces[trace].size())
			traces[trace][tsettings[trace].pos]=p;
		else
			traces[trace].push_back(p);
		tsettings[trace].pos++;
		if(tsettings[trace].pos>=max_points)
			tsettings[trace].pos=0;
	}
	else
	{
		traces[trace].push_back(p);
	}
}

inline void MULTIPLOT::set_pointsize(unsigned int trace, float psize)
{
	if(tsettings.size()<=trace)
	{
		int n=trace-tsettings.size()+1;
		for(int a=0;a<n;a++)
		{
			tsettings.push_back(TRACE());			
		}
	}
	tsettings[trace].point_size=psize;
}

inline void MULTIPLOT::set_linewidth(unsigned int trace, float width)
{
	if(tsettings.size()<=trace)
	{
		int n=trace-tsettings.size()+1;
		for(int a=0;a<n;a++)
		{
			tsettings.push_back(TRACE());			
		}
	}
	tsettings[trace].line_width=width;
}

inline void MULTIPLOT::clear()
{
	for(unsigned int a=0;a<traces.size();a++)
	{
		traces[a].clear();
		tsettings[a].pos=0;
	}	
}

inline void MULTIPLOT::clear(unsigned int trace)
{
		traces[trace].clear();
		tsettings[trace].pos=0;
}


inline void MULTIPLOT::set_max_points(int mx)
{
	max_points=mx;
}

inline void MULTIPLOT::set_grid_color(float r, float g, float b)
{
	grid_color.r=r;
	grid_color.g=g;
	grid_color.b=b;
}

inline void MULTIPLOT::set_bg_color(float r, float g, float b)
{
	bg_color.r=r;
	bg_color.g=g;
	bg_color.b=b;
	glClearColor(bg_color.r, bg_color.g, bg_color.b, 1);		// Set The background color
}


inline void MULTIPLOT::set_scrolling(int max_points_to_plot)
{
	if(max_points_to_plot<=0)
	{
		scroll=false;
		return;
	}

	scroll=true;
	max_points=max_points_to_plot;
	for(unsigned int t=0;t<traces.size();t++)
		if(traces[t].capacity()<max_points)
			traces[t].reserve(max_points);
	
}

inline void MULTIPLOT::set_grid(int ggridx, int ggridy, float ggridx_step, float ggridy_step, float w)
{
	gridx=ggridx;
	gridy=ggridy;
	gridx_step=ggridx_step;
	gridy_step=ggridy_step;
	grid_linewidth=w;
}


inline void MULTIPLOT::initgl()
{
	glViewport(0 , 0,w() ,h());	// Set Up A Viewport
	glMatrixMode(GL_PROJECTION);								// Select The Projection Matrix
	glLoadIdentity();											// Reset The Projection Matrix
	glOrtho( 0, w(), 0, h(), -1, 1 );							// Select Ortho Mode (640x480)
	//gluPerspective(50, (float)w()/(float)h(), 5,  2000); // Set Our Perspective
	glMatrixMode(GL_MODELVIEW);									// Select The Modelview Matrix
	glLoadIdentity();											// Reset The Modelview Matrix
	glDisable(GL_DEPTH_TEST);									// Enable Depth Testing
	glDisable(GL_LIGHTING);
	glShadeModel(GL_SMOOTH);									// Select Smooth Shading
	glClearColor(bg_color.r, bg_color.g, bg_color.b, 1);		// Set The background color
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);			// Clear The Screen And Depth Buffer
}	



inline PLOT_POINT MULTIPLOT::draw_grid()
{
	double xstep=0;
	double ystep=0;
	if(gridx==MP_LINEAR_GRID)
	{
		double diff=maximum.x - minimum.x;
		if(diff==0)return PLOT_POINT(0,0);
		double exp=floor(log10(fabs(diff)));
		double shiftscale=pow(10,exp);
		// get the starting point for the grid
		double startx=shiftscale*floor(minimum.x / shiftscale);

		if(gridx_step>0)
		{
			xstep=gridx_step;
		}
		else	// auto grid size
		{			
			xstep=shiftscale*1.0;
			if(diff/xstep < 4) // draw more lines
				xstep*=0.5;

			/*
			xstep=floor(maximum.x / shiftscale) - floor(minimum.x / shiftscale);
			xstep=2*floor(0.5*xstep);// if uneven, make even, because uneven stepsizes will cause uneven distributed lines around the (0,0) koordinate frame axis
			xstep=shiftscale*xstep / 4.0;
			*/
			
		}
		double x=startx;
		
		glLineWidth(grid_linewidth);
		glColor3f(grid_color.r, grid_color.g, grid_color.b);
		glBegin(GL_LINES);
		int bailout=0; // bailout is a safety to avoid endless recursions caused maybe through numerical errors..
		while(x<maximum.x && bailout<100)
		{
			x+=xstep;
			bailout++;
			glVertex2f((GLfloat)(x-offset.x)*scale.x	,0.0f  );
			glVertex2f((GLfloat)(x-offset.x)*scale.x	,(float)h());
		}
		glEnd();
	}
	if(gridy==MP_LINEAR_GRID)
	{
		double diff=maximum.y - minimum.y;
		if(diff==0)return PLOT_POINT(0,0);
		double exp=floor(log10(fabs(diff)));
		double shiftscale=pow(10,exp);
		// get the starting point for the grid
		double starty=shiftscale*floor(minimum.y / shiftscale);

		if(gridy_step>0)
		{
			ystep=gridy_step;
		}
		else	// auto grid size
		{
			ystep=shiftscale*1.0;
			if(diff/ystep < 4) // draw more lines
				ystep*=0.5;

			/*
			ystep=floor(maximum.y / shiftscale) - floor(minimum.y / shiftscale);
			ystep=2*floor(0.5*ystep); // if uneven, make even, because uneven stepsizes will cause uneven distributed lines around the (0,0) koordinate frame axis
			ystep=shiftscale*ystep / 4.0;
			*/

		}
		double y=starty;
		glLineWidth(1.0);
		glColor3f(grid_color.r, grid_color.g, grid_color.b);
		glBegin(GL_LINES);
		int bailout=0; // bailout is a safety to avoid endless recursions caused maybe through numerical errors..
		while(y<maximum.y && bailout<100)
		{
			y+=ystep;
			bailout++;
			glVertex2f(0.0f					,(GLfloat)(y-offset.y)*scale.y);
			glVertex2f((float)w()			,(GLfloat)(y-offset.y)*scale.y);
		}
		glEnd();
	}

	return PLOT_POINT((float)xstep,(float)ystep);
}

inline void MULTIPLOT::draw()
{
	if(!valid())
	{
		initgl();
		valid(1);
	}
	
	glClear(GL_COLOR_BUFFER_BIT);// | GL_DEPTH_BUFFER_BIT);			// Clear The Screen And Depth Buffer



	// draw the grid
	PLOT_POINT grid_spacing=draw_grid();

	// draw the coordinate cross with center (0,0)
	glLineWidth(2.0f*grid_linewidth);
	glColor3f(grid_color.r, grid_color.g, grid_color.b);
	glBegin(GL_LINES);
		glVertex2f(0.0f					,0-offset.y*scale.y);
		glVertex2f((float)w()			,0-offset.y*scale.y);
		glVertex2f(0-offset.x*scale.x	,0.0f  );
		glVertex2f(0-offset.x*scale.x	,(float)h());		
	glEnd();



	maximum.x=maximum.y=-FLT_MAX;
	minimum.x=minimum.y=FLT_MAX;


	if(scroll)
	{
		PLOT_POINT p;
		for(unsigned int t=0;t<traces.size();t++)
		{
			glLineWidth(tsettings[t].line_width);
			glBegin(GL_LINE_STRIP);
			unsigned int ps=tsettings[t].pos;
			unsigned int n=(unsigned int)traces[t].size();
			if(max_points<n)n=max_points;
			for(unsigned int a=0;a<n;a++)
			{
				if(ps>=traces[t].size())
					ps=0;
				p=traces[t][ps];
				if(tsettings[t].line_width>0)
				{
					glColor3f(p.r,p.g,p.b);
					glVertex2f((p.x-offset.x)*scale.x,(p.y-offset.y)*scale.y);
				}
				if(p.x>maximum.x)maximum.x=p.x;
				if(p.x<minimum.x)minimum.x=p.x;
				if(p.y>maximum.y)maximum.y=p.y;
				if(p.y<minimum.y)minimum.y=p.y;		
	
				ps++;
			}
			glEnd();
			glLineWidth(1.0);

			
			ps=tsettings[t].pos;
			n=(unsigned int)traces[t].size();
			glPointSize(tsettings[t].point_size);
			glBegin(GL_POINTS);
			for(unsigned int a=0;a<n;a++)
			{
				if(ps>=traces[t].size())
					ps=0;
				p=traces[t][ps];
				if(tsettings[t].point_size>0)
				{
					glColor3f(p.r,p.g,p.b);
					glVertex2f((p.x-offset.x)*scale.x,(p.y-offset.y)*scale.y);
				}
				ps++;
			}
			glEnd();
		}

	}
	else
	{
		PLOT_POINT p;
		// for speedup skip points - this is useful if we have many many points
		float step=1;
	
		for(unsigned int t=0;t<traces.size();t++)
		{
			if(traces[t].size()>max_points)
				step=traces[t].size()/float(max_points);

			glLineWidth(tsettings[t].line_width);
			glBegin(GL_LINE_STRIP);
			for(float a=0;a<traces[t].size();a+=step)
			{
				p=traces[t][int(a)];
				if(tsettings[t].line_width>0)
				{
					glColor3f(p.r,p.g,p.b);
					glVertex2f((p.x-offset.x)*scale.x,(p.y-offset.y)*scale.y);
				}
				if(p.x>maximum.x)maximum.x=p.x;
				if(p.x<minimum.x)minimum.x=p.x;
				if(p.y>maximum.y)maximum.y=p.y;
				if(p.y<minimum.y)minimum.y=p.y;		
			}
			glEnd();
			
			if(tsettings[t].point_size>0)
			{
				glPointSize(tsettings[t].point_size);
				glBegin(GL_POINTS);
				for(float a=0;a<traces[t].size();a+=step)
				{
					p=traces[t][int(a)];
					glColor3f(p.r,p.g,p.b);
					glVertex2f((p.x-offset.x)*scale.x,(p.y-offset.y)*scale.y);
				}
				glEnd();
			}
		}
	}

	float diff;	
	diff=maximum.x-minimum.x;	if(diff!=0.0)scale.x=w()/diff;
	diff=maximum.y-minimum.y;	if(diff!=0.0)scale.y=h()/diff;
	offset.x=minimum.x;
	offset.y=minimum.y;

	
	stringstream ss;
	ss << "range: x=[" <<  minimum.x << "; " << maximum.x << "] ";
	ss << "y=[" <<  minimum.y << "; " << maximum.y << "] ";	
	ss << "grid: w=[" <<  grid_spacing.x << "] ";
	ss << "h=[" <<  grid_spacing.y << "] ";
	this->label((ss.str() + title).c_str());

		
}

#endif
