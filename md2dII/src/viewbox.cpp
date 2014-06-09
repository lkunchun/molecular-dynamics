//
// "$Id: viewbox.cxx,v 1.1 2002-12-18 16:39:34-08 kunchun Exp kunchun $"
//

#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Scroll.H>
#include <FL/Fl_Hor_Value_Slider.H>
#include <FL/fl_draw.H>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <unistd.h>
#include <math.h>
#include "multiplot.h"

//#include "drawing_area.h"

class ViewBox : public Fl_Widget {
  
  double *pos;  //in units of sigma
  double *box;  //in units of sigma
  double *vel;  //use for color, maximum = 3kT
  double *size;
  double scale;  //scale varies from 1 to --
  int natom;
  int center[2];  //varies from 0 to 1
  int scenter[2];
  int sx,sy,sw,sh; //selection box
  Fl_Slider *scale_slider;
  

  #define ViewBox_DIM 2
  #define ViewBox_DEFAULT_MIN_VIEW 5

  void erase_box() {
     window()->make_current();
     fl_overlay_clear();
  }

  int handle(int event) {
     static int ix, iy;
     static int dragged=0;
     static int button;
     int x2,y2;
     switch (event) {
     case FL_PUSH:
       erase_box();
       ix = Fl::event_x(); if (ix<x()) ix=x(); if (ix>=x()+w()) ix=x()+w()-1;
       iy = Fl::event_y(); if (iy<y()) iy=y(); if (iy>=y()+h()) iy=y()+h()-1;
       dragged = 0;
       button = Fl::event_button();
       //fprintf(stderr, "button is %d \n", button);
       //cerr << "hello" << endl;
       return 1;
     case FL_DRAG:
       dragged = 1;
       erase_box();
       x2 = Fl::event_x(); if (x2<x()) x2=x(); if (x2>=x()+w()) x2=x()+w()-1;
       y2 = Fl::event_y(); if (y2<y()) y2=y(); if (y2>=y()+h()) y2=y()+h()-1;
       if (button != 1) {ix = x2; iy = y2; return 1;}
       if (ix < x2) {sx = ix; sw = x2-ix;} else {sx = x2; sw = ix-x2;}
       if (iy < y2) {sy = iy; sh = y2-iy;} else {sy = y2; sh = iy-y2;}
       window()->make_current();
       if (sw>sh) {
          sh = sw;
       } else {
          sw = sh;
       }
       fl_overlay_rect(sx,sy,sw,sh);
       return 1;
     case FL_RELEASE:
       erase_box();
       if (dragged == 1) {
          scale = w()/sw;    
          scenter[0] = ((sx+sw/2));
          scenter[1] = ((sy+sw/2));
          fprintf(stderr, "center: %d, %d ixy: %d, %d sw: %d\n", center[0], center[1],ix,iy,sw);
       } else {
          scale = 1;
          scenter[0] = x()+w()/2;
          scenter[1] = y()+h()/2;
       }
       draw();
       return 1;
     }
     return 0;
  }

  void draw() {
    int i;
    double temp_pos[2];
    double temp_vel2;
    double temp_size;
    fl_clip(x(),y(),w(),h());
    fl_color(FL_BLACK);
    fl_rectf(x(),y(),w(),h());
    fl_push_matrix();

    //fl_translate(x(),y());
    //fl_translate(x()+w()/2.0, y()+h()/2.0);
    fl_translate(center[0], center[1]);
    //fl_translate(x()+w()/2.0, y()+h()/2.0);
    fl_scale(w()*scale,h()*scale);

/* for rotation
    if (args[5]) {
      fl_translate(x()+w()/2.0, y()+h()/2.0);
      fl_rotate(args[5]);
      fl_translate(-(x()+w()/2.0), -(y()+h()/2.0));
    }
*/

    fl_color(FL_GREEN);

    if (pos != 0) {
       if(scale > box[0]/ViewBox_DEFAULT_MIN_VIEW) {
          scale = box[0]/ViewBox_DEFAULT_MIN_VIEW;
       }
    }

    for(i=0;i<natom;i++) {
       temp_pos[0] = pos[2*i];
       temp_pos[1] = pos[2*i+1];
       temp_size = size[i];
       temp_pos[0] /= box[0];
       temp_pos[0] = temp_pos[0] - rint(temp_pos[0]);
       temp_pos[0] += ((float) center[0]-x())/w() - ((float) scenter[0]-x())/w();
       //temp_pos[0] *= scale*((float ) (w()));
       //temp_pos[0] *= ((float ) (w()));
       //temp_pos[0] += center[0]*scale*((float ) (w()));
       temp_pos[1] /= box[1];
       temp_pos[1] = temp_pos[1] - rint(temp_pos[1]);
       temp_pos[1] += ((float) center[1]-x())/h() - ((float) scenter[1]-y())/h();
       //temp_pos[1] *= scale*((float ) (h()));
       //temp_pos[1] *= ((float ) (h()));
       //temp_pos[1] += center[1]*scale*((float ) (w()));
       temp_size /= box[0];
       //temp_size *= scale*((float ) (w()-x()));
       //temp_size *= ((float ) (w()));

       temp_vel2 = vel[2*i]*vel[2*i];
       temp_vel2 += vel[2*i+1]*vel[2*i+1];
       temp_vel2 = temp_vel2*255/5.0;
       //fl_color((uchar) ((int) temp_vel2), 0, 255);
       fl_color(255, 0, 0);
       //view box only look at the left bottom corner

       //if(fabs(temp_pos[0]) < 0.5 && fabs(temp_pos[1]) < 0.5) { 
          fl_begin_polygon();
          fl_circle((float) temp_pos[0],(float) temp_pos[1],(float) temp_size/2.0);
          fl_end_polygon();

//          fl_color(FL_BLUE);
       //} else {
       //   fprintf(stderr, "bummer\n");
       //}

       
    }

    //fl_translate(-(x()+w()/2.0), -(y()+h()/2.0));
    //fl_translate(center[0], center[1]);
    //fl_scale(w()*scale,h()*scale);

    fl_pop_matrix();

//    fl_push_matrix();
//       fl_translate(center[0],center[1]);
//       fl_scale(w()*scale,h()*scale);
//    fl_pop_matrix();

    fl_pop_clip();

  }

//  void slider_cb(Fl_Widget* o, void* v) {
//     Fl_Slider *temp_slider = (Fl_Slider*) o;
 //    scale = temp_slider->value();
     //temp_vb->redraw();
  //}

public:

  ViewBox(int X,int Y,int W,int H) : Fl_Widget(X,Y,W,H) { pos = 0;
     size=0; scale = 1; box=0; vel=0;natom=0; scenter[0]=X+W/2;
     scenter[1]=Y+H/2; center[0] = scenter[0]; center[1]=scenter[1];
//     scale_slider = new Fl_Hor_Value_Slider(50,310,240,25,"Scale"); 
//     scale_slider->minimum(1.0); s->maximum(ViewBox_DEFAULT_MIN_VIEW);
//     scale_slider->callback(slider_cb,(void*) 1);
//     scale_slider->align(FL_ALIGN_LEFT);
  }

  void setup(int n,double* p, double *s, double *b, double *v, int *c,double sc) {
     natom = n; pos = p; size =s; box=b;vel=v;center[0]=*c;center[1]=*(c+1);scale=sc;
  }

  void vb_rescale(double s) {
     if(s >= 1.0) scale = s;
  }

  double boxL() {
     return box[0];
  }

  void vb_recenter(int x, int val) {
     center[x] = val;
  }

};

extern "C"
{
   void ballview_(int*,double*,double*,double*,double*);
}


static int c[] = {150,150};
void ballview_(int* natom, double* pos, double* size, double* bx, double *vel) {

   static Fl_Double_Window *window=0;
   static ViewBox *vbox=0;
   static Fl_Slider* s = 0;
   static int waiting_time = 10;
   static int waiting_counter = 0;
   if(window==0) {
      window = new Fl_Double_Window(310,360);
      vbox = new ViewBox(10,10,280,280);
      vbox->setup( *natom, pos, size,bx,vel,c,1.0);
      window->end();
      window->show();
   }

   if(Fl::check()) {
      if((waiting_counter%=waiting_time)==0) {
         vbox->redraw();
      }
      //Fl::wait(0.2);
      usleep(100);
      waiting_counter++;
   } else {
      fprintf(stderr,"error\n");
   }

   return;
 
}



/*

static double p[] = {0.0,0.0,4.0,4.0, -4.0,-4.0};
static double b[] = {10.0,10.0};
static double s[] = {1.0,1.0,1.0};
static int c[] = {150,150};
static long nat = 3;

void center_h_cb(Fl_Widget* o, void* v) {
   ((ViewBox*) v) -> vb_recenter(0,(((Fl_Scrollbar*) o) -> value()));
   ((ViewBox*) v) -> redraw();
}

void center_v_cb(Fl_Widget* o, void* v) {
   ((ViewBox*) v) -> vb_recenter(1,(((Fl_Scrollbar*) o) -> value()));
   ((ViewBox*) v) -> redraw();
}

int main(int argc, char**argv) {

  Fl_Double_Window window(310,360);

  ViewBox vb (10,10,280,280);
  vb.setup(nat, p,s,b,c,1.0);


*
  Fl_Scrollbar* centering_h = new Fl_Scrollbar(10,290,280,15);
  centering_h->type(FL_HORIZONTAL); 
  centering_h->minimum(10); centering_h->maximum(280);
  Fl_Scrollbar* centering_v = new Fl_Scrollbar(290,10,15,280);
  centering_v->minimum(10); centering_v->maximum(280);
  centering_h -> callback(center_h_cb, (void*) &vb);
  centering_h -> linesize(1); 
  centering_v -> callback(center_v_cb, (void*) &vb);
  centering_v -> linesize(1);
*


  Fl_Slider* s = new Fl_Hor_Value_Slider(50,310,240,25,"Scale"); 
  s->minimum(1.0); s->maximum(vb.boxL()/ViewBox_DEFAULT_MIN_VIEW);
  s->callback(slider_cb,(void*) &vb);
  s->align(FL_ALIGN_LEFT);
  
//  window.resizable(vb);
//  window.resizable(s);
  window.end();
  window.show(argc,argv);

  return Fl::run();


}

*/

