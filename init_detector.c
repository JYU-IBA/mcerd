#include "general.h"
#include "prototypes.h"

void init_detector(Global *global, Detector *det) {
    Point p1, p2, p3;
    double d, tmax1, tmax2;
    int i;

    for (i = 0; i < det->nfoils; i++) {
        switch (det->foil[i].type) {
            case FOIL_CIRC:
                p1.x = sin(det->angle) * det->foil[i].dist;
                p1.y = 0.0;
                p1.z = cos(det->angle) * det->foil[i].dist;
                p2 = p1;
                p2.y += det->foil[i].size[0];
                p3 = p1;
                p3.x += cos(det->angle) * det->foil[i].size[0];
                p3.z -= sin(det->angle) * det->foil[i].size[0];
                break;
            case FOIL_RECT:
                p1.x = sin(det->angle) * det->foil[i].dist;
                p1.y = 0.0;
                p1.z = cos(det->angle) * det->foil[i].dist;
                p2 = p1;
                p2.y += det->foil[i].size[1];
                p3 = p1;
                p3.x += cos(det->angle) * det->foil[i].size[0];
                p3.z -= sin(det->angle) * det->foil[i].size[0];
                break;
                break;
            default:
                break;
        }
        det->foil[i].center = p1;
        det->foil[i].plane = get_plane_params(&p1, &p2, &p3);
        det->foil[i].virtual = FALSE;
    }

    if (global->virtualdet) {
        det->vfoil = det->foil[0];
        if (det->vfoil.type == FOIL_CIRC)
            det->vfoil.size[1] = det->vfoil.size[0] * det->vsize[1];
        else
            det->vfoil.size[1] *= det->vsize[1];
        det->vfoil.size[0] *= det->vsize[0];
        det->vfoil.virtual = TRUE;
    }

    if (det->foil[0].type == FOIL_CIRC) {
        p1.x = global->bspot.x;
        p1.y = global->bspot.y;
        p1.z = p1.x / tan(C_PI / 2.0 - global->beamangle);
        p2.x = det->foil[0].center.x - cos(det->angle) * det->foil[0].size[0];
        p2.y = det->foil[0].size[0];
        p2.z = det->foil[0].center.z + sin(det->angle) * det->foil[0].size[0];
        p3 = coord_transform(p1, det->angle, C_PI, p2, FORW);
        d = sqrt(ipow2(p3.x) + ipow2(p3.y) + ipow2(p3.z));
        tmax1 = acos(p3.z / d);
        p1.x = -p1.x;
        p1.y = -p1.y;
        p1.z = -p1.z;
        p2.x = det->foil[0].center.x + cos(det->angle) * det->foil[0].size[0];
        p2.y = -det->foil[0].size[0];
        p2.z = det->foil[0].center.z - sin(det->angle) * det->foil[0].size[0];
        p3 = coord_transform(p1, det->angle, C_PI, p2, FORW);
        d = sqrt(ipow2(p3.x) + ipow2(p3.y) + ipow2(p3.z));
        tmax2 = acos(p3.z / d);
        det->thetamax = max(tmax1, tmax2);
    } else if (det->foil[0].type == FOIL_RECT) {
        p1.x = global->bspot.x;
        p1.y = global->bspot.y;
        p1.z = p1.x / tan(C_PI / 2.0 - global->beamangle);
        p2.x = det->foil[0].center.x - cos(det->angle) * det->foil[0].size[0];
        p2.y = det->foil[0].size[1];
        p2.z = det->foil[0].center.z + sin(det->angle) * det->foil[0].size[0];
        p3 = coord_transform(p1, det->angle, C_PI, p2, FORW);
        d = sqrt(ipow2(p3.x) + ipow2(p3.y) + ipow2(p3.z));
        tmax1 = acos(p3.z / d);
        p1.x = -p1.x;
        p1.y = -p1.y;
        p1.z = -p1.z;
        p2.x = det->foil[0].center.x + cos(det->angle) * det->foil[0].size[0];
        p2.y = -det->foil[0].size[1];
        p2.z = det->foil[0].center.z - sin(det->angle) * det->foil[0].size[0];
        p3 = coord_transform(p1, det->angle, C_PI, p2, FORW);
        d = sqrt(ipow2(p3.x) + ipow2(p3.y) + ipow2(p3.z));
        tmax2 = acos(p3.z / d);
        det->thetamax = max(tmax1, tmax2);
    }

}

Plane get_plane_params(Point *p1, Point *p2, Point *p3) {
    Plane p;

    if (p3->x == p1->x || p2->y == p1->y) {
        fprintf(stderr, "Geometry is not supported\n");
        exit(3);
    }

    p.type = GENERAL_PLANE;

    p.a = (p3->z - p1->z) / (p3->x - p1->x);
    p.b = (p2->z - p1->z) / (p2->y - p1->y);
    p.c = p1->z - p.a * p1->x - p.b * p1->y;

    return (p);
}

#if 0
/*
   det->thetamin = PI;
   det->thetamax = 0.0;
   det->fiimin = 2.0*PI;
   det->fiimax = 0.0;

   det->brect = get_beam_rect(global->beamangle,det->bsize,det->bcenter);

   det->aerect = get_mask_rect(det->adist,det->dangle,det->asize,
                              det->aheigh,det->enlarge_x,det->enlarge_y);

   det->arect = get_mask_rect(det->adist,det->dangle,det->asize,
                              det->aheigh,1.0,1.0);

   det->absrect = get_mask_rect(det->absdist,det->dangle,det->asize,
                                det->aheigh,1.0,1.0);

   det->drect = get_mask_rect(det->ddist,det->dangle,det->dsize,
                              det->dheigh,1.0,1.0);

   det->aplane = get_plane_params(&(det->arect));

   det->absplane = get_plane_params(&(det->absrect));

   det->dplane = get_plane_params(&(det->drect));      

*/
/*
   for(i=0;i<det->nfoils;i++){
      det->fcirc[i].r = det->diam[i];
      det->fcirc[i].p.x = det->dist[i]*sin(det->angle);
      det->fcirc[i].p.y = 0.0;
      det->fcirc[i].p.z = det->dist[i]*cos(det->angle);
      det->fplane[i] = get_circ_plane(&(det->fcirc[i]),det->angle);
   }
*/
/*   
   for(i=0;i<4;i++){
      for(j=0;j<4;j++){
         dir = get_lab_angle(&(det->brect.p[i]),&(det->aerect.p[j]));
         if(dir.theta < det->thetamin)
            det->thetamin = dir.theta;
         if(dir.theta > det->thetamax)
            det->thetamax = dir.theta;
         if(dir.fii < det->fiimin)
            det->fiimin = dir.fii;
         if(dir.fii > det->fiimax)
            det->fiimax = dir.fii;
      }
   }
*/      

Rect get_beam_rect(double angle,Point2 size,Point2 center)
{
   Rect r;
   double dx,dy,dz;
   
   dx = 0.5*size.x;
   dy = 0.5*size.y;
   dz = 0.5*size.x/cos(angle);
   
   r.p[0].x = center.x - dx;
   r.p[0].y = center.y - dy;
   r.p[0].z = -dz;

   r.p[1].x = r.p[0].x;
   r.p[1].y = r.p[0].y + 2.0*dy;
   r.p[1].z = r.p[0].z;

   r.p[2].x = r.p[1].x + 2.0*dx;
   r.p[2].y = r.p[1].y;
   r.p[2].z = r.p[1].z + 2.0*dz;

   r.p[3].x = r.p[2].x;
   r.p[3].y = r.p[2].y - 2.0*dy;
   r.p[3].z = r.p[2].z;

   return(r);

}
Vector get_lab_angle(Point *p1,Point *p2)
{
   Vector d;
   Point dp;
   double r;
   
   dp.x = p2->x - p1->x;
   dp.y = p2->y - p1->y;
   dp.z = p2->z - p1->z;      

   r = sqrt(ipow2(dp.x) + ipow2(dp.y) + ipow2(dp.z));

   d.theta = acos(dp.z/r);
   d.fii = atan2(dp.y,dp.x);
   
   return(d);
   
}
Rect get_mask_rect(double dist,double angle,Point2 size,double
                   heigh,double enlarge_x,double enlarge_y)
{
   Rect r;
   Point center;
   double dx,dy,dz;

   center.x = dist*sin(angle);
   center.y = heigh;
   center.z = dist*cos(angle);

   dx = 0.5*enlarge_x*size.x*cos(angle);
   dy = 0.5*enlarge_y*size.y;
   dz = 0.5*enlarge_x*size.x*sin(angle);

   r.p[0].x = center.x - dx;
   r.p[0].y = center.y - dy;
   r.p[0].z = center.z + dz;

   r.p[1].x = center.x - dx;
   r.p[1].y = center.y + dy;
   r.p[1].z = center.z + dz;

   r.p[2].x = center.x + dx;
   r.p[2].y = center.y + dy;
   r.p[2].z = center.z - dz;

   r.p[3].x = center.x + dx;
   r.p[3].y = center.y - dy;
   r.p[3].z = center.z - dz;


   return(r);

}

Plane get_circ_plane(Circ *circ,double angle)
{
   Plane p;

   p.type = GENERAL_PLANE;

   p.a = -tan(angle);
   p.b = 0.0;
   p.c = circ->p.z - p.a*circ->p.x - p.b*circ->p.y;

   return(p);

}

#endif
