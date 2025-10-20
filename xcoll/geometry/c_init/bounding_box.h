// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_GEOM_BOUNDING_BOX_H


/*gpufun*/
void BoundingBox_sync(BoundingBox box){
    double rC     = BoundingBox_get_rC(box);
    double sin_tb = BoundingBox_get_sin_tb(box);
    double cos_tb = BoundingBox_get_cos_tb(box);
    double sin_tC = BoundingBox_get_sin_tC(box);
    double cos_tC = BoundingBox_get_cos_tC(box);
    double proj_l = rC * (cos_tb*cos_tC + sin_tb*sin_tC); // projection of position vector on length
    BoundingBox_set_proj_l(box, proj_l);
    double proj_w = rC * (cos_tb*sin_tC - sin_tb*cos_tC); // projection of position vector on width
    BoundingBox_set_proj_w(box, proj_w);
}


/*gpufun*/
void BoundingBox_set(BoundingBox box, double rC, double sin_tC, double cos_tC,
                     double l, double w, double sin_tb, double cos_tb){
    BoundingBox_set_rC(box, rC);
    BoundingBox_set_sin_tC(box, sin_tC);
    BoundingBox_set_cos_tC(box, cos_tC);
    BoundingBox_set_l(box, l);
    BoundingBox_set_w(box, w);
    BoundingBox_set_sin_tb(box, sin_tb);
    BoundingBox_set_cos_tb(box, cos_tb);
    BoundingBox_sync(box);
}


/*gpufun*/
int8_t BoundingBox_overlaps(BoundingBox b1, BoundingBox b2){ // double overlaps[8]){
    // v1-v4 are the four vertices of the first box in counterclockwise order
    // w1-w4 are the four vertices of the second box in counterclockwise order
    // e1-e2 are the two axes of the first box
    // f1-f2 are the two axes of the second box
    // 1 : overlap, 0: no overlap
    // overlaps = [min_e1, max_e1, min_e2, max_e2, min_f1, max_f1, min_f2, max_f2]

    double rC1     = BoundingBox_get_rC(b1);
    double l1      = BoundingBox_get_l(b1);
    double w1      = BoundingBox_get_w(b1);
    double proj_l1 = BoundingBox_get_proj_l(b1);
    double proj_w1 = BoundingBox_get_proj_w(b1);
    double sin_tb1 = BoundingBox_get_sin_tb(b1);
    double cos_tb1 = BoundingBox_get_cos_tb(b1);
    double sin_tC1 = BoundingBox_get_sin_tC(b1);
    double cos_tC1 = BoundingBox_get_cos_tC(b1);

    double rC2     = BoundingBox_get_rC(b2);
    double l2      = BoundingBox_get_l(b2);
    double w2      = BoundingBox_get_w(b2);
    double proj_l2 = BoundingBox_get_proj_l(b2);
    double proj_w2 = BoundingBox_get_proj_w(b2);
    double sin_tb2 = BoundingBox_get_sin_tb(b2);
    double cos_tb2 = BoundingBox_get_cos_tb(b2);
    double sin_tC2 = BoundingBox_get_sin_tC(b2);
    double cos_tC2 = BoundingBox_get_cos_tC(b2);

    double projs[4];
    double cos_tb1_tb2 = cos_tb1 * cos_tb2 + sin_tb1 * sin_tb2;  // cos (tb1 - tb2)
    double sin_tb1_tb2 = sin_tb1 * cos_tb2 - cos_tb1 * sin_tb2;  // sin (tb1 - tb2)
    // The length of the projection of vector W on E is given by: |W| cos(theta) where theta is the angle between W and E
    // projection of vertices of box 2 on the length axis of box 1 (e1)
    projs[0] = rC2 * (cos_tb1 * cos_tC2 + sin_tb1 * sin_tC2);  // first vertex w1:  |w1| cos (tb1 - tc2)
    projs[1] = projs[0] + l2 * cos_tb1_tb2; // second vertex w2 = w1 + f1
    projs[2] = projs[1] + w2 * sin_tb1_tb2; // third vertex w3 = w1 + f1 + f2
    projs[3] = projs[0] + w2 * sin_tb1_tb2; // fourth vertex w4 = w1 + f2
    sort_array_of_4_double(projs);
    if (!INTERVALS_OVERLAP(proj_l1, proj_l1 + l1, projs[0], projs[3])){
        return 0;
    }
    // length of projection of vertices of box 2 on the width axis of box 1 (e2)
    projs[0] = rC2 * (cos_tb1 * sin_tC2 - sin_tb1 * cos_tC2);  // first vertex w1:  |w1| cos(tb1 + pi/2 - tc2)
    projs[1] = projs[0] - l2 * sin_tb1_tb2; // second vertex w2 = w1 + f1 TODO: check sign
    projs[2] = projs[1] + w2 * cos_tb1_tb2; // third vertex w3 = w1 + f1 + f2
    projs[3] = projs[0] + w2 * cos_tb1_tb2; // fourth vertex w4 = w1 + f2
    sort_array_of_4_double(projs);
    if (!INTERVALS_OVERLAP(proj_w1, proj_w1 + w1, projs[0], projs[3])){
        return 0;
    }
    // length of projection of vertices of box 1 on the length axis of box 2 (f1)
    sin_tb1_tb2 = - sin_tb1_tb2; // due to sin being asymmetric
    projs[0] = rC1 * (cos_tb2 * cos_tC1 + sin_tb2 * sin_tC1);  // first vertex v1:  |v1| cos (tb1 - tc2)
    projs[1] = projs[0] + l1 * cos_tb1_tb2; // second vertex v2 = v1 + e1
    projs[2] = projs[1] + w1 * sin_tb1_tb2; // third vertex v3 = v1 + e1 + e2
    projs[3] = projs[0] + w1 * sin_tb1_tb2; // fourth vertex v4 = v1 + e2
    sort_array_of_4_double(projs);
    if (!INTERVALS_OVERLAP(proj_l2, proj_l2 + l2, projs[0], projs[3])){
        return 0;
    }
    // length of projection of vertices of box 1 on the width axis of box 2 (f2)
    projs[0] = rC1 * (cos_tb2 * sin_tC1 - sin_tb2 * cos_tC1);  // first vertex v1:  |v1| cos(tb1 + pi/2 - tc2)
    projs[1] = projs[0] - l1 * sin_tb1_tb2; // second vertex v2 = v1 + e1
    projs[2] = projs[1] + w1 * cos_tb1_tb2; // third vertex v3 = v1 + e1 + e2
    projs[3] = projs[0] + w1 * cos_tb1_tb2; // fourth vertex v4 = v1 + e2
    sort_array_of_4_double(projs);
    if (!INTERVALS_OVERLAP(proj_w2, proj_w2 + w2, projs[0], projs[3])){
        return 0;
    }
    return 1;
}

#endif /* XCOLL_GEOM_BOUNDING_BOX_H */