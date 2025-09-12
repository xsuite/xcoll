// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #


/*gpufun*/
void BoundingBox_set_params(BoundingBox box, double rC, double sin_tC, double cos_tC,
                            double l, double w, double sin_tb, double cos_tb){
    BoundingBox_set_rC(box, rC); // length of position vector to first vertex
    BoundingBox_set_sin_tC(box, sin_tC); // angle of position vector to first vertex
    BoundingBox_set_cos_tC(box, cos_tC);
    BoundingBox_set_l(box, l); // length of the box
    BoundingBox_set_w(box, w); // width of the box
    BoundingBox_set_sin_tb(box, sin_tb);  // orientation of the box (angle of length wrt horizontal)
    BoundingBox_set_cos_tb(box, cos_tb);
    BoundingBox_set_proj_l(box, rC * (cos_tb*cos_tC + sin_tb*sin_tC)); // projection of position vector on length: rC * (cos_t*cos_tC + sin_t*sin_tC)
    BoundingBox_set_proj_w(box, rC * (cos_tb*sin_tC - sin_tb*cos_tC)); // projection of position vector on width:  rC * (cos_t*sin_tC - sin_t*cos_tC)
}

/*gpufun*/
int8_t BoundingBox_overlaps(BoundingBox b1, BoundingBox b2){ // double overlaps[8]){
    // v1-v4 are the four vertices of the first box in counterclockwise order
    // w1-w4 are the four vertices of the second box in counterclockwise order
    // e1-e2 are the two axes of the first box
    // f1-f2 are the two axes of the second box
    // 1 : overlap, 0: no overlap
    // overlaps = [min_e1, max_e1, min_e2, max_e2, min_f1, max_f1, min_f2, max_f2]
    double cos_tb_b1 = BoundingBox_get_cos_tb(b1);
    double sin_tb_b1 = BoundingBox_get_sin_tb(b1);
    double cos_tC_b1 = BoundingBox_get_cos_tC(b1);
    double sin_tC_b1 = BoundingBox_get_sin_tC(b1);
    double l_b1      = BoundingBox_get_l(b1);
    double w_b1      = BoundingBox_get_w(b1);
    double rC_b1     = BoundingBox_get_rC(b1);
    double proj_l_b1 = BoundingBox_get_proj_l(b1);
    double proj_w_b1 = BoundingBox_get_proj_w(b1);

    double cos_tb_b2 = BoundingBox_get_cos_tb(b2);
    double sin_tb_b2 = BoundingBox_get_sin_tb(b2);
    double cos_tC_b2 = BoundingBox_get_cos_tC(b2);
    double sin_tC_b2 = BoundingBox_get_sin_tC(b2);
    double l_b2      = BoundingBox_get_l(b2);
    double w_b2      = BoundingBox_get_w(b2);
    double rC_b2     = BoundingBox_get_rC(b2);
    double proj_l_b2 = BoundingBox_get_proj_l(b2);
    double proj_w_b2 = BoundingBox_get_proj_w(b2);

    double projs[4];
    double cos_tb1_tb2 = cos_tb_b1 * cos_tb_b2 + sin_tb_b1 * sin_tb_b2;  // cos (tb1 - tb2)
    double sin_tb1_tb2 = sin_tb_b1 * cos_tb_b2 - cos_tb_b1 * sin_tb_b2;  // sin (tb1 - tb2)
    // The length of the projection of vector W on E is given by: |W| cos(theta) where theta is the angle between W and E
    // projection of vertices of box 2 on the length axis of box 1 (e1)
    projs[0] = rC_b2 * (cos_tb_b1 * cos_tC_b2 + sin_tb_b1 * sin_tC_b2);  // first vertex w1:  |w1| cos (tb1 - tc2)
    projs[1] = projs[0] + l_b2 * cos_tb1_tb2; // second vertex w2 = w1 + f1
    projs[2] = projs[1] + w_b2 * sin_tb1_tb2; // third vertex w3 = w1 + f1 + f2
    projs[3] = projs[0] + w_b2 * sin_tb1_tb2; // fourth vertex w4 = w1 + f2
    sort_array_of_4_double(projs);
    if (!INTERVALS_OVERLAP(proj_l_b1, proj_l_b1 + l_b1, projs[0], projs[3])){
        return 0;
    }
    // length of projection of vertices of box 2 on the width axis of box 1 (e2)
    projs[0] = rC_b2 * (cos_tb_b1 * sin_tC_b2 - sin_tb_b1 * cos_tC_b2);  // first vertex w1:  |w1| cos(tb1 + pi/2 - tc2)
    projs[1] = projs[0] - l_b2 * sin_tb1_tb2; // second vertex w2 = w1 + f1 TODO: check sign
    projs[2] = projs[1] + w_b2 * cos_tb1_tb2; // third vertex w3 = w1 + f1 + f2
    projs[3] = projs[0] + w_b2 * cos_tb1_tb2; // fourth vertex w4 = w1 + f2
    sort_array_of_4_double(projs);
    if (!INTERVALS_OVERLAP(proj_w_b1, proj_w_b1 + w_b1, projs[0], projs[3])){ 
        return 0;
    }
    // length of projection of vertices of box 1 on the length axis of box 2 (f1)
    sin_tb1_tb2 = - sin_tb1_tb2; // due to sin being asymmetric
    projs[0] = rC_b1 * (cos_tb_b2 * cos_tC_b1 + sin_tb_b2 * sin_tC_b1);  // first vertex v1:  |v1| cos (tb1 - tc2)
    projs[1] = projs[0] + l_b1 * cos_tb1_tb2; // second vertex v2 = v1 + e1
    projs[2] = projs[1] + w_b1 * sin_tb1_tb2; // third vertex v3 = v1 + e1 + e2
    projs[3] = projs[0] + w_b1 * sin_tb1_tb2; // fourth vertex v4 = v1 + e2
    sort_array_of_4_double(projs);
    if (!INTERVALS_OVERLAP(proj_l_b2, proj_l_b2 + l_b2, projs[0], projs[3])){ 
        return 0;
    }
    // length of projection of vertices of box 1 on the width axis of box 2 (f2)
    projs[0] = rC_b1 * (cos_tb_b2 * sin_tC_b1 - sin_tb_b2 * cos_tC_b1);  // first vertex v1:  |v1| cos(tb1 + pi/2 - tc2)
    projs[1] = projs[0] - l_b1 * sin_tb1_tb2; // second vertex v2 = v1 + e1
    projs[2] = projs[1] + w_b1 * cos_tb1_tb2; // third vertex v3 = v1 + e1 + e2
    projs[3] = projs[0] + w_b1 * cos_tb1_tb2; // fourth vertex v4 = v1 + e2 
    sort_array_of_4_double(projs);
    if (!INTERVALS_OVERLAP(proj_w_b2, proj_w_b2 + w_b2, projs[0], projs[3])){ 
        return 0;
    }
    return 1;
}
