// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #


typedef struct BoundingBox_s {
    double rC;        // length of position vector to first vertex
    double sin_tC;    // angle of position vector to first vertex, [radians]
    double cos_tC;
    double l;         // length of the box
    double w;         // width of the box
    double sin_tb;    // orientation of the box (angle of length wrt horizontal), [radians]
    double cos_tb;
    double proj_l;    // projection of position vector on length: rC * (cos_t*cos_tC + sin_t*sin_tC)
    double proj_w;    // projection of position vector on width:  rC * (cos_t*sin_tC - sin_t*cos_tC)
} BoundingBox_s;
typedef BoundingBox_s* BoundingBox;


/*gpufun*/
void BoundingBox_sync(BoundingBox box){
    // projection of position vector on length: rC * (cos_t*cos_tC + sin_t*sin_tC)
    box->proj_l = box->rC * (box->cos_tb*box->cos_tC + box->sin_tb*box->sin_tC);
    // projection of position vector on width:  rC * (cos_t*sin_tC - sin_t*cos_tC)
    box->proj_w = box->rC * (box->cos_tb*box->sin_tC - box->sin_tb*box->cos_tC);
}

/*gpufun*/
int8_t BoundingBox_overlaps(BoundingBox b1, BoundingBox b2){ // double overlaps[8]){
    // v1-v4 are the four vertices of the first box in counterclockwise order
    // w1-w4 are the four vertices of the second box in counterclockwise order
    // e1-e2 are the two axes of the first box
    // f1-f2 are the two axes of the second box
    // 1 : overlap, 0: no overlap
    // overlaps = [min_e1, max_e1, min_e2, max_e2, min_f1, max_f1, min_f2, max_f2]

    double projs[4];
    double cos_tb1_tb2 = b1->cos_tb * b2->cos_tb + b1->sin_tb * b2->sin_tb;  // cos (tb1 - tb2)
    double sin_tb1_tb2 = b1->sin_tb * b2->cos_tb - b1->cos_tb * b2->sin_tb;  // sin (tb1 - tb2)
    // The length of the projection of vector W on E is given by: |W| cos(theta) where theta is the angle between W and E
    // projection of vertices of box 2 on the length axis of box 1 (e1)
    projs[0] = b2->rC * (b1->cos_tb * b2->cos_tC + b1->sin_tb * b2->sin_tC);  // first vertex w1:  |w1| cos (tb1 - tc2)
    projs[1] = projs[0] + b2->l * cos_tb1_tb2; // second vertex w2 = w1 + f1
    projs[2] = projs[1] + b2->w * sin_tb1_tb2; // third vertex w3 = w1 + f1 + f2
    projs[3] = projs[0] + b2->w * sin_tb1_tb2; // fourth vertex w4 = w1 + f2
    sort_array_of_4_double(projs);
    if (!INTERVALS_OVERLAP(b1->proj_l, b1->proj_l + b1->l, projs[0], projs[3])){
        return 0;
    }
    // length of projection of vertices of box 2 on the width axis of box 1 (e2)
    projs[0] = b2->rC * (b1->cos_tb * b2->sin_tC - b1->sin_tb * b2->cos_tC);  // first vertex w1:  |w1| cos(tb1 + pi/2 - tc2)
    projs[1] = projs[0] - b2->l * sin_tb1_tb2; // second vertex w2 = w1 + f1 TODO: check sign
    projs[2] = projs[1] + b2->w * cos_tb1_tb2; // third vertex w3 = w1 + f1 + f2
    projs[3] = projs[0] + b2->w * cos_tb1_tb2; // fourth vertex w4 = w1 + f2
    sort_array_of_4_double(projs);
    if (!INTERVALS_OVERLAP(b1->proj_w, b1->proj_w + b1->w, projs[0], projs[3])){
        return 0;
    }
    // length of projection of vertices of box 1 on the length axis of box 2 (f1)
    sin_tb1_tb2 = - sin_tb1_tb2; // due to sin being asymmetric
    projs[0] = b1->rC * (b2->cos_tb * b1->cos_tC + b2->sin_tb * b1->sin_tC);  // first vertex v1:  |v1| cos (tb1 - tc2)
    projs[1] = projs[0] + b1->l * cos_tb1_tb2; // second vertex v2 = v1 + e1
    projs[2] = projs[1] + b1->w * sin_tb1_tb2; // third vertex v3 = v1 + e1 + e2
    projs[3] = projs[0] + b1->w * sin_tb1_tb2; // fourth vertex v4 = v1 + e2
    sort_array_of_4_double(projs);
    if (!INTERVALS_OVERLAP(b2->proj_l, b2->proj_l + b2->l, projs[0], projs[3])){
        return 0;
    }
    // length of projection of vertices of box 1 on the width axis of box 2 (f2)
    projs[0] = b1->rC * (b2->cos_tb * b1->sin_tC - b2->sin_tb * b1->cos_tC);  // first vertex v1:  |v1| cos(tb1 + pi/2 - tc2)
    projs[1] = projs[0] - b1->l * sin_tb1_tb2; // second vertex v2 = v1 + e1
    projs[2] = projs[1] + b1->w * cos_tb1_tb2; // third vertex v3 = v1 + e1 + e2
    projs[3] = projs[0] + b1->w * cos_tb1_tb2; // fourth vertex v4 = v1 + e2
    sort_array_of_4_double(projs);
    if (!INTERVALS_OVERLAP(b2->proj_w, b2->proj_w + b2->w, projs[0], projs[3])){
        return 0;
    }
    return 1;
}


// Expose functions to Xobject test interface
// ------------------------------------------


/*gpufun*/
void BoundingBoxTest_to_BoundingBox(BoundingBoxTest box_in, BoundingBox box_out){
    box_out->cos_tb = BoundingBoxTest_get_cos_tb(box_in);
    box_out->sin_tb = BoundingBoxTest_get_sin_tb(box_in);
    box_out->l = BoundingBoxTest_get_l(box_in);
    box_out->w = BoundingBoxTest_get_w(box_in);
    box_out->rC = BoundingBoxTest_get_rC(box_in);
    box_out->cos_tC = BoundingBoxTest_get_cos_tC(box_in);
    box_out->sin_tC = BoundingBoxTest_get_sin_tC(box_in);
    box_out->proj_l = BoundingBoxTest_get_proj_l(box_in);
    box_out->proj_w = BoundingBoxTest_get_proj_w(box_in);
}

/*gpufun*/
void BoundingBox_to_BoundingBoxTest(BoundingBox box_in, BoundingBoxTest box_out){
    BoundingBoxTest_set_cos_tb(box_out, box_in->cos_tb);
    BoundingBoxTest_set_sin_tb(box_out, box_in->sin_tb);
    BoundingBoxTest_set_l(box_out, box_in->l);
    BoundingBoxTest_set_w(box_out, box_in->w);
    BoundingBoxTest_set_rC(box_out, box_in->rC);
    BoundingBoxTest_set_cos_tC(box_out, box_in->cos_tC);
    BoundingBoxTest_set_sin_tC(box_out, box_in->sin_tC);
    BoundingBoxTest_set_proj_l(box_out, box_in->proj_l);
    BoundingBoxTest_set_proj_w(box_out, box_in->proj_w);
}

/*gpufun*/
void BoundingBoxTest_set_params(BoundingBoxTest box, double rC, double sin_tC, double cos_tC,
                                double l, double w, double sin_tb, double cos_tb){
    BoundingBox_s box1;
    box1.rC = rC;
    box1.sin_tC = sin_tC;
    box1.cos_tC = cos_tC;
    box1.l = l;
    box1.w = w;
    box1.sin_tb = sin_tb;
    box1.cos_tb = cos_tb;
    BoundingBox_sync(&box1);
    BoundingBox_to_BoundingBoxTest(&box1, box);
}

/*gpufun*/
int8_t BoundingBoxTest_overlaps(BoundingBoxTest b1, BoundingBoxTest b2){
    BoundingBox_s box1, box2;
    BoundingBoxTest_to_BoundingBox(b1, &box1);
    BoundingBoxTest_to_BoundingBox(b2, &box2);
    return BoundingBox_overlaps(&box1, &box2);
}
