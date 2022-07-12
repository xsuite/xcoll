#!/bin/bash

# python _make_json_from_initial.py

for beam in B1 B2
do
    cd SixTrack_$beam
    mkdir -p output
    while read coll
    do
        echo "Collimator "$coll
        mkdir temp
        # Add dump markers
        sed 's/%COLL%/'${coll}'/g' fort.3_template > temp/fort.3
        # (replace second occurence of marker)
        sub='mt_'${coll}'_mken'
        awk '/'${sub}'/{c+=1}{if(c==2){sub("'${sub}'","GO '${sub}'",$0)};print}' fort.2 > temp/fort.2
        cd temp
            ln -fns ../CollDB-Testing.dat .
            ln -fns ../fort3.limi .
            ln -fns ../fort.8 .
            ln -fns ../initial.dat .
            ln -fns ../sixtrack .
            ./sixtrack | tee ../output/six_${coll}.out 2>&1
            mv dump* ../output/
            json=../../Collimators/${coll}.json
            grep -i $coll -A20 colltrack.out > tempcoll
            material=$( head -2 tempcoll | tail -1 | awk '{print $2;}' )
            length=$( head -3 tempcoll | tail -1 | awk '{print $3;}' )
            angle=$( head -4 tempcoll | tail -1 | awk '{print $3;}' )
            offset=$( head -5 tempcoll | tail -1 | awk '{print $3;}' )
            jaw=$( head -17 tempcoll | tail -1 | awk '{print $4;}' )
            gap=$( grep -i 'NSIG_FAM CRY7' CollDB-Testing.dat | awk '{print $3;}' )
            cryinfo=$( grep -i 'CRYSTAL  '$coll CollDB-Testing.dat )
            bend=$( echo $cryinfo | awk '{print $3;}')
            xdim=$( echo $cryinfo | awk '{print $4;}')
            ydim=$( echo $cryinfo | awk '{print $5;}')
            thick=$( echo $cryinfo | awk '{print $6;}')
            crytilt=$( echo $cryinfo | awk '{print $7;}')
            miscut=$( echo $cryinfo | awk '{print $8;}')
            orient=$( echo $cryinfo | awk '{print $9;}')
            if [ "${orient}" == "" ]; then orient=0; fi
            echo '{"__class__": "K2Crystal", "n_alloc": 50000, "random_generator_seed": 6574654,' > $json
            echo '    "active_length": '${length}', "inactive_front": 0, "inactive_back": 0, "angle": '${angle}',' >> $json
            echo '    "jaw_F_L": '${jaw}', "jaw_F_R": -'${jaw}', "jaw_B_L": '${jaw}', "jaw_B_R": -'${jaw}', ' >> $json
            grep -i ' '$coll linopt_coupled.dat | tail -1 | awk '{print "    \"dx\": "$30", \"dpx\": "$31", \"dy\": "$32", \"dpy\": "$33",";}' >> $json
            echo '    "offset": '${offset}', "tilt": [0.0, 0.0], "is_active": 1, "record_impacts": 0,' >> $json
            echo '    "bend": '${bend}', "xdim": '${xdim}', "ydim": '${ydim}', "thick": '${thick}', "gap": '${gap}', ' >> $json
            echo '    "crytilt": '${crytilt}', "miscut": '${miscut}', "orient": '${orient}', "emittance": [3.5e-6,3.5e-6],' >> $json
            echo '    "onesided": 0, "material": "'${material}'"}' >> $json
        cd ..
   #    grep -i tilt six_${coll}.out > tilt_${coll}.out
        rm -r temp
    done <../collimators_$beam
    cd ..
done
python _correct_json_collimators.py
python _make_json_references.py
