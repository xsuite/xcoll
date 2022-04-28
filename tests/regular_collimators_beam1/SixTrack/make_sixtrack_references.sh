while read coll
do
    echo "Collimator "$coll
    mkdir temp
    # Add dump markers
    sed 's/%DUMP%/'${coll}'/g' fort.3_template > temp/fort.3
    # Add GO statement in front of collimator in the structure list
    # (replace 5th occurence because of the 1) extraction markers and 2) single element list:
    # SINGLE ELEMENTS ( .. mt_coll_mken coll mt_coll_mkex .. ) STRUCTURE ( .. mt_coll_mken GO coll .. )
    #awk '/'${coll}'/{c+=1}{if(c==5){sub("'${coll}'","GO '${coll}'",$0)};print}' fort.2_template > temp/fort.2
    # (replace second occurence of marker)
    sub='mt_'${coll}'_mken'
    awk '/'${sub}'/{c+=1}{if(c==2){sub("'${sub}'","GO '${sub}'",$0)};print}' fort.2_template > temp/fort.2
    cd temp
        ln -fns ../CollDB-Testing.dat .
        ln -fns ../fort3.limi .
        ln -fns ../fort.8 .
        ln -fns ../initial.dat .
        ln -fns ../sixtrack .
        ./sixtrack | tee ../output/six_${coll}.out 2>&1
        mv dump* ../output/
    cd ..
#    grep -i tilt six_${coll}.out > tilt_${coll}.out
    rm -r temp
done <collimators

