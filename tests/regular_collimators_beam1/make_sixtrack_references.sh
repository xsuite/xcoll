    cp fort.3_template fort.3
    sed -i 's/%DUMP%/${coll}/g' fort.3
    ./sixtrack > six_${coll}.out


