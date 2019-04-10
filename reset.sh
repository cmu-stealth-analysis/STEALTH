for pycfile in *.pyc; do
    if [ "${pycfile}" == "*.pyc" ]; then
        echo "No pyc files to clean."
    else
        rm -v ${pycfile}
    fi
done

git fetch my-stealth && git reset --hard my-stealth/tanmay-devel
cd eventSelection
make
cd ..
cd getMCSystematics
make
cd ..
