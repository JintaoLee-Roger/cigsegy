if [[ "$OSTYPE" == "linux-gnu"* ]]; then
  echo "build in linux"
  echo "In docker"
  # for cpyroot in `ls -d /opt/_internal/cpython-3.*`;                                      
  # do
  #   result=$(echo $cpyroot | grep "3.5")
  #   if [ "$result" = "" ] ; then
  #   cpy=$cpyroot/bin/python3
  #   pypi=$cpyroot/bin/pip3
  #   $pypi install pybind11
  #   $pypi wheel .
  #   fi
  # done
  
  cp -r /share/thridpart/fmt_2_24/* /usr/local/
  rm -rf /share/build
  mkdir /share/build
  rm -rf /share/wheels
  mkdir /share/wheels
  cd /share/build
  
  for cpyroot in `ls -d /opt/_internal/cpython-3.*`;
  do
    version="${cpyroot#/opt/_internal/cpython-}"
    version_min="${version%.*}"
    if [ "$version_min" != "3.5" ] ; then
      cpy=$cpyroot/bin/python3
      cpylib=$cpyroot/lib
      pypi=$cpyroot/bin/pip3
      $pypi install pybind11
      cpybind=${cpyroot}/lib/python${version_min}/site-packages/pybind11
      cmake .. -DCMAKE_INSTALL_PREFIX=/share/dist_${version_min}/ -Dpybind11_ROOT=${cpybind} -DPYTHON_EXECUTABLE=${cpy} -DPYTHON_LIBRARIES=${cpylib}
      make -j8 
      make install
      cd /share/dist_$version_min/python 
      $pypi wheel .
      whl=`ls *.whl`
      prefix="${whl%linux_x86_64.whl}"
      cp $whl /share/wheels/${prefix}manylinux_2_24_x86_64.whl
      # auditwheel repair $whl -w /share/wheelhouse/
      cd /share/build
      rm CMakeCache.txt
    fi
  done
  
  rm -rf /share/build
  rm -rf /share/dist_*
  
  # for dist_path in `ls -d /share/dist_3.*`;
  # do
  #   cd $dist_path/python
  #   pip wheel .
  #   whl=`ls *.whl`
  #   auditwheel repair $whl -w /share/wheelhouse/
  # done
  
  
  rm -rf /share/thridpart/cigsegy/cigsegy*
  rm -rf /share/thridpart/cigsegy/dist
  rm -rf /share/thridpart/cigsegy/src
  rm -rf /share/thridpart/cigsegy/python
  
  cp -r /share/src /share/thridpart/cigsegy/
  cp -r /share/python /share/thridpart/cigsegy/
  cp -r /share/setup.py /share/thridpart/cigsegy/
  
  sed -i "s/fmt_root = ''/fmt_root = '.'/" /share/thridpart/cigsegy/setup.py
  
  cd /share/thridpart/cigsegy/
  python=/opt/_internal/cpython-3.9.15/bin/python3
  ${python} setup.py sdist
  mv /share/thridpart/cigsegy/dist/*.tar.gz /share/wheels/

elif [[ "$OSTYPE" == "darwin"* ]]; then
  echo "build on MacOS"

  rm -rf build
  rm -rf dist*
  mkdir build 
  mkdir wheels
  cd  build 

  for env in py38 py39 base py311;
  do
    source activate $env
    python -m pip install --upgrade pip
    cmake .. -DCMAKE_INSTALL_PREFIX=../dist_$env
    make 
    make install 
    rm CMakeCache.txt
    cd ../dist_${env}/python 
    pip wheel .
    mv *.whl ../../wheels/
    cd ../../build/
  done

  cd ../
  source activate base

  rm -rf ./thridpart/cigsegy/cigsegy*
  rm -rf ./thridpart/cigsegy/dist
  rm -rf ./thridpart/cigsegy/src
  rm -rf ./thridpart/cigsegy/python
  rm ./thridpart/cigsegy/setup.py
  
  cp -r ./src ./thridpart/cigsegy/
  cp -r ./python ./thridpart/cigsegy/
  cp -r ./setup.py ./thridpart/cigsegy/
  
  sed -i "" "s/fmt_root = ''/fmt_root = '.'/" ./thridpart/cigsegy/setup.py

  cd ./thridpart/cigsegy/
  python setup.py sdist
  mv ./dist/*.tar.gz ../../wheels/


else
  echo "Unkown system"
fi

