# MEPP2 installation report for OSX 10.14 (Mojave)
````
$ defaults read loginwindow SystemVersionStampAsString
   ---> 10.14.4
````

### Installation of the core and filters
````
$ brew install cgal               # 4.14
   ---> pulls boost                 1.69.0
              eigen                 3.3.7
$ brew install open-mesh          # 7.1
$ brew install cmake              # 3.14.4
$
$ git clone https://github.com/MEPP-team/MEPP2.git
$ cd MEPP2 && mkdir Build && cd Build
$ cmake -DBUILD_USE_GUI=OFF -DCMAKE_BUILD_TYPE=Release ..
$ make
$ ctest
````

### Building of the GUI
````
$ brew install qt                 # 5.12.3
$ brew install open-scene-graph   # 3.6.3
$ brew cask info xquartz          # 2.7.11
$
$ ccmake -DBUILD_USE_GUI=ON -DCMAKE_BUILD_TYPE=Release -DBUILD_USE_QT5=ON -DQT5_DIR=/usr/local/Cellar/qt/5.12.3 ..
$ make
$ ctest
````

### Building the documentation
````
$ brew install doxygen            # 1.8.15
$ brew install graphviz           # 2.40.1
$ ccmake -DBUILD_USE_GUI=OFF -DBUILD_DOCUMENTATION=ON -DCMAKE_BUILD_TYPE=Release ..
$ make
````
