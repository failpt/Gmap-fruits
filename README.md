# Gmap-fruits

Running the code requires [CMake](https://cmake.org/download/) 3.20 or newer, a C++20 compiler, and [CGAL 6](https://doc.cgal.org/latest/Manual/installation.html) with [Qt6](https://doc.qt.io/qt-6/get-and-install-qt.html).

## Usage

Build the project with:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

Run with:

```bash
./build/Gmap_fruits [fruit] [step] [cut]
```
*Replace `./build/Gmap_fruits` with `build\Release\Gmap_fruits.exe` on Windows.*

The project currently has grammars for three fruits: tomato, cherry and lime. The fruit is left uncut by default, and adding a `cut` argument cuts off its upper half. 
If unspecified, the last step of the grammar is visualized, which with the current settings is 11, 21 and 7 for tomato, lime and cherry respectively (changing precision of the first two would yield a different number of steps).

## Credits

The methodology is taken from:

> Evans Bohl, Olivier Terraz and Djamchid Ghazanfarpour.
> 
> **[Modeling fruits and their internal structure using parametric 3Gmap L-systems.](https://link.springer.com/article/10.1007/s00371-015-1108-9)**
> 
> *The Visual Computer* Volume 31, pages 819–829 (2015).

Grammar 1 (tomato) and Grammar 2 (cherry) are transcribed from the paper, with the order of gluing and adding layers being slightly modified in Grammar 2 for stylistic purposes.

CGAL provides the generalized maps, linear cell complex and viewer.

## License

Apache License 2.0, see [LICENSE](LICENSE).
