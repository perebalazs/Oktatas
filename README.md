# 2D-s statikai és dinamikai feladat megoldása végeselem módszerrel  Julia+gmsh segítségével

* 2D-s *síkfeszültség*, *síkalakváltozás* és *3D-s test*  feladat megoldása.
  * Felületi és térfogati terhelés tetszőleges peremen és térfogaton:
    * 2D-s: vonal menti és felületi terhelés,
    * 3D-s: felületi és térfogati terhelés.
  * Megfogás tetszőleges peremen (pont?, vonal, felület) $x$, $y$ és/vagy $z$ irányban.
  * Tetszőleges elemekkel:
    * 2D: háromszög és/vagy négyszög elemekkel (akár vegyesen is),
    * 3D: tetraéder, hexaéder, ék és/vagy piramis elemekkel (akár vegyesen is).
  * Közelítés fokszáma egytől tizedfokig (nagyon magas fokszámok esetén nem biztos hogy működik - miért?).
  * Végeselemek: csak peremen lévő csomópontok, vagy belső csomópontok is (magasabb fokszámnál).
  * Elmozdulásmező szemléltetése.
  * Feszültségmező szemléltetése (szakadásos vagy simított mező).
  * Dinamikai feladatnál (explicit) animáció készítése bármelyik elmozdulás, sebesség vagy feszültség mezőről.
  * Tetszőleges vonal (spline, körív, stb.) mentén tetszőleges eredmény (elmozulás, sebesség vagy feszültség) ábrázolása grafikonon.
  * Geometria rajzolása és hálózás a *gmsh* programban, mentés *.geo fájlba (hálózó parancsokkal együtt).
* Példák az *examples* könyvtárban.
* Minden fent említett funciót a *FEM.jl* fájl számol (jelenleg 889 sor).
* A használathoz telepíteni kell a [julia](https://julialang.org/) és [gmsh](https://gmsh.info) programokat. Utóbbiból az *SDK*-t.

Megköszönök minden észrevételt és hibajelentést.
