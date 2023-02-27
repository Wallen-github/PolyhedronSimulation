# PolyhedronSimulation
This Project uses LMGC90 to simulate Polyhendron's behavior under gravitation. The LMGC90 is one of the discrete element methods which can handle detailed contact dynamics.

There are a few papers about the code and the method (by [F. Dubois](https://scholar.google.com/citations?user=boV9fugAAAAJ&hl=en&oi=ao) and others) that you can find for yourself.  However, this is a [paper](https://www.sciencedirect.com/science/article/pii/S0019103521001238) we wrote about the implementation of self-gravity.  The code also has a [wiki](https://git-xen.lmgc.univ-montp2.fr/lmgc90/lmgc90_user/-/wikis/home) page and a very comprehensive [documentation page](http://www.lmgc.univ-montp2.fr/~mozul/LMGC90_USER/UserDoc/docs_2019/#) that you will find very useful once you start using it.

## General Structure

- Preprocessor: Construct a numerical model
- Computation
- Visualize and analyze results

## Preprocessor

Functions in preprocessing are categorized into three domains: granular, masonry, and mesh manipulation. Any functionality filed under any category can be used or combined in any other context.

## wiki of pykdgrav
https://github.com/mikegrudic/pytreegrav/issues/1

# version20230207

bug: this code could change the raw data in gen_sample.py and Computation.py.

tries: I tried to input my raw data before 'poly1 = rigidPolyhedron()' in line 51 of gen_sample.py. I also tried to input the raw data in BODIES.DAT in DATBOX directory.

The raw data and configuration which I want to create are shown in below

<table>
    <tr>
        <td ><center><img src="./version20230207/ResultPic/Rawdata.png" > </center></td>
        <td ><center><img src="./version20230207/ResultPic/RawConfig.jpeg" ></center></td>
    </tr>
</table>

But, the results got from LMGC90 code are
    
<table>
    <tr>
        <td ><center><img src="./version20230207/ResultPic/Pic1.png" > </center></td>
        <td ><center><img src="./version20230207/ResultPic/Pic2.png" ></center></td>
    </tr>
</table>

Comparing these results and raw data, we can find two considerable differences. 1. The shape of a concave polyhedron is revised. 2. The positions of two polyhedrons are changed too. 

In 'rigidContactor3D.py', some lines from 668 change the positions of vertices, but it can't explain the different values in BODIES.DAT. These lines compute the barycenter and inertial momentum matrix, and vertices positions are updated w.r.t the barycenter.

# version20230214

In this version, we can create any polyhedron we want from 'gen_sample.py', whether convex or concave. But the concave polyhedron will meet errors in 'Computation.py', such as 'Error DiscreteGeometry::build_HE_Hdl: Humm contour not closed impossible. Error: impossible to create the HE structure'

<table>
    <tr>
        <td ><center><img src="./version20230212/ResultPic/genPic.png" > </center></td>
        <td ><center><img src="./version20230212/ResultPic/genPic2.png" ></center></td>
    </tr>
</table>

The left is the concave polyhedron, and the right one is the convex polyhedron. For the convex configuration, an animation can be generated, shown in following

https://user-images.githubusercontent.com/38872598/221667287-256bf48b-62ab-4483-93cd-5e270db99a5a.mp4
