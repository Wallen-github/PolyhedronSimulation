# PolyhedronSimulation
This Project uses LMGC90 to simulate Polyhendron's behavior under gravitation. The LMGC90 is one of discrete element methonds which can handle detailed contact dynamics.

## version20230207

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

## version20230214

In this version, we can create any polyhedron we want from 'gen_sample.py', whether convex or concave. But the concave polyhedron will meet errors in 'Computation.py', such as 'Error DiscreteGeometry::build_HE_Hdl: Humm contour not closed impossible. Error: impossible to create the HE structure'

<table>
    <tr>
        <td ><center><img src="./version20230212/ResultPic/genPic.png" > </center></td>
        <td ><center><img src="./version20230212/ResultPic/genPic2.png" ></center></td>
    </tr>
</table>

The left is the concave polyhedron, and the right one is the convex polyhedron. For the convex configuration, an animation can be generated, shown in following

https://user-images.githubusercontent.com/38872598/221667287-256bf48b-62ab-4483-93cd-5e270db99a5a.mp4
