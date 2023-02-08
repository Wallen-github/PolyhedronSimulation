# PolyhedronSimulation
This Project uses LMGC90 to simulate Polyhendron's behaviour under gravitation. The LMGC90 is one of discrete element methonds which can handle detailed contact dynamics.

## version20230207
bug: this code could change the raw data in gen_sample.py and Computation.py. 
tries: I tried to input my raw data before 'poly1 = rigidPolyhedron()' in line 51 of gen_sample.py. I also tried to input the raw data in BODIES.DAT in DATBOX directory.

The raw data and configuration which I want to create are shown in below
<table>
    <tr>
        <td ><center><img src="./version20230207/ResultPic/Rawdata.png" > </center></td>
        <td ><center><img src="./ResultPic/RawConfig.jpeg" ></center></td>
    </tr>
</table>
But, the results got from LMGC90 code are
<table>
    <tr>
        <td ><center><img src="./version20230207/ResultPic/Pic1.png" > </center></td>
        <td ><center><img src="./version20230207/ResultPic/Pic2.png" ></center></td>
    </tr>
</table>

In the 'rigidContactor3D.py', there are some lines changing the vertices positions from line 668. These lines compute the barycenter and inertial momentum. Vertices positions are changed w.r.t the barycenter, but it cant explain the values in BODIES.DAT.
