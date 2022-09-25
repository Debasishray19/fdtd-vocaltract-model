# A three-dimensional (3D) articulatory model
:mailbox: <b>Summary: </b>
<br>The implemented 3D articulatory model synthesizes speech acoustic for static vocal tract shapes. The model couples a lumped-element vocal fold model [1] *(i.e., two-mass model)* with the 3D vocal tract model to generate synthetic audio output for various vocal tract shapes. During the simulation, the model uses the 1D area functions dataset [2] (i.e., change in vocal tract cross-sectional shape from glottis to lips) to generate the 3D vocal tract contour. The synthesizer then uses the finite-difference time-domain numerical scheme [3][4] to discretize acoustic components (pressure and velocity) on a staggered othrogonal grid and computes the acoustic pressure wave propagation. The implemented model is motivated from Takemoto's 3D vocal tract model [3] and "Aerophones In Flatland" [5] research articles. It also faciliates the transfer function analysis of static vocal tracts having different cross-sectional shapes (i.e., circular, elliptical and square).

<img src="img/rotating_tract.gif" width="500" height="400">
<img src="img/3d_geometries.JPG" width="500" height="300">

:video_game: <b>Start the simulator:</b>
<br> - To simulate the articulatory model, use the main.m file in the MATLAB environemnt.
<br> - For transfer function analysis, use the plotFrequencyPhase.m file.
<br> - To geneate audio samples, use the generateAudio.m file.

:books: <b>References:</b>
<br>If you use the code for your research work, please cite the following papers -
<br>[1] <a href ="https://www.isca-speech.org/archive/interspeech_2022/mohapatra22_interspeech.html">"Three-dimensional finite-difference time-domain acoustic analysis of simplified vocal tract shapes"</a>  by Mohapatra et al.
```
@inproceedings{mohapatra22_interspeech,
  author={Debasish Mohapatra and Mario Fleischer and Victor Zappi and Peter Birkholz and Sidney Fels},
  title={{Three-dimensional finite-difference time-domain acoustic analysis of simplified vocal tract shapes}},
  year=2022,
  booktitle={Proc. Interspeech 2022},
  pages={764--768},
  doi={10.21437/Interspeech.2022-10649}
}
```
<br>[2] <a href ="https://www.isca-speech.org/archive/interspeech_2019/mohapatra19_interspeech.html">"An Extended Two-Dimensional Vocal Tract Model for Fast Acoustic Simulation of Single-Axis Symmetric Three-Dimensional Tubes"</a>  by Mohapatra et al.
```
@inproceedings{mohapatra19_interspeech,
  author={Debasish Ray Mohapatra and Victor Zappi and Sidney Fels},
  title={{An Extended Two-Dimensional Vocal Tract Model for Fast Acoustic Simulation of Single-Axis Symmetric Three-Dimensional Tubes}},
  year=2019,
  booktitle={Proc. Interspeech 2019},
  pages={3760--3764},
  doi={10.21437/Interspeech.2019-1764}
}
```

:golf: <b>Future work: </b>
<br> - Implementation of bent vocal tract geometries having simplified cross-sections.
<br> - Studying the effect of side branches/cavities.
<br> - Simulation of realistic vocal tract models.
<br> - Implementation of mouth radiation.
<br> - Implementation of an accelerated 3D FDTD vocal tract model using modern GPUs [CUDA/OpenCL environment]

:warning: <b>Note:</b>
<br> For bugs/suggestion/collaboration, please contact: debasishray@ece.ubc.ca
