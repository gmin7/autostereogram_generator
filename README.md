## Random-Dot Autostereogram (RDA) Generator

Using the paper from https://www.cs.waikato.ac.nz/~ihw/papers/94-HWT-SI-IHW-SIRDS-paper.pdf to turn my old ray caster into an autostereogram generator.  
Perlin implementation adapted from source code found: http://eastfarthing.com/blog/2015-04-21-noise/  
Can you spot the bunny?

Grayscale:  
<img src="/images/autostereogram_bunny_GS.png"/>  
RGB:  
<img src="/images/autostereogram_bunny_RGB.png"/>  
Depth map(Spoiler!):  
<img src="/images/depth_map_bunny.png"/>  

### Perlin Noise
I wanted to see what it would be like using Perlin Noise to generate the source image prior to generating an RDA. My implementation works heavily off of the JS source code used in the very helpful demo found here: http://eastfarthing.com/blog/2015-04-21-noise/

However, as this version of Perlin Noise is not periodic, I must keep it at a high granularity to minimize the artifacts which arise:  
<img src="/images/autostereogram_bunny_RGBPerlin_spaces4.png"/>  
