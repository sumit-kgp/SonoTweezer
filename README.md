# SonoTweezer

---page under development---

# Introduction
Our device to focus ultrasound waves to trap particles underwater based on the twin trap principle, details can be found in the article: https://research.utwente.nl/en/publications/sonotweezer-an-acoustically-powered-end-effector-for-underwater-m
![image](https://github.com/user-attachments/assets/837599d8-3706-4596-bab6-e9326519e023)

We present a matlab implementation of the original concept on airborne acoustic levitation kit, Ultraino, developed by Asier Marzo (https://github.com/asiermarzo/Ultraino). The given code allows simulation of acoustic fields generated by an array of symmetrically-placed transducers in air (Ultraino) and water (SonoTweezer) based on piston-source model. This model allows to simulate different mechanisms based on vortex and twin trap configurations. The device in the published article, SonoTweezer, uses six MHz-range immersible ultrasonic transducers to generate a twin trap for manipulating sub-millimeter particles.

# Implementation

Please use MATLAB v2019b onwards to use this code. 
![image](https://github.com/user-attachments/assets/6210154c-f65a-477e-ac6d-f0eda30e2480)

The code is divided into main programs:
- Ultraino_original: script to generate acoustic fields in air using both twin trap and vortex trap configurations with multiple rings of distributed transducers (parameters defined in the script)
- SonoTweezer_original: script to generate acoustic fields in water using twin trap configuration with a single ring of distributed transducers (parameters based on a model version from Imasonic SAS, can be tailored to preference)
- SonoTweezer_parametric_variation: script to tailor different configurations of acoustic fields in water besides SonoTweeezer configuration
- Directivity: program to plot acoustic pressure distribution around a single transducer 
Some parameters used in this study can be found in the table below:

![image](https://github.com/user-attachments/assets/0957e911-cbd8-4992-b865-51f3766affb6)

These codes use the following different MATLAB functions:
- transducerloc: function that defines distribution of transducers in a space
- complex_pressure: computing pressure at a spatial point (Equation 4 in the article)
- directfunc: computing directivity function of a transducer based on defined properties (wavelength of sound in the medium, density and size of transducer)
- distance_and_angle_xy: allows computation of acoustic fields in XY-plane i.e. lateral plane (shown in fig. below)
- distance_and_angle_yz: allows computation of acoustic fields in YZ-plane i.e. longitudinal plane (shown in fig. below)
- potential field: allows computation of Gor'kov potential based on acoustic fields (Equation 2 in the article). Same function can be used with acoustic fields in all the possible planes.
  
![image](https://github.com/user-attachments/assets/4cd0fa12-6028-4867-a0c0-a323ace6f76f)

# Acknowledgement

This program is an outcome of master's thesis work of Robbert-Jan Fidder under the supervision of Sumit Mohanty and Sarthak Misra at Surgical Robotics Lab, University of Twente. Valuable theoretical insights for modeling acoustic fields was provided by Nathan Blanken at Physics of Fluids group, University of Twente. Other authors who contributed to the experimental realization of the concept are Christoff Heunis (CEO, Flux Robotics), Pedro Matos (CWI) and Mert Kaya (University of Twente) in the Netherlands.

This work is supported by funds from the Netherlands Organization for Scientific Research (Innovational Research Incentives Scheme– VIDI: SAMURAI project # 14855)
