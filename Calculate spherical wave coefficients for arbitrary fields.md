## Error definition
We need to find the constant vector C (Or we can use k), which minimizes the residual e defined as:
<center> 
    $ er= \Delta Er \cdot \Delta Er^{\dagger} $
 </center>
 Where the error vector,
 <center> 
    $ 
    Er(k)=
\begin{bmatrix}
    E_{des}-E^k_{sph} \\
    H_{des}-H^k_{sph}\\
\end{bmatrix} $
 </center>
 Where subscripts and refer to waves reconstructed from spherical waves and desired waves, respectively.
 
 ## Plane wave notation and polarizations 
 
### Plane wave:
 <center> 
    $ \vec{E}(\vec{r})=\vec{E_{0}}e^{i\vec{k_e}\cdot\vec{r}} $
 </center>
 and 
  <center> 
    $ \vec{H}(\vec{r})=\sqrt{\frac{\varepsilon}{\mu}}\vec{e_k} \times \vec{E}(r)  $
 </center>
 ### Basis constructions: 
From the lab basis {$e_x$,$e_y$,$e_z$}, we can constract a set of {$e_k$,$e_a$,$e_b$} by two angles (a$\phi$, b$\theta$) for the following, 

 <center> 
    $ e_k=sin(b)cos(a)\vec{e_x} +sin(b)sin(a)\vec{e_y}+cos(b)\vec{e_z}\\
    e_b=cos(b)sin(a)\vec{e_x} +cos(b)sin(a)\vec{e_y}-sin(b)\vec{e_z}\\
     e_a=-sin(a)\vec{e_x} +cos(a)\vec{e_y}\\$
    
 </center>

By considering a plane wave with arbitrary isotropic media, direction and polarizations:
 where the E-field can be expressed as the following: 
 <center> 
    $ \vec{E}(\vec{r})=\vec{E_{0}}e^{i\vec{k_e}\cdot\vec{r}} $
 </center>
 Where 
 <center> 
    $ \vec{E}(\vec{r})=E_{a}\vec{e_a}+E_{b}\vec{e_b} $
 </center>
 Then, respects to the basis, any polarization can be expressed as two angles $\psi$ and $\chi$ in a vibration ellipse[1].$\\[0.5cm]$  
 
 ![Screen%20Shot%202020-03-30%20at%2012.49.13%20PM.png](attachment:Screen%20Shot%202020-03-30%20at%2012.49.13%20PM.png)
  <center> 
    $ E_{a}=E_{0}(cos\chi sin\psi-isin\chi cos\psi) $
 </center>
 and 
   <center> 
    $ E_{b}=E_{0}(cos\chi sin\psi+isin\chi cos\psi) $
 </center>
 #### Stokes parameters and degree of polarization
 For quasi-monochromatic light, the amplitude of the electric field fluctuate in time, which reduced the degree of polarization. 
  <center> 
    $ I_{e}=E_0^2 $
 </center>
 $\\[0.5cm]$
   <center> 
    $ Q_{e}=-E_0^2cos2\chi cos2\psi $
 </center>
 $\\[0.5cm]$
  <center> 
    $ U_{e}=-E_0^2cos2\chi sin2\psi $
 </center>
 $\\[0.5cm]$
  <center> 
    $ V_{e}=-E_0^2sin2\chi $
 </center>
  $\\[0.5cm]$
  <center> 
    $ P=\frac{\sqrt{Q_e^2+U_e^2+V_e^2}}{I_e} $
 </center>
 
 ### Associated Legendre functions and angular functions 
 The angular function can relate to the associated Legendre functions by the following
   <center> 
    $ \pi^m_n(\theta)=\sqrt{\frac{2n+1}{2}\frac{(n-m)!}{(n+m)!}}\frac{P^m_n(cos(\theta))}{sin(\theta)} $
 </center>
and 
   <center> 
    $ \tau^m_n(\theta)=\sqrt{\frac{2n+1}{2}\frac{(n-m)!}{(n+m)!}}\frac{d}{d\theta}P^m_n(cos(\theta)) $
 </center>
 Sphere harmonics
 
 ![NumberedEquation6.gif](attachment:NumberedEquation6.gif)
 
 ### vector spherical waves expansion
<center> 
    $ \vec{E_e(r)}=\sum_{n=1}^{\infty} \sum_{m=-n}^{n} a_{mn}\vec{M_{mn}^1}(kr)+ b_{mn}\vec{N_{mn}^1}(kr)$
 </center>
 Where, with two angles (a$\phi$, b$\theta$), we have: 
    <center> 
    $ a_{mn}=-\frac{4i^n}{\sqrt{2n(n+1)}} \vec{e_{polar}} \cdot \{im\pi^{|m|}_n(b)\vec{e_b}+\tau^{|m|}_n(b)\vec{e_a}\}e^{-ima}   $
 </center>
 and 
     <center> 
    $ b_{mn}=-\frac{4i^{n+1}}{\sqrt{2n(n+1)}} \vec{e_{polar}} \cdot \{-im\tau^{|m|}_n(b)\vec{e_b}+\pi^{|m|}_n(b)\vec{e_a}\}e^{-ima}   $
 </center>
 The harmonics give: 
  <center> 
    $ \vec{m}_{mn}=\frac{1}{\sqrt{2n(n+1)}}\{-\tau^{|m|}_n(b)\vec{e_b}+im\pi^{|m|}_n(b)\vec{e_a}\}e^{-ima}   $
 </center>
 and 
   <center> 
    $ \vec{n}_{mn}=\frac{1}{\sqrt{2n(n+1)}}\{\tau^{|m|}_n(b)\vec{e_b}+im\pi^{|m|}_n(b)\vec{e_a}\}e^{-ima}   $
 </center>
 Thus, 
     <center> 
    $ a_{mn}=4i^n \vec{e_{polar}} \cdot \vec{m}_{mn}^*  $
 </center>
 and 
     <center> 
    $ b_{mn}=-4i^{n+1} \vec{e_{polar}} \cdot \vec{n}_{mn}^*  $
 </center>
 The vector wave functions are:
      <center> 
    $ M^{1,3}_{mn}=z^{1,3}_{n} \vec{m}_{mn}  $
 </center>
 and 
     <center> 
    $  N^{1,3}_{mn}=\sqrt{\frac{n(n+1)}{2}}\frac{z^{1,3}_{n}}{kr}Y_{mn}\vec{e_r}+(\frac{\frac{d}{d(kr)}krz^{1,3}_{n}}{kr})\vec{n}_{mn} $
 </center>
 where,
<center> 
    $ z^{1}_{n}=j_{n}\\
     z^3_n=h^{(1)}_n   $
    
 </center>


```python

```
