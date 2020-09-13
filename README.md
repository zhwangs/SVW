# Arbitrary polarization plane wave reconstruction in Spherical vector waves basis


## Objective
Decompose an arbitrary polarization plane wave into spherical vector waves (SVW) basis. This script is closely following the book, 'Scattering, Absorption, and Emission of Light by Small Particles', by (M. I. MISHCHENKO, et al). The explicit form of this expansion can be found in MISHCHENKO Appendix C (page 377, C56). Plane wave generation is closely following the book, 'Light Scattering by Systems of Particles', by (Adrian Doicu, et al). 


##  The function script  
The following script contains all functions that are used in plane wave reconstruction process. SphEM-Mesh_angle/test/Exp/WD_fuc.py 

### General functions
The following functions are general conversion functions.

#### <li> ```nth_derivative(f, wrt, n) ```
Taking n'th derivative for a flatten torch tensor </li>

#### <li> ```unit2cart(nu,phi,detach=True) -> (r_hat, nu_hat, phi_hat)```
A conversion function that generates an unit torch tensor in cartesian coordinates. When detach is true, the output will be a numpy array. Otherwise, it will be a torch tensor.  </li>


### WignerD functions
The following functions are related to WignerD functions. 



#### <li> ```nu2x(nu), x2nu(x)```
Convert a torch tensor argument (nu) into a torch tensor and its inverse (MISHCHENKO, B8)</li>

#### <li> ```A_s_mn(s,m,n)```
WignerD functions Coefficient with $s>=0$, $-s<=m$, $n<=s$(MISHCHENKO, B10) </li>
    
#### <li> ```WignerD_fuc(x,s,m,n)```
WignerD functions with torch tensor input/output (MISHCHENKO, B12) </li>
 


### Angular functions
The following functions are related to the spherical vector waves (SVW) angular parts. 

#### <li>``` pi_mn(nu,m,n)```
Angular scalar PI functions (order n,m) with torch tensor input/output (MISHCHENKO, C19). </li>
 
#### <li> ```tau_mn(nu,m,n)```
Angular scalar TAU functions (order n,m) with torch tensor input/output (MISHCHENKO, C20). </li>

#### <li> ```CBP_mn(nu,phi,m,n,CBP_theta=False)-> (C_mn,B_mn,P_mn)```
Angular vector functions (order n,m) with torch tensor input/output (MISHCHENKO, C19,C20,C21). The CBP_theta represents the Phi dependence. The default value is CBP_theta=False. When CBP_theta=True, Angular vector functions will have the variable phi dependence (MISHCHENKO, C16,C17,C18). </li>

#### <li> ```gamma_mn(m,n)```
Angular vector functions Coefficient (MISHCHENKO, C22). </li>


### Radial functions
The following functions are related to the spherical vector waves (SVW) radial parts. 

#### <li> ```bessel(x,n,kind=1,derivative=0)```
Bessel functions of two kinds: first(kind=1)/second(kind=2). Nth derivative (derivative=N). A torch tensor input (x) is required. </li>

#### <li> ```hankel(x,n,kind=1,derivative=0) ```
Hankel functions of two kinds: first(kind=1)/second(kind=2). Nth derivative (derivative=N). A torch tensor input (x) is required. </li>
 
 
### Spherical Vector Waves (SVW)
The following functions are related to the spherical vector waves (SVW). 


#### <li> ```MN_mn(r,nu,phi,e_k,m,n,RG=True) -> (M_mn,N_mn)```
M wave and N waves are defined in (MISHCHENKO, C14,C15): r,nu,phi are tensor objects. e_k defines the wavevector.When RG=True, the output gives the regular M and N waves.   </li>

### Coefficient Matrix of Expansion
The following functions are related to the coefficient Matrix for the plane wave reconstruction in Spherical vector waves basis.


#### <li>  ```ab_matrix(n,m,e_polar,e_k) -> (a_mn,b_mn)```
a_mn and b_mn are defined in (MISHCHENKO, C57,C58,C59): e_k defines the wavevector and e_polar defines the polarization of the plane wave.   </li>


### Plane wave Reconstruction
The plane wave reconstruction in Spherical vector waves basis.  

#### <li>  ```rec_plane_wave(r,nu,phi,e_k,e_polar,n_max) -> (E_total_x, E_total_y,E_total_z)```
e_k defines the wavevector and e_polar defines the polarization of the plane wave. n_max defines expansion length. E_total_x, E_total_y,E_total_z are the reconstruction plane wave in cartisian cordinates. (MISHCHENKO, C56)</li>


### Plane Wave generation
The following functions generate the plane wave in arbitrary polarization. 


#### <li> ```plane_wave_dir(polar_psi,polar_chi,E0,k,direction=False (True) ) -> e_polar (e_a,e_b)```
polar_psi (orientation angle) defines the orientation of the vibration ellipse (respects to initial vibration basis). polar_chi (ellipticity angle) defines the ellipticity of the vibration ellipse. E0 defines the scalar field magnitude. k defines the wavevector. When direction=True, it returns the initial vibration basis, e_a,e_b. (Adrian Doicu, 1.16)</li>


#### <li> ```angle_from_ek(ek) ->(nu, phi)```
for a given numpy array, it returns the polar angles. </li>

#### <li> ```plane_wave_grid(x,y,z,e_polar,K,omega=0,t0=0) ->(E_field_x,E_field_y,E_field_z)```
for a given numpy 2D grid,x,y,z, with polarization e_polar, wavevector K and angular frequency Omega, it returns the vector field values at the time t=t0.  </li>


### Error Map 
The following functions compared the ralative error of the reconstruction results agaist the ground truth. 


#### <li>```Emap_K(Geo_mesh,geo_mag,Geo_mesh_polar,K_mesh,K_mag,n_max=6,polar_psi=0,polar_chi=0,E0=1)-> (Error_x_array,Error_y_array,Error_z_array)```
    
This error map calculates the ralative error for all given wavevector directions with a fixed geometry size and wavevector magnitude. Geo_mesh (Geo_mesh_polar, in polar form) defines the input geometry mesh, and geo_mag defines the input geometry size. K_mesh defines the input wavevector mesh, and K_mag defines the input wavevector magnitude. Others are defined above. It returns the relative error in three components. One can visualize it by using 3D geometric heatmap.</li>


#### <li>```error_map(n_max, polar_psi, polar_chi,E0,K, K_mag_range,theta, phi, radius,radius_mag_range)-> (Error_map_x,Error_map_y,Error_map_z)```
    
This error map calculates the ralative error in a span of geometry sizes and wavevector magnitude, with a fixed wavevector direction. K_mag_range defines the range of wavevector magnitude, and radius_mag_range defines the range of geometry sizes. Others are defined above.It returns the relative error in three components. One can visualize it by using 2D plane heatmap. </li>


##  The geometry mesh generation script  
The following script contains all functions that are used in the mesh generation process. SphEM-Mesh_angle/src/Mesh_inter.py

### Class: ```Mesh_inter```

#### ```__init__(self,gmsh_file):```
process the gmesh file into numpy meshes.

### Class Method: ```Mesh_inter```


#### ```mesh_interpo(self,mesh_density) -> (N/A)```
Interpolate the triangular mesh of the given mesh file into a desired mesh size (defined by mesh_density). 

#### ```L_grid(self,N) -> (e_r, Data_grid_align,Radius_grid_NxN,Theta,Phi)```
Equally spaced NxN (grid) construction. e_r defines the unit direction vector shape(3,N**2), Data_grid_align defines the location of the desired location. Radius_grid_NxN,Theta,Phi, are the polar components equal space grid.  


#### ```mesh_save_rec(self,path) -> (N/A)```
Save the equal spaced mesh file into matlab variable at the desired path.

#### ```generate_rec_mesh(self,Theta,Phi,R) -> (x,y,z)```
generate the cartisian cordinate mesh from polar mesh. 
