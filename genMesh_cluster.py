'''This is a script to generate and save triangular finite element meshes with polydisperse spherical exclusions'''

import dolfin as df
import os
import shutil
import scipy.io as sio
import mshr
import sys
import time


# Global parameters
nMeshes = 10000
nElements = 256  # PDE discretization
foldername1 = './data/meshSize=' + str(nElements)


# Parameters only for 'circles' mode
nExclusionsDist='logn'
nExclusionParams = [7.8, 0.2]
coordinateDist = 'GP'
# to avoid circles on boundaries. Min. distance of circle centers to (lo., r., u., le.) boundary
# negative margin means no margin
margins = [0.003, 0.003, 0.003, 0.003]
origin_margin = .03
substractCorners = False     # Substracts circles from domain corners s.t. flow cannot pass
radiiDist = 'lognGP'
r_params = [-5.53, 0.2]
# for x~gauss
coordinate_cov = [[0.55, -0.45], [-0.45, 0.55]]
coordinate_mu = [.8, .8]
# for x~GP
covFun = 'squaredExponential'
cov_l = 0.08
sig_scale = 1.2
sigmaGP_r = 0.01
lengthScale_r = .05


# we need to remove '.0' for correct path names
if nExclusionParams[0] % 1 == 0:
    nExclusionParams[0] = int(nExclusionParams[0])
if nExclusionParams[1] % 1 == 0:
    nExclusionParams[1] = int(nExclusionParams[1])
if r_params[0] % 1 == 0:
    r_params[0] = int(r_params[0])
if r_params[1] % 1 == 0:
    r_params[1] = int(r_params[1])
if cov_l % 1 == 0:
    cov_l = int(cov_l)
if sig_scale % 1 == 0:
    sig_scale = int(sig_scale)


print('Generating mesh with non-overlapping circular exclusions...')
foldername = foldername1 + '/nonOverlappingDisks/margins=' + str(margins[0]) + '_' + str(margins[1]) + '_' + \
             str(margins[2]) + '_' + str(margins[3]) + '/N~logn/mu=' + str(nExclusionParams[0]) + '/sigma=' +\
             str(nExclusionParams[1]) + '/x~' + coordinateDist

if coordinateDist == 'gauss':
    foldername += '/mu=' + str(coordinate_mu[0]) + '_' + str(coordinate_mu[1]) + '/cov=' + \
                  str(coordinate_cov[0][0]) + '_' + str(coordinate_cov[0][1]) + '_' + str(coordinate_cov[1][1]) +\
                  '/'
elif coordinateDist == 'gauss_randmu':
    foldername += '/mu=rand' + '/cov=' + \
                  str(coordinate_cov[0][0]) + '_' + str(coordinate_cov[0][1]) + '_' + str(coordinate_cov[1][1]) + \
                  '/'
elif coordinateDist == 'GP':
    foldername += '/cov=' + covFun + '/l=' + str(cov_l) + '/sig_scale=' + str(sig_scale) + '/'

elif coordinateDist == 'engineered' or coordinateDist == 'tiles':
    foldername += '/'

foldername += 'r~' + radiiDist
if radiiDist == 'lognGP':
    foldername += '/mu=' + str(r_params[0]) + '/sigma=' + str(r_params[1]) +\
                  '/sigmaGP_r=' + str(sigmaGP_r) + '/l=' + str(lengthScale_r)
else:
    foldername += '/mu=' + str(r_params[0]) + '/sigma=' + str(r_params[1])
if origin_margin:
    foldername += '/origin_rejection=' + str(origin_margin)
if not os.path.exists(foldername):
    os.makedirs(foldername)


'''first copy 'microstructureInformation_nomesh' to 'microstructureInformation', to give signal that mesh is
generated, so that no other job is taking the same microstructure to generate a mesh'''
mesh_name_iter = 0
meshfile = foldername + '/mesh' + str(mesh_name_iter) + '.mat'
nomeshInfoFile = foldername + '/microstructureInformation_nomesh' + str(mesh_name_iter) + '.mat'
while mesh_name_iter < nMeshes:
    # create computation_started.txt if not existent
    if not os.path.isfile(foldername + '/computation_started.txt'):
        started_file = open(foldername + '/computation_started.txt', 'w')
        started_file.close()

    started_file = open(foldername + '/computation_started.txt', 'r')
    started_computations = started_file.readlines()
    started_file.close()
    print('started_computations == ', started_computations)

    while ((not os.path.isfile(nomeshInfoFile)) or os.path.isfile(meshfile)
           or ((str(mesh_name_iter) + '\n') in started_computations)):
        mesh_name_iter += 1
        meshfile = foldername + '/mesh' + str(mesh_name_iter) + '.mat'
        nomeshInfoFile = foldername + '/microstructureInformation_nomesh' + str(mesh_name_iter) + '.mat'

    with open(foldername + '/computation_started.txt', 'a') as started_file:
        started_file.write(str(mesh_name_iter) + '\n')
    print('Generating mesh number ', mesh_name_iter)

    print('Loading microstructural data...')
    matfile = sio.loadmat(nomeshInfoFile)
    diskCenters = matfile['diskCenters']
    diskRadii = matfile['diskRadii']
    diskRadii = diskRadii.flatten()
    print('... microstructural data loaded.')

    # Generate domain object
    print('Generating domain object...')
    domain = mshr.Rectangle(df.Point(0.0, 0.0), df.Point(1.0, 1.0))
    nCircles = diskCenters.shape[0]
    for blob in range(0, nCircles):
        c = df.Point(diskCenters[blob, 0], diskCenters[blob, 1])
        domain -= mshr.Circle(c, diskRadii[blob])
    print('...domain object generated.')

    # Final step: generate mesh using mshr
    print('generating FE mesh...')
    sys.stdout.flush()
    mesh = mshr.generate_mesh(domain, nElements)
    print('...FE mesh generated.')

    # save to local directory first
    print('Saving mesh to NFS...')
    sys.stdout.flush()
    t0 = time.time()
    sio.savemat(meshfile, {'x': mesh.coordinates(), 'cells': mesh.cells() + 1}, do_compression=True)
    t1 = time.time()
    print('...done. Time: ', t1 - t0)
    sys.stdout.flush()

    # save vertex coordinates and cell connectivity to mat file for easy read-in to matlab
    # move microstructureInformation file, this is the signal that a job is already generating a mesh
    shutil.move(foldername + '/microstructureInformation_nomesh' + str(mesh_name_iter) + '.mat',
             foldername + '/microstructureInformation' + str(mesh_name_iter) + '.mat')
    print('... ./microstructureInformation_nomesh' + str(mesh_name_iter) + '.mat renamed.')
    sys.stdout.flush()


