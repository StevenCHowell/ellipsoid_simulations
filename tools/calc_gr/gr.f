c         1         2         3         4         5         6         7
c123456789012345678901234567890123456789012345678901234567890123456789012

        subroutine update_gr(coor,boxl,gr,dr,nbins,natoms)
        double precision boxl, rcutsq
        double precision coor(natoms,3)
        integer natoms,i,j,nbins,ri
        double precision gr(nbins)
        double precision xij,yij,zij,rijsq,dr

cf2py	intent(in):: coor,boxl,dr,nbins,natoms
cf2py	intent(in,out):: gr
cf2py	intent(hide):: i,j,ri
cf2py   depend(nbins):: gr
cf2py	intent(hide):: xij,yij,zij,rijsq,rcutsq

        rcutsq=(boxl/2.0)*(boxl/2.0)

        do 200,i=1,natoms-1
                do 300,j=i+1,natoms
                        xij=coor(i,1)-coor(j,1)
                        yij=coor(i,1)-coor(j,2)
                        zij=coor(i,1)-coor(j,3)

                        xij=xij-boxl*(ANINT(xij/boxl))
                        yij=yij-boxl*(ANINT(yij/boxl))
                        zij=zij-boxl*(ANINT(zij/boxl))

                        rijsq=xij*xij+yij*yij+zij*zij

                        if (rijsq .lt. rcutsq) then
                                ri = int(dsqrt(rijsq)/dr)
                                gr(ri) = gr(ri) + 2.0
                        endif

  300   continue
  200   continue

        end

c         1         2         3         4         5         6         7
c123456789012345678901234567890123456789012345678901234567890123456789012
c
c     to setup and incorporate into python:
c
c     python setup_gr.py build
c     cp build/lib*/distance.so ./
c
c     to call this from python:
c
c     import sys ; sys.path.append('./')
c     from distance import distance
c
c     dist = distance(coor, dist)
c
c         1         2         3         4         5         6         7
c123456789012345678901234567890123456789012345678901234567890123456789012
