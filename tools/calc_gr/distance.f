        subroutine distance(coor,dist,n)
        double precision coor(n,3),dist(n,n)
        integer n,i,j
        double precision xij,yij,zij

cf2py	intent(in):: coor,n
cf2py	intent(in,out):: dist
cf2py	intent(hide):: xij,yij,zij,

c         1         2         3         4         5         6         7
c123456789012345678901234567890123456789012345678901234567890123456789012
c
c     to setup and incorporate into python:
c
c     python setup_distance.py build
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

        do 200,i=1,n-1
            do 300,j=i+1,n
                xij=coor(i,1)-coor(j,1)
                yij=coor(i,2)-coor(j,2)
                zij=coor(i,3)-coor(j,3)

                dist(i,j)=sqrt(xij*xij+yij*yij+zij*zij)

  300   continue
  200   continue

        end

c         1         2         3         4         5         6         7
c123456789012345678901234567890123456789012345678901234567890123456789012
