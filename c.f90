program rcb
implicit none
external  RefineRule

integer, parameter :: nvmax = 100000, ntmax = 2 * nvmax, nbmax = 50000

integer nv, nt, nb

real*8 vrt(2,nvmax)
integer tri(3,ntmax), labelT(ntmax)

integer bnd(2,nbmax), labelB(nbmax)

integer, parameter :: maxlevel = 150

integer, parameter :: MaxWi = 5500000
integer  iW(MaxWi)

logical history(maxlevel*ntmax)
integer  nlevel
integer  ilevel, iERR

nlevel = 10

nv = 4
vrt(1,1) = 0d0
vrt(2,1) = 0d0
vrt(1,2) = 0d0
vrt(2,2) = 1d0
vrt(1,3) = 1d0
vrt(2,3) = 1d0
vrt(1,4) = 1d0
vrt(2,4) = 0d0

nt = 2
tri(1,1) = 1
tri(2,1) = 2
tri(3,1) = 3
labelT(1) = 1

tri(1,2) = 3
tri(2,2) = 4
tri(3,2) = 1
labelT(2) = 2

nb = 4
bnd(1,1) = 1
bnd(2,1) = 2
labelB(1) = 1

bnd(1,2) = 2
bnd(2,2) = 3
labelB(2) = 1

bnd(1,3) = 3
bnd(2,3) = 4
labelB(3) = 2

bnd(1,4) = 4
bnd(2,4) = 1
labelB(4) = 1

write(*, *) 'Initial mesh:   numbers of nodes and triangles:',nv, nt

call graph_demo(nv,vrt, nt,tri, 'mesh_initial.ps','Initial mesh')
      
call InitializeRCB(nt, ntmax, vrt, tri, MaxWi, iW, iERR)

do ilevel = 1, nlevel
    call LocalRefine(nv, nvmax, nb, nbmax, nt, ntmax, vrt, tri, bnd, labelB, &
    labelT, RefineRule, ilevel, maxlevel, history, MaxWi, iW, iERR)
end do
write(*, *) 'Refined mesh:   numbers of nodes and triangles:',nv, nt

call graph_demo(nv,vrt, nt,tri, 'mesh_final1.ps', 'Locally refined mesh')
      
end program rcb


real*8 function func(x)
implicit none
real*8 :: x
real*8, parameter :: t = 0.1
func = 0.5 + 0.1 * sin(4 * 3.14 * (x - t))
return
end function func


Subroutine RefineRule (nE, IPE, XYP, verf, ilevel)
implicit none
      
integer nE
integer IPE(3,*)
real*8 XYP(2,*)
integer verf(*)
integer ilevel
integer i, j
real*8 x1, y1, x2, y2, x3, y3, cx, cy
integer up, down
integer, parameter :: parts = 10
real*8 func
external func

if (ilevel .le. 0) then
    do i = 1, nE
        verf(i) =  0
    end do
else
    do i = 1, nE
        up = 0
        down = 0
        x1 = XYP(1, IPE(1, i))
        y1 = XYP(2, IPE(1, i))
        x2 = XYP(1, IPE(2, i))
        y2 = XYP(2, IPE(2, i))
        x3 = XYP(1, IPE(3, i))
        y3 = XYP(2, IPE(3, i))

        do j = 0, parts
            cx = x1 + (x2 - x1) / parts * j
            cy = y1 + (y2 - y1) / parts * j
            if (func(cx) .gt. cy) then
                up = 1
            else
                down = 1
            end if
        end do
        do j = 0, parts
            cx = x1 + (x3 - x1) / parts * j
            cy = y1 + (y3 - y1) / parts * j
            if (func(cx) .gt. cy) then
                up = 1
            else
                down = 1
            end if
        end do
        do j = 0, parts
            cx = x3 + (x2 - x3) / parts * j
            cy = y3 + (y2 - y3) / parts * j
            if (func(cx) .gt. cy) then
                up = 1
            else
                down = 1
            end if
        end do
        
        if (up .eq. 1 .and. down .eq. 1) then
            verf(i) =  2
        else
            verf(i) = 0
        end if
    end do
end if
end subroutine
