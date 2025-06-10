! Selecting loops and remove them, generating a new CA file: "newca.input"
!
! Note:A box boundary information is necessary, because the atom coordinates
!      in CA files are not in periodic conditions. In order to satisfy 
!      Cluster Analysis after nodes-counting, coordinates of nodes should
!      satisfy PPP boundary conditions.
!
! Inputs: 1. some positions near to removing loops
!         2. dump.0 (using its box size informations for periodic boundary)
!         3. disline.input (CA file), removing loops and create a new CA file
!  
!lwh 
!2023.10.25  
!2024.6.21 modified

program main
implicit none
!------------------------------state the variables-------------------
!added variables
integer :: i , j , k !omit
integer :: lines,l0
character(len=20) :: char1

!----------- remove loops information
integer :: nloops,flag_vertex_line
double precision,allocatable,dimension(:,:) :: loops
double precision :: dist,cutoff

!-----------read auxiliary datas
integer,allocatable,dimension(:) :: LineNum , StruType ,NumLinepos !the order-number & 
!of dislocation lines , structure type , vertex quantity
double precision,allocatable,dimension(:,:) :: LineBergers    !atom positions , &
!-----------read vertex positions datas
double precision,allocatable,dimension(:,:,:) :: pos    !atom positions , &
!the bergers of dislocation lines 
integer :: vertexnum,newNUM
integer :: Num , Nmax  !the quantity of lines , the quantity of positions
character(len=20) :: char0    !omit
!-------------------------------- part 2 varias : box.dat
integer :: n0
double precision :: xlo,xhi,ylo,yhi,zlo,zhi

!---------------------------------initialization-----------------------

! Removing loops dislocations in disline.input CA file by
! defining a position which near to the loop, see cutoff
! parameter. Note that 1) do not set the position near to 
! other dislocations in case of deleting original dislocation 
! networks; 2) the dislocation number in new CA file should
! change manually, be attention! 3) To assure one loop atom
! is in the cutoff-radius sphere is enough, we don't have 
! to wrap the whole line!!!
!
! input1: disline.input, original CA file
! input2: dump.0, using its x/y/zlo and x/y/zhi, only, for perodic
! output1: newca.input, new CA file but dislocation number should
!          be changed manually.
! output2: nodes.dump, print nodes positions.

nloops = 1 ! how many loops are about to remove
allocate(loops(1:nloops,1:3))
loops(1,1:3) = (101.096, 96.4252, -99.2811)
! loops(2,1:3) = (/97.991, 90.7169, -67.9176/) ! positions
cutoff = 15.0 ! cutoff distance 


!read in words.dat
open(1,file= 'words.dat' ,status="old") !access how many lines in this file
read(1,*)lines
close(1)
write(*,*)"-->Total lines:",lines
!which line is start line
! 1-st Reading "disline.input" to identify first-reading-line
open(1,file= './datafiles_inputs/disline.input' ,status="old") !dislocations start at line-i
do i = 1 , lines
        read(1,*)char1
        if(char1=='DISLOCATIONS')then
                l0 = i !line starts
                write(*,*)char1,'start at line',l0
        endif
enddo
close(1)
write(*,*)"-->Start line:",l0

!------------------------------read data from input file,except Pos------------
! 2-nd Reading "disline.input" to identify how many atoms at each disline, so to
!       identify the maximum number of atoms in single disline and to allocate 
!       the position matrix
open(10,file='./datafiles_inputs/disline.input',status='old')
do i = 1 , l0-1 !omit first 'l0-1' lines
        read(10,*)
enddo
read(10,*)char0,Num ! how many dislocation lines
allocate(LineNum(1:Num))
allocate(StruType(1:Num))
allocate(NumLinepos(1:Num))
allocate(LineBergers(1:Num,1:3))

do i = 1,Num
  read(10,*)LineNum(i) ! id of line, from 0 to n !!!
  read(10,*)LineBergers(i,1:3)
  read(10,*)StruType(i)
  read(10,*)NumLinepos(i) ! num of vertexes of the line
  do j = 1,NumLinepos(i)
    read(10,*)
  enddo
enddo

close(10)

!-----------------------allocate and load and moving Pos----------------- 
!-----------------------and select nodes sites--------------------------- 
! Reading "dump.0" to using its box size infomations
open(1,file='./datafiles_inputs/dump.0',status='old') ! DUMP, get box boundary
read(1,*)
read(1,*)
read(1,*)
read(1,*)n0 !number-of-atoms
read(1,*)
read(1,*)xlo,xhi
read(1,*)ylo,yhi
read(1,*)zlo,zhi
close(1)
write(*,*)"-->The num of atoms in dump.0:",n0

!! allocate vertex matrix, which contains all vertex of all 
!! dislocation lines. "pos(i,j,k)". Total vertex counting,
!! the "vertexnum".
Nmax = maxval(NumLinepos(1:Num)) ! maximum atoms in single line
allocate(pos(Num,Nmax,1:3)) ! (i,j,k), i-th disline,j-th atom of i-th disline
vertexnum = 0

! 3-rd Reading "disline.input" to read-in positions of vertexes
open(10,file='./datafiles_inputs/disline.input',status='old')
do i = 1 , l0-1 !omit first 'l0-1' lines
        read(10,*)
enddo
read(10,*)
do i = 1,Num
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  do j = 1,NumLinepos(i)
        vertexnum = vertexnum + 1 ! counting total vertex
        read(10,*)pos(i,j,1:3)
        if(pos(i,j,1)<xlo)  pos(i,j,1)=pos(i,j,1)+xhi-xlo
        if(pos(i,j,1)>=xhi) pos(i,j,1)=pos(i,j,1)-xhi+xlo
        if(pos(i,j,2)<ylo)  pos(i,j,2)=pos(i,j,2)+yhi-ylo
        if(pos(i,j,2)>=yhi) pos(i,j,2)=pos(i,j,2)-yhi+ylo
        if(pos(i,j,3)<zlo)  pos(i,j,3)=pos(i,j,3)+zhi-zlo
        if(pos(i,j,3)>=zhi) pos(i,j,3)=pos(i,j,3)-zhi+zlo
  enddo
enddo
close(10)
write(*,*)"Total vertex ", vertexnum


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Working on vertexs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! new ca file and dump vertex position file
newNUM = Num-nloops ! Num, originally number of dislocation lines; newNUM,after removing target loops
write(*,*)"@@@@ newNUM calculated:",newNUM
open(1,file='newca.input',status='replace') ! output file
write(1,*)"New CA file but Removed Loop Dislines"
write(1,*)"DISLOCATIONS ", newNUM ! new dislocation number, and using the guessed dislocation number for the moment

open(2,file='nodes.dump',status='replace') ! output file
write(2,'(A)') "ITEM: TIMESTEP"
write(2,'(A)') "0"
write(2,'(A)') "ITEM: NUMBER OF ATOMS"
write(2,*) vertexnum ! using total vertex for the moment
write(2,'(A)') "ITEM: BOX BOUNDS pp pp pp"
write(2,*)  xlo,xhi 
write(2,*)  ylo,yhi
write(2,*)  zlo,zhi
write(2,'(A)') "ITEM: ATOMS id type x y z"

newNUM = 0 ! as re-counting, since we don't assure its quantity
vertexnum = 0 ! recounting the remained vertex as a pointer, and changing in dump file.
loop1:do i = 1,Num ! i-th vertex-line, pos(i,j,1:3)
        flag_vertex_line = 0 ! 1:near to target loop
        loop2:do j = 1,NumLinepos(i)  ! j-th vertex of i-th line
                loop3:do k = 1 , nloops ! k-th target loop position, loops(k.1:3)
                        dist = sqrt((pos(i,j,1)-loops(k,1))**2+(pos(i,j,2)-loops(k,2))**2+(pos(i,j,3)-loops(k,3))**2)
                        if (dist<=10) then
                                flag_vertex_line = 1 ! near to target loop, remove
                                exit loop2
                        endif
                enddo loop3
        enddo loop2
        if(flag_vertex_line==0)then ! away from target loop, save
                newNUM = newNUM + 1
                write(1,*) LineNum(i) , 'NOT the real disline number !!!'
                write(1,*) LineBergers(i,1:3)
                write(1,*) StruType(i)
                write(1,*) NumLinepos(i)
                do j = 1,NumLinepos(i)
                        vertexnum = vertexnum + 1
                        write(1,*) pos(i,j,1:3)
                        write(2,*) vertexnum,"1",pos(i,j,1:3)
                enddo
        endif
enddo loop1
write(*,*)"@@@@ newNUM counted:",newNUM
write(*,*)"@@@@ If the two newNUM is not the same, using the second 'counted' one!!!"

close(1)
close(2)
write(*,*)"!!! Changing This Value To 'nodes.dump' file :"
write(*,*)"Non-loop vertex ", vertexnum

!----------------------end-----------------------------------
end
