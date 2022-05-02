program pegarconf_hb
implicit none

integer*8                   :: i, j, k, l, nsoluto, nwaters, nconf, npc, conta, ag, aa
integer                     :: stat
real*8, allocatable         :: mol_qm(:,:), mol_mm(:,:,:,:)
real*8                      :: cargasp(3)
integer*8, allocatable      :: n(:)
character(100), allocatable :: title(:)
character(5), allocatable   :: atomname(:)
character(100)              :: xyzfile, fgaussian
character(5)                :: dummy
character(9)                :: junk

!!!gaussian cabecalho
character(150)              :: cabecalho(9)

cargasp(1)   = -0.8476d0
cargasp(2:3) =  0.4238d0 

read(5,'(a9,a)',iostat=stat) junk, xyzfile  !!tem que ser o aquivo _all.xyz
call read_stat(xyzfile,stat)
call read_flag(junk,'XYZFILE =')

read(5,'(a9,a)',iostat=stat) junk, fgaussian
call read_stat(fgaussian,stat)
call read_flag(junk,'GAUFILE =')

read(5,'(a9,i8)',iostat=stat) junk, nsoluto
call read_stat('nsoluto',stat)
call read_flag(junk,'NSOLUTO =')

read(5,'(a9,i8)',iostat=stat) junk, nwaters
call read_stat('nwaters',stat)
call read_flag(junk,'NWATERS =')

read(5,'(a9,i8)',iostat=stat) junk, nconf
call read_stat('nconf',stat)
call read_flag(junk,'NCONNFS =')

read(5,'(a9,i8)',iostat=stat) junk, npc
call read_stat('npc',stat)
call read_flag(junk,'NPCHARG =')

allocate(title(nconf))
allocate(n(nconf))

xyzfile   = trim(adjustl(xyzfile))
fgaussian = trim(adjustl(fgaussian))

open(1,file=xyzfile,iostat=stat)
if(stat .ne. 0) stop'Abort, problem opening the xyzfile!'

open(2,file=fgaussian,iostat=stat)
if(stat .ne. 0) stop'Abort, problem opening the cabecalho do gaussian!'
do l=1,9
  read(2,'(a150)') cabecalho(l)
end do

open(3,file='confs_hb.xyz')
if(stat .ne. 0) stop'Abort, problem opening the file confs_hb.xyz!'

open(4,file='confs_hb_tudoqm.com')
if(stat .ne. 0) stop'Abort, problem opening the file confs_hb_tudoqm.com!'


if( mod(npc,3) .ne. 0 ) stop'Abort, numero de moleculas de água é quebrado!!'

conta = 0

do i=1,nconf
  read(1,*,iostat=stat) n(i)
  call read_stat('n(i)',stat) 
  if( (n(i)-npc) .eq. (nsoluto+3*nwaters) ) then
      conta = conta + 1
      read(1,'(a)',iostat=stat) title(i)
      call read_stat('title(i)',stat)
      allocate( mol_qm(n(i),3), mol_mm(n(i),npc/3,3,3) )
      allocate( atomname(n(i)-npc) )
      write(3,'(i0)') n(i)
      write(3,'(a)') title(i)
      do l=1,6
        write(6,'(a)') trim(adjustl(cabecalho(l)))
        write(4,'(a)') trim(adjustl(cabecalho(l)))
      end do
      write(6,'(a,i0)') trim(adjustl(cabecalho(7)))//', Config = '//trim(adjustl(title(i)(24:34)))//', conta = ', conta
      write(6,'(a)') trim(adjustl(cabecalho(8)))
      write(6,'(a)') trim(adjustl(cabecalho(9)))

      write(4,'(a,i0)') trim(adjustl(cabecalho(7)))//', Config = '//trim(adjustl(title(i)(24:34)))//', conta = ', conta
      write(4,'(a)') trim(adjustl(cabecalho(8)))
      write(4,'(a)') trim(adjustl(cabecalho(9)))

      do j=1,(n(i)-npc)
        read(1,*,iostat=stat) atomname(j), (mol_qm(j,k),k=1,3)
        call read_stat('(mol_qm(j,k),k=1,3)',stat)
        write(6,'(a,4x,3(f13.9,2x))') atomname(j), (mol_qm(j,k),k=1,3)
        write(4,'(a,4x,3(f13.9,2x))') atomname(j), (mol_qm(j,k),k=1,3)
        write(3,'(a,4x,3(f13.9,2x))') atomname(j), (mol_qm(j,k),k=1,3)
      enddo
      do ag=1,npc/3
        do aa=1,3
          read(1,*,iostat=stat) dummy, (mol_mm(i,ag,aa,k),k=1,3)
          call read_stat('(mol_mm(j,k),k=1,3)',stat)
          write(3,'(a,4x,3(f13.9,2x))') dummy, (mol_mm(i,ag,aa,k),k=1,3)
        end do
      enddo
      write(6,*)
      do ag=1,npc/3
        do aa=1,3
          write(6,'(4(f13.9,2x))') (mol_mm(i,ag,aa,k),k=1,3), cargasp(aa)
          if(aa.eq.1) write(4,'(a,2x,3(f13.9,2x))') 'O', (mol_mm(i,ag,aa,k),k=1,3)   
          if(aa.eq.2 .or. aa.eq.3) write(4,'(a,2x,3(f13.9,2x))') 'H', (mol_mm(i,ag,aa,k),k=1,3) 
        end do
      enddo
      write(6,*)
      write(4,*)
      if(i .ne. nconf) write(6,'(a)') '--link1--'
      if(i .ne. nconf) write(4,'(a)') '--link1--'
      deallocate(mol_qm, mol_mm, atomname)
  else if ( (n(i)-npc) .ne. (nsoluto+3*nwaters) .and. i .ne. nconf ) then
    do l=1,(n(i)+1)
      read(1,*)
    end do
  end if
end do

contains

subroutine read_stat(cvar,stat)
implicit none

integer, intent(in) :: stat
character(len=*), intent(in) :: cvar

if(stat .ne. 0) then
  write(6,'(2a)') 'Problem reading: ', cvar
  stop'Abort, see the mesage in the output file!'
end if
end subroutine read_stat

subroutine read_flag(junk,flag)
implicit none

character(len=*), intent(in) :: junk, flag

if(trim(adjustl(junk)) .ne. flag) then
  write(6,'(5a)') 'Abort: ',junk,' is not equal to ',flag,'!!'
  stop'Abort, see the mesage in the output file!'
end if
end subroutine read_flag


end program pegarconf_hb
