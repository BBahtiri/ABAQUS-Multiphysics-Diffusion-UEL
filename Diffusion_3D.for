! User element subroutine for the diffusion of water in epoxy
! Units: s, mm, MPa
      

      module kvisual
      implicit none
      real*8 UserVar(70000,16,8)
      integer nelem
      save
      end module
      
      subroutine uel(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,
     1 props,nprops,coords,mcrd,nnode,u,du,v,a,jtype,time,dtime,
     2 kstep,kinc,jelem,params,ndload,jdltyp,adlmag,predef,npredf,
     3 lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,njpro,period)

      use kvisual
      include 'aba_param.inc' !implicit real(a-h o-z)
      
      dimension rhs(mlvarx,*),amatrx(ndofel,ndofel),props(*),svars(*),
     1 energy(*),coords(mcrd,nnode),u(ndofel),du(mlvarx,*),v(ndofel),
     2 a(ndofel),time(2),params(*),jdltyp(mdload,*),adlmag(mdload,*),
     3 ddlmag(mdload,*),predef(2,npredf,nnode),lflags(*),jprops(*)

      parameter(ndim=3,ntens=6,ninpt=8,nsvint=16)
      
      dimension wght(ninpt),dN(1,nnode),dNdz(ndim,nnode),dNS(ninpt),
     2 dNdx(ndim,nnode),b(ntens,nnode*ndim),ddsdde(ntens,ntens),
     3 stress(ntens),stran(ntens),bC(ndim,nnode),xm(nnode,nnode),
     4 xk(nnode,nnode),BB(nnode,nnode),SHa(nnode,1),coord320(ndim,nnode),
     5 dstran(ntens),statevLocal(nsvint),COFACTOR(ndim,ndim)
      
      data wght /1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0/
        
!     initialising
      do k1=1,ndofel
       rhs(k1,1)=0.d0
      end do
      amatrx=0.d0
      
!     find number of elements          
      if (dtime.eq.0.d0) then
       if (jelem.eq.1) then
        nelem=jelem
       else
        if (jelem.gt.nelem) nelem=jelem 
       endif 
      endif      
      
!     reading parameters
	  if(kstep.eq.1) then
		D=props(3) ! Diffusion coefficient
	  endif
	  if(kstep.eq.2) then
	    D=props(4)
	  endif
      Vh=8000.d0 ! Molar volume of H
      T=300.d0 ! Temperature
      R=8314.5d0 ! Gas constant

!     compute the hydrostatic stress
      SHa=0.d0
      coord320=0.d0
      coord320(1,1)=-1.d0
      coord320(2,1)=-1.d0
	  coord320(3,1)=-1.d0
      coord320(1,2)=1.d0
      coord320(2,2)=-1.d0
      coord320(3,2)=-1.d0
      coord320(1,3)=1.d0  
      coord320(2,3)=1.d0
      coord320(3,3)=-1.d0
      coord320(1,4)=-1.d0
      coord320(2,4)=1.d0
      coord320(3,4)=-1.d0  
      coord320(1,5)=-1.d0        
      coord320(2,5)=-1.d0
	  coord320(3,5)=1.d0
	  coord320(1,6)=1.d0
	  coord320(2,6)=-1.d0
	  coord320(3,6)=1.d0
	  coord320(1,7)=1.d0
	  coord320(2,7)=1.d0
	  coord320(3,7)=1.d0
	  coord320(1,8)=-1.d0
	  coord320(2,8)=1.d0
	  coord320(3,8)=1.d0
	  
	  coord320(2,9)=-1.d0
	  coord320(3,9)=-1.d0
	  coord320(2,11)=1.d0
	  coord320(3,11)=-1.d0
	  coord320(2,13)=-1.d0
	  coord320(3,13)=1.d0
	  coord320(2,15)=1.d0
	  coord320(3,15)=1.d0
	  
	  coord320(1,10)=1.d0
	  coord320(3,10)=-1.d0
	  coord320(1,12)=-1.d0
	  coord320(3,12)=-1.d0
	  coord320(1,14)=1.d0
	  coord320(3,14)=1.d0
	  coord320(1,16)=-1.d0
	  coord320(3,16)=1.d0
	  
	  coord320(1,17)=-1.d0
	  coord320(2,17)=-1.d0
	  coord320(1,18)=1.d0
	  coord320(2,18)=-1.d0
	  coord320(1,19)=1.d0
	  coord320(2,19)=1.d0
	  coord320(1,20)=-1.d0
	  coord320(2,20)=1.d0
	  
	  
      do inod=1,nnode
       g=dsqrt(3.d0)*coord320(1,inod)
       h=dsqrt(3.d0)*coord320(2,inod)
	   l=dsqrt(3.d0)*coord320(3,inod)
       dNS(1)=(1.d0-g)*(1.d0-h)*(1.d0-l)/8.d0
       dNS(2)=(1.d0+g)*(1.d0-h)*(1.d0-l)/8.d0
       dNS(3)=(1.d0-g)*(1.d0+h)*(1.d0-l)/8.d0
       dNS(4)=(1.d0+g)*(1.d0+h)*(1.d0-l)/8.d0
	   dNs(5)=(1.d0-g)*(1.d0-h)*(1.d0+l)/8.d0
	   dNs(6)=(1.d0+g)*(1.d0-h)*(1.d0+l)/8.d0
	   dNs(7)=(1.d0-g)*(1.d0+h)*(1.d0+l)/8.d0
	   dNs(8)=(1.d0+g)*(1.d0+h)*(1.d0+l)/8.d0
       do i=1,ninpt
        isvinc=(i-1)*nsvint
        SHa(inod,1)=SHa(inod,1)+dNS(i)*svars(isvinc+10)
       end do
      end do      
      
      do kintk=1,ninpt
!     evaluate shape functions and derivatives
       call kshapefcn(kintk,ninpt,nnode,ndim,dN,dNdz)      
       call kjacobian(jelem,ndim,nnode,coords,dNdz,djac,dNdx,mcrd,COFACTOR)
       dvol=wght(kintk)*djac
       
!     form B-matrix
       b=0.d0
	   do inod=1,nnode
	    bC(1,inod)=dNdx(1,inod)
		bC(2,inod)=dNdx(2,inod)
		bC(3,inod)=dNdx(3,inod)
!		b(1,inod)=dNdx(1,inod)
!		b(2,inod)=dNdx(2,inod)
!		b(3,inod)=dNdx(3,inod)
	   end do	
	   do k=1,nnode
	     b(1,1+ndim*(k-1))=dNdx(1,k)
		 b(2,2+ndim*(k-1))=dNdx(2,k)
		 b(3,3+ndim*(k-1))=dNdx(3,k)
		 b(4,1+ndim*(k-1))=dNdx(2,k)
		 b(4,2+ndim*(k-1))=dNdx(1,k)
		 b(5,2+ndim*(k-1))=dNdx(3,k)
		 b(5,3+ndim*(k-1))=dNdx(2,k)
		 b(6,1+ndim*(k-1))=dNdx(3,k)
		 b(6,3+ndim*(k-1))=dNdx(1,k)
	   end do
	   
!	   b(1,1)=dndx(1,1)
!	   b(2,2)=dndx(2,1)
!	   b(3,3)=dndx(3,1)
!	   b(4,1)=dndx(2,1)
!	   b(4,2)=dndx(1,1)
!	   b(5,2)=dndx(3,1)
!	   b(5,3)=dndx(2,1)
!	   b(6,1)=dndx(3,1)
!	   b(6,3)=dndx(1,1)
!      do inod=1,nnode-1
!       b(1,3*inod+1)=dndx(1,inod+1)
!       b(2,3*inod+2)=dndx(2,inod+1)
!		b(3,3*inod+3)=dndx(3,inod+1)
!		
!       b(4,3*inod+1)=dNdx(2,inod+1)
!		b(4,3*inod+2)=dNdx(1,inod+1)
!		b(5,3*inod+2)=dNdx(3,inod+1)
!		b(5,3*inod+3)=dNdx(2,inod+1)
!		b(6,3*inod+1)=dNdx(3,inod+1)
!		b(6,3*inod+3)=dNdx(1,inod+1)
!      end do                     
       
!     compute from nodal values

       cL=0.d0
       do inod=1,nnode
        cL=cL+dN(1,inod)*u(3*nnode+inod)
       end do   
       
           
!     compute the increment of strain 
       dstran=matmul(b,du(1:ndim*nnode,1))
       call kstatevar(kintk,nsvint,svars,statevLocal,1)
       stress=statevLocal(1:ntens)
       stran(1:ntens)=statevLocal((ntens+1):(2*ntens))
       
       
!     call umat to obtain stresses and constitutive matrix 
       call kumat(props,ddsdde,stress,dstran,ntens,statevLocal,ndim)
       stran=stran+dstran
       
       
       statevLocal(1:ntens)=stress(1:ntens)
       statevLocal((ntens+1):(2*ntens))=stran(1:ntens)
       statevLocal(2*ntens+2)=(stress(1)+stress(2)+stress(3))/3.d0 !SH
       
       call kstatevar(kintk,nsvint,svars,statevLocal,0)

!	amatrx for stress contribution
 
       amatrx(1:60,1:60)=amatrx(1:60,1:60)+dvol*
     1 matmul(matmul(transpose(b),ddsdde),b)
        
       rhs(1:60,1)=rhs(1:60,1)-
     1 dvol*(matmul(transpose(b),stress))       
           
       
  
!	Diffusion process
        
       xm=matmul(transpose(dN),dN)/D
       BB=matmul(transpose(bC),bC)
       xk=BB-Vh/(R*T)*matmul(BB,matmul(SHa,dN))
        
		
       amatrx(61:80,61:80)=amatrx(61:80,61:80)+dvol*(xm/dtime+xk)
        
       rhs(61:80,1)=rhs(61:80,1)-dvol*(matmul(xk,u(61:80))+
     1 matmul(xm,du(61:80,1))/dtime)
                   
! output
       UserVar(jelem,1:6,kintk)=statevLocal(1:6)
       UserVar(jelem,7:14,kintk)=statevLocal((ntens+1):(2*ntens+2))
       UserVar(jelem,(2*ntens+3),kintk)=cL
               write(*,*) '-------------------------------------------------'
			   write(*,*) 'DU: ',du(1:80,1)
			   write(*,*) 'Stress:', stress
			   write(*,*) '-------------------------------------------------'	   
		  
      end do       ! end loop on material integration points
      RETURN
      END
      
      subroutine kshapefcn(kintk,ninpt,nnode,ndim,dN,dNdz)
c
      include 'aba_param.inc'
c
      parameter (gaussCoord=0.577350269d0)
      dimension dN(1,nnode),dNdz(ndim,*),coord38(3,8)
      
      data  coord38 /-1.d0, -1.d0, -1.d0,
     2                1.d0, -1.d0, -1.d0,
     3               -1.d0,  1.d0, -1.d0,
     4                1.d0,  1.d0, -1.d0,
     5               -1.d0, -1.d0,  1.d0,
     6                1.d0, -1.d0,  1.d0,
     7               -1.d0,  1.d0,  1.d0,
     8                1.d0,  1.d0,  1.d0/	 
!     2D 4-nodes

!     determine (g,h,r)
      g=coord38(1,kintk)*gaussCoord
      h=coord38(2,kintk)*gaussCoord
	  r=coord38(3,kintk)*gaussCoord

!     shape functions 
      dN(1,1)=0.125d0*(1.d0-g)*(1.d0-h)*(1.d0-r)*(-g-h-r-2.d0)
	  dN(1,2)=0.125d0*(1.d0+g)*(1.d0-h)*(1.d0-r)*(g-h-r-2.d0)
	  dN(1,3)=0.125d0*(1.d0+g)*(1.d0+h)*(1.d0-r)*(g+h-r-2.d0)
	  dN(1,4)=0.125d0*(1.d0-g)*(1.d0+h)*(1.d0-r)*(-g+h-r-2.d0)
	  dN(1,5)=0.125d0*(1.d0-g)*(1.d0-h)*(1.d0+r)*(-g-h+r-2.d0)
	  dN(1,6)=0.125d0*(1.d0+g)*(1.d0-h)*(1.d0+r)*(g-h+r-2.d0)
	  dN(1,7)=0.125d0*(1.d0+g)*(1.d0+h)*(1.d0+r)*(g+h+r-2.d0)
	  dN(1,8)=0.125d0*(1.d0-g)*(1.d0+h)*(1.d0+r)*(-g+h+r-2.d0)
	  
	  dN(1,9)=0.25d0*(1.d0-g*g)*(1.d0-h)*(1.d0-r)
	  dN(1,11)=0.25d0*(1.d0-g*g)*(1.d0+h)*(1.d0-r)
	  dN(1,13)=0.25d0*(1.d0-g*g)*(1.d0-h)*(1.d0+r)
	  dN(1,15)=0.25d0*(1.d0-g*g)*(1.d0+h)*(1.d0+r)
	  
	  dN(1,10)=0.25d0*(1.d0-h*h)*(1.d0+g)*(1.d0-r)
	  dN(1,12)=0.25d0*(1.d0-h*h)*(1.d0-g)*(1.d0-r)
	  dN(1,14)=0.25d0*(1.d0-h*h)*(1.d0+g)*(1.d0+r)
	  dN(1,16)=0.25d0*(1.d0-h*h)*(1.d0-g)*(1.d0+r)
	  
	  dN(1,17)=0.25d0*(1.d0-r*r)*(1.d0-g)*(1.d0-h)
	  dN(1,18)=0.25d0*(1.d0-r*r)*(1.d0+g)*(1.d0-h)
	  dN(1,19)=0.25d0*(1.d0-r*r)*(1.d0+g)*(1.d0+h)
	  dN(1,20)=0.25d0*(1.d0-r*r)*(1.d0-g)*(1.d0+h)
	  
	  
!     derivative d(Ni)/d(g)
      dNdz(1,1)=0.125d0*(h-1.d0)*(r-1.d0)*(2.d0*g+r+h+1.d0)
	  dNdz(1,2)=0.125d0*(h-1.d0)*(r-1.d0)*(2.d0*g-r-h-1.d0)
	  dNdz(1,3)=-0.125d0*(h+1.d0)*(r-1.d0)*(2.d0*g-r+h-1.d0)
	  dNdz(1,4)=-0.125d0*(h+1.d0)*(r-1.d0)*(2.d0*g+r-h+1.d0)
	  dNdz(1,5)=-0.125d0*(h-1.d0)*(r+1.d0)*(2.d0*g-r+h+1.d0)
	  dNdz(1,6)=-0.125d0*(h-1.d0)*(r+1.d0)*(2.d0*g+r-h-1.d0)
	  dNdz(1,7)=0.125d0*(h+1.d0)*(r+1.d0)*(2.d0*g+r+h-1.d0)
	  dNdz(1,8)=0.125d0*(h+1.d0)*(r+1.d0)*(2.d0*g-r-h+1.d0)
	  
	  dNdz(1,9)=-0.5d0*(h-1.d0)*(r-1.d0)*g
	  dNdz(1,11)=0.5d0*(h+1.d0)*(r-1.d0)*g
	  dNdz(1,13)=0.5d0*(h-1.d0)*(r+1.d0)*g
	  dNdz(1,15)=-0.5d0*(h+1.d0)*(r+1.d0)*g
	  
	  dNdz(1,10)=0.25d0*(h*h-1.d0)*(r-1.d0)
	  dNdz(1,12)=-0.25d0*(h*h-1.d0)*(r-1.d0)
	  dNdz(1,14)=0.25d0*(1.d0-h*h)*(r+1.d0)
	  dNdz(1,16)=0.25d0*(h*h-1.d0)*(r+1.d0)
	  
	  dNdz(1,17)=-0.25d0*(h-1.d0)*(r*r-1.d0)
	  dNdz(1,18)=0.25d0*(h-1.d0)*(r*r-1.d0)
	  dNdz(1,19)=0.25d0*(h+1.d0)*(1.d0-r*r)
	  dNdz(1,20)=0.25d0*(h+1.d0)*(r*r-1.d0)
    

!     derivative d(Ni)/d(h)
      dNdz(2,1)=0.125d0*(g-1.d0)*(r-1.d0)*(2.d0*h+r+g+1.d0)
	  dNdz(2,2)=-0.125d0*(g+1.d0)*(r-1.d0)*(2.d0*h+r-g+1.d0)
	  dNdz(2,3)=-0.125d0*(g+1.d0)*(r-1.d0)*(2.d0*h-r+g-1.d0)
	  dNdz(2,4)=0.125d0*(g-1.d0)*(r-1.d0)*(2.d0*h-r-g-1.d0)
	  dNdz(2,5)=-0.125d0*(g-1.d0)*(r+1.d0)*(2.d0*h-r+g+1.d0)
	  dNdz(2,6)=0.125d0*(g+1.d0)*(r+1.d0)*(2.d0*h-r-g+1.d0)
	  dNdz(2,7)=0.125d0*(g+1.d0)*(r+1.d0)*(2.d0*h+r+g-1.d0)
	  dNdz(2,8)=-0.125d0*(g-1.d0)*(r+1.d0)*(2.d0*h+r-g-1.d0)
		   
	  dNdz(2,9)=-0.25d0*(g*g-1.d0)*(r-1.d0)
	  dNdz(2,11)=0.25d0*(g*g-1.d0)*(r-1.d0)
	  dNdz(2,13)=0.25d0*(g*g-1.d0)*(r+1.d0)
	  dNdz(2,15)=0.25d0*(1.d0-g*g)*(r+1.d0)
		   
	  dNdz(2,10)=0.5d0*(g+1.d0)*(r-1.d0)*h
	  dNdz(2,12)=-0.5d0*(g-1.d0)*(r-1.d0)*h
	  dNdz(2,14)=-0.5d0*(g+1.d0)*(r+1.d0)*h
	  dNdz(2,16)=0.5d0*(g-1.d0)*(r+1.d0)*h
		   
	  dNdz(2,17)=-0.25d0*(g-1.d0)*(r*r-1.d0)
	  dNdz(2,18)=0.25d0*(g+1.d0)*(r*r-1.d0)
	  dNdz(2,19)=0.25d0*(g+1.d0)*(1.d0-r*r)
	  dNdz(2,20)=0.25d0*(g-1.d0)*(r*r-1.d0)
	  
	  
!		derivative d(Ni)/d(r)
      dNdz(3,1)=0.125d0*(g-1.d0)*(h-1.d0)*(2.d0*r+h+g+1.d0)
	  dNdz(3,2)=-0.125d0*(g+1.d0)*(h-1.d0)*(2.d0*r+h-g+1.d0)
	  dNdz(3,3)=0.125d0*(g+1.d0)*(h+1.d0)*(2.d0*r-h-g+1.d0)
	  dNdz(3,4)=-0.125d0*(g-1.d0)*(h+1.d0)*(2.d0*r-h+g+1.d0)
	  dNdz(3,5)=0.125d0*(g-1.d0)*(h-1.d0)*(2.d0*r-h-g-1.d0)
	  dNdz(3,6)=-0.125d0*(g+1.d0)*(h-1.d0)*(2.d0*r-h+g-1.d0)
	  dNdz(3,7)=0.125d0*(g+1.d0)*(h+1.d0)*(2.d0*r+h+g-1.d0)
	  dNdz(3,8)=-0.125d0*(g-1.d0)*(h+1.d0)*(2.d0*r+h-g-1.d0)
		   
	  dNdz(3,9)=-0.25d0*(g*g-1.d0)*(h-1.d0)
	  dNdz(3,11)=0.25d0*(g*g-1.d0)*(h+1.d0)
	  dNdz(3,13)=0.25d0*(g*g-1.d0)*(h-1.d0)
	  dNdz(3,15)=0.25d0*(1.d0-g*g)*(h+1.d0)
		   
	  dNdz(3,10)=0.25d0*(g+1.d0)*(h*h-1.d0)
	  dNdz(3,12)=-0.25d0*(g-1.d0)*(h*h-1.d0)
	  dNdz(3,14)=0.25d0*(g+1.d0)*(1.d0-h*h)
	  dNdz(3,16)=0.25d0*(g-1.d0)*(h*h-1.d0)
		   
	  dNdz(3,17)=-0.5d0*(g-1.d0)*(h-1.d0)*r
	  dNdz(3,18)=0.5d0*(g+1.d0)*(h-1.d0)*r
	  dNdz(3,19)=-0.5d0*(g+1.d0)*(h+1.d0)*r
	  dNdz(3,20)=0.5d0*(g-1.d0)*(h+1.d0)*r
      
      return
      end 

      subroutine kjacobian(jelem,ndim,nnode,coords,dNdz,djac,dNdx,mcrd,COFACTOR)
!     Notation: djac - Jac determinant; xjaci - inverse of Jac matrix 
!     dNdx - shape functions derivatives w.r.t. global coordinates
      include 'aba_param.inc'

      dimension xjac(ndim,ndim),xjaci(ndim,ndim),coords(mcrd,nnode),
     1 dNdz(ndim,nnode),dNdx(ndim,nnode),COFACTOR(ndim,ndim)

      xjac=0.d0
	  djac_inv=0.d0

      do inod=1,nnode
       do idim=1,ndim
        do jdim=1,ndim
         xjac(jdim,idim)=xjac(jdim,idim)+
     1        dNdz(jdim,inod)*coords(idim,inod)      
        end do
       end do 
      end do

      djac=xjac(1,1)*(xjac(2,2)*xjac(3,3)-xjac(3,2)*xjac(2,3))-
     & xjac(2,1)*(xjac(1,2)*xjac(3,3)-xjac(3,2)*xjac(1,3))+
     & xjac(3,1)*(xjac(1,2)*xjac(2,3)-xjac(2,2)*xjac(1,3))
	  
	  
	  
!      djac=xjac(1,1)*xjac(2,2)*xjac(3,3)-xjac(1,1)*xjac(2,3)*xjac(3,2)-
!     & xjac(1,2)*xjac(2,1)*xjac(3,3)+xjac(1,2)*xjac(2,3)*xjac(3,1)+
!     & xjac(1,3)*xjac(2,1)*xjac(3,2)-xjac(1,3)*xjac(2,2)*xjac(3,1)
      if (djac.gt.0.d0) then ! jacobian is positive - o.k.
	  djac_inv=1.d0/djac
	  
	  xjaci(1,1)=djac_inv*(xjac(2,2)*xjac(3,3)-xjac(3,2)*xjac(2,3))
	  xjaci(1,2)=djac_inv*(xjac(3,2)*xjac(1,3)-xjac(1,2)*xjac(3,3))
	  xjaci(1,3)=djac_inv*(xjac(1,2)*xjac(2,3)-xjac(2,2)*xjac(1,3))
	  xjaci(2,1)=djac_inv*(xjac(3,1)*xjac(2,3)-xjac(2,1)*xjac(3,3))
	  xjaci(2,2)=djac_inv*(xjac(1,1)*xjac(3,3)-xjac(3,1)*xjac(1,3))
	  xjaci(2,3)=djac_inv*(xjac(2,1)*xjac(1,3)-xjac(1,1)*xjac(2,3))
	  xjaci(3,1)=djac_inv*(xjac(2,1)*xjac(3,2)-xjac(3,1)*xjac(2,2))
	  xjaci(3,2)=djac_inv*(xjac(3,1)*xjac(1,2)-xjac(1,1)*xjac(3,2))
	  xjaci(3,3)=djac_inv*(xjac(1,1)*xjac(2,2)-xjac(2,1)*xjac(1,2))
	  
	  
!	  COFACTOR(1,1) = +(xjac(2,2)*xjac(3,3)-xjac(2,3)*xjac(3,2))
!	  COFACTOR(1,2) = -(xjac(2,1)*xjac(3,3)-xjac(2,3)*xjac(3,1))
!	  COFACTOR(1,3) = +(xjac(2,1)*xjac(3,2)-xjac(2,2)*xjac(3,1))
!	  COFACTOR(2,1) = -(xjac(1,2)*xjac(3,3)-xjac(1,3)*xjac(3,2))
!	  COFACTOR(2,2) = +(xjac(1,1)*xjac(3,3)-xjac(1,3)*xjac(3,1))
!	  COFACTOR(2,3) = -(xjac(1,1)*xjac(3,2)-xjac(1,2)*xjac(3,1))
!	  COFACTOR(3,1) = +(xjac(1,2)*xjac(2,3)-xjac(1,3)*xjac(2,2))
!	  COFACTOR(3,2) = -(xjac(1,1)*xjac(2,3)-xjac(1,3)*xjac(2,1))
!	  COFACTOR(3,3) = +(xjac(1,1)*xjac(2,2)-xjac(1,2)*xjac(2,1))
!	  xjaci = transpose(COFACTOR)/djac
      else ! negative or zero jacobian
       write(7,*)'WARNING: element',jelem,'has neg. Jacobian'
      endif
	  
      dNdx=matmul(xjaci,dNdz) 
      return
      end

c*****************************************************************
      subroutine kstatevar(npt,nsvint,statev,statev_ip,icopy)
c
c     Transfer data to/from element-level state variable array from/to
c     material-point level state variable array.
c
      include 'aba_param.inc'

      dimension statev(*),statev_ip(*)

      isvinc=(npt-1)*nsvint     ! integration point increment

      if (icopy.eq.1) then ! Prepare arrays for entry into umat
       do i=1,nsvint
        statev_ip(i)=statev(i+isvinc)
       enddo
      else ! Update element state variables upon return from umat
       do i=1,nsvint
        statev(i+isvinc)=statev_ip(i)
       enddo
      end if

      return
      end

c*****************************************************************
      subroutine kumat(props,ddsdde,stress,dstran,ntens,statev,ndim)
c
c     Subroutine with the material model
c
      include 'aba_param.inc' !implicit real(a-h o-z)
      
      dimension props(*),ddsdde(ntens,ntens),stress(ntens),statev(*),
     + dstran(ntens)

!     Initialization
      ddsdde=0.d0
      E=props(1) ! Young's modulus
      xnu=props(2) ! Poisson's ratio
      K = E/(3.d0*(1-2*xnu))
!     Build stiffness matrix
      EBULK3 = E/(1.d0-2.d0*xnu)
	  EBULK = EBULK3/3.d0
      eg2=E/(1.d0+xnu)
	  eg = eg2/2.d0
      elam=(E/(1.d0-2.d0*xnu)-eg2)/3.d0
	  
!	  In 3D :

      !Building material matrix
      DO K1=1, 3
          DO K2=1, 3
              ddsdde(K2, K1)=elam
          END DO
          ddsdde(K1, K1)=eg2+elam
      END DO
      DO K1=ndim+1, ntens
          ddsdde(K1, K1)=eg
      END DO
	  
      
!     Plane strain
    !  do k1=1,3
    !   do k2=1,3
    !    ddsdde(k2,k1)=elam
    !   end do
    !   ddsdde(k1,k1)=eg2+elam
    !  end do
    !  ddsdde(4,4)=eg2/2.d0
      
!     NOTE: ABAQUS uses engineering shear strains,
!     i.e. stran(ndi+1) = 2*e_12, etc...

	



      stress=stress+matmul(ddsdde,dstran)   

      return
      end

c*****************************************************************
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,rpl,ddsddt,
     1 drplde,drpldt,stran,dstran,time,dtime,temp2,dtemp,predef,dpred,
     2 cmname,ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     3 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,jstep,kinc)

      use kvisual
      include 'aba_param.inc' !implicit real(a-h o-z)

      character*8 cmname
      dimension stress(ntens),statev(nstatv),ddsdde(ntens,ntens),
     1 ddsddt(ntens),drplde(ntens),stran(ntens),dstran(ntens),
     2 time(2),predef(1),dpred(1),props(nprops),coords(3),drot(3,3),
     3 dfgrd0(3,3),dfgrd1(3,3),jstep(4)

      ddsdde=0.0d0
      noffset=noel-nelem
      statev(1:nstatv)=UserVar(noffset,1:nstatv,npt)
     
      return
      end
