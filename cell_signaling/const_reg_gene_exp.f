
      program const_reg_gene_exp

c     integrates NEQU first order explicit ODE's using RK4 algorithm
c     The equations to be integrated are to be Defined in subroutine DERIVADES
c     RK4 is implemented in a subroutine and is independent of the specific ODEs
  
      implicit none
      integer NEQU
      parameter (NEQU=4)
c     t: independent variable
c     h: integration step
c     y: vector of NEQU variables at independent variable value t
c     yout: vector of NEQU variables at independent variable value t+h
c     tf: final value of t to integrate (initial value is assumed to be 0)
    
      double precision y(nequ),yout(nequ)
      double precision t,h,tf,E1_ini,E2_ini

c     auxiliary variables for loops and writing 
      integer npas,i,j,nwrite

c     parameters 
      double precision k1,k2,k3,k4,k5,k6,E1,E2

      common /param/ k1,k2,k3,k4,k5,k6,E1,E2
      
      open (unit=20,file="expression.dat",status="unknown")



ccccccccccccccccccc parameter values cccccccccccccccccc

cccc  parameters to be changed   ccccccccc
 
      k1 = 1.d0/150.d0
      k2 = 1.d0
      k3 = 1.d0
      k4 = 1.d0/150.d0
      k5 = 1.d0
      k6 = 1.d0

ccccccccccccccccccc initial conditions ccccccccccccccccc      
      t=0.0d0
      y(1) = 3.d0
      y(2) = 0.d0
      y(3) = 0.d0
      y(4) = 0.d0
      E1_ini = 5.d0
      E2_ini = 0.3d0
      E1 = E1_ini
      E2 = E2_ini
      write(20,200) t,y(1),y(2),y(3),y(4)
cccccccccccccccccc integration over time cccccccccccccccc               
c until final time tf with time step h, saving data every nwrite cccccc      
   
      tf=10000.0d0
      h=0.01d0
      npas=int(tf/h)
      nwrite=50

      do j=1,npas
         
         call RK4(nequ,t,h,y,yout)
      
ccccc    update independent variable
         t=dble(j)*h

ccccc    write data 
         if(mod(j,nwrite).eq.0)then
             write(20,200)t,yout(1),yout(2),yout(3),yout(4)
         endif
     
ccccc    update variables for next step
         do i=1,nequ
            y(i)=yout(i)
         enddo

         E1 = E1_ini - y(2)
         E2 = E2_ini - y(4)
      enddo
    
200   format(8(e20.8,1x))

         
         end 






      subroutine DFUNCTIONS(NEQU,t,y,dy)
      implicit none
      double precision y(nequ),dy(nequ)
      double precision t
      integer nequ

c     parameters 
      double precision k1,k2,k3,k4,k5,k6,E1,E2

      common /param/ k1,k2,k3,k4,k5,k6,E1,E2

        dy(1)= -k1*y(1)*E1+k2*y(2)+k6*y(4)
        dy(2)= k1*(y(1)*E1)-(k2+k3)*y(2)
        dy(3)= -k4*(y(3)*E2)+k5*y(4)+k3*y(2) 
        dy(4)= k4*(y(3)*E2)-(k5+k6)*y(4)

      return
      end



      subroutine RK4(nequ,x,h,y,yout)
      implicit none 
      double precision y(nequ),yout(nequ)
      double precision dy1(nequ),dy2(nequ),dy3(nequ),dy4(nequ)
      double precision x1,x3
      double precision y1(nequ),y2(nequ),y3(nequ)
      double precision x,h
      integer i,nequ

         call DFUNCTIONS(NEQU,x,y,dy1)
         x1=x+h*0.5d0
         y1=y+h*0.5d0*dy1  
         call DFUNCTIONS(NEQU,x1,y1,dy2)
         y2=y+h*0.5d0*dy2  
         call DFUNCTIONS(NEQU,x1,y2,dy3)
         x3=x+h
         y3=y+h*dy3
         call DFUNCTIONS(NEQU,x3,y3,dy4)

         yout=y+h*(dy1+2.0d0*(dy2+dy3)+dy4)/6.0d0
      return
      end


      
