! * libdescrip - Vertical grid descriptor library for FORTRAN programming
! * Copyright (C) 2016  Direction du developpement des previsions nationales
! *                     Centre meteorologique canadien
! *
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation,
! * version 2.1 of the License.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! *
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library; if not, write to the
! * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! * Boston, MA 02111-1307, USA.
program print_vcode
   
   use vGrid_Descriptors, only:vgd_print,VGD_ERROR
   use Unit_Testing, only: ut_report   

   implicit none
   
   logical :: ok=.true.

   ! Print all known Vcode
   if(vgd_print(-1).eq.VGD_ERROR)ok=.false.

   ! print specific Vcode
   if(vgd_print(1001).eq.VGD_ERROR)ok=.false.

   if(vgd_print(1002).eq.VGD_ERROR)ok=.false.

   if(vgd_print(1003).eq.VGD_ERROR)ok=.false.

   if(vgd_print(2001).eq.VGD_ERROR)ok=.false.

   if(vgd_print(5001).eq.VGD_ERROR)ok=.false.

   if(vgd_print(5002).eq.VGD_ERROR)ok=.false.

   if(vgd_print(5003).eq.VGD_ERROR)ok=.false.

   if(vgd_print(5004).eq.VGD_ERROR)ok=.false.

   if(vgd_print(5005).eq.VGD_ERROR)ok=.false.

   if(vgd_print(5100).eq.VGD_ERROR)ok=.false.

   if(vgd_print(5999).eq.VGD_ERROR)ok=.false.

   if(vgd_print(21001).eq.VGD_ERROR)ok=.false.

   if(vgd_print(21002).eq.VGD_ERROR)ok=.false.

   call ut_report(ok,'print_vcode')

end program print_vcode
