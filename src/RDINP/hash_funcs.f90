!
!       hash_funcs.f -- a library of hash table management routines
!
!                                      by
!
!                              Herbert J. Bernstein
!                                Bernstein + Sons
!                    5 Brewster Lane, Bellport, NY 11713-0177, USA
!                       email: yaya@bernstein-plus-sons.com
!
!       work on these routines done in part at Brookhaven National
!       Laboratory, under contract to the U.S. Department of Energy
!
!       work on these routines as been done in part at Dowling College 
!       under contract to the International Union of Crystallography 
!       and under grants from the National Science Foundation and 
!       the U.S. Department of Energy.
!
!       Copyright (C) Herbert J. Bernstein 1997 -- 2006
!
!       hash_funcs.f is free software; you can redistribute this software 
!       and/or modify this software under the terms of the GNU General 
!       Public License as published by the Free Software Foundation; 
!       either version 2 of the License, or (at your option) any later version.
!
!       Alternatively you may reditribute and/or modify hash_funcs.f
!       under the terms of the GNU Lesser General Public 
!       License as published by the Free Software Foundation; either 
!       version 2.1 of the License, or (at your option) any later version.
!
!       This software is distributed in the hope that it will be useful,
!       but WITHOUT ANY WARRANTY; without even the implied warranty of
!       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!       GNU General Public License for more details.
!
!       You should have received a copy of the GNU General Public License
!       along with this software; if not, write to the Free Software
!       Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!       You should have received a copy of the GNU Lesser General Public License
!       along with this software; if not, write to the Free Software
!       Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!-------------------------------------------------------------------------------
!
!       Routines
!
!       hash_init          Initializes a hash table controlled list
!                          call hash_init(data_structure_args)
!
!       hash_find          Searches for a string in a list
!                          call hash_find(name,data_structure_args,ifind)
!
!       hash_fnext         Searches for the next matching string in a list
!                          call hash_fnext(name,data_structure_args,ifind,icurr)
!
!       hash_store         Inserts a new string in a list
!                          call hash_store(name,data_structure_args,ifind)
!
!       hash_snext         Inserts a string in a list, allowing duplicates
!                          call hash_store(name,data_structure_args,ifind,icurr)
!
!       hash_value         Integer function returns index into hash_list
!                          ih = hash_value(name,hash_length)
!
!       The necessary data_structure_args for these routines are
!          name_list   -- an array of character strings
!                         character*(*) name_list(list_length)
!          chain_list  -- chain pointers for searches
!                         integer chain_list(list_length)
!          list_length -- the size of the list arrays
!                         integer list_length
!          num_list    -- number of entries in the list
!                         integer num_list
!          hash_table  -- the initial hashed pointers
!                         integer hash_table
!          hash_length -- the size of the hash table
!                         integer hash_length
!
!
!       The three remaining arguments are
!          name        -- string to search for
!                         character*(*) name
!          ifind       -- return value, 0 for not found (hash_find)
!                         or list full (hash_store), otherwise
!                         the index in name_list of the entry
!          icurr       -- the prior matching index
!
!       The relationship among the arrays used is:
!
!       hash_table is an array (preferably of a modest prime
!       dimension) which starts containing all zeros, which are
!       replaced by pointers to entries in name_list, based
!       values returned by hash_value ranging from 1 to hash_length.
!       Each name is placed in name_list.  A initial zero is placed
!       in the matching entry in chain_list, when the first entry
!       is made.  When a new entry with the same hash_value must be
!       placed a pointer is inserted into chain_list to hook the
!       values together.
!
        subroutine hash_init(name_list,chain_list,list_length,num_list, &
     &                       hash_table,hash_length)
!
!       initialization routine for a hash table controlled list
!          name_list   -- a list of character strings
!          chain_list  -- chain pointers for searches
!          list_length -- the size of the list arrays
!          num_list    -- number of entries in the list
!          hash_table  -- the initial hashed pointers
!          hash_length -- the size of the hash table
!
           integer list_length
           character*(*) name_list(list_length)
           integer hash_length,num_list,i
           integer chain_list(list_length)
           integer hash_table(hash_length)
           num_list=0
           do i = 1,hash_length
           hash_table(i)=0
           enddo
           return
           end          

        subroutine                                                      &
     &  hash_find(name,name_list,chain_list,list_length,num_list,       &
     &                       hash_table,hash_length,ifind)
!
!       search routine for a hash table controlled list
!          name        -- string to find
!          name_list   -- a list of character strings
!          chain_list  -- chain pointers for searches
!          list_length -- the size of the list arrays
!          num_list    -- number of entries in the list
!          hash_table  -- the initial hashed pointers
!          hash_length -- the size of the hash table
!          ifind       -- returned index or 0
!
           character*(*) name
           integer list_length, hash_length
           integer lenn
           character*(*) name_list(list_length)
           integer chain_list(list_length)
           integer hash_table(hash_length)
           integer hash_value
           integer ifind,num_list,ih,ip
           integer lastnb
           ifind=0
           ih=hash_value(name,hash_length)
           ip=hash_table(ih)
           lenn = lastnb(name)
 100       if (ip.eq.0) return
           if (name_list(ip).eq.name(1:lenn)) then
             ifind=ip
             return
           else
             ip=chain_list(ip)
             go to 100
           endif
           end
           
        subroutine                                                      &
     &  hash_fnext(name,name_list,chain_list,list_length,num_list,      &
     &                       hash_table,hash_length,ifind,icurr)
!
!       search routine for a hash table controlled list
!          name        -- string to find
!          name_list   -- a list of character strings
!          chain_list  -- chain pointers for searches
!          list_length -- the size of the list arrays
!          num_list    -- number of entries in the list
!          hash_table  -- the initial hashed pointers
!          hash_length -- the size of the hash table
!          ifind       -- returned index or 0
!          icurr       -- current match or 0
!
           character*(*) name
           integer hash_length
           integer list_length
           integer lenn
           character*(*) name_list(list_length)
           integer chain_list(list_length)
           integer hash_table(hash_length)
           integer hash_value
           integer ifind,num_list,ih,ip
           integer lastnb
           integer icurr
           
           ifind=0
           if (icurr.eq.0) then
             ih = hash_value(name,hash_length)
             ip = hash_table(ih)
           else
             ip = chain_list(icurr)
           endif
           lenn = lastnb(name)
 100       if (ip.eq.0) return
           if (name_list(ip).eq.name(1:lenn)) then
             ifind=ip
             return
           else
             ip=chain_list(ip)
             go to 100
           endif
           end

        subroutine                                                      &
     &  hash_store(name,name_list,chain_list,list_length,num_list,      &
     &                       hash_table,hash_length,ifind)
!
!       store routine for a hash table controlled list
!          name        -- string to find
!          name_list   -- a list of character strings
!          chain_list  -- chain pointers for searches
!          list_length -- the size of the list arrays
!          num_list    -- number of entries in list
!          hash_table  -- the initial hashed pointers
!          hash_length -- the size of the hash table
!          ifind       -- index of entry or 0 (table full)
!
           integer list_length
           character*(*) name
           character*(*) name_list(list_length)
           integer hash_length
           integer lenn
           integer chain_list(list_length)
           integer hash_table(hash_length)
           integer hash_value
           integer ifind,num_list,ih,ip,iq
           integer lastnb

           ifind=0
           ih = hash_value(name,hash_length)
           ip=hash_table(ih)
           iq=0
           lenn = lastnb(name)
 100       if (ip.eq.0) go to 200
           if (name_list(ip).eq.name(1:lenn)) then
             ifind=ip
             return
           else
             iq=ip
             ip=chain_list(ip)
             go to 100
           endif
 200       if (num_list.lt.list_length) then
             num_list=num_list+1
             name_list(num_list)=name(1:lenn)
             chain_list(num_list)=0
             if (iq.eq.0) then
               hash_table(ih)=num_list
             else
               chain_list(iq)=num_list
             endif
             ifind=num_list
             return
           else
             ifind = 0
             return
           endif
           end


        subroutine                                                      &
     &  hash_snext(name,name_list,chain_list,list_length,num_list,      &
     &                       hash_table,hash_length,ifind,icurr)
!
!       store routine for a hash table controlled list
!          name        -- string to find
!          name_list   -- a list of character strings
!          chain_list  -- chain pointers for searches
!          list_length -- the size of the list arrays
!          num_list    -- number of entries in list
!          hash_table  -- the initial hashed pointers
!          hash_length -- the size of the hash table
!          ifind       -- index of entry or 0 (table full)
!          icurr       -- current match or 0
!
           integer list_length
           character*(*) name
           character*(*) name_list(list_length)
           integer hash_length
           integer lenn
           integer chain_list(list_length)
           integer hash_table(hash_length)
           integer hash_value
           integer ifind,num_list,ih,ip,iq
           integer lastnb
           integer icurr

           ifind=0
           ih = 0
           if (icurr.eq.0) then
             ih = hash_value(name,hash_length)
             ip=hash_table(ih)
             iq = 0
           else
             ip = chain_list(icurr)
             iq = icurr
           endif
           lenn = lastnb(name)
 100       if (ip.eq.0) go to 200
           if (name_list(ip).eq.name(1:lenn)) then
             ifind=ip
             return
           else
             iq=ip
             ip=chain_list(ip)
             go to 100
           endif
 200       if (num_list.lt.list_length) then
             num_list=num_list+1
             name_list(num_list)=name(1:lenn)
             chain_list(num_list)=0
             if (iq.eq.0) then
               hash_table(ih)=num_list
             else
               chain_list(iq)=num_list
             endif
             ifind=num_list
             return
           else
             ifind = 0
             return
           endif
           end

      integer function hash_value(name,hash_length)
!
!     function to return a hash value of string name to fit
!     a hash table of length hash_length
      character*(*) name
      integer lastnb
      
      integer hash_length,id,ii,i,ic,lenn
      lenn = lastnb(name)
      hash_value=1
      id = 0
      do ii = 1,lenn
        i = 1+lenn-ii
        ic = ichar(name(i:i))
        if (ic.ge.65) then
          hash_value=mod(hash_value*(ic-64),hash_length)+1
          id = id+1
          if (id.gt.3) return
        endif
      enddo
      return
      end
        
