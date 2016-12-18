C**********************************************
C    Umfpack method   2012年03月13日 星期二 17时57分08秒 
C**********************************************
      subroutine umfdglr(nm,f,ida,jda,da,x,Ax,NAp,NAi,
     +             index,kmax,irep)
      implicit none
      integer nm,irep,kmax
      double precision f(nm),da(kmax),Ax(kmax),x(nm),start,finish
      integer ida(kmax),jda(kmax)
      integer NAi(kmax), NAp(nm+1), index(nm)
C     irep=1,同样的方程组系数矩阵，已求解过，矩阵LU分解信息保存在文件 s64.unf 或 s32.umf 中 
C     nm  总的自由度数（未知量数），f 右端向量， da 包含所有单刚信息
C     da 第一行储存的方程组矩阵行号，第二行储存的方程组矩阵列号， 第三行对应行列的值
C     kmax 所有保存的单刚信息数     
C     非0元素的计算先列后行
C!!!  da 中元素的先后顺序:同一行元素顺序无关，同一列元素一定要按行号从小到大的顺序在da中出现，当然无需连续出现
C     这是因为矩阵按向量存储是先列后行，且每一列第一个非0元素已在NAp中确定，所以列的先后顺序在Ax中一定满足，但同一列
C     中行的顺序必须事先满足，否则无法调用 UMFPACK;已经加了排序的代码，测试通过 2011年09月07日 星期三 08时25分26秒 

      integer symbolic,numeric,filenum,status,sys
      double precision Axtemp(1024),control(20),info(90)
      integer NAitemp(1024),k,i,j,inde,jcol,jNAp,nflag,ktemp

C      print*,nm,kmax,irep,size(control)
      do k=1,nm
       index(k)=0
       x(k)=0.d0    !解赋0值
      end do
      if(irep.eq.1) goto 100  !减少重复计算

      do k=1,kmax
       Ax(k)=0.0d0
       NAi(k)=0
       j=jda(k) ! da 的第二行自由度编号，即总刚的第j列
       index(j)=index(j)+1  ! index(j) 总刚的第j列非0元素个数，用于生成 NAp
      enddo

C-------generate Ap()---------------------  
      NAp(1)=0
      do j=1,nm
       if(index(j).lt.1.or.index(j).gt.nm)print*,index(j),j
       NAp(j+1)=NAp(j)+index(j) ! index(j) 总刚的第j列非0元素个数，NAp（j+1）第j（包括j）列之前的非0元素个数
       index(j)=0  !清0留它用
      enddo

C-------generate Ax()---Ai----------------  
C
       do k=1,kmax
        i=ida(k)  !行标
        j=jda(k)  !列标
        index(j)=index(j)+1   !第j列非0元素个数记录
        ktemp=NAp(j)+index(j)  !确定当前元素为第ktemp个非0元素
        Ax(ktemp)=da(k) !刚度矩阵向量存储方式
        NAi(ktemp)=i-1  !记录行标，注意标号从0开始，即第一行为标号0
        if(i.lt.1.or.i.gt.nm)print*,i
        if(j.lt.1.or.j.gt.nm)print*,j
       end do
        inde=NAp(nm+1) !非0元素总数
              call cpu_time(start)
c       set default parameters
        symbolic=0
        call umf4def (control)

c       print control parameters.  set control (1) to 1 to print
c       error messages only
c        control (1) = 2
c        call umf4pcon (control)

c       pre-order and symbolic analysis
        call umf4sym(nm,nm,NAp,NAi,Ax,symbolic,control,info)
c       check umf4sym error condition
        if (info (1) .lt. 0) then
            print *, 'Error occurred in umf4sym: ', info (1)
            stop
        endif
	filenum = 64
C        write(*,*),nm,inde
        call umf4ssym (symbolic, filenum, status)
        if (status .lt. 0) then
            print *, 'Error occurred in umf4ssym: ', status
            stop
        endif
        goto 200

c       load the symbolic factorization back in (filename: s32.umf)
100        filenum=64
        call umf4lsym (symbolic, filenum, status)
        if (status .lt. 0) then
            print *, 'Error occurred in umf4lsym: ', status
            stop
        endif

c       numeric factorization

200        call umf4num (NAp,NAi,Ax,symbolic,numeric,control,info)
c       check umf4num error condition
        if (info (1) .lt. 0) then
            print *, 'Error occurred in umf4num: ', info (1)
            stop
        endif
c       free the symbolic analysis
        call umf4fsym (symbolic)


c       solve Ax=b, without iterative refinement
        sys = 0
        call umf4sol (sys, x, f, numeric,control,info)
        if (info (1) .lt. 0) then
            print *, 'Error occurred in umf4sol: ', info (1)
            stop
        endif

c       solve Ax=b, with iterative refinement
c        sys = 0
c        call umf4solr (sys, NAp,NAi,Ax, x, f, numeric,control,info)
c        if (info (1) .lt. 0) then
c            print *, 'Error occurred in umf4solr: ', info (1)
c            stop
c        endif

c       free the numeric factorization
        call umf4fnum (numeric)
      call cpu_time(finish)
      print '(" Time of solving the algebra system = ",
     +          f6.3," seconds.")',finish-start
        return
       end
