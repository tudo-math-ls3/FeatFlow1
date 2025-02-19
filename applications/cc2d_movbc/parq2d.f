      DOUBLE PRECISION FUNCTION PARX(T,IBCT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(MAXI=23,MAXJ=7,K=2*MAXI+2*(MAXJ-2)) 
      DIMENSION X(1:K)
      DIMENSION TPAR(1:K+1) 
      DATA (X(I),I=1,K)/
     1     0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5,
     1     0.55,0.6, 0.7,  0.8, 0.9,  1.0, 1.1,  1.2, 1.4,  1.6,      
     2     1.8, 2.0, 2.2,    
     2     2.2, 2.2, 2.2, 2.2, 2.2,   
     3     2.2, 2.0, 1.8, 1.6, 1.4, 1.2, 1.1, 1.0, 0.9, 0.8,
     3     0.7, 0.6, 0.55,0.5, 0.45,0.4, 0.35,0.3, 0.25,0.2,
     3     0.15, 0.1, 0.0,
     4     0.0, 0.0, 0.0, 0.0, 0.0/     
     
        DATA (TPAR(J),J=1,K+1)/

     1       0.0,  1.0,  2.0,  3.0,  4.0,  5.0, 6.0, 7.0, 8.0, 9.0,
     1       10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0,
     2       18.0, 19.0, 20.0,       
     3       21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0,        
     3       30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0,       
     4       39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0,
     5       48.0, 49.0, 50.0, 51.0, 52.0, 53.0, 54.0, 55.0,56.0      /

      IF (T.GE.55D0) GOTO 55
      IF (T.GE.54D0) GOTO 54
      IF (T.GE.53D0) GOTO 53
      IF (T.GE.52D0) GOTO 52
      IF (T.GE.51D0) GOTO 51
      IF (T.GE.50D0) GOTO 50
      IF (T.GE.49D0) GOTO 49
      IF (T.GE.48D0) GOTO 48
      IF (T.GE.47D0) GOTO 47
      IF (T.GE.46D0) GOTO 46
      IF (T.GE.45D0) GOTO 45
      IF (T.GE.44D0) GOTO 44
      IF (T.GE.43D0) GOTO 43
      IF (T.GE.42D0) GOTO 42
      IF (T.GE.41D0) GOTO 41
      IF (T.GE.40D0) GOTO 40
      IF (T.GE.39D0) GOTO 39
      IF (T.GE.38D0) GOTO 38
      IF (T.GE.37D0) GOTO 37
      IF (T.GE.36D0) GOTO 36
      IF (T.GE.35D0) GOTO 35
      IF (T.GE.34D0) GOTO 34
      IF (T.GE.33D0) GOTO 33
      IF (T.GE.32D0) GOTO 32
      IF (T.GE.31D0) GOTO 31
      IF (T.GE.30D0) GOTO 30
      IF (T.GE.29D0) GOTO 29 
      IF (T.GE.28D0) GOTO 28
      IF (T.GE.27D0) GOTO 27
      IF (T.GE.26D0) GOTO 26
      IF (T.GE.25D0) GOTO 25
      IF (T.GE.24D0) GOTO 24
      IF (T.GE.23D0) GOTO 23
      IF (T.GE.22D0) GOTO 22
      IF (T.GE.21D0) GOTO 21
      IF (T.GE.20D0) GOTO 20
      IF (T.GE.19D0) GOTO 19
      IF (T.GE.18D0) GOTO 18
      IF (T.GE.17D0) GOTO 17 
      IF (T.GE.16D0) GOTO 16
      IF (T.GE.15D0) GOTO 15
      IF (T.GE.14D0) GOTO 14
      IF (T.GE.13D0) GOTO 13
      IF (T.GE.12D0) GOTO 12
      IF (T.GE.11D0) GOTO 11
      IF (T.GE.10D0) GOTO 10
      IF (T.GE.9D0) GOTO 9
      IF (T.GE.8D0) GOTO 8 
      IF (T.GE.7D0) GOTO 7
      IF (T.GE.6D0) GOTO 6
      IF (T.GE.5D0) GOTO 5
      IF (T.GE.4D0) GOTO 4
      IF (T.GE.3D0) GOTO 3
      IF (T.GE.2D0) GOTO 2
      IF (T.GE.1D0) GOTO 1
   


        PARX=X(1)+(X(2)-X(1))/(TPAR(2)-TPAR(1))*(T-TPAR(1)) 
        GOTO 99999

 1      PARX=X(2)+(X(3)-X(2))/(TPAR(3)-TPAR(2))*(T-TPAR(2)) 
        GOTO 99999

 2      PARX=X(3)+(X(4)-X(3))/(TPAR(4)-TPAR(3))*(T-TPAR(3)) 
        GOTO 99999

 3      PARX=X(4)+(X(5)-X(4))/(TPAR(5)-TPAR(4))*(T-TPAR(4)) 
        GOTO 99999

 4      PARX=X(5)+(X(6)-X(5))/(TPAR(6)-TPAR(5))*(T-TPAR(5)) 
        GOTO 99999

 5      PARX=X(6)+(X(7)-X(6))/(TPAR(7)-TPAR(6))*(T-TPAR(6)) 
        GOTO 99999

 6      PARX=X(7)+(X(8)-X(7))/(TPAR(8)-TPAR(7))*(T-TPAR(7)) 
        GOTO 99999

 7      PARX=X(8)+(X(9)-X(8))/(TPAR(9)-TPAR(8))*(T-TPAR(8)) 
        GOTO 99999

 8      PARX=X(9)+(X(10)-X(9))/(TPAR(10)-TPAR(9))*(T-TPAR(9)) 
        GOTO 99999

 9      PARX=X(10)+(X(11)-X(10))/(TPAR(11)-TPAR(10))*(T-TPAR(10)) 
        GOTO 99999

 10     PARX=X(11)+(X(12)-X(11))/(TPAR(12)-TPAR(11))*(T-TPAR(11)) 
        GOTO 99999

 11     PARX=X(12)+(X(13)-X(12))/(TPAR(13)-TPAR(12))*(T-TPAR(12)) 
        GOTO 99999

 12     PARX=X(13)+(X(14)-X(13))/(TPAR(14)-TPAR(13))*(T-TPAR(13)) 
        GOTO 99999

 13     PARX=X(14)+(X(15)-X(14))/(TPAR(15)-TPAR(14))*(T-TPAR(14)) 
        GOTO 99999

 14     PARX=X(15)+(X(16)-X(15))/(TPAR(16)-TPAR(15))*(T-TPAR(15)) 
        GOTO 99999

 15     PARX=X(16)+(X(17)-X(16))/(TPAR(17)-TPAR(16))*(T-TPAR(16)) 
        GOTO 99999

 16     PARX=X(17)+(X(18)-X(17))/(TPAR(18)-TPAR(17))*(T-TPAR(17)) 
        GOTO 99999

 17     PARX=X(18)+(X(19)-X(18))/(TPAR(19)-TPAR(18))*(T-TPAR(18)) 
        GOTO 99999

 18     PARX=X(19)+(X(20)-X(19))/(TPAR(20)-TPAR(19))*(T-TPAR(19)) 
        GOTO 99999

 19     PARX=X(20)+(X(21)-X(20))/(TPAR(21)-TPAR(20))*(T-TPAR(20)) 
        GOTO 99999

 20     PARX=X(21)+(X(22)-X(21))/(TPAR(22)-TPAR(21))*(T-TPAR(21)) 
        GOTO 99999
 
 21     PARX=X(22)+(X(23)-X(22))/(TPAR(23)-TPAR(22))*(T-TPAR(22)) 
        GOTO 99999

 22     PARX=X(23)+(X(24)-X(23))/(TPAR(24)-TPAR(23))*(T-TPAR(23)) 
        GOTO 99999

 23     PARX=X(24)+(X(25)-X(24))/(TPAR(25)-TPAR(24))*(T-TPAR(24)) 
        GOTO 99999

 24     PARX=X(25)+(X(26)-X(25))/(TPAR(26)-TPAR(25))*(T-TPAR(25)) 
        GOTO 99999

 25     PARX=X(26)+(X(27)-X(26))/(TPAR(27)-TPAR(26))*(T-TPAR(26)) 
        GOTO 99999

 26     PARX=X(27)+(X(28)-X(27))/(TPAR(28)-TPAR(27))*(T-TPAR(27)) 
        GOTO 99999

 27     PARX=X(28)+(X(29)-X(28))/(TPAR(29)-TPAR(28))*(T-TPAR(28)) 
        GOTO 99999

 28     PARX=X(29)+(X(30)-X(29))/(TPAR(30)-TPAR(29))*(T-TPAR(29)) 
        GOTO 99999

 29     PARX=X(30)+(X(31)-X(30))/(TPAR(31)-TPAR(30))*(T-TPAR(30)) 
        GOTO 99999

 30     PARX=X(31)+(X(32)-X(31))/(TPAR(32)-TPAR(31))*(T-TPAR(31)) 
        GOTO 99999

 31     PARX=X(32)+(X(33)-X(32))/(TPAR(33)-TPAR(32))*(T-TPAR(32)) 
        GOTO 99999

 32     PARX=X(33)+(X(34)-X(33))/(TPAR(34)-TPAR(33))*(T-TPAR(33)) 
        GOTO 99999

 33     PARX=X(34)+(X(35)-X(34))/(TPAR(35)-TPAR(34))*(T-TPAR(34)) 
        GOTO 99999

 34     PARX=X(35)+(X(36)-X(35))/(TPAR(36)-TPAR(35))*(T-TPAR(35)) 
        GOTO 99999

 35     PARX=X(36)+(X(37)-X(36))/(TPAR(37)-TPAR(36))*(T-TPAR(36)) 
       GOTO 99999
       
 36    PARX=X(37)+(X(38)-X(37))/(TPAR(38)-TPAR(37))*(T-TPAR(37)) 
        GOTO 99999

 37     PARX=X(38)+(X(39)-X(38))/(TPAR(39)-TPAR(38))*(T-TPAR(38)) 
        GOTO 99999


 38     PARX=X(39)+(X(40)-X(39))/(TPAR(40)-TPAR(39))*(T-TPAR(39)) 
        GOTO 99999


 39     PARX=X(40)+(X(41)-X(40))/(TPAR(41)-TPAR(40))*(T-TPAR(40)) 
        GOTO 99999

 40     PARX=X(41)+(X(42)-X(41))/(TPAR(42)-TPAR(41))*(T-TPAR(41)) 
        GOTO 99999

 41     PARX=X(42)+(X(43)-X(42))/(TPAR(43)-TPAR(42))*(T-TPAR(42)) 
        GOTO 99999

 42     PARX=X(43)+(X(44)-X(43))/(TPAR(44)-TPAR(43))*(T-TPAR(43)) 
        GOTO 99999

 43     PARX=X(44)+(X(45)-X(44))/(TPAR(45)-TPAR(44))*(T-TPAR(44)) 
        GOTO 99999

 44     PARX=X(45)+(X(46)-X(45))/(TPAR(46)-TPAR(45))*(T-TPAR(45)) 
        GOTO 99999

 45     PARX=X(46)+(X(47)-X(46))/(TPAR(47)-TPAR(46))*(T-TPAR(46)) 
        GOTO 99999

 46     PARX=X(47)+(X(48)-X(47))/(TPAR(48)-TPAR(47))*(T-TPAR(47)) 
        GOTO 99999 
 
 47     PARX=X(48)+(X(49)-X(48))/(TPAR(49)-TPAR(48))*(T-TPAR(48)) 
        GOTO 99999
  
 48     PARX=X(49)+(X(50)-X(49))/(TPAR(50)-TPAR(49))*(T-TPAR(49)) 
        GOTO 99999  

 49     PARX=X(50)+(X(51)-X(50))/(TPAR(51)-TPAR(50))*(T-TPAR(50)) 
        GOTO 99999
  
 50     PARX=X(51)+(X(52)-X(51))/(TPAR(52)-TPAR(51))*(T-TPAR(51)) 
        GOTO 99999
  
 51     PARX=X(52)+(X(53)-X(52))/(TPAR(53)-TPAR(52))*(T-TPAR(52)) 
        GOTO 99999 
 
 52     PARX=X(53)+(X(54)-X(53))/(TPAR(54)-TPAR(53))*(T-TPAR(53)) 
        GOTO 99999 
 
 53     PARX=X(54)+(X(55)-X(54))/(TPAR(55)-TPAR(54))*(T-TPAR(54)) 
        GOTO 99999 
 
 54     PARX=X(55)+(X(56)-X(55))/(TPAR(56)-TPAR(55))*(T-TPAR(55)) 
        GOTO 99999 
 
 55     PARX=X(56)+(X(1)-X(56))/(TPAR(57)-TPAR(56))*(T-TPAR(56)) 
        GOTO 99999  


C
99999 END
C

       DOUBLE PRECISION FUNCTION PARY(T,IBCT)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       PARAMETER(MAXI=23,MAXJ=7,K=2*MAXI+2*(MAXJ-2)) 
       DIMENSION Y(1:K)
       DIMENSION TPAR(1:K+1) 
        DATA (Y(J),J=1,K)/

     1       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     1       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     1       0.0, 0.0, 0.0, 
     2       0.1, 0.15,0.2, 0.25,0.3,  
     3       0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,
     3       0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,
     3       0.41,0.41,0.41,
     4       0.3, 0.25,0.2, 0.15,0.1/

    
        DATA (TPAR(J),J=1,K+1)/

     1       0.0,  1.0,  2.0,  3.0,  4.0,  5.0, 6.0, 7.0, 8.0, 9.0,
     1       10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0,
     2       18.0, 19.0, 20.0,       
     3       21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0,        
     3       30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0,       
     4       39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0,
     5       48.0, 49.0, 50.0, 51.0, 52.0, 53.0, 54.0, 55.0,56.0      /
 
      IF (T.GE.55D0) GOTO 55
      IF (T.GE.54D0) GOTO 54
      IF (T.GE.53D0) GOTO 53
      IF (T.GE.52D0) GOTO 52
      IF (T.GE.51D0) GOTO 51
      IF (T.GE.50D0) GOTO 50
      IF (T.GE.49D0) GOTO 49
      IF (T.GE.48D0) GOTO 48
      IF (T.GE.47D0) GOTO 47
      IF (T.GE.46D0) GOTO 46
      IF (T.GE.45D0) GOTO 45
      IF (T.GE.44D0) GOTO 44
      IF (T.GE.43D0) GOTO 43
      IF (T.GE.42D0) GOTO 42
      IF (T.GE.41D0) GOTO 41
      IF (T.GE.40D0) GOTO 40
      IF (T.GE.39D0) GOTO 39
      IF (T.GE.38D0) GOTO 38
      IF (T.GE.37D0) GOTO 37
      IF (T.GE.36D0) GOTO 36
      IF (T.GE.35D0) GOTO 35
      IF (T.GE.34D0) GOTO 34
      IF (T.GE.33D0) GOTO 33
      IF (T.GE.32D0) GOTO 32
      IF (T.GE.31D0) GOTO 31
      IF (T.GE.30D0) GOTO 30
      IF (T.GE.29D0) GOTO 29 
      IF (T.GE.28D0) GOTO 28
      IF (T.GE.27D0) GOTO 27
      IF (T.GE.26D0) GOTO 26
      IF (T.GE.25D0) GOTO 25
      IF (T.GE.24D0) GOTO 24
      IF (T.GE.23D0) GOTO 23
      IF (T.GE.22D0) GOTO 22
      IF (T.GE.21D0) GOTO 21
      IF (T.GE.20D0) GOTO 20
      IF (T.GE.19D0) GOTO 19
      IF (T.GE.18D0) GOTO 18
      IF (T.GE.17D0) GOTO 17 
      IF (T.GE.16D0) GOTO 16
      IF (T.GE.15D0) GOTO 15
      IF (T.GE.14D0) GOTO 14
      IF (T.GE.13D0) GOTO 13
      IF (T.GE.12D0) GOTO 12
      IF (T.GE.11D0) GOTO 11
      IF (T.GE.10D0) GOTO 10
      IF (T.GE.9D0) GOTO 9
      IF (T.GE.8D0) GOTO 8 
      IF (T.GE.7D0) GOTO 7
      IF (T.GE.6D0) GOTO 6
      IF (T.GE.5D0) GOTO 5
      IF (T.GE.4D0) GOTO 4
      IF (T.GE.3D0) GOTO 3
      IF (T.GE.2D0) GOTO 2
      IF (T.GE.1D0) GOTO 1
   
        PARY=Y(1)+(Y(2)-Y(1))/(TPAR(2)-TPAR(1))*(T-TPAR(1)) 
        GOTO 99999

 1      PARY=Y(2)+(Y(3)-Y(2))/(TPAR(3)-TPAR(2))*(T-TPAR(2)) 
        GOTO 99999

 2      PARY=Y(3)+(Y(4)-Y(3))/(TPAR(4)-TPAR(3))*(T-TPAR(3)) 
        GOTO 99999

 3      PARY=Y(4)+(Y(5)-Y(4))/(TPAR(5)-TPAR(4))*(T-TPAR(4)) 
        GOTO 99999

 4      PARY=Y(5)+(Y(6)-Y(5))/(TPAR(6)-TPAR(5))*(T-TPAR(5)) 
        GOTO 99999

 5      PARY=Y(6)+(Y(7)-Y(6))/(TPAR(7)-TPAR(6))*(T-TPAR(6)) 
        GOTO 99999

 6      PARY=Y(7)+(Y(8)-Y(7))/(TPAR(8)-TPAR(7))*(T-TPAR(7)) 
        GOTO 99999

 7      PARY=Y(8)+(Y(9)-Y(8))/(TPAR(9)-TPAR(8))*(T-TPAR(8)) 
        GOTO 99999

 8      PARY=Y(9)+(Y(10)-Y(9))/(TPAR(10)-TPAR(9))*(T-TPAR(9)) 
        GOTO 99999

 9      PARY=Y(10)+(Y(11)-Y(10))/(TPAR(11)-TPAR(10))*(T-TPAR(10)) 
        GOTO 99999

 10     PARY=Y(11)+(Y(12)-Y(11))/(TPAR(12)-TPAR(11))*(T-TPAR(11)) 
        GOTO 99999

 11     PARY=Y(12)+(Y(13)-Y(12))/(TPAR(13)-TPAR(12))*(T-TPAR(12)) 
        GOTO 99999

 12     PARY=Y(13)+(Y(14)-Y(13))/(TPAR(14)-TPAR(13))*(T-TPAR(13)) 
        GOTO 99999

 13     PARY=Y(14)+(Y(15)-Y(14))/(TPAR(15)-TPAR(14))*(T-TPAR(14)) 
        GOTO 99999

 14     PARY=Y(15)+(Y(16)-Y(15))/(TPAR(16)-TPAR(15))*(T-TPAR(15)) 
        GOTO 99999

 15     PARY=Y(16)+(Y(17)-Y(16))/(TPAR(17)-TPAR(16))*(T-TPAR(16)) 
        GOTO 99999

 16     PARY=Y(17)+(Y(18)-Y(17))/(TPAR(18)-TPAR(17))*(T-TPAR(17)) 
        GOTO 99999

 17     PARY=Y(18)+(Y(19)-Y(18))/(TPAR(19)-TPAR(18))*(T-TPAR(18)) 
        GOTO 99999

 18     PARY=Y(19)+(Y(20)-Y(19))/(TPAR(20)-TPAR(19))*(T-TPAR(19)) 
        GOTO 99999

 19     PARY=Y(20)+(Y(21)-Y(20))/(TPAR(21)-TPAR(20))*(T-TPAR(20)) 
        GOTO 99999

 20     PARY=Y(21)+(Y(22)-Y(21))/(TPAR(22)-TPAR(21))*(T-TPAR(21)) 
        GOTO 99999
 
 21     PARY=Y(22)+(Y(23)-Y(22))/(TPAR(23)-TPAR(22))*(T-TPAR(22)) 
        GOTO 99999

 22     PARY=Y(23)+(Y(24)-Y(23))/(TPAR(24)-TPAR(23))*(T-TPAR(23)) 
        GOTO 99999

 23     PARY=Y(24)+(Y(25)-Y(24))/(TPAR(25)-TPAR(24))*(T-TPAR(24)) 
        GOTO 99999

 24     PARY=Y(25)+(Y(26)-Y(25))/(TPAR(26)-TPAR(25))*(T-TPAR(25)) 
        GOTO 99999

 25     PARY=Y(26)+(Y(27)-Y(26))/(TPAR(27)-TPAR(26))*(T-TPAR(26)) 
        GOTO 99999

 26     PARY=Y(27)+(Y(28)-Y(27))/(TPAR(28)-TPAR(27))*(T-TPAR(27)) 
        GOTO 99999

 27     PARY=Y(28)+(Y(29)-Y(28))/(TPAR(29)-TPAR(28))*(T-TPAR(28)) 
        GOTO 99999

 28     PARY=Y(29)+(Y(30)-Y(29))/(TPAR(30)-TPAR(29))*(T-TPAR(29)) 
        GOTO 99999

 29     PARY=Y(30)+(Y(31)-Y(30))/(TPAR(31)-TPAR(30))*(T-TPAR(30)) 
        GOTO 99999

 30     PARY=Y(31)+(Y(32)-Y(31))/(TPAR(32)-TPAR(31))*(T-TPAR(31)) 
        GOTO 99999

 31     PARY=Y(32)+(Y(33)-Y(32))/(TPAR(33)-TPAR(32))*(T-TPAR(32)) 
        GOTO 99999

 32     PARY=Y(33)+(Y(34)-Y(33))/(TPAR(34)-TPAR(33))*(T-TPAR(33)) 
        GOTO 99999

 33     PARY=Y(34)+(Y(35)-Y(34))/(TPAR(35)-TPAR(34))*(T-TPAR(34)) 
        GOTO 99999

 34     PARY=Y(35)+(Y(36)-Y(35))/(TPAR(36)-TPAR(35))*(T-TPAR(35)) 
        GOTO 99999

 35     PARY=Y(36)+(Y(37)-Y(36))/(TPAR(37)-TPAR(36))*(T-TPAR(36)) 
        GOTO 99999

 36     PARY=Y(37)+(Y(38)-Y(37))/(TPAR(38)-TPAR(37))*(T-TPAR(37)) 
        GOTO 99999

 37     PARY=Y(38)+(Y(39)-Y(38))/(TPAR(39)-TPAR(38))*(T-TPAR(38)) 
        GOTO 99999

 38     PARY=Y(39)+(Y(40)-Y(39))/(TPAR(40)-TPAR(39))*(T-TPAR(39)) 
        GOTO 99999

 39     PARY=Y(40)+(Y(41)-Y(40))/(TPAR(41)-TPAR(40))*(T-TPAR(40)) 
        GOTO 99999

 40     PARY=Y(41)+(Y(42)-Y(41))/(TPAR(42)-TPAR(41))*(T-TPAR(41)) 
        GOTO 99999

 41     PARY=Y(42)+(Y(43)-Y(42))/(TPAR(43)-TPAR(42))*(T-TPAR(42)) 
        GOTO 99999

 42     PARY=Y(43)+(Y(44)-Y(43))/(TPAR(44)-TPAR(43))*(T-TPAR(43)) 
        GOTO 99999

 43     PARY=Y(44)+(Y(45)-Y(44))/(TPAR(45)-TPAR(44))*(T-TPAR(44)) 
        GOTO 99999

 44     PARY=Y(45)+(Y(46)-Y(45))/(TPAR(46)-TPAR(45))*(T-TPAR(45)) 
        GOTO 99999

 45     PARY=Y(46)+(Y(47)-Y(46))/(TPAR(47)-TPAR(46))*(T-TPAR(46)) 
        GOTO 99999

 46     PARY=Y(47)+(Y(48)-Y(47))/(TPAR(48)-TPAR(47))*(T-TPAR(47)) 
        GOTO 99999 
 
 47     PARY=Y(48)+(Y(49)-Y(48))/(TPAR(49)-TPAR(48))*(T-TPAR(48)) 
        GOTO 99999
  
 48     PARY=Y(49)+(Y(50)-Y(49))/(TPAR(50)-TPAR(49))*(T-TPAR(49)) 
        GOTO 99999  

 49     PARY=Y(50)+(Y(51)-Y(50))/(TPAR(51)-TPAR(50))*(T-TPAR(50)) 
        GOTO 99999
  
 50     PARY=Y(51)+(Y(52)-Y(51))/(TPAR(52)-TPAR(51))*(T-TPAR(51)) 
        GOTO 99999
  
 51     PARY=Y(52)+(Y(53)-Y(52))/(TPAR(53)-TPAR(52))*(T-TPAR(52)) 
        GOTO 99999 
 
 52     PARY=Y(53)+(Y(54)-Y(53))/(TPAR(54)-TPAR(53))*(T-TPAR(53)) 
        GOTO 99999 
 
 53     PARY=Y(54)+(Y(55)-Y(54))/(TPAR(55)-TPAR(54))*(T-TPAR(54)) 
        GOTO 99999 
 
 54     PARY=Y(55)+(Y(56)-Y(55))/(TPAR(56)-TPAR(55))*(T-TPAR(55)) 
        GOTO 99999 
 
 55     PARY=Y(56)+(Y(1)-Y(56))/(TPAR(57)-TPAR(56))*(T-TPAR(56)) 
        GOTO 99999  
C
99999 END
C
      DOUBLE PRECISION FUNCTION TMAX(IBCT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       TMAX=56D0

      END
 











