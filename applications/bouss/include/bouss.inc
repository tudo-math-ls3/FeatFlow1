      DOUBLE PRECISION VOLEXK
      COMMON /BOUSSD/ ALFA,OMEGAT,VOLEXK,DMPTMP,g(2)
      COMMON /BOUSSI/ LDT,LTOLD,KT,KFT,KAUXT,NUPT,KNUPT(NNLEV),
     *                ILMINT,ILMAXT,IINTT,ISMT,ISLT,NSMT,NSLT
      COMMON /TMPTIM/ TIMTAL,TITLIN,TIMTMG,TIMTPO
C
C     TIMTAL: TIMe Temperature ALl: Gesamtzeit fuer die Temperaturberechnung
C     ueber alle Zeitschritte. 
C
