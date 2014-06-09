      parameter(MAX_BUFFER=80,MAX_DATA=12)
      parameter(PARSER_ERROR=-2, QUITTING=-1, RUN=0, ATOM_INPUT=1,
     1          TEMPERATURE_COEFF=2, BOX_COEFF=3, TIME_COEFF=4,  
     2          ANDERSEN_COEFF=5, VELOCITY_SCALE_COEFF=6, 
     3          CONFIG_OUTFILE=7, COMMENT=8, SET_SEED=9,
     4          RESTART_OUTFILE=10, THERMO_OUTFILE=11,
     5          GROW_CRYSTAL=12)

      character*80 string_buffer
      double precision double_buffer(MAX_DATA)
      integer int_buffer
      
      common /str_buf/ string_buffer 
      common /int_buf/ int_buffer
      common /double_buf/  double_buffer
