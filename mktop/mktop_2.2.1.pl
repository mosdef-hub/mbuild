#!/usr/bin/perl

$gromacs_dir="/usr/local/gromacs/share/gromacs/top";


## The force field parameters must be in directory $gromacs_dir (usually /usr/local/gromacs/share/gromacs/top)


###########################################################################################################################
###############################       MKTOP - Topology Builder       ######################################################
######## This is free software, distributed under the GNU license.


## This program was developed by Andre A. S. T. Ribeiro*, Bruno A. C. Horta and Ricardo B. de Alencastro, at the Federal University of Rio de Janeiro, Brazil

#contact email: aastr@iq.ufrj.br (Andre)

## Cite this work as:
# Ribeiro, A.A.S.T.; Horta, B.A.C.; de Alencastro, R.B.  J. Braz. Chem. Soc., Vol. 19, No. 7, 1433-1435, 2008.

# Brief version history:

# 2.1 Fixed bug with improper dihedrals (Thanks to Pramod Akula-Bala)
# 2.0 Amber ff included
# 1.0 OPLS-AA only


##Default bond lengths

## These are essential for the correct bond list determination

$length[1][6]=1.1;  ## H-C bond  (angstroms!!!)
$length[1][7]=1.0;
$length[1][16]=1.2;
$length[6][6]=1.4;
$length[6][7]=1.3;
$length[6][8]=1.4;
$length[6][9]=1.3;
$length[6][15]=1.8;
$length[6][17]=1.8;
$length[6][36]=1.9;
$length[6][53]=2.1;
$length[7][7]=1.3;
$length[7][8]=1.3;
$length[8][1]=1.0;
$length[8][15]=1.6;
$length[8][17]=1.8;
$length[16][16]=2.0;


$deltab=0.3;  # atoms are considered bonded if the distance is within $deltab of the $length value

$deltaX[1] = 0.2;  ## same as before, but now only for bonds containing H (Z=1). Fell free to expand this. 

## Feel free to change these values if you have problems with the bond list (check the topology!!!) 

##Maximum coordination allowed:
$valM[1] = 1;
$valM[6] = 4;
$valM[7] = 4; ##so we can treat a protonated nitrogen
$valM[8] = 2; ##wont do hydronium

##-----------------------Do not change anything below this line unless you know what you are doing!!!---------------


$forcefield = "opls";
$conect = -1;

for ($i = 0; $i <= @ARGV; $i++) {
 if ($ARGV[$i] eq "-o") {
  $topology = "$ARGV[$i+1]";
 }
 if ($ARGV[$i] eq "-i") {
  $input = "$ARGV[$i+1]";
 }
 if ($ARGV[$i] eq "-ff") {
  $forcefield = "$ARGV[$i+1]";
 }
 if ($ARGV[$i] eq "-c") {
  $charges = "$ARGV[$i+1]";
 }
 if ($ARGV[$i] eq "-conect") {
  if ($ARGV[$i+1] eq "yes")
  {
   $conect = 1;
  }
  elsif ($ARGV[$i+1] eq "no")
  {
   $conect = 0;
  }
  else
  {
   print "Unknown value. Use yes or no. Aborting.\n";
   exit 1;
  }
 }
 if ($ARGV[$i] eq "-h") {
  print "MKTOP - General Molecular Topology Builder\n";
  print "Commands:\n";
  print "-i PDB input\n";
  print "-o MKTOP GROMACS Topology output\n";
  print "-c File containing the charges\n";
  print "-ff Force field (opls/amber)\n";
  print "-connect Use Connect PDB Information (yes/no)\n";
  exit 1;
 }
}
unless (-e $input) {
 print "Fatal error: the input file does not exist\n";
 exit 1;
}
if ($conect == -1)
{
 print "\nYou have to set the -conect option. Aborting.\n";
 exit 1;
}

if($forcefield ne "opls" and $forcefield ne "amber")
{
 printf("Forcefield must be either opls or amber\n");
 exit 1;
}
my @bond;

##defaults


##Including AMBER 94 FF ###

$amber{'opls_135'} = "CT";
$amber{'opls_136'} = "CT";
$amber{'opls_137'} = "CT";
$amber{'opls_138'} = "CT";
$amber{'opls_139'} = "CT";
$amber{'opls_140'} = "HC";
$amber{'opls_141'} = "CM";
$amber{'opls_142'} = "CM";
$amber{'opls_143'} = "CM";
$amber{'opls_144'} = "HC";
$amber{'opls_145'} = "CA";
$amber{'opls_145B'} = "CA";
$amber{'opls_146'} = "HA";
$amber{'opls_147'} = "CA";
$amber{'opls_148'} = "CT";
$amber{'opls_149'} = "CT";
$amber{'opls_150'} = "CM";
$amber{'opls_151'} = "Cl";
$amber{'opls_152'} = "CT";
$amber{'opls_153'} = "H1";
$amber{'opls_154'} = "OH";
$amber{'opls_155'} = "HO";
$amber{'opls_156'} = "H1";
$amber{'opls_157'} = "CT";
$amber{'opls_158'} = "CT";
$amber{'opls_159'} = "CT";
$amber{'opls_160'} = "CT";
$amber{'opls_161'} = "CT";
$amber{'opls_162'} = "OH";
$amber{'opls_163'} = "HO";
$amber{'opls_164'} = "F";
$amber{'opls_165'} = "H3";
$amber{'opls_166'} = "CT";
$amber{'opls_167'} = "OH";
$amber{'opls_168'} = "HO";
$amber{'opls_169'} = "OH";
$amber{'opls_170'} = "HO";
$amber{'opls_171'} = "OH";
$amber{'opls_172'} = "HO";
$amber{'opls_173'} = "CT";
$amber{'opls_174'} = "CT";
$amber{'opls_175'} = "CT";
$amber{'opls_176'} = "H1";
$amber{'opls_178'} = "CM";
$amber{'opls_179'} = "OS";
$amber{'opls_180'} = "OS";
$amber{'opls_181'} = "CT";
$amber{'opls_182'} = "CT";
$amber{'opls_183'} = "CT";
$amber{'opls_184'} = "CT";
$amber{'opls_185'} = "H1";
$amber{'opls_186'} = "OS";
$amber{'opls_187'} = "OH";
$amber{'opls_188'} = "HO";
$amber{'opls_189'} = "CT";
$amber{'opls_190'} = "H2";
$amber{'opls_191'} = "CT";
$amber{'opls_192'} = "H2";
$amber{'opls_193'} = "CT";
$amber{'opls_194'} = "H2";
$amber{'opls_195'} = "CT";
$amber{'opls_196'} = "H2";
$amber{'opls_197'} = "CT";
$amber{'opls_198'} = "CT";
$amber{'opls_200'} = "SH";
$amber{'opls_201'} = "SH";
$amber{'opls_202'} = "S";
$amber{'opls_203'} = "S";
$amber{'opls_204'} = "HS";
$amber{'opls_205'} = "HS";
$amber{'opls_206'} = "CT";
$amber{'opls_207'} = "CT";
$amber{'opls_208'} = "CT";
$amber{'opls_209'} = "CT";
$amber{'opls_210'} = "CT";
$amber{'opls_211'} = "CT";
$amber{'opls_212'} = "CT";
$amber{'opls_213'} = "CT";
$amber{'opls_214'} = "CT";
$amber{'opls_215'} = "CT";
$amber{'opls_216'} = "CT";
$amber{'opls_217'} = "CT";
$amber{'opls_218'} = "CT";
$amber{'opls_219'} = "CT";
$amber{'opls_220'} = "CT";
$amber{'opls_221'} = "CT";
$amber{'opls_222'} = "S";
$amber{'opls_223'} = "CT";
$amber{'opls_223B'} = "CT";
$amber{'opls_224'} = "CT";
$amber{'opls_224B'} = "CT";
$amber{'opls_225'} = "CT";
$amber{'opls_225B'} = "CT";
$amber{'opls_226'} = "Cl";
$amber{'opls_227'} = "CM";
$amber{'opls_228'} = "CT";
$amber{'opls_229'} = "CT";
$amber{'opls_230'} = "CT";
$amber{'opls_231'} = "C";
$amber{'opls_232'} = "C";
$amber{'opls_233'} = "C";
$amber{'opls_234'} = "C";
$amber{'opls_235'} = "C";
$amber{'opls_236'} = "O";
$amber{'opls_237'} = "N";
$amber{'opls_238'} = "N";
$amber{'opls_239'} = "C";
$amber{'opls_240'} = "H";
$amber{'opls_241'} = "H";
$amber{'opls_242'} = "CT";
$amber{'opls_243'} = "CT";
$amber{'opls_244'} = "CT";
$amber{'opls_245'} = "CT";
$amber{'opls_246'} = "CT";
$amber{'opls_247'} = "C";
$amber{'opls_248'} = "O";
$amber{'opls_249'} = "N";
$amber{'opls_250'} = "H";
$amber{'opls_251'} = "N";
$amber{'opls_252'} = "C";
$amber{'opls_253'} = "O";
$amber{'opls_254'} = "H";
$amber{'opls_255'} = "H2";
$amber{'opls_256'} = "CT";
$amber{'opls_257'} = "CT";
$amber{'opls_258'} = "CT";
$amber{'opls_259'} = "CT";
$amber{'opls_260'} = "Error: no sp in amber";
$amber{'opls_261'} = "Error: no sp in amber";
$amber{'opls_262'} = "Error: no sp in amber";
$amber{'opls_263'} = "CA";
$amber{'opls_264'} = "Cl";
$amber{'opls_265'} = "N";
$amber{'opls_266'} = "CA";
$amber{'opls_267'} = "C";
$amber{'opls_268'} = "OH";
$amber{'opls_269'} = "O";
$amber{'opls_270'} = "HO";
$amber{'opls_271'} = "C";
$amber{'opls_272'} = "O2";
$amber{'opls_273'} = "CT";
$amber{'opls_274'} = "CT";
$amber{'opls_275'} = "CT";
$amber{'opls_276'} = "CT";
$amber{'opls_277'} = "C";
$amber{'opls_278'} = "O";
$amber{'opls_278'} = "O";
$amber{'opls_279'} = "H1";
$amber{'opls_280'} = "C";
$amber{'opls_281'} = "O";
$amber{'opls_282'} = "HC";
$amber{'opls_283'} = "CT";
$amber{'opls_284'} = "CT";
$amber{'opls_285'} = "CT";
$amber{'opls_286'} = "N3";
$amber{'opls_287'} = "N3";
$amber{'opls_288'} = "N3";
$amber{'opls_289'} = "N";
$amber{'opls_290'} = "N";
$amber{'opls_291'} = "CT";
$amber{'opls_292'} = "CT";
$amber{'opls_292B'} = "CT";
$amber{'opls_293'} = "CT";
$amber{'opls_293B'} = "CT";
$amber{'opls_294'} = "CT";
$amber{'opls_295'} = "CT";
$amber{'opls_296'} = "CT";
$amber{'opls_297'} = "CT";
$amber{'opls_298'} = "CT";
$amber{'opls_299'} = "CT";
$amber{'opls_300'} = "N2";
$amber{'opls_301'} = "H";
$amber{'opls_302'} = "CM";
$amber{'opls_303'} = "N2";
$amber{'opls_304'} = "H";
$amber{'opls_305'} = "CT";
$amber{'opls_306'} = "CT";
$amber{'opls_307'} = "CT";
$amber{'opls_308'} = "CT";
$amber{'opls_309'} = "N3";
$amber{'opls_310'} = "H";
$amber{'opls_311'} = "NC";
$amber{'opls_312'} = "CM";
$amber{'opls_313'} = "N2";
$amber{'opls_314'} = "H";
$amber{'opls_315'} = "CM";
$amber{'opls_316'} = "H4";
$amber{'opls_317'} = "CM";
$amber{'opls_318'} = "H4";
$amber{'opls_319'} = "N";
$amber{'opls_319B'} = "N";
$amber{'opls_320'} = "C";
$amber{'opls_321'} = "N";
$amber{'opls_322'} = "C";
$amber{'opls_323'} = "CM";
$amber{'opls_324'} = "CM";
$amber{'opls_325'} = "H";
$amber{'opls_326'} = "O";
$amber{'opls_327'} = "H";
$amber{'opls_328'} = "O";
$amber{'opls_329'} = "H1";
$amber{'opls_330'} = "H1";
$amber{'opls_331'} = "CT";
$amber{'opls_332'} = "HC";
$amber{'opls_333'} = "NA";
$amber{'opls_333B'} = "NA";
$amber{'opls_334'} = "C";
$amber{'opls_335'} = "NC";
$amber{'opls_336'} = "CA";
$amber{'opls_337'} = "CA";
$amber{'opls_338'} = "CA";
$amber{'opls_339'} = "H";
$amber{'opls_340'} = "O";
$amber{'opls_341'} = "N2";
$amber{'opls_342'} = "H";
$amber{'opls_343'} = "H";
$amber{'opls_344'} = "H4";
$amber{'opls_345'} = "H4";
$amber{'opls_346'} = "NC";
$amber{'opls_347'} = "CQ";
$amber{'opls_348'} = "NC";
$amber{'opls_349'} = "CB";
$amber{'opls_350'} = "CB";
$amber{'opls_351'} = "CA";
$amber{'opls_352'} = "NB";
$amber{'opls_353'} = "CK";
$amber{'opls_354'} = "NA";
$amber{'opls_354B'} = "N*";
$amber{'opls_355'} = "H5";
$amber{'opls_356'} = "N2";
$amber{'opls_357'} = "H";
$amber{'opls_358'} = "H";
$amber{'opls_359'} = "H5";
$amber{'opls_360'} = "H";
$amber{'opls_361'} = "NA";
$amber{'opls_362'} = "CA";
$amber{'opls_363'} = "NC";
$amber{'opls_364'} = "CB";
$amber{'opls_365'} = "CB";
$amber{'opls_366'} = "C";
$amber{'opls_367'} = "H";
$amber{'opls_368'} = "N2";
$amber{'opls_369'} = "H";
$amber{'opls_370'} = "O";

$amber{'opls_520'} = "NC";
$amber{'opls_521'} = "CA";
$amber{'opls_522'} = "CA";
$amber{'opls_523'} = "CA";
$amber{'opls_524'} = "H4";
$amber{'opls_525'} = "HA";
$amber{'opls_526'} = "HA";
$amber{'opls_527'} = "NC";
$amber{'opls_528'} = "CA";
$amber{'opls_529'} = "H4";
$amber{'opls_530'} = "NC";
$amber{'opls_531'} = "CQ";
$amber{'opls_532'} = "CA";
$amber{'opls_533'} = "CA";
$amber{'opls_534'} = "H5";
$amber{'opls_535'} = "H4";
$amber{'opls_536'} = "HA";
$amber{'opls_537'} = "NC";
$amber{'opls_538'} = "CA";
$amber{'opls_539'} = "CA";
$amber{'opls_540'} = "H4";
$amber{'opls_541'} = "HA";
$amber{'opls_542'} = "NA";
$amber{'opls_543'} = "CW";
$amber{'opls_544'} = "C*";
$amber{'opls_545'} = "H";
$amber{'opls_546'} = "H4";
$amber{'opls_547'} = "HA";
$amber{'opls_548'} = "NA";
$amber{'opls_549'} = "NB";
$amber{'opls_550'} = "CV";
$amber{'opls_551'} = "C*";
$amber{'opls_552'} = "CW";
$amber{'opls_553'} = "H";
$amber{'opls_554'} = "H4";
$amber{'opls_555'} = "HA";
$amber{'opls_556'} = "H4";
$amber{'opls_557'} = "NA";
$amber{'opls_558'} = "CR";
$amber{'opls_559'} = "NB";
$amber{'opls_560'} = "CV";
$amber{'opls_561'} = "CW";
$amber{'opls_562'} = "H";
$amber{'opls_563'} = "H5";
$amber{'opls_564'} = "H4";
$amber{'opls_565'} = "H4";
$amber{'opls_566'} = "No suitable atom type";
$amber{'opls_567'} = "No suitable atom type";
$amber{'opls_568'} = "No suitable atom type";
$amber{'opls_569'} = "No suitable atom type";
$amber{'opls_570'} = "No suitable atom type";
$amber{'opls_571'} = "No suitable atom type";
$amber{'opls_572'} = "No suitable atom type";
$amber{'opls_573'} = "No suitable atom type";
$amber{'opls_574'} = "No suitable atom type";
$amber{'opls_575'} = "No suitable atom type";
$amber{'opls_576'} = "No suitable atom type";
$amber{'opls_577'} = "No suitable atom type";
$amber{'opls_578'} = "No suitable atom type";
$amber{'opls_579'} = "No suitable atom type";
$amber{'opls_580'} = "No suitable atom type";
$amber{'opls_581'} = "No suitable atom type";
$amber{'opls_582'} = "No suitable atom type";
$amber{'opls_583'} = "No suitable atom type";
$amber{'opls_584'} = "No suitable atom type";
$amber{'opls_585'} = "No suitable atom type";
$amber{'opls_586'} = "No suitable atom type";
$amber{'opls_587'} = "NA";
$amber{'opls_588'} = "CW";
$amber{'opls_589'} = "C*";
$amber{'opls_590'} = "CA";
$amber{'opls_591'} = "CA";
$amber{'opls_592'} = "CA";
$amber{'opls_593'} = "CA";
$amber{'opls_594'} = "CB";
$amber{'opls_595'} = "CB";
$amber{'opls_596'} = "H";
$amber{'opls_597'} = "H2";
$amber{'opls_598'} = "HA";
$amber{'opls_599'} = "HA";
$amber{'opls_600'} = "HA";
$amber{'opls_601'} = "HA";
$amber{'opls_602'} = "HA";
$amber{'opls_603'} = "NC";
$amber{'opls_604'} = "CA";
$amber{'opls_605'} = "CA";
$amber{'opls_606'} = "CA";
$amber{'opls_607'} = "CA";
$amber{'opls_608'} = "CA";
$amber{'opls_609'} = "CA";
$amber{'opls_610'} = "CA";
$amber{'opls_611'} = "CA";
$amber{'opls_612'} = "CA";
$amber{'opls_613'} = "H4";
$amber{'opls_614'} = "HA";
$amber{'opls_615'} = "HA";
$amber{'opls_616'} = "HA";
$amber{'opls_617'} = "HA";
$amber{'opls_618'} = "HA";
$amber{'opls_619'} = "HA";
$amber{'opls_620'} = "NC";
$amber{'opls_621'} = "CQ";
$amber{'opls_622'} = "NC";
$amber{'opls_623'} = "CB";
$amber{'opls_624'} = "CB";
$amber{'opls_625'} = "CA";
$amber{'opls_626'} = "NB";
$amber{'opls_627'} = "CK";
$amber{'opls_628'} = "NA";
$amber{'opls_629'} = "H5";
$amber{'opls_630'} = "H4";
$amber{'opls_631'} = "H5";
$amber{'opls_632'} = "H";
$amber{'opls_633'} = "No suitable atom type";
$amber{'opls_634'} = "No suitable atom type";
$amber{'opls_635'} = "No suitable atom type";
$amber{'opls_636'} = "No suitable atom type";
$amber{'opls_637'} = "No suitable atom type";
$amber{'opls_638'} = "No suitable atom type";
$amber{'opls_639'} = "No suitable atom type";
$amber{'opls_640'} = "No suitable atom type";


$charmm{'opls_135'} = "CT3x";
$charmm{'opls_136'} = "CT2x";
$charmm{'opls_137'} = "CT1x";
$charmm{'opls_138'} = "CT3x";
$charmm{'opls_139'} = "CT";
$charmm{'opls_140'} = "HA";
$charmm{'opls_141'} = "CE1";
$charmm{'opls_142'} = "CE1";
$charmm{'opls_143'} = "CE2";
$charmm{'opls_144'} = "HE1";
$charmm{'opls_145'} = "CA";
$charmm{'opls_145B'} = "CA";
$charmm{'opls_146'} = "HP";
$charmm{'opls_147'} = "CA";
$charmm{'opls_148'} = "CT3x";
$charmm{'opls_149'} = "CT2x";
$charmm{'opls_150'} = "CE1";
$charmm{'opls_151'} = "CLAL";
$charmm{'opls_152'} = "CT2x";
$charmm{'opls_153'} = "HA";
$charmm{'opls_154'} = "OH1";
$charmm{'opls_155'} = "H";
$charmm{'opls_156'} = "HA";
$charmm{'opls_157'} = "CT3x";
$charmm{'opls_158'} = "CT1x";
$charmm{'opls_159'} = "CT";
$charmm{'opls_160'} = "CT2x";
$charmm{'opls_161'} = "CF3";
$charmm{'opls_162'} = "OH1";
$charmm{'opls_163'} = "H";
$charmm{'opls_164'} = "F3";
$charmm{'opls_165'} = "HA";
$charmm{'opls_166'} = "CA";
$charmm{'opls_167'} = "OH1";
$charmm{'opls_168'} = "H";
$charmm{'opls_169'} = "OH1";
$charmm{'opls_170'} = "H";
$charmm{'opls_171'} = "OH1";
$charmm{'opls_172'} = "H";
$charmm{'opls_173'} = "CT2x";
$charmm{'opls_174'} = "CT1x";
$charmm{'opls_175'} = "CT";
$charmm{'opls_176'} = "HA";
$charmm{'opls_178'} = "HE1";
$charmm{'opls_179'} = "OS";
$charmm{'opls_180'} = "OS";
$charmm{'opls_181'} = "CT3x";
$charmm{'opls_182'} = "CT2x";
$charmm{'opls_183'} = "CT1x";
$charmm{'opls_184'} = "CT";
$charmm{'opls_185'} = "HA";
$charmm{'opls_186'} = "OS";
$charmm{'opls_187'} = "OH1";
$charmm{'opls_188'} = "H";
$charmm{'opls_189'} = "CT2x";
$charmm{'opls_190'} = "HA";
$charmm{'opls_191'} = "CT2x";
$charmm{'opls_192'} = "HA";
$charmm{'opls_193'} = "CT1x";
$charmm{'opls_194'} = "HA";
$charmm{'opls_195'} = "CT1x";
$charmm{'opls_196'} = "HA";
$charmm{'opls_197'} = "CT";
$charmm{'opls_198'} = "CT";
$charmm{'opls_199'} = "CT3x";
$charmm{'opls_200'} = "S";
$charmm{'opls_201'} = "S";
$charmm{'opls_202'} = "S";
$charmm{'opls_203'} = "SM";
$charmm{'opls_204'} = "H";
$charmm{'opls_205'} = "H";
$charmm{'opls_206'} = "CT2x";
$charmm{'opls_207'} = "CT1x";
$charmm{'opls_208'} = "CT";
$charmm{'opls_209'} = "CT3x";
$charmm{'opls_210'} = "CT2x";
$charmm{'opls_211'} = "CT1x";
$charmm{'opls_212'} = "CT";
$charmm{'opls_213'} = "CT3x";
$charmm{'opls_214'} = "CT2x";
$charmm{'opls_215'} = "CT1x";
$charmm{'opls_216'} = "CT";
$charmm{'opls_217'} = "CT3x";
$charmm{'opls_218'} = "CT2x";
$charmm{'opls_219'} = "CT1x";
$charmm{'opls_220'} = "CT";
$charmm{'opls_221'} = "CT2x";
$charmm{'opls_222'} = "S";
$charmm{'opls_223'} = "CT2x";
$charmm{'opls_223B'} = "CT2x";
$charmm{'opls_224'} = "CT1x";
$charmm{'opls_224B'} = "CT1x";
$charmm{'opls_225'} = "CT";
$charmm{'opls_225B'} = "CT";
$charmm{'opls_226'} = "CLAL";
$charmm{'opls_227'} = "CE1";
$charmm{'opls_228'} = "CT3x";
$charmm{'opls_229'} = "CT1x";
$charmm{'opls_230'} = "CT";
$charmm{'opls_231'} = "C";
$charmm{'opls_232'} = "C";
$charmm{'opls_233'} = "C";
$charmm{'opls_234'} = "C";
$charmm{'opls_235'} = "C";
$charmm{'opls_236'} = "O";
$charmm{'opls_237'} = "NH2";
$charmm{'opls_238'} = "NH1";
$charmm{'opls_239'} = "N";
$charmm{'opls_240'} = "H";
$charmm{'opls_241'} = "H";
$charmm{'opls_242'} = "CT3x";
$charmm{'opls_243'} = "CT3x";
$charmm{'opls_244'} = "CT2x";
$charmm{'opls_245'} = "CT2x";
$charmm{'opls_246'} = "CT2x";
$charmm{'opls_247'} = "C";
$charmm{'opls_248'} = "O";
$charmm{'opls_249'} = "NH2";
$charmm{'opls_250'} = "H";
$charmm{'opls_251'} = "N";
$charmm{'opls_252'} = "C";
$charmm{'opls_253'} = "O";
$charmm{'opls_254'} = "H";
$charmm{'opls_255'} = "HA";
$charmm{'opls_256'} = "CT3x";
$charmm{'opls_257'} = "CT2x";
$charmm{'opls_258'} = "CT1x";
$charmm{'opls_259'} = "CT";
$charmm{'opls_260'} = "CA";
$charmm{'opls_261'} = "CN";
$charmm{'opls_262'} = "NC";
$charmm{'opls_263'} = "CA";
$charmm{'opls_264'} = "CLAL";
$charmm{'opls_265'} = "NH1";
$charmm{'opls_266'} = "CA";
$charmm{'opls_267'} = "C";
$charmm{'opls_268'} = "OH1";
$charmm{'opls_269'} = "O";
$charmm{'opls_270'} = "H";

if($forcefield eq "amber")
{
 $ffdir = "amber03.ff";
}
if($forcefield eq "opls")
{
 $ffdir = "oplsaa.ff";
}

$nbpath = "$gromacs_dir/$ffdir/ffnonbonded.itp";
$bonpath = "$gromacs_dir/$ffdir/ffbonded.itp";

unless(-e "$nbpath") {
 print("ERROR: File $nbpath\n");
 print("Please check the gromacs_dir variable at the beginning of the script\n");
 exit 1;
}

open (FF, "<$nbpath");
$PARAMETERS = 0;
  LINE: while (<FF>) {
       ( /^#/)      and next LINE;
       ( /^\;/)      and next LINE;
       ( /^\[/)      and next LINE;
       (/^\s*$/)   and next LINE;
       chomp;
       @data = &split_blank($_);
        $PARAMETERS++;
        @a = split(/\t/, $data[1]); 
        $opls{$data[0]}=$a[0];
  }
  close FF;

for ($j = 1; $j <= 100; $j++) {
 if ($valM[$i] == 0) {
  $valM[$i] == 6;
 }
 for ($k = 1; $k <= 100; $k++) {
  if ($length[$k][$j] =~ /\d/) {
   $length[$j][$k]=$length[$k][$j];
  }
 }
}

for ($j = 2; $j <= 18; $j++) {
 for ($k = 2; $k <= 18; $k++) {
  if ($length[$k][$j] !~ /\d/) {  ## if there is not anything already defined use this... (may not work)
   $length[$k][$j]=1.6;
   $length[$j][$k]=1.6;
  }
 }
}

for ($j = 2; $j <= 50; $j++) {
 for ($k = 19; $k <= 50; $k++) { ## same thing, different value...
  if ($length[$k][$j] !~ /\d/) {
   $length[$k][$j]=2.0;
  }
  if ($length[$j][$k] !~ /\d/) {
   $length[$j][$k]=2.0;
  }
 }
}

$mol[1]=' H';
$mol[6]=' C';
$mol[7]=' N';
$mol[8]=' O';
$mol[9]=' F';
$mol[14]='SI';
$mol[15]=' P';
$mol[16]=' S';
$mol[17]='CL';
$mol[26]='FE';
$mol[36]='BR';
$mol[53]=' I';

$mass[1]=1.00800;
$mass[6]=12.01100;
$mass[7]=14.00670;
$mass[8]=15.99940;
$mass[9]=18.99840;
$mass[14]=28.08000;
$mass[15]=30.97376;
$mass[16]=32.06000;
$mass[17]=35.45300;
$mass[26]=55.84700;
$mass[36]=79.90400;
$mass[53]=126.90450;


print "\n**************MKTOP 2.0 - General Molecular Topology Builder**************\n";
printf "\nCite this work as:\n";
printf "Ribeiro, A.A.S.T.; Horta, B.A.C.; de Alencastro, R.B.  J. Braz. Chem. Soc., Vol. 19, No. 7, 1433-1435, 2008.\n";

print "\nThe PDB input is: $input\n";

$linha=1;
$line=1;
    open (DB, "<$input");
    LINE: while (<DB>) {
       $line++;
       ( /^#/)      and next LINE;
       (/^\s*$/)   and next LINE;
       chomp;
       @data = &split_blank($_);
       if ($data[0] eq 'HETATM' or $data[0] eq 'ATOM') {
        $atom[$linha]=substr($_,12,4);
        $atom[$linha] =~ s/^\s+//;
        $atom[$linha] =~ s/^\d+//;
        $atom[$linha] =~ s/\s+$//;
        $atom[$linha] =~ s/\d+$//;
        $res[$linha]=substr($_,17,3);
        $x[$linha]=substr($_,30,8);
        $y[$linha]=substr($_,38,8);
        $z[$linha]=substr($_,46,8);
        $symbol[$linha]=substr($_,76,2);
        if(!$symbol[$linha])
        {
          printf("\nYour PDB lacks element symbol information (columns 77-78, right-justified). Aborting\n");
          exit 1;
        }

        for ($j = 1; $j <= 100; $j++) {
         if ($symbol[$linha] eq $mol[$j]) {
          $Z[$linha]=$j;
          $count[$j]=1;
          last;
         }
        }
        if($j eq 101)
        {
         print "\nWe could not determine atom type. Check your input (especially columns 77-78). Aborting\n";
         printf("The problem was on line %d.\n",$line-1);
         exit 1;
        }
           # last LINE;
        $linha=($linha+1);
       }
       if ($data[0] eq 'CONECT' and $conect == 1) {
        $i = $data[1];
         $bond_list[$i] = "$data[2]";
         $bond[$i][$Z[$data[2]]]++;
        for($j=3;$j<=$#data;$j++)
        {
         $bond_list[$i] = "$bond_list[$i]:$data[$j]";
         $bond[$i][$Z[$data[$j]]]++;
        }
       }
      }
    close DB;
$linha--;
 if ($linha == 0) {
  print "There is something wrong with the PDB input.\n";
  print "MKTOP did not find any atoms present in the system.\n";
  print "Atoms should be specified as either 'ATOM' or 'HETATM'\n";
  exit 1;
 }

print "\nEntering Topology Builder...\n\n";
if ($topology !~ /\w/ and $topology !~ /\d/) {
 print "You must specify an output file with the -o option\n";
 exit 1;
}

print "\nEntering Topology Builder...\n\n";

$read_charge=0;

if ($charges =~ /\d/ or $charges =~ /\w/) {
 $read_charge=1;
 print "Molecular Charges: $charges\n";
 unless (-e $charges) {
  print "Fatal error: $charges not found\n";
  exit 1;
 }
}
else {
 print "Charges will not be included in the topology because they were not provided!\n\n";
}

if (-e $topology) {
 print "Fatal error: the topology output file already exists\n";
 exit 1;
}

 print "Number of atoms: $linha\n";

 for ( $i=1; $i <= $linha; $i++) {
         $a=$Z[$i];
         $control=0;
         while (@list) {
          pop(@list);
         }
         $kill=0;
         for ( $j = $linha; $j > $i; $j--) {
             $b=$Z[$j];
             $deltabond = $deltab;
             for ($X = 1; $X <= 100; $X++) {
              if ($Z[$i] == $X) {
               if ($deltaX[$X] > 0) {
                $deltabond = $deltaX[$X];
               }
               if ($valence[$i] == $valM[$X]) {
                last;  
               }
              }   
             }
             $minbond=($length[$a][$b]-$deltabond);
             $maxbond=($length[$a][$b]+$deltabond);
              $difx[$i][$j]=($x[$j]-$x[$i]);
              $difx[$j][$i]=($x[$i]-$x[$j]);
              $dify[$i][$j]=($y[$j]-$y[$i]);
              $dify[$j][$i]=($y[$i]-$y[$j]);
              $difz[$i][$j]=($z[$j]-$z[$i]);
              $difz[$j][$i]=($z[$i]-$z[$j]);
              $dist[$i][$j]=sqrt(($difx[$i][$j]**2) + ($dify[$i][$j]**2) + ($difz[$i][$j]**2));
              $dist[$j][$i]=$dist[$i][$j];
             if ($dist[$i][$j] <= $maxbond and $dist[$i][$j] >= $minbond and $i != $j) {
              $valence[$i]++;
              $valence[$j]++;
              for ($k = 0; $k <= 10; $k++) {
               if ($array[$i][$k] == 0) {
                $array[$i][$k] = $j;
                last;
               }
              }
              for ($k = 0; $k <= 10; $k++) {
               if ($array[$j][$k] == 0) {
                $array[$j][$k] = $i;
                last;
               }
              } 
             }
        }
      for ($k = 0; $k <= 10; $k++) {
       if ($array[$i][$k] > 0) {
        push(@list,$array[$i][$k]);
       }
      }
  if($conect == 0)
  {
       $bond_list[$i] = join (':', @list);
      foreach $atom (@list) {
       if ($bond[$i][$Z[$atom]] =~ /\d/) {
        $bond[$i][$Z[$atom]]++;
       }
       else {
       $bond[$i][$Z[$atom]]=1;
       }
      }
  }
  else
  {
   if(!$bond_list[$i])
   {
    print "You said we could use PDB Conect information (default choice) but we did not find information for atom $i.\nAborting\n";
    exit 1;
   }
  }
 }
open (OUT, ">>$topology");

if(!$conect) 
{
 print "\nWe followed your instructions and did not use PDB Conect information\nThis means CHECK OUT YOUR TOPOLOGY!!!!!\n\n";
}


print OUT qq|;
; This topology was generated by MKTOP ;

; Please cite this work.

; Ribeiro, A.A.S.T.; Horta, B.A.C.; de Alencastro, R.B.  J. Braz. Chem. Soc., Vol. 19, No. 7, 1433-1435, 2008.

; MKTOP does not generate charges!!! You will need to take care of that yourself!!

; Check your topology!! Make sure the atom types make sense.

#include "$gromacs_dir/$ffdir/forcefield.itp"


[ moleculetype ]
; Name            nrexcl
MOL             3

[ atoms ]
;   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB      massB
|;

$GROUP=1;

 if ($read_charge == 1) {
  open (RESP, "<$charges");
   LINE: while (<RESP>) {
       ( /^#/)      and next LINE;
       (/^\s*$/)   and next LINE;
       chomp;
       @data = &split_blank($_);
       $i = $data[0];
       $q[$i]=$data[1];
   }
  close RESP;
 }

   ## now it begins!!!!!
 for ($i = 1; $i <= $linha; $i++) {
    
    $Z = $Z[$i];    
   
    @array = split(/:/, $bond_list[$i]);
    if ($Z[$i] == 6) {
     $group[$i]=$GROUP;
     $GROUP++;
    }
    if ($Z == 1) {
     if ($bond[$i][6] == 1) { 
      $type[$i]="opls_140";
     }
    }
    if ($Z == 6) {
     if (&bondif($i,7,1,8,1,6,1) and $valence[$i] == 3) {
      $type[$i] = "opls_235";
     }
     if (&bondif($i,7,1,8,1,1,1) and $valence[$i] == 3) {
      $type[$i] = "opls_235";
     }
     if (&bondif($i,7,2,8,1) and $valence[$i] == 3) {
      $type[$i] = "opls_247";
     }
     if (&bondif($i,8,2) and $valence[$i] == 4) {
      $diol=0;
      foreach $atom (@array) {
       if ($Z[$atom] == 8) {
        @array2 = split(/:/, $bond_list[$atom]);
        foreach $atom2 (@array2) {
         if ($Z[$atom2] == 1) {
          $DIOL[$diol] = $atom;
          $diol++;
         }
        }
       }
      }
      if ($diol == 2) {
       $type[$DIOL[0]] = "opls_169";
       $type[$DIOL[1]] = "opls_169";
      }
     }
     if (&bondif($i,8,3) and $valence[$i] == 4) {
      $triol=0;
      foreach $atom (@array) {
       if ($Z[$atom] == 8) {
        @array2 = split(/:/, $bond_list[$atom]);
        foreach $atom2 (@array2) {
         if ($Z[$atom2] == 1) {
          $TRIOL[$triol] = $atom;
          $triol++;
         }
        }
       }
      }
      if ($triol == 3) {
       $type[$TRIOL[0]] = "opls_171";
       $type[$TRIOL[1]] = "opls_171";
       $type[$TRIOL[2]] = "opls_171";
      }
     }
     if (&bondif($i,6,1,1,3)) {
      $type[$i] = "opls_135";  ##isobutane
     }
     if (&bondif($i,6,2,1,2)) {
      $type[$i] = "opls_136";
     }
     if (&bondif($i,6,3,1,1)) {
      $type[$i] = "opls_137";
     }
     if (&bondif($i,6,4)) {
      $type[$i] = "opls_139";
     }
     if (&bondif($i,9,3)) {
      $type[$i] = "opls_961";
     }
     if (&bondif($i,9,2)) {
      $type[$i] = "opls_962";
     }
     if ($valence[$i] == 2) {
      if (&bondif($i,6,2)) {
       $gone=0;
       foreach $atom (@array) {
        if ($gone == 0) {
         if (&bondif($atom,1,3) or &bondif($atom,1,2)) {
          $type[$i] = "opls_927";
         }
         elsif (&bondif($atom,1,1)) {
          $type[$i] = "opls_928";
         }
         else {
          $type[$i] = "opls_929";
         }
        }
        if (&bondif($atom,6,2) and $valence[$atom] == 2) {
         $gone = 1;
         $type[$i] = "opls_931";
         $type[$atom] = "opls_931";
        }
        if (&bondif($atom,1,1) and $valence[$atom] == 2) {
         $type[$atom] = "opls_925";
        }
       }
      }
     }
     if (&bondif($i,7,1,6,1) and $valence[$i] == 2) { ##nitrile
      $type[$i] = "opls_754";
     }
    }
    if ($Z == 7) {
     if ($valence[$i] == 4) {
      if (&bondif($i,1,0,6,4)) {
       $type[$i] = "opls_288"; 
      }
      if (&bondif($i,1,1,6,3)) {
       $type[$i] = "opls_940";
      } 
      if (&bondif($i,1,3,6,1) or &bondif($i,1,2,6,2)) {
       $type[$i] = "opls_287";
      } 
     }
     else {
      $amide=0;
      $sulfonamide=0; 
      @array = split(/:/, $bond_list[$i]);
      foreach $atom (@array) {
       if ($Z[$atom] == 6 and $valence[$atom] == 3) {
        @array2 = split (/:/, $bond_list[$atom]);
        foreach $atom2 (@array2) {
         if ($Z[$atom2] == 8) {
          $amide=1;
         }
        }
       }
       if ($Z[$atom] == 16 and &bondif($atom,8,2)) {
        $sulfonamide=1;
       }
      }
      if (&bondif($i,1,1)) {
       $type[$i] = "opls_901"; ## secondary amine
       if ($amide == 1) {
        $type[$i] = "opls_238"; ## peptide bond / secondary amide
       }
       if ($sulfonamide == 1) {
        $type[$i] = "opls_480";
       }
      }
      if (&bondif($i,1,2)) {
       if ($amide == 1) {
        $type[$i] = "opls_237";  ##primary amide
       }
       elsif ($sulfonamide == 1) {
        $type[$i] = "opls_478";
       }
       else {
        $type[$i] = "opls_900"; ## primary amine
       }
      }
      if (&bondif($i,1,0)) {
       if ($amide == 1) {
        if (&ring(0,5,$i,4)) {
         $type[$i] = "opls_238"; #proline
         if (&bondif($array_r[1],1,1)) {
          $type[$array_r[1]] = "opls_246";
          $type[$array_r[4]] = "opls_245"; 
         }
         if (&bondif($array_r[1],1,2)) {
          $type[$array_r[1]] = "opls_245";
          $type[$array_r[4]] = "opls_246";
         }
        }
        else { 
         $type[$i] = "opls_239";  ##tertiary amide
        }
       }
       elsif ($sulfonamide == 1) {
        $type[$i] = "opls_480";
       }
       else {
        $type[$i] = "opls_902"; ## tertiary amine
       }
      }
      if (&bondif($i,8,2)) {
       $type[$i] = "opls_760";
      }
      if (&bondif($i,6,1) and $valence[$i] == 1) { ##nitrile
       $type[$i] = "opls_753"; 
      }
     }
    }
    if ($Z == 8) {
     @array = split(/:/, $bond_list[$i]);   
     if ($bond[$i][6] == 1 and $valence[$i] == 1) { ##carbonyl
      foreach $atom (@array) {
       if ($Z[$atom] == 6) {
        @array2 = split (/:/, $bond_list[$atom]);
        foreach $atom2 (@array2) {
         if ($Z[$atom2] == 7) {
          $type[$i]="opls_236"; ## amide
         }
         if ($Z[$atom2] == 1 and $type[$i] !~ /opls/) {
          $type[$i] = "opls_278"; ## aldehyde
         }
        }
        if ($type[$i] !~ /opls/) {
         $type[$i] = "opls_281"; ## ketone
        }
       }
      }
     }
     if ($bond[$i][1] == 1 and $valence[$i] == 2) { ##alcohol
      $type[$i] = "opls_154";
      foreach $atom (@array) {
       if ($Z[$atom] == 6) {
        if (&bondif($atom,1,3)) {
         $type[$atom] = "opls_157"; ## methyl alcohol
        }
        if (&bondif($atom,1,2)) {
         $type[$atom] = "opls_157"; ## ethyl
        }
        if (&bondif($atom,1,1)) {
         $type[$atom] = "opls_158"; ## i-propyl
        }
        if (&bondif($atom,1,0)) {
         $type[$atom] = "opls_159";  ## t-butyl
        }
       }
      }
     }
     if (&bondif($i,6,2) and $valence[$i] == 2) {
      $type[$i] = "opls_180";  ## ether
      foreach $atom (@array) {
       if (&bondif($atom,1,3)) {
        $type[$atom] = "opls_181"; ## methyl ether
       }
       if (&bondif($atom,1,2)) {
        $type[$atom] = "opls_182"; ## ethyl
       }
       if (&bondif($atom,1,1)) {
        $type[$atom] = "opls_183"; ## i-propyl
       }
       if (&bondif($atom,1,0)) {
        $type[$atom] = "opls_182";  ## t-butyl
       }
      }
     }
    }
    if ($Z == 15) {
     if (&bondif($i,6,4)) {
      $type[$i] = "opls_781";
     } 
    }
    if ($Z == 16) {
     if ($valence[$i] == 2) {
      if (&bondif($i,6,1,1,1)) {
       $type[$i] = "opls_200"; ##thiol
      }
      if (&bondif($i,6,2)) {
       $type[$i] = "opls_202"; #sulfide 
      }
      if (&bondif($i,16,1)) {
       $type[$i] = "opls_203"; ##disulfide
      }
     }
     if ($valence[$i] == 3) {
      if (&bondif($i,8,1)) {
       $type[$i] = "opls_496";
      }
     }
     if ($valence[$i] == 4) {
      if (&bondif($i,8,2,6,2) or &bondif($i,8,3,6,1)) { ##sulfone
       $type[$i] = "opls_493";
      }
      if (&bondif($i,8,2,7,1)) {
       $type[$i] = "opls_474";
      }
     }
    }
    if ($Z == 17) {
     $type[$i] = "opls_151";
    }
    if ($Z == 36) {
     $type[$i] = "opls_722";
    }
   }


for ($i = 1; $i <= $linha; $i++) {
 if ($Z[$i] == 15) {
  if (&bondif($i,8,4) and $valence[$i] == 4) { ##phosphate
   @list = split(/:/, $bond_list[$i]);
   $c=0;
   $phosp[0][1] = "opls_445";
   $phosp[0][2] = "opls_440";
   $phosp[0][3] = "opls_440";
   $phosp[1][1] = "opls_446";
   $phosp[1][2] = "opls_441";
   $phosp[1][3] = $phosp[1][2];
   $phosp[2][1] = "opls_447";
   $phosp[2][2] = "opls_442";
   $phosp[2][3] = $phosp[2][2];
   foreach $atom (@list) {
    if (&bondif($atom,6,1)) {
     $c++;
    }
   }
   foreach $atom (@list) {
    if (&bondif($atom,6,1)) {
     $type[$atom] = $phosp[2][$c];
    }
    else {
     $type[$atom] = $phosp[1][$c];
    }
    if (&bondif($atom,1,1)) {
     $type[$atom] = "opls_154";
    }
   }
   $type[$i] = $phosp[0][$c];
  }
  if (&bondif($i,8,3,6,1) and $valence[$i] == 4) {
   $type[$i] = "opls_450";
   @list = split(/:/, $bond_list[$i]);
   foreach $atom (@list) {
    if (&bondif($atom,6,1) and $Z[$atom] == 8) {
     $type[$atom] = "opls_452";
    }
    elsif ($Z[$atom] == 6) {
     $type[$atom] = "opls_455";
    }
    else {
     $type[$atom] = "opls_451";
    }
   }
  }
 }
}
#updating the carbon types...

for ($i = 1; $i <= $linha; $i++) {
 @list = split(/:/, $bond_list[$i]);
 if ($Z[$i] == 6) {
  if ($valence[$i] == 4) {
   foreach $atom (@list) {
    if ($type[$atom] eq "opls_238" and $type[$i] !~ /opls/) {  ## secondary amide
     foreach $atom2 (@list) {
      if ($type[$atom2] eq "opls_235") {
       if (&bondif($i,6,1,1,2)) {
        $type[$i] = "opls_223B"; # glycine CAlpha
       }
       if (&bondif($i,6,2,1,1) or &bondif($i,6,3)) {
        $type[$i] = "opls_224B"; ## AA Calpha
       }
      }
     }
     if ($type[$i] !~ /opls/) {
      if (&bondif($i,1,3)) {
       $type[$i] = "opls_242";
      } 
      if (&bondif($i,1,2,6,1) or &bondif($i,1,1,6,2) or &bondif($i,6,3)) {
       $type[$i] = "opls_244";
      }
     }
    }
    if ($type[$atom] eq "opls_239") { ## tertiary amide
     if (&bondif($i,1,3)) {
      $type[$i] = "opls_243";
     }
     if (&bondif($i,6,1,1,2)) {
      $type[$i] = "opls_245";
     }
     if (&bondif($i,6,2,1,1) or &bondif($i,6,3)) {
      $type[$i] = "opls_246";
     }
    }
    if ($type[$atom] eq "opls_900") {  ## primary amine
     if (&bondif($i,6,1,1,2)) {
      $type[$i] = "opls_906";
     }
     if (&bondif($i,1,3)) {
      $type[$i] = "opls_903";
     }
    }
    if ($type[$atom] eq "opls_901") { ## secondary amine
     if (&bondif($i,6,1,1,2)) {
      $type[$i] = "opls_907";
     }
     if (&bondif($i,1,3)) {
      $type[$i] = "opls_904";
     }
    }
    if ($type[$atom] eq "opls_902") { ## tertiary amine
     if (&bondif($i,6,1,1,2)) {
      $type[$i] = "opls_908";
     }
     if (&bondif($i,1,3)) {
      $type[$i] = "opls_905";
     }
    }
    if ($type[$atom] eq "opls_940") {
     if (&bondif($i,6,3)) {
      $type[$i] = "opls_945";
     }
     if (&bondif($i,6,2,1,1)) {
      $type[$i] = "opls_944";
     }
     if (&bondif($i,6,1,1,2)) {
      $type[$i] = "opls_943";
     }
     if (&bondif($i,1,3)) {
      $type[$i] = "opls_942";
     }
    }
    if ($type[$atom] eq "opls_442") {
     $type[$i] = "opls_443";
    }
    if ($type[$atom] eq "opls_447") {
     $type[$i] = "opls_448";
    }
    if ($type[$atom] eq "opls_760") {
     if (&bondif($i,1,3)) {
      $type[$i] = "opls_762";
     }
     if (&bondif($i,1,2)) {
      $type[$i] = "opls_765";
     }
     if (&bondif($i,1,1)) {
      $type[$i] = "opls_766";
     }
    }
    if ($type[$atom] eq "opls_474") {
     $type[$i] = "opls_476"; ## C - S sulfonamide
    }
    if ($type[$atom] eq "opls_480") {
     if (&bondif($i,1,3)) {
      $type[$i] = "opls_482";
     }
     else {
      $type[$i] = "opls_484";
     }
    }
    if ($type[$atom] eq "opls_754") {
     if (&bondif($i,6,4)) {
      $type[$i] = "opls_758";
     }
     if (&bondif($i,6,3,1,1)) {
      $type[$i] = "opls_757";
     }
     if (&bondif($i,6,2,1,2)) {
      $type[$i] = "opls_756";
     }
     if (&bondif($i,6,1,1,3)) {
      $type[$i] = "opls_755";
     }
    }
   }
  }
  if ($valence[$i] == 3) {

   ##alkenes
   if (&bondif($i,6,3)) {
    $type[$i] = "opls_141"; 
   }
   for ($k = 9; $k <= 50; $k++) {
    if (&bondif($i,$k,1) and $valence[$k] > 1) {
     if (&bondif($i,1,1)) {
      $type[$i] = "opls_142";
     }
     else {
      $type[$i] = "opls_141";
     }
    }
   }
   if (&bondif($i,6,1,1,1,17,1)) {
    $type[$i] = "opls_227";
   }
   if (&bondif($i,6,2,1,1)) {
    $type[$i] = "opls_142";
   }
   if (&bondif($i,6,1,1,2)) {
    $type[$i] = "opls_143";
   }
   
   ## anything with a carboxyl

   $acid=0;
   
   if (&bondif($i,8,2,6,1) or &bondif($i,8,3)) {
    $acid=1;
   }
   foreach $atom (@list) {
    if ($type[$atom] eq "opls_281" and $type[$i] ne "opls_267" and $type[$i] ne "opls_465") {
     $type[$i] = "opls_280";
    }
    if ($type[$atom] eq "opls_278") {
     $type[$i] = "opls_277";
    }
    if ($acid == 1) {
     if ($valence[$atom] == 2 and $Z[$atom] == 8) {
      $acid++;
      @array = split(/:/, $bond_list[$atom]);
      foreach $atom2 (@array) {
       if ($atom2 != $i) {
        if ($Z[$atom2] == 1) {
         $type[$i] = "opls_267";
         $type[$atom] = "opls_268";
         $type[$atom2] = "opls_270";
        }
        if ($Z[$atom2] == 6) {
         $type[$i] = "opls_465";
         $type[$atom] = "opls_467";
         $type[$atom2] = "opls_468";
        }
       }
      }
     }
    }
   }
   if ($acid == 1) {
    $type[$i] = "opls_271";
   }
  }
 }
}


##acetals

for ($i = 1; $i <= $linha; $i++) {
 if ($Z[$i] == 6 and $valence[$i] == 4) {
  $k=0;
  $H=0;
  @array = &splitb($i);
  foreach $atom (@array) {
   if ($Z[$atom] == 8 and $valence[$atom] == 2) {
    $k++;
   }
   if ($Z[$atom] == 1) {
    $H++;
   }
  }
  if ($k == 2) {
   $hemi=0;
   foreach $atom (@array) {
    if ($Z[$atom] == 8) {
     @array2 = &splitb($atom);
     foreach $atom2 (@array2) {
      if ($Z[$atom2] == 1) {
       $hemi = 1;
       $type[$atom] = "opls_187";
       $type[$atom2] = "opls_188"; 
       if ($H == 0) {
        $type[$i] = "opls_198";
       }
       if ($H == 1) {
        $type[$i] = "opls_195";
       }
       if ($H == 2) {
        $type[$i] = "opls_191";
       } 
      }
      if ($Z[$atom2] == 6 and $atom2 != $i) {
       $type[$atom] = "opls_186";
      }
     }
    }
   }
   if ($hemi == 0) {
    if ($H == 0) {
     $type[$i] = "opls_197";
    }
    if ($H == 0) {
     $type[$i] = "opls_193";
    }
    if ($H == 0) {
     $type[$i] = "opls_189";
    }
   }
   foreach $atom (@array) {
    if ($H == 1 and $Z[$atom] == 1) {
     if ($hemi == 1) {
      $type[$atom] = "opls_196";
     }
     else {
      $type[$atom] = "opls_194";
     }
    }
    if ($H == 2 and $Z[$atom] == 1) {
     if ($hemi == 1) {
      $type[$atom] = "opls_192";
     }
     else {
      $type[$atom] = "opls_190";
     }
    }
   }
  }
 }
}

$A=1;
$D=1;
$P=1;
$B=1;

    ## Now comes the hard part, ring determination!!!

for ($i = 1; $i <= $linha; $i++) {
 if ($Z[$i] != 1) { 
  if ($Z[$i] == 6) {
   if ($valence[$i] == 3 and $type[$i] ne "opls_145" and $type[$i] ne "opls_145B") {  
    if (&ring(0,6,$i,3)) { ##benzene
     for ($k = 0; $k < 6; $k++) {
      $array2[$k] = $array_r[$k];
     } 
     for ($k = 0; $k < 6; $k++) {
      $type[$array2[$k]] = "opls_145";
      @list = split(/:/, $bond_list[$array2[$k]]);
      foreach $atom (@list) {
       $ok=1;
       for ($j = 0; $j < 6; $j++) {
        if ($atom == $array2[$j]) {
         $ok=0;
        }
       }
       if ($ok == 1 and &ring(0,6,$atom,3)) {
        $type[$array2[$k]] = "opls_145B";
       } 
      }
     }
    }
    if (&ring(0,10,$i,3)) {     ## naphthalene
     for ($k = 0; $k < 10; $k++) {
      $type[$array_r[$k]] = "opls_145";
     }
     for ($k = 0; $k < 10; $k++) {
      $ring_bond=0;
      @list = split(/:/, $bond_list[$array_r[$k]]);
      foreach $atom (@list) {
       for ($j = 0; $j < 10; $j++) {
        if ($atom == $array_r[$j]) {
         $ring_bond++;
         if ($ring_bond == 3) {
          $type[$array_r[$k]] = "opls_147";
         }
        }
       }
      }
     }
    }
    if (&ring(0,14,$i,3)) {
     for ($k = 0; $k < 14; $k++) {
      $type[$array_r[$k]] = "opls_145";
     }
    }
   }
   if ($valence[$i] == 4) {
    if (&ring(0,3,$i,4)) { ## cyclopropane
     for ($k = 0; $k < 3; $k++) {
      if ($bond[$array_r[$k]][1] == 2) {
       $type[$array_r[$k]] = "opls_711";
      }
      if ($bond[$array_r[$k]][1] == 1) {
       $type[$array_r[$k]] = "opls_712";
      }
      if ($bond[$array_r[$k]][1] == 0) {
       $type[$array_r[$k]] = "opls_713";
      }
     }
    }
    else {
     $ring[$i][3]=1;
    }
   }
  }
  if ($Z[$i] == 7 and $purine[$i] != 1) { 
   if (&ring(0,6,$i,3) and $valence[$i] == 2) { ## pyridine
    $type[$i] = "opls_520";
    $type[$array_r[1]] = "opls_521";
    $type[$array_r[2]] = "opls_522";
    $type[$array_r[3]] = "opls_523";
    $type[$array_r[4]] = "opls_522";
    $type[$array_r[5]] = "opls_521";
   }
   if (&ring(0,6,$i,3,3,7)) { ## pyrazine
    $type[$i] = "opls_527";
    $type[$array_r[1]] = "opls_528";
    $type[$array_r[2]] = "opls_528";
    $type[$array_r[3]] = "opls_527";
    $type[$array_r[4]] = "opls_528";
    $type[$array_r[5]] = "opls_528";
   }
   if (&ring(0,5,$i,3,2,8)) { ##oxazole
    $type[$i] = "opls_573";
    $type[$array_r[1]] = "opls_572";
    $type[$array_r[2]] = "opls_571";
    $type[$array_r[3]] = "opls_575";
    $type[$array_r[4]] = "opls_574";
   }
   if (&ring(0,5,$i,3,1,8)) { ##isoxazole
    $type[$i] = "opls_580";
    $type[$array_r[1]] = "opls_579";
    $type[$array_r[2]] = "opls_583";
    $type[$array_r[3]] = "opls_582";
    $type[$array_r[4]] = "opls_581";
   }
   if ($valence[$i] == 3) {
    if (&ring(0,5,$i,0,1,7)) { ##pyrazole
     $type[$i] = "opls_548";
     $type[$array_r[1]] = "opls_549"; ## we use the parameters for the Ns if the ring is modified
    }
    if (&ring(0,5,$i,3,1,7)) { 
     $type[$i] = "opls_548";
     $type[$array_r[1]] = "opls_549";
     $type[$array_r[2]] = "opls_550"; ## and for all the atoms if everything is OK
     $type[$array_r[3]] = "opls_551";
     $type[$array_r[4]] = "opls_552";
    }
    if (&ring(0,5,$i,3,2,7) and $type[$i] ne "opls_354" and $type[$i] ne "opls_628") {  ##imidazole
     $type[$i] = "opls_557";
     $type[$array_r[1]] = "opls_558";
     $type[$array_r[2]] = "opls_559"; ## and for all the atoms if everything is OK
     $type[$array_r[3]] = "opls_560";
     $type[$array_r[4]] = "opls_561";
    }
    if (&ring(0,5,$i,3)) { ## pyrrole
     $type[$i] = "opls_542";
     $type[$array_r[1]] = "opls_543";
     $type[$array_r[2]] = "opls_544";
     $type[$array_r[3]] = "opls_544";
     $type[$array_r[4]] = "opls_543";
    }
   }
   if (&ring(0,6,$i,3,2,7) and $type[$i] ne "opls_321" and $type[$i] ne "opls_335") { ## pyrimidine
    $type[$i] = "opls_530";
    $type[$array_r[1]] = "opls_531";
    $type[$array_r[2]] = "opls_530";
    $type[$array_r[3]] = "opls_532";
    $type[$array_r[4]] = "opls_533";
    $type[$array_r[5]] = "opls_532";
    ##the derivatives now...
    if (&bondif($array_r[2],1,1) and &bondif($array_r[3],8,1) and &bondif($array_r[1],8,1)) { ##uracil and thymine
     $type[$i] = "opls_319";
     $type[$array_r[1]] = "opls_320";
     $type[$array_r[2]] = "opls_321";
     $type[$array_r[3]] = "opls_322";
     $type[$array_r[4]] = "opls_323";
     $type[$array_r[5]] = "opls_324";
    }
    if ($valence[$array_r[2]] == 2 and &bondif($array_r[3],7,2) and &bondif($array_r[1],8,1)) { ##cytosine
     $type[$i] = "opls_333";
     $type[$array_r[1]] = "opls_334";
     $type[$array_r[2]] = "opls_335";
     $type[$array_r[3]] = "opls_336";
     $type[$array_r[4]] = "opls_337";
     $type[$array_r[5]] = "opls_338";
    }
   }
   if (&ring(0,9,$i,3,2,7,4,7,6,7)) {
    if ($valence[$array_r[4]] == 3) {
     for ($k = 0; $k <= 9; $k++) {
      $purine[$array_r[$k]] = 1;
     }
     if (&bondif($array_r[8],1,1)) { ##purine
      $type[$i] = "opls_620";
      $type[$array_r[1]] = "opls_621";
      $type[$array_r[2]] = "opls_622";
      $type[$array_r[3]] = "opls_623";
      $type[$array_r[4]] = "opls_628";
      $type[$array_r[5]] = "opls_627";
      $type[$array_r[6]] = "opls_626";
      $type[$array_r[7]] = "opls_624";
      $type[$array_r[8]] = "opls_625";
     }
     if (&bondif($array_r[8],8,1)) { ##guanine
      $type[$i] = "opls_361";
      $type[$array_r[1]] = "opls_362";
      $type[$array_r[2]] = "opls_363";
      $type[$array_r[3]] = "opls_364";
      $type[$array_r[4]] = "opls_354";
      $type[$array_r[5]] = "opls_353";
      $type[$array_r[6]] = "opls_352";
      $type[$array_r[7]] = "opls_365";
      $type[$array_r[8]] = "opls_366";
     }
     if (&bondif($array_r[8],7,2)) { ##adenine
      $type[$i] = "opls_346";
      $type[$array_r[1]] = "opls_347";
      $type[$array_r[2]] = "opls_348";
      $type[$array_r[3]] = "opls_349";
      $type[$array_r[4]] = "opls_354";
      $type[$array_r[5]] = "opls_353";
      $type[$array_r[6]] = "opls_352";
      $type[$array_r[7]] = "opls_350";
      $type[$array_r[8]] = "opls_351";
     }
    }
   }
  }
  if ($Z[$i] == 16) { ## thiazole
   if (&ring(0,5,$i,3,2,7)) {
    $type[$i] = "opls_633";
    $type[$array_r[1]] = "opls_634";
    $type[$array_r[2]] = "opls_635";
    $type[$array_r[3]] = "opls_636";
    $type[$array_r[4]] = "opls_637";
   }
  }
 }
}

for ($i = 1; $i <= $linha; $i++) { ##quinoline
 if ($Z[$i] == 7 and $valence[$i] == 2) {
  $array_r[0]=$i;
  if (&ring(0,10,$i,3)) {
   $type[$i] = "opls_603";

   for ($k = 0; $k <= 10; $k++) {
    $list[$k]=$array_r[$k];
   }

   $array_r[0]=$list[1]; 

   if(&ring(0,6,$array_r[0],3)) { ## wrong order
    for ($k = 0; $k <= 10; $k++) {
     $array_r[0+$k]=$list[10-$k];
    }
   }
   else { ## correct order
    for ($k = 0; $k <= 10; $k++) {
     $array_r[$k]=$list[$k];
    }
   }

   $type[$array_r[1]] = "opls_604";
   $type[$array_r[2]] = "opls_605";
   $type[$array_r[3]] = "opls_606";
   $type[$array_r[4]] = "opls_612";
   $type[$array_r[5]] = "opls_607"; 
   $type[$array_r[6]] = "opls_608"; 
   $type[$array_r[7]] = "opls_609";
   $type[$array_r[8]] = "opls_610";
   $type[$array_r[9]] = "opls_611";
  }
 }
}


##updating hydrogens and some of the heavy atoms...

for ($i = 1; $i <= $linha; $i++) {
 $j = $bond_list[$i];
 if ($Z[$i] == 1) {
  if($forcefield eq "opls")
  {
   if ($type[$j] eq "opls_145") {
    $type[$i] = "opls_146"; ## aromatic
   }
   if ($type[$j] eq "opls_142" or $type[$j] eq "opls_143") {
    $type[$i]="opls_144";
   }
   if ($type[$j] eq "opls_154") {
    $type[$i] = "opls_155"; ##  not methanol!
   }
   if ($type[$j] eq "opls_287") {
    $type[$i] = "opls_290"; ##  not methanol!
   }
   if ($type[$j] eq "opls_940") {
    $type[$i] = "opls_941";
   }
   if ($type[$j] eq "opls_309") {
    $type[$i] = "opls_310";
   } 
   if ($type[$j] eq "opls_238") {
    $type[$i] = "opls_241";
   }
   if ($type[$j] eq "opls_237") {
    $type[$i] = "opls_240";
   }
   if ($type[$j] eq "opls_900") {
    $type[$i] = "opls_909";
   }
   if ($type[$j] eq "opls_901") {
    $type[$i] = "opls_910";
   }
   if ($type[$j] eq "opls_542") {
    $type[$i] = "opls_545";
   }
   if ($type[$j] eq "opls_543") {
    $type[$i] = "opls_546";
   }
   if ($type[$j] eq "opls_544") {
    $type[$i] = "opls_547";
   }
   if ($type[$j] eq "opls_604") {
    $type[$i] = "opls_613";
   }
   if ($type[$j] eq "opls_605") {
    $type[$i] = "opls_614";
   }
   if ($type[$j] eq "opls_606") {
    $type[$i] = "opls_615";
   }
   if ($type[$j] eq "opls_607") {
    $type[$i] = "opls_616";
   }
   if ($type[$j] eq "opls_608") {
    $type[$i] = "opls_617";
   }
   if ($type[$j] eq "opls_609") {
    $type[$i] = "opls_618";
   }
   if ($type[$j] eq "opls_610") {
    $type[$i] = "opls_619";
   }
   if ($type[$j] eq "opls169") {
    $type[$i] = "opls_170";
   }
   if ($type[$j] eq "opls_171") {
    $type[$i] = "opls_172";
   }
   if ($type[$j] eq "opls_249") {  ## urea
    $type[$i] = "opls_250";
   }
   if ($type[$j] eq "opls_443") {
    $type[$i] = "opls_444";
   }
   if ($type[$j] eq "opls_277") {
    $type[$i] = "opls_279";
   }
   if ($type[$j] eq "opls_762" or $type[$j] eq "opls_765" or $type[$j] eq "opls_766") {
    $type[$i] = "opls_763";
   }
   if ($type[$j] eq "opls_581") {
    $type[$i] = "opls_584";
   }
   if ($type[$j] eq "opls_582") {
    $type[$i] = "opls_585";
   }
   if ($type[$j] eq "opls_583") {
    $type[$i] = "opls_586";
   }
   if ($type[$j] eq "opls_476") {
    $type[$i] = "opls_477";
   }
   if ($type[$j] eq "opls_478") {
    $type[$i] = "opls_479";
   }
   if ($type[$j] eq "opls_480") {
    $type[$i] = "opls_481";
   }
   if ($type[$j] eq "opls_482") {
    $type[$i] = "opls_483";
   }
   if ($type[$j] eq "opls_484") {
    $type[$i] = "opls_485";
   }
   if ($type[$j] eq "opls_548") {
    $type[$i] = "opls_553";
   }
   if ($type[$j] eq "opls_550") {
    $type[$i] = "opls_554";
   }
   if ($type[$j] eq "opls_551") {
    $type[$i] = "opls_555";
   }
   if ($type[$j] eq "opls_552") {
    $type[$i] = "opls_556";
   }
   if ($type[$j] eq "opls_521") {
    $type[$i] = "opls_524";
   }
   if ($type[$j] eq "opls_522") {
    $type[$i] = "opls_525";
   }
   if ($type[$j] eq "opls_523") {
    $type[$i] = "opls_526";
   }
   if ($type[$j] eq "opls_634") {
    $type[$i] = "opls_638";
   }
   if ($type[$j] eq "opls_636") {
    $type[$i] = "opls_639";
   }
   if ($type[$j] eq "opls_637") {
    $type[$i] = "opls_640";
   }
   if ($type[$j] eq "opls_523") {
    $type[$i] = "opls_526";
   }
   if ($type[$j] eq "opls_528") {
    $type[$i] = "opls_529";
   }
   if ($type[$j] eq "opls_531") {
    $type[$i] = "opls_534";
   }
   if ($type[$j] eq "opls_532") {
    $type[$i] = "opls_535";
   }
   if ($type[$j] eq "opls_533") {
    $type[$i] = "opls_536";
   }
   if ($type[$j] eq "opls_157") {
    $type[$i] = "opls_156";
   }
   if ($type[$j] eq "opls_200") {
    $type[$i] = "opls_204";
   }
   if ($type[$j] eq "opls_572") {
    $type[$i] = "opls_576";
   }
   if ($type[$j] eq "opls_574") {
    $type[$i] = "opls_577";
   }
   if ($type[$j] eq "opls_575") {
    $type[$i] = "opls_578";
   }
   if ($type[$j] eq "opls_321") {
    $type[$i] = "opls_327";
   }
   if ($type[$j] eq "opls_323") {
    $type[$i] = "opls_239";
   }
   if ($type[$j] eq "opls_324") {
    $type[$i] = "opls_330";
   }
   if ($type[$j] eq "opls_319") {
    $type[$i] = "opls_325";
   }
   if ($type[$j] eq "opls_333") {
    $type[$i] = "opls_339";
   }
   if ($type[$j] eq "opls_337") {
    $type[$i] = "opls_344";
   }
   if ($type[$j] eq "opls_338") {
    $type[$i] = "opls_345";
   }
   if ($type[$j] eq "opls_557") {
    $type[$i] = "opls_562";
   }
   if ($type[$j] eq "opls_558") {
    $type[$i] = "opls_563";
   }
   if ($type[$j] eq "opls_560") {
    $type[$i] = "opls_564";
   }
   if ($type[$j] eq "opls_561") {
    $type[$i] = "opls_565";
   }
   if ($type[$j] eq "opls_361") {
    $type[$i] = "opls_367";
   }
   if ($type[$j] eq "opls_347") {
    $type[$i] = "opls_355";
   }
   if ($type[$j] eq "opls_353") {
    $type[$i] = "opls_359";
   }
   if ($type[$j] eq "opls_354") {
    $type[$i] = "opls_360";
   }
   if ($type[$j] eq "opls_621") {
    $type[$i] = "opls_629";
   }
   if ($type[$j] eq "opls_625") {
    $type[$i] = "opls_630";
   }
   if ($type[$j] eq "opls_627") {
    $type[$i] = "opls_631";
   }
   if ($type[$j] eq "opls_628") {
    $type[$i] = "opls_632";
   }
  }
  if($forcefield eq "amber")
  {
   #With amber we need to pay special attention to hydrogens
   if($Z[$j] == 6)
   {
    @array = split(/:/, $bond_list[$j]);
    $nel=0;
    foreach $atom (@array)
    {
     if(&is_el($Z[$atom]))
     {
      $nel++;
     }
    }
    if($amber{$type[$j]} eq "CT" or $amber{$type[$j]} eq "C" or $amber{$type[$j]} eq "CM")
    {
     $type[$i] = "HC";
     @array = split(/:/, $bond_list[$j]);
     if($nel == 1)
     {
      $type[$i] = "H1";
     }
     if($nel == 2)
     {
      $type[$i] = "H2";
     }
     if($nel == 3)
     {
      $type[$i] = "H3";
     }
    }
    else
    {
     $type[$i] = "HA";
     if($nel == 1)
     {
      $type[$i] = "H4";
     }
     if($nel == 2)
     {
      $type[$i] = "H5";
     }
    }
   }  
   if($Z[$j] == 8)
   {
    $type[$i] = "HO";
   }
   if($Z[$j] == 7)
   {
    $type[$i] = "H";
   }
   if($Z[$j] == 16)
   {
    $type[$i] = "HS";
   }
  }
 }
 if ($Z[$i] == 17) {
  if ($type[$j] eq "opls_145") {
   $type[$i] = "opls_264";
   $type[$j] = "opls_263";
  }
  if ($type[$j] eq "opls_227") {
   $type[$i] = "opls_226";
  }
  if ($Z[$j] == 6 and $valence[$j] == 3 and $type[$i] !~ /opls/) {
   $type[$i] = "opls_264";
  }
 }
 if ($Z[$i] == 9) {
  if ($type[$j] eq "opls_145") { ##fluorobenzene
   $type[$j] = "opls_718";
   $type[$i] = "opls_719";
  }
  if ($type[$j] eq "opls_961" or $type[$j] eq "opls_962") {
   $type[$i] = "opls_965" 
  }
  if ($type[$i] !~ /opls/ and $Z[$j] == 6) {
   @array = split(/:/, $bond_list[$j]);
   $type[$i] = "opls_956";
   $h=0;
   foreach $atom (@array) {
    if ($Z[$atom] == 1) {
     $h++;
     $type[$atom] = "opls_958";
    }
   }
   if ($h == 0) {
    $type[$j] = "opls_960";
   }
   if ($h == 1) {
    $type[$j] = "opls_959";
   }
   if ($h == 2) {
    $type[$j] = "opls_957";
   }
  }
 }
 if ($Z[$i] == 17) {
  if ($type[$j] eq "opls_145") { ##chlorobenzene
   $type[$i] = "opls_264";
   $type[$j] = "opls_263";
  }
 }
 if ($Z[$i] == 36) {
  if ($type[$j] eq "opls_145") {
   $type[$i] = "opls_730";
   $type[$j] = "opls_729";
  }
 }
 if ($Z[$i] == 53) {
  if ($type[$j] eq "opls_145") {
   $type[$i] = "opls_732";
   $type[$j] = "opls_731";
  }
 }

 @array = &splitb($i);
 foreach $atom (@array) {
  if ($Z[$i] == 6) {
   if ($type[$atom] eq "opls_145") {
    if (&bondif($i,6,1,1,3) and $type[$i] eq "opls_135") {
     $type[$i] = "opls_148";
    }
    if (&bondif($i,6,2,1,2) and $type[$i] eq "opls_136") {
     $type[$i] = "opls_149";
    }
    if (&bondif($i,6,3,1,1) and $type[$i] eq "opls_137") {
     $type[$i] = "opls_515";
    }
    if (&bondif($i,6,4) and $type[$i] eq "opls_139") {
     $type[$i] = "opls_515";
    }
    if ($type[$i] eq "opls_267") {
     $type[$i] = "opls_470";
    }
    if ($type[$i] eq "opls_754") {
     $type[$i] = "opls_261";
     $type[$atom] = "opls_260";
    }
    foreach $atom2 (@array) {
     if ($type[$atom2] eq "opls_154") {
      if (&bondif($i,1,2)) {
       $type[$i] = "opls_218";
      }
      if (&bondif($i,1,1)) {
       $type[$i] = "opls_219";
      }
      if (&bondif($i,1,0)) {
       $type[$i] = "opls_220";
      }
     }
     if ($type[$atom2] eq "opls_753") {
      $type[$atom2] = "opls_262";
     }
    }
   }
   if ($type[$atom] eq "opls_277" or $type[$atom] eq "opls_280") {
    foreach $atom2 (@array) {
     if ($Z[$atom2] == 1) {
      $type[$atom2] = "opls_282";
     }
    }
   }
   if ($type[$atom] eq "opls_474" and $type[$i] eq "opls_145") {
    $type[$i] = "opls_488";  #aromatic C - S sulfonamide
   }
   if ($type[$atom] eq "opls_200") {
    if (&bondif($i,1,3)) {
     $type[$i] = "opls_217";
    }
    if (&bondif($i,1,2)) {
     $type[$i] = "opls_206";
    }
    if (&bondif($i,1,1)) {
     $type[$i] = "opls_207";
    }
    if (&bondif($i,1,0)) {
     $type[$i] = "opls_208";
    }
   }
   if ($type[$atom] eq "opls_202" and $valence[$i] == 4) {
    if (&bondif($i,1,3)) {
     $type[$i] = "opls_209";
    }
    if (&bondif($i,1,2)) {
     $type[$i] = "opls_210";
    }
    if (&bondif($i,1,1)) {
     $type[$i] = "opls_211";
    }
    if (&bondif($i,1,0)) {
     $type[$i] = "opls_212";
    }
   }
   if ($type[$atom] eq "opls_203") {
    if (&bondif($i,1,3)) {
     $type[$i] = "opls_213";
    }
    if (&bondif($i,1,2)) {
     $type[$i] = "opls_214";
    }
    if (&bondif($i,1,1)) {
     $type[$i] = "opls_215";
    }
    if (&bondif($i,1,0)) {
     $type[$i] = "opls_216";
    }
   }
   if ($type[$atom] eq "opls_323") {
    $type[$i] = "opls_331";
    foreach $atom2 (@array) {
     if ($Z[$atom2] == 1) {
      $type[$atom2] = "opls_332";
    
     }
    }
   }
   if ($type[$atom] eq "opls_151" and $valence[$i] == 4) {
    $type[$i] = "opls_152";
    foreach $atom2 (@array) {
     if ($Z[$atom2] == 1) {
      $type[$atom2] = "opls_153";
     }
    }
   }
   if ($type[$atom] eq "opls_271") {
    if (&bondif($i,1,2,6,1)) {
     foreach $atom2 (@array) {
      if ($type[$atom2] eq "opls_238") {
       $type[$i] = "opls_284"; #GLY C-terminal
      }
     }
    }
    if (&bondif($i,1,1,6,2)) {
     foreach $atom2 (@array) {
      if ($type[$atom2] eq "opls_238") {
       $type[$i] = "opls_283";
       if (&ring(0,$i,5,4,1,7)) {
        $type[$i] = "opls_285"; # proline
       }
      }
     }
    }
   }
   if ($type[$atom] eq "opls_287") {
    if (&bondif($i,1,2,6,1)) {
     $type[$i] = "opls_292";
     foreach $atom2 (@array) {
      if ($type[$atom2] eq "opls_271") {
       $type[$i] = "opls_298"; #GLY zwitterion
      }
      if ($type[$atom2] eq "opls_235") {
       $type[$i] = "opls_292B"; #GLY N-terminal
      }
     }
    }
    if (&bondif($i,1,1,6,2)) {
     $type[$i] = "opls_293";
     foreach $atom2 (@array) {
      if ($type[$atom2] eq "opls_271") {
       $type[$i] = "opls_299";
      }
      if ($type[$atom2] eq "opls_235") {
       $type[$i] = "opls_293B";
       if (&ring(0,$i,5,4,1,7)) {
        $type[$i] = "opls_295"; # proline
        $type[$array_r[2]] = "opls_296";
       }
      }
     }
    }
    if (&bondif($i,1,3)) {
     $type[$i] = "opls_291";
    }
   }
  }
  if ($Z[$i] == 7) {
   if ($type[$atom] eq "opls_247") {
    $type[$i] = "opls_249";
    foreach $atom2 (@array) {
     if ($Z[$$atom2] == 1) {
      $type[$atom2] = "opls_250";
     }
    }
   }
   if ($type[$atom] eq "opls_145" and $type[$i] eq "opls_760") { #nitrobenzene
    $type[$i] = "opls_767";
    $type[$atom] = "opls_768";
   }
   if ($type[$atom] eq "opls_336" and $type[$i] ne "opls_335") {
    $type[$i] = "opls_341";
    @array_atom = &splitb($atom);
    foreach $atom2 (@array_atom) {
     if ($Z[$atom2] == 6) {
      $C5 = $atom2;
     }
     if ($Z[$atom2] == 7 and $atom2 != $i) {
      $N3 = $atom2;
     }
    }
    foreach $atom2 (@array) {
     if ($Z[$atom2] == 1) {
      if ($dist[$atom2][$C5] <= $dist[$atom2][$N3]) {
       $type[$atom2] = "opls_343";
      }
      else {
       $type[$atom2] = "opls_342";
      }
     }
    }
   }
   if ($type[$atom] eq "opls_362" and $purine[$i] != 1) {
    $type[$i] = "opls_368";
    foreach $atom2 (@array) {
     if ($Z[$atom2] == 1) {
      $type[$atom2] = "opls_369";
     }
    }
   } 
   if ($type[$atom] eq "opls_351" and $purine[$i] != 1) {
    $type[$i] = "opls_356";
    @array_atom = &splitb($atom);
    foreach $atom2 (@array_atom) {
     if ($Z[$atom2] == 6) {
      $C5 = $atom2;
     }
     if ($Z[$atom2] == 7 and $atom2 != $i) {
      $N1 = $atom2;
     }
    }
    foreach $atom2 (@array) {
     if ($Z[$atom2] == 1) {
      if ($dist[$atom2][$C5] <= $dist[$atom2][$N1]) {
       $type[$atom2] = "opls_358";
      }
      else {
       $type[$atom2] = "opls_357";
      }
     }
    }
   }
  }
  if ($Z[$i] == 8) {
   if ($type[$atom] eq "opls_145") {  
    if ($type[$i] eq "opls_154") {
     $type[$i] = "opls_167";
     foreach $atom2 (@array) {
      if ($Z[$atom2] == 6) {
       $type[$atom2] = "opls_166";
      }
      if ($Z[$atom2] == 1) {
       $type[$atom2] = "opls_168";
      }
     }
    }
    if (&bondif($i,6,1,1,0)) {  ##anisole
     $type[$i] = "opls_179";
    }
   }
   if ($type[$atom] eq "opls_247") {
    $type[$i] = "opls_248";
   }
   if ($type[$atom] eq "opls_271") { ##carboxylate
    $type[$i] = "opls_272";
   }
   if ($type[$atom] eq "opls_267" and $valence[$i] == 1) { # 0=C in carboxyl
    $type[$i] = "opls_269";
   }
   if ($type[$atom] eq "opls_465" and $valence[$i] == 1) { # O=C in ester
    $type[$i] = "opls_466";
   }
   if ($type[$atom] eq "opls_493" and $valence[$i] == 1) { ## O=S in sulfone
    $type[$i] = "opls_494";
   }
   if ($type[$atom] eq "opls_474") {
    $type[$i] = "opls_475";
   }
   if ($type[$atom] eq "opls_767" or $type[$atom] eq "opls_760") { #nitro
    $type[$i] = "opls_761";
   }
   if ($type[$atom] eq "opls_496") {
    $type[$i] = "opls_497";
    @array2 = split (/:/, $bond_list[$atom]);
    foreach $atom2 (@array2) {
     if ($Z[$atom2] == 6) {
      if (&bondif($atom2,1,3)) {
       $type[$atom2] = "opls_498";
      }
      if (&bondif($atom2,1,2,6,1)) {
       $type[$atom2] = "opls_499";
      }
     }
    }
   }
   if ($type[$atom] eq "opls_320") {
    $type[$i] = "opls_326";
   }
   if ($type[$atom] eq "opls_322") {
    $type[$i] = "opls_328";
   }
   if ($type[$atom] eq "opls_334") {
    $type[$i] = "opls_340";
   }

  }
  if ($Z[$i] == 16) {
   if ($type[$atom] eq "opls_145") {
    if ($valence[$i] == 2) {
     if (&bondif($i,6,2)) {
      $type[$i] = "opls_222";
      $type[$atom] = "opls_228";
     }
     if (&bondif($i,6,1,1,1)) {
      $type[$i] = "opls_734";
      $type[$atom] = "opls_735";
      foreach $atom2 (@array) {
       if ($Z[$atom2] == 1) {
        $type[$atom2] = "opls_204";
       }
      }
     }
    }
    if ($type[$i] eq "opls_474") {
     $type[$atom] = "opls_488";
    }
   }
  }
 }
}

##dienes
print "Defining bond orders for alkene carbons...\n";

for ($i = 1; $i <= $linha; $i++) {
 $BO_ene[$i] = 0;
 $ene[$i] = 0;
 if ($type[$i] eq "opls_141") {
  $ene[$i] = 3;
 }
 if ($type[$i] eq "opls_142") {
  $ene[$i] = 2;
 }
 if ($type[$i] eq "opls_143") {
  $ene[$i] = 1;
 }
}
for ($i = 1; $i <= $linha; $i++) {
 if ($ene[$i] >= 1) {
  @array = &splitb($i);
  foreach $atom (@array) {
   if ($ene[$atom] == 0) {
    $BO[$i][$atom] = 1;
    $BO[$atom][$i] = 1;
   }
   if ($ene[$atom] >= 1) {
    $BO_ene[$i]++;
   }
  }
  if ($BO_ene[$i] == 1) {
   foreach $atom (@array) {
    if ($ene[$atom] >= 1) {
     $BO[$i][$atom] = 2;
     $BO[$atom][$i] = 2;
    }
   }
  }
 }
}
$kill = 0;
$k = 0;
while ($kill == 0) {
 $ok = 1;
 $k++;
 for ($i = 1; $i <= $linha; $i++) {
  $assign=0;
  if ($BO_ene[$i] >= 2) {
   @array = &splitb($i);
   if ($BO_ene[$i] == 2) {
    foreach $atom (@array) {
     if ($ene[$atom] >= 1 and $BO[$i][$atom] == 2) {
      $assign=1;
     }
     if ($ene[$atom] >= 1 and $BO[$i][$atom] == 1) {
      $assign=2;
     }
    }
   }
   if ($BO_ene[$i] == 3) {
    foreach $atom (@array) {
     if ($ene[$atom] >= 1 and $BO[$i][$atom] == 2) {
      $assign = 1;
     }
    }
   }
   if ($assign) {
    foreach $atom (@array) {
     if ($ene[$atom] >= 1 and $BO[$i][$atom] == 0) {
      if ($assign == 1) {
       $BO[$i][$atom] = 1;
       $BO[$atom][$i] = 1;
       if ($ene[$i] == 2 and $ene[$atom] == 2) {
        $type[$i] = "opls_150";
        $type[$atom] = "opls_150";
       }
       if ($ene[$i] == 3 and $ene[$atom] == 3) {
        $type[$i] = "opls_178";
        $type[$atom] = "opls_178";
       }
      }
      if ($assign == 2) {
       $BO[$i][$atom] = 2;
       $BO[$atom][$i] = 2;
      }
     }
    }
   }
   $sum=0;
   foreach $atom (@array) {
    $sum = ($sum + $BO[$i][$atom]);
   }
   if ($sum != 4) {
    $ok = 0;
   }
  }
 }
 if ($ok == 1 or $k > 50) {
  $kill = 1;
 }
}


for ($i = 1; $i <= $linha; $i++) {  ## this is the if-everything-else-fails section
 if ($type[$i] !~ /./) {
  if ($Z[$i] == 1) {
   $type[$i] = "opls_140";
  }
  if ($Z[$i] == 6) {
   if ($valence[$i] == 4) {
    $type[$i] = "opls_135";
   }
   if ($valence[$i] == 3) {
    $type[$i] = "opls_145";     
   }
  }
  if ($Z[$i] == 7) {
   $type[$i] = "opls_900";
  }
  if ($Z[$i] == 8) {
   if (&bondif($i,1,1)) {
    $type[$i] = "opls_154";
   }
   if (&bondif($i,1,0,6,1)) {
    $type[$i] = "opls_180";
   }
   if (&bondif($i,15,1) or &bondif($i,15,2)) {
    $type[$i] = "opls_441";
   }
   if (&bondif($i,16,1)) {
    $type[$i] = "opls_496";
   }
  }
  if ($Z[$i] == 14) {
   $type[$i] = "SI";
  }
  if ($Z[$i] == 15) {
   $type[$i] = "opls_440";
  }
  if ($Z[$i] == 16) {
   if (&bondif($i,1,1)) {
    $type[$i] = "opls_200";
   }
   elsif (&bondif($i,8,1) or &bondif($i,8,2) or &bondif($i,8,3) or &bondif($i,8,4)) {
    $type[$i] = "opls_497";
   }
   else {
    $type[$i] = "opls_202";
   } 
  } 
 }
}

@defgroup = ('opls_474','opls_493','opls_496','opls_393','opls_440','opls_445','opls_450','opls_781');
for ($i = 1; $i <= $linha; $i++) {
 foreach $group (@defgroup) {
  if ($type[$i] eq $group) {
   $makeg[$i] = 1;
   $group[$i] = $GROUP;
   $GROUP++;
  }
 }
}


for ($i = 1; $i <= $linha; $i++) {
 if ($Z[$i] != 6 and $makeg[$i] != 1) {
  $valenceg=10;
  foreach (&splitb($i)) {
   if ($makeg[$_] == 1) {
    $group[$i] = $group[$_];
    last;
   }
   if ($Z[$_] == 6 and $valence[$_] < $valenceg) {
    $valenceg = $valence[$_];
    $group[$i] = $group[$_];
   }
  }
 }
}

for ($i = 1; $i <= $linha; $i++) {
 if ($group[$i] !~ /\d/) {
  foreach (&splitb($i)) {
   $group[$i] = $group[$_];
  }
 }
}

for ($i = 1; $i <= $linha; $i++) 
{
 if($type[$i] =~ /opls/)
 {
  $btype[$i] = $$forcefield{$type[$i]};
  if($forcefield ne "opls")
  {
   $type[$i] = $$forcefield{$type[$i]};
  }
 }
 else
 {
  $btype[$i] = $type[$i];
 }
}

 unless($read_charge or 1 eq 1) # this needs further testing...
 {
  $i=0;
  open (FF, "<$nbpath");
    LINE: while (<FF>) {
        ( /^#/)      and next LINE;
  #     (/^\s*$/)   and next LINE;
       chomp;
       @data = &split_blank($_);
       if ($data[3] =~ /\d/) {
         for ($j = 1; $j <= $linha; $j++) {
          if ($data[0] eq $type[$j]) {
           $q[$j] = $data[3];
          }
         }
        }
    }
    close FF;
   }
 $i=0;
 open (FF, "<$bonpath");
   LINE: while (<FF>) {
       ( /^#/)      and next LINE;
       (/^\s*$/)   and next LINE;
       chomp;
       @data = &split_blank($_);
       if ($_  =~ /constrainttypes/) {
        last LINE;
       }
       if ($data[4] =~ /\d/ and $data[2] == 1) {
        for ($k = 0; $k <= 1; $k++) {
         $ok=0;
         for ($j = 1; $j <= $linha; $j++) {
          if ($data[$k] eq $btype[$j]) {
           $j = ($linha+1);
           $ok = 1;
          }
         }
         if ($ok == 0) {
          $k = 3;
         }
        }
        if ($ok == 1) {
         $bond_ff[$i] = $_;
         $i++;
        }
       }
   }
   close FF;
  $i=0;
  $read=0;
  open (FF, "<$bonpath");
      LINE: while (<FF>) {
       ( /^#/)      and next LINE;
       (/^\s*$/)   and next LINE;
       chomp;
       if ($_  =~ /dihedraltypes/) {
        last LINE;
       }
       if ($_ =~ /angletypes/) {
        $read = 1;
       }
       if ($read == 1) {
        @data = &split_blank($_);
        for ($k = 0; $k <= 2; $k++) {
         $ok=0;
         for ($j = 1; $j <= $linha; $j++) {
          if ($data[$k] eq $btype[$j]) {
           $j = ($linha+1);
           $ok = 1; 
          }
         }
         if ($ok == 0) {
          $k = 3;
         }
        }
        if ($ok == 1) {
         $angle_ff[$i] = $_;
         $i++;
        }
       }
  }
  $i = 0;
print "Defining bonded terms...\n\n";
for ( $i=1; $i <= $linha; $i++) {
 @array = &splitb($i);
 foreach $atom (@array) {
  if ($atom > $i) {
   $ok=0;
   if ($type_b{$btype[$i]}{$btype[$atom]} > 0) {
    $BOND[$B] = "$i $atom 1 $dist[$i][$atom] $type_b{$btype[$i]}{$btype[$atom]}";
    $ok=1;
    $B++;
   }
   else {
    foreach $entry (@bond_ff) {
     @data = &split_blank($entry);
     if (($data[0] eq "$btype[$i]" and $data[1] eq "$btype[$atom]") or ($data[1] eq "$btype[$i]" and $data[0] eq "$btype[$atom]")) {
      $BOND[$B] = "$i $atom 1 $dist[$i][$atom] $data[4]";
      $type_b{$btype[$i]}{$btype[$atom]}=$data[4];
      $type_b{$btype[$atom]}{$btype[$i]}=$data[4];
      $ok=1;
      $B++;
      last;
     }
    }
   }
   if ($ok == 0) {
    $fix_ff=1;
    $BOND[$B] = "$i $atom 1 $dist[$i][$atom] 0";
    $B++;
   }
  }
  @array2 = &splitb($atom);
  foreach $atom2 (@array2) {
   if ($atom2 > $i) {
    $cos = ((($difx[$atom][$i]*$difx[$atom][$atom2])+($dify[$atom][$i]*$dify[$atom][$atom2])+($difz[$atom][$i]*$difz[$atom][$atom2]))/($dist[$atom][$i]*$dist[$atom][$atom2]));
    $angle = &arccos($cos);
    
    $ok=0;
    if ($type_a{$btype[$i]}{$btype[$atom]}{$btype[$atom2]} > 0) {
     $angle[$A] = "$i $atom $atom2 1 $angle $type_a{$btype[$i]}{$btype[$atom]}{$btype[$atom2]}";
     $ok=1;
     $A++;
    }
    else {
     foreach $entry (@angle_ff) {
      @data = &split_blank($entry);
      if ($data[1] ne $btype[$atom]) {
       next;
      }
      elsif (($data[0] eq $btype[$i] and $data[2] eq $btype[$atom2]) or ($data[2] eq $btype[$i] and $data[0] eq $btype[$atom2])) {
       $angle[$A] = "$i $atom $atom2 1 $angle $data[5]";
       $type_a{$btype[$i]}{$btype[$atom]}{$btype[$atom2]}=$data[5];
       $type_a{$btype[$atom2]}{$btype[$atom]}{$btype[$i]}=$data[5];
       $ok=1;
       $A++;
       last;
      }
     }
    }
    if ($ok == 0) {
     $fix_ff=1;
     $angle[$A] = "$i $atom $atom2 1 $angle 0";
     $A++;
    }
   }
   foreach $atom3 (@array) {
     $c=0;
     foreach $var ($i,$atom,$atom2,$atom3) {
      foreach $var2 ($i,$atom,$atom2,$atom3) {
       if ($var == $var2) {
        $c++;
       }
      }
     }
     if ($c == 4 and $dihedral_gone[$atom3][$i][$atom][$atom2] != 1) {
      $dihedral_gone[$atom3][$i][$atom][$atom2]=1;
      $dihedral_gone[$atom2][$atom][$i][$atom3]=1;
      $ok=0;
      $a = $type[$atom3];
      $b = $type[$i];
      $c = $type[$atom];
      $d = $type[$atom2];
      @data = &split_blank($entry);
      if($forcefield eq "opls")
      {
       $dihedral[$D] = "$atom3 $i $atom $atom2 3";
      }
      if($forcefield eq "amber")
      {
       $dihedral[$D] = "$atom3 $i $atom $atom2 9";
      }
      $D++;
      $pair[$P]="$atom2 $atom3 1";
      $P++;
     }
   }
  }
 }
}


##find improper dihedrals
print "Defining improper dihedrals...\n\n";


$I=1;
for ($i = 1; $i <= $linha; $i++) {
 if ($btype[$i] eq "C" or $btype[$i] eq "C_2" or $btype[$i] eq "C_3") {
  @array = split (/:/, $bond_list[$i]);
  $a = 0;
  $j = 0;
  foreach $atom (@array) {
   if ($btype[$atom] eq "O" or $btype[$atom] eq "O_2") {
    $j = $atom; 
   }
   else {
    $k[$a]=$atom;
    $a++;
   }
  }
  if($forcefield eq "opls")
  {
   $improper[$I]="$k[0] $k[1] $i $j 1    improper_O_C_X_Y";
  }
  if($forcefield eq "amber")
  {
   $improper[$I]="$k[0] $k[1] $i $j 4";
  }
  $I++;
 }
 if ($btype[$i] eq "CM" or $btype[$i] eq "C=") {
  @array = split (/:/, $bond_list[$i]);
  if($forcefield eq "opls")
  {
   $improper[$I]="$array[0] $array[1] $i $array[2] 1   improper_Z_CM_X_Y";
  }
  if($forcefield eq "amber")
  {
   $improper[$I]="$array[0] $array[1] $i $array[2] 4";
  }
  $I++;
 }

 if ($btype[$i] eq "N") {
  @array = split (/:/, $bond_list[$i]);
  if($forcefield eq "opls")
  {
   $improper[$I]="$array[0] $array[1] $i $array[2] 1   improper_Z_N_X_Y";
  }
  if($forcefield eq "amber")
  {
   $improper[$I]="$array[0] $array[1] $i $array[2] 4";
  }
  $I++;
 }
 if ($btype[$i] eq "NO") {
  @array = split (/:/, $bond_list[$i]);
  $Y=0;
  foreach $atom (@array) {
   if ($btype[$atom] ne "ON") {
    $X=$atom;
   }
   elsif($Y==0) {
    $Y=$atom;
   }
   else {
    $Z=$atom;
   }
  }
  if($forcefield eq "opls")
  {
   $improper[$I]="$X $Y $i $Z 1   improper_X_NO_ON_NO";
  }
  if($forcefield eq "amber")
  {
   $improper[$I]="$X $Y $i $Z 4";
  }
  $I++;
 }

 @list = ("CA","CB","CC","CN","CV","CW","CR","CK","CQ","CS","C*");
 foreach $var (@list) {
  if ($btype[$i] eq $var) {
   @array = split (/:/, $bond_list[$i]);
   if($forcefield eq "opls")
   {
    $improper[$I]="$array[0] $array[1] $i $array[2] 1   improper_Z_CA_X_Y";
   }
   if($forcefield eq "amber")
   {
    $improper[$I]="$array[0] $array[1] $i $array[2] 4";
   }
   $I++;
  }
 }

}




#This is for impropers along a ring (Thanks for D van der Spoel and M Hong)

#for ($i=1; $i<$D; $i++) {
# @array = split(/ /,$dihedral[$i]);
# $k = 0;
# for ($j=0;$j<4;$j++)  {
#  foreach $var (@list) {
#   if ($opls{$type[$array[$j]]} eq $var) {
#    $k++;
#   }
#  }
# }
# if($k==4) {
#  $improper[$I] = "$array[0] $array[1] $array[2] $array[3] 1 improper_Z_CA_X_Y";
#  $I++;
# }
#}

## define some special parameters

if ($fix_ff == 1) {
 print "The following section lists bonded terms for which there is no parameterization in the force field and some adaptations were made:\n";
 print "You should refine them...\n\n";

 for ($i=1; $i < $B; $i++) {
  @list = split (/ /, $BOND[$i]);
  if ($list[4] == 0) {
   $a = $type[$list[0]];
   $b = $type[$list[1]];
   $fix="200000.0 kJ mol(-1) nm(-2)";
   $k_fix="200000.0";
   if ($Z[$list[0]] == 6 and $Z[$list[1]] == 6) {
    $fix = "CT and CT";
    $k_fix = "224262.4";
   }
   if (($Z[$list[0]] == 1 and $Z[$list[1]] == 6) or ($Z[$list[0]] == 6 and $Z[$list[1]] == 1)) {
    $fix = "CT and HC";
    $k_fix = "284512.0";
   }
   if (($opls{$a} eq "CM" and $opls{$b} eq "CW") or ($opls{$a} eq "CW" and $opls{$b} eq "CM")) {
    $fix = "CM and CA";
    $k_fix = "357313.6";
   }
   print "bond $list[0] ($btype[$list[0]]) and $list[1] ($btype[$list[1]]) - will use $fix\n";
   $BOND[$i] = "$list[0] $list[1] $list[2] $list[3] $k_fix";
  }
 }
 
 for ($i=1; $i < $A; $i++) {
  @list = split(/ /, $angle[$i]);
  if ($list[5] == 0) {
   $a = $type[$list[0]];
   $b = $type[$list[1]];
   $c = $type[$list[2]];
   $fix = "200.000 kJ mol(-1) rad(-2)";
   $k_fix="200.000";
   if ($Z[$list[0]] == 6 and $Z[$list[1]] == 6 and $Z[$list[2]] == 6) {
    $fix = "CT  - C - CT";
    $k_fix = "585.760";
   }
   foreach $d (0,1,2) {
    if ($Z[$list[$d]] == 1) {
     $fix = "CT - CT - HC";
     if($forcefield = "opls")
     {
      $k_fix = "513.800";
     }
     if($forcefield = "amber")
     {
      $k_fix = "418.400";
     }
     foreach $e (0,1,2) {
      if ($Z[$list[$e]] == 1 and $e != $d) {
       $fix = "HC - CT - HC";
       if($forcefield = "opls")
       {
        $k_fix = "276.144";
       }
       if($forcefield = "amber")
       {
        $k_fix = "292.880";
       }
      }
     } 
    }
    if ($Z[$list[$d]] == 7) {
     $fix = "CT - CT - N";
     $k_fix = "669.440";
    }
   }
   if ($btype[$list[1]] eq "CM") {
    if ($Z[$list[0]] == 1 or $Z[$list[2]] == 1) {
     if($forcefield = "opls")
     {
      $fix = "CM  - CM - HC";
      $k_fix = "292.880";
     }
     if($forcefield = "amber")
     {
      $fix = "CM  - CM - HA";
      $k_fix = "418.400";
     }
    }
    else {
     if($forcefield = "opls")
     {
      $fix = "C  - CM - CM";
      $k_fix = "711.280";
     }
     if($forcefield = "amber")
     {
      $fix = "CT  - CM - CM";
      $k_fix = "585.760";
     }
    }
   }
   if ($opls{$b} eq "S") {
    $fix = "CT  - S - CT";
    $k_fix = "518.816";
   }
   print "angle $list[0]($btype[$list[0]]) - $list[1]($btype[$list[1]]) - $list[2]($btype[$list[2]]) - will use $fix\n";
   $angle[$i] = "$list[0] $list[1] $list[2] $list[3] $list[4] $k_fix";
  }
 }
 
 print "******END OF SPECIAL PARAMETERS LIST******\n";

}

for ($i = 1; $i <= 100; $i++) {
 $count[$i]=1;
}


for ($i = 1; $i <= $linha; $i++) {
 $Z=$Z[$i];
    printf OUT "      %3d", $i;
    print OUT "  $type[$i]   1   $res[$i]   ";
    $atom=qq| $mol[$Z]$count[$Z]|;
    printf OUT " %4s", $atom;
    print OUT "    $group[$i]  ";
    printf OUT " %6.5f", $q[$i];
    print OUT "     $mass[$Z]   \n";
       $count[$Z]++;
   }


print "\nWriting all the stuff...\n";

print OUT "\n\n";

 print OUT "\n\n[ bonds ]\n";
 for ($i=1; $i < $B; $i++) {
  @list = split (/ /, $BOND[$i]);
  print OUT"$list[0] $list[1] $list[2]  ";
  printf OUT "%6.3f  ", ($list[3]/10);
  print OUT "$list[4]\n";
 }

 print OUT "\n\n[ angles ]\n";
 for ($i=1; $i < $A; $i++) {
  @list = split (/ /, $angle[$i]);
  print OUT"$list[0] $list[1] $list[2] $list[3]  ";
  printf OUT "%6.3f  ", ($list[4]);
  print OUT "$list[5]\n";
 }

 print OUT "\n\n[ dihedrals ]\n";
 for ($i=1; $i < $D; $i++) {
  @list = split(/ /, $dihedral[$i]);
 # if ($list[5] != 0 or $list[6] != 0 or $list[7] != 0 or $list[8] != 0) {
   print OUT "$dihedral[$i]\n";
 # }
 }
 print OUT "\n\n[ dihedrals ]\n";
 for ($i=1; $i < $I; $i++) {
  print OUT "$improper[$i]\n";
 }
 print OUT "\n\n[ pairs ]\n";
 for ($i=1; $i < $P; $i++) {
  print OUT "$pair[$i]\n";
 }


print OUT qq|\n\n[ system ]
; Name
MKTOP

[ molecules ]
; Compound        #mols
MOL             1|;




close OUT;

print "The GROMACS Topology output was written at: $topology\n";
sub split_decode {
# --------------------------------------------------------
# Takes one line of the database as input and returns an
# array of all the values. It replaces special mark up that
# join_encode makes such as replacing the '``' symbol with a
# newline and the '~~' symbol with a database delimeter.

    my ($input) = shift;
    my (@array) = split (/\Q$db_delim\E/o, $input, $#db_cols+1);
    foreach (@array) {
        s/~~/$db_delim/g;   # Retrieve Delimiter..
        s/``/\n/g;          # Change '' back to newlines..
    }
    return @array;
}

sub split_blank {
#split entry
    my ($input) = shift;
    my @data;
    my (@array) = split (/ /, $input);
    my $k =0;
    foreach $var (@array) {
     if ($var =~ /\d/ or $var =~ /\w/ ) {
      $data[$k]=$var;
      $k++
     }
    }

    return @data;
}


sub bondif {
 my $i = $_[0];
 my $ok=1;
 my $j;
 for ($j = 1; $j <= @_; $j++) {
  if (int($j/2) != ($j/2)) {
   $a=$_[$j];
   $b=$_[$j+1];
   if ($bond[$i][$a] != $b and $b =~ /\d/) {
    $ok = 0;
   }
  }
 }
 return $ok;
}

sub arccos {
 $pi = 3.1415926535;
 my $x = shift;
 my $var = ($pi/2);
 my $i;
 my $ok;
 my $delta = 0.00001;
 $ok = 0;
 $i = 0;
 if (abs($x) > 0.98) {
  if ($x < 0) {
   $i = $pi;
  }
  while ($ok == 0) {
   if ($x > 0) {
    if (cos($i) < $x) {
     $var = $i;
     $ok = 1;
    }
    else {
     $i = ($i + $delta);
    }
   }
   else {
    if (cos($i) > $x) {
     $var = $i;
     $ok = 1;
    }
    else {
     $i = ($i - $delta);
    }
   }
  }
 }
 else {
  for ($i = 0; $i <= 50; $i++) {
   $var = ($var - ((&factorial(2*$i)*($x**((2*$i)+1)))/((&factorial($i)**2)*(2**(2*$i))*((2*$i)+1))));
  }
 }
 $var=($var*180/$pi);
 return ($var);
}


sub factorial {
 my $x = shift;
 my $i;
 my $fact = 1;
 for ($i = $x; $i > 0; $i--) {
  $fact = ($fact*$i);
 }
 return $fact; 
}




sub ring {
# This routine takes at least 4 arguments
# 1st should be always zero, as this sets how many times has the routine called itself(i.e. which atom are we counting)
# 2nd is the total number of atoms in the ring we are trying to find
# 3rd is the initial atom (this can be carbon or non-carbon atom)
# 4th is the number of bonded atoms the all the atoms along the ring need to have
# now come the optional arguments
# if there is one or more non-carbon atoms in the ring, these need to be specified.
# examples:
# &ring(0,6,$i,3) will look for a 6-atom ring starting at $i, so if $i is carbon we would be looking for a benzene-type ring
# &ring(0,6,$i,3,3,7) will do the same as above, except the atom in position 3 needs to be nitrogen. 

 my $j = $_[0];
 my $R = $_[1];
 my $i = $_[2];
 my $valence = $_[3];
 my $m = $_[4];
 my $n = $_[5];
 my $o = $_[6];
 my $p = $_[7];
 my $q = $_[8];
 my $r = $_[9];
 my $k;
 my $gone=0;
 my @atoms;
 my @F;
 my $a;
 my $return;
 my $var;
 my @list;
 my $atom;
 if ($j == 0) {
  $array_r[0] = $i;
 }
 if ($j == $R and $array_r[$R] == $i) {
  return 1;
 }

 for ($k = $R; $k > $j; $k--) {
  $array_r[$k]=0;
 }
 $go_j=0;
 if ($n > 0) {
  if ($j == $m) {
   if ($Z[$array_r[$j]] != $n) {
    return 0;
   }
   else {
    $go_j = 1;
   }
  }
 }
 if ($o > 0) {
  if ($j == $o) {
   if ($Z[$array_r[$j]] != $p) {
    return 0;
   }
   else {
    $go_j = 1;
   }
  }
 }
 if ($q > 0) {
  if ($j == $q) {
   if ($Z[$array_r[$j]] != $r) {
    return 0;
   }
   else {
    $go_j = 1;
   }
  }
 }


 if (($Z[$array_r[$j]] == 6 and ($valence[$array_r[$j]] == $valence or $valence == 0)) or $j == 0 or $go_j == 1) {
  @list = split(/:/, $bond_list[$array_r[$j]]);
  $k=0;
  foreach $atom (@list) {
   unless ($Z[$atom] == 1 or $ring[$atom][$R] == 1) {
    $atoms[$k][$j] = $atom;
    $k++;
   }
  }
  $F[$j]=($k-1);
  for ($k = 0; $k <= $F[$j]; $k++) {
   $ok=1;
   for ($a = 1; $a < $j; $a++) {
    if ($array_r[$a] == $atoms[$k][$j]) {
     $ok=0;
    }
   }
   if ($ok == 1) {
    $array_r[$j+1]=$atoms[$k][$j];
    $gone++;
    $var = &ring(($j+1),$R,$i,$valence,$m,$n,$o,$p,$q,$r);
    if ($var == 1) {
     return $var;
     last;
    }
   }
  }
 }
 if ($gone == 0) {
  return 0;
 }
}

sub ring2 {
# This routine takes at least 4 arguments
# 1st should be always zero, as this sets how many times has the routine called itself(i.e. which atom are we counting)
# 2nd is the total number of atoms in the ring we are trying to find
# 3rd is the initial atom (this can be carbon or non-carbon atom)
# 4th is the number of bonded atoms the all the atoms along the ring need to have
# now come the optional arguments
# if there is one or more non-carbon atoms in the ring, these need to be specified.
# examples:
# &ring(0,6,$i,3) will look for a 6-atom ring starting at $i, so if $i is carbon we would be looking for a benzene-type ring
# &ring(0,6,$i,3,3,7) will do the same as above, except the atom in position 3 needs to be nitrogen. 

 my $j = $_[0];
 my $R = $_[1];
 my $i = $_[2];
 my $valence = $_[3];
 my $m = $_[4];
 my $n = $_[5];
 my $o = $_[6];
 my $p = $_[7];
 my $q = $_[8];
 my $r = $_[9];
 my $k;
 my $gone=0;
 my @atoms;
 my @F;
 my $a;
 my $return;
 my $var;
 my @list;
 my $atom;

 printf("$j $R $i $gone $array_r[$j]\n");

 if ($j == 0) {
  $array_r[0] = $i;
 }
 if ($j == $R and $array_r[$R] == $i) {
  return 1;
 }

 for ($k = $R; $k > $j; $k--) {
  $array_r[$k]=0;
 }
 $go_j=0;
 if ($n > 0) {
  if ($j == $m) {
 printf("la1 $j $R $i $gone $array_r[$j]\n");
   if ($Z[$array_r[$j]] != $n) {
    return 0;
   }
   else {
    $go_j = 1;
   }
  }
 }
 if ($o > 0) {
  if ($j == $o) {
 printf("la 2 $j $R $i $gone $array_r[$j]\n");
   if ($Z[$array_r[$j]] != $p) {
    return 0;
   }
   else {
    $go_j = 1;
   }
  }
 }
 if ($q > 0) {
  if ($j == $q) {
 printf("la3 $j $R $i $gone $array_r[$j]\n");
   if ($Z[$array_r[$j]] != $r) {
    return 0;
   }
   else {
    $go_j = 1;
   }
  }
 }

    if($array_r[$j] == 11)
    {
     printf("onze $go_j\n");
    }

 if (($Z[$array_r[$j]] == 6 and ($valence[$array_r[$j]] == $valence or $valence == 0)) or $j == 0 or $go_j == 1) {
    if($array_r[$j] == 11)
    {
     printf("onze\n");
    }
  @list = split(/:/, $bond_list[$array_r[$j]]);
  $k=0;
  foreach $atom (@list) {
   unless ($Z[$atom] == 1 or $ring[$atom][$R] == 1) {
    $atoms[$k][$j] = $atom;
    $k++;
   }
  }
  $F[$j]=($k-1);
  for ($k = 0; $k <= $F[$j]; $k++) {
   $ok=1;
    if($array_r[$j] == 11)
    {
     printf("onze $atoms[$k][$j]\n");
    }
   for ($a = 1; $a < $j; $a++) {
    if ($array_r[$a] == $atoms[$k][$j]) {
     $ok=0;
    }
   }
   if ($ok == 1) {
    if($array_r[$j] == 11)
    {
     printf("onze $atoms[$k][$j]\n");
    }
    $array_r[$j+1]=$atoms[$k][$j];
    $gone++;
    $var = &ring2(($j+1),$R,$i,$valence,$m,$n,$o,$p,$q,$r);
    if ($var == 1) {
     return $var;
     last;
    }
   }
  }
 }
 if ($gone == 0) {
  return 0;
 }
}
sub splitb {
 my $i = $_[0];
 my  @list = split(/:/, $bond_list[$i]);
 return @list;
}

sub is_el {
my $i = $_[0];
 if($i == 7 or $i == 8 or $i == 9 or $i == 15 or $i == 16 or $i == 17)
 {
  return 1;
 }
 else
 {
  return 0;
 }

}

