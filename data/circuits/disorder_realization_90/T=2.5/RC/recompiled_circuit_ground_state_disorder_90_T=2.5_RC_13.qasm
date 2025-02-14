OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6381792) q[0];
sx q[0];
rz(-1.8621651) q[0];
sx q[0];
rz(2.3655565) q[0];
rz(-2.4322721) q[1];
sx q[1];
rz(-1.5172989) q[1];
sx q[1];
rz(0.59901839) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31604345) q[0];
sx q[0];
rz(-0.49015309) q[0];
sx q[0];
rz(-1.8328299) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.232595) q[2];
sx q[2];
rz(-2.4233492) q[2];
sx q[2];
rz(-2.3423549) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.73246519) q[1];
sx q[1];
rz(-2.5149087) q[1];
sx q[1];
rz(0.62925689) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.78035155) q[3];
sx q[3];
rz(-1.2557185) q[3];
sx q[3];
rz(-1.9596069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0029995) q[2];
sx q[2];
rz(-2.004576) q[2];
sx q[2];
rz(-0.19621672) q[2];
rz(1.044322) q[3];
sx q[3];
rz(-0.82572562) q[3];
sx q[3];
rz(-0.31061068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1021378) q[0];
sx q[0];
rz(-0.65334833) q[0];
sx q[0];
rz(0.85773221) q[0];
rz(2.6013382) q[1];
sx q[1];
rz(-1.0539571) q[1];
sx q[1];
rz(-2.3589755) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7359813) q[0];
sx q[0];
rz(-2.1850039) q[0];
sx q[0];
rz(3.0354397) q[0];
rz(-2.1925794) q[2];
sx q[2];
rz(-1.3781386) q[2];
sx q[2];
rz(-2.4198143) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5021584) q[1];
sx q[1];
rz(-2.5089988) q[1];
sx q[1];
rz(0.32986395) q[1];
x q[2];
rz(2.2563062) q[3];
sx q[3];
rz(-1.1210223) q[3];
sx q[3];
rz(1.3441844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9767849) q[2];
sx q[2];
rz(-2.3471577) q[2];
sx q[2];
rz(0.59107333) q[2];
rz(2.1137386) q[3];
sx q[3];
rz(-0.63575345) q[3];
sx q[3];
rz(-3.0052321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9117821) q[0];
sx q[0];
rz(-2.6214143) q[0];
sx q[0];
rz(-1.0562563) q[0];
rz(2.1504869) q[1];
sx q[1];
rz(-1.3737563) q[1];
sx q[1];
rz(0.80449218) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5493889) q[0];
sx q[0];
rz(-0.69800663) q[0];
sx q[0];
rz(2.7207359) q[0];
rz(-pi) q[1];
rz(0.095011906) q[2];
sx q[2];
rz(-2.1974051) q[2];
sx q[2];
rz(-2.1206405) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.04136297) q[1];
sx q[1];
rz(-1.3303262) q[1];
sx q[1];
rz(-2.9036456) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25193416) q[3];
sx q[3];
rz(-2.478699) q[3];
sx q[3];
rz(2.5980662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6911917) q[2];
sx q[2];
rz(-0.85019008) q[2];
sx q[2];
rz(2.5737393) q[2];
rz(1.9495226) q[3];
sx q[3];
rz(-0.53693938) q[3];
sx q[3];
rz(0.15538628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6464226) q[0];
sx q[0];
rz(-1.3986873) q[0];
sx q[0];
rz(-1.1871185) q[0];
rz(0.25539708) q[1];
sx q[1];
rz(-1.4030158) q[1];
sx q[1];
rz(-0.38937169) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1215243) q[0];
sx q[0];
rz(-0.25857718) q[0];
sx q[0];
rz(0.54988396) q[0];
rz(-1.537048) q[2];
sx q[2];
rz(-0.72740388) q[2];
sx q[2];
rz(-0.45798618) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7070605) q[1];
sx q[1];
rz(-1.7599306) q[1];
sx q[1];
rz(-0.906859) q[1];
rz(-pi) q[2];
x q[2];
rz(0.98020966) q[3];
sx q[3];
rz(-0.23543963) q[3];
sx q[3];
rz(-0.93577945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4046459) q[2];
sx q[2];
rz(-0.32361042) q[2];
sx q[2];
rz(-1.3673937) q[2];
rz(2.6020452) q[3];
sx q[3];
rz(-1.5893785) q[3];
sx q[3];
rz(-1.3246983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0747727) q[0];
sx q[0];
rz(-0.19817752) q[0];
sx q[0];
rz(-2.5812126) q[0];
rz(0.21891521) q[1];
sx q[1];
rz(-1.0812662) q[1];
sx q[1];
rz(-2.970649) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8290299) q[0];
sx q[0];
rz(-1.1033022) q[0];
sx q[0];
rz(-2.3947021) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4291228) q[2];
sx q[2];
rz(-1.3686841) q[2];
sx q[2];
rz(-1.0712445) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.39300181) q[1];
sx q[1];
rz(-2.1868863) q[1];
sx q[1];
rz(-0.93008091) q[1];
rz(-pi) q[2];
x q[2];
rz(0.070329026) q[3];
sx q[3];
rz(-0.95169765) q[3];
sx q[3];
rz(0.17507565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7818452) q[2];
sx q[2];
rz(-1.1283987) q[2];
sx q[2];
rz(2.183059) q[2];
rz(0.16252276) q[3];
sx q[3];
rz(-0.84238094) q[3];
sx q[3];
rz(-1.5490279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.573134) q[0];
sx q[0];
rz(-0.064366654) q[0];
sx q[0];
rz(-2.3934613) q[0];
rz(2.023078) q[1];
sx q[1];
rz(-2.7225814) q[1];
sx q[1];
rz(0.31707877) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7359809) q[0];
sx q[0];
rz(-1.5026341) q[0];
sx q[0];
rz(2.8951485) q[0];
rz(-pi) q[1];
rz(-0.19642475) q[2];
sx q[2];
rz(-1.5228009) q[2];
sx q[2];
rz(0.53507198) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2851781) q[1];
sx q[1];
rz(-0.76907255) q[1];
sx q[1];
rz(-1.673322) q[1];
x q[2];
rz(0.77692399) q[3];
sx q[3];
rz(-1.5424929) q[3];
sx q[3];
rz(1.6203209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.76724425) q[2];
sx q[2];
rz(-2.2129462) q[2];
sx q[2];
rz(-1.413215) q[2];
rz(0.080538571) q[3];
sx q[3];
rz(-1.7929411) q[3];
sx q[3];
rz(-0.22182375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11248511) q[0];
sx q[0];
rz(-1.2487829) q[0];
sx q[0];
rz(1.2741733) q[0];
rz(1.796465) q[1];
sx q[1];
rz(-0.74256623) q[1];
sx q[1];
rz(1.3412195) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9233415) q[0];
sx q[0];
rz(-0.78235432) q[0];
sx q[0];
rz(-0.73114354) q[0];
rz(-1.8597827) q[2];
sx q[2];
rz(-1.1062876) q[2];
sx q[2];
rz(1.4204841) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0772895) q[1];
sx q[1];
rz(-1.3796564) q[1];
sx q[1];
rz(-0.15644507) q[1];
rz(1.7246095) q[3];
sx q[3];
rz(-0.65585583) q[3];
sx q[3];
rz(-2.3668309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.207927) q[2];
sx q[2];
rz(-1.9542481) q[2];
sx q[2];
rz(-2.6742317) q[2];
rz(-2.3815239) q[3];
sx q[3];
rz(-1.4773388) q[3];
sx q[3];
rz(1.2490341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6811328) q[0];
sx q[0];
rz(-1.0204027) q[0];
sx q[0];
rz(1.6946174) q[0];
rz(2.815822) q[1];
sx q[1];
rz(-0.21369801) q[1];
sx q[1];
rz(-1.5209341) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5012799) q[0];
sx q[0];
rz(-0.58459832) q[0];
sx q[0];
rz(-0.49638293) q[0];
x q[1];
rz(-0.09163945) q[2];
sx q[2];
rz(-1.5704201) q[2];
sx q[2];
rz(-0.74972744) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.28315) q[1];
sx q[1];
rz(-1.9236739) q[1];
sx q[1];
rz(2.9550885) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8160024) q[3];
sx q[3];
rz(-2.4225261) q[3];
sx q[3];
rz(1.9089501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1763566) q[2];
sx q[2];
rz(-0.74702817) q[2];
sx q[2];
rz(-2.526324) q[2];
rz(-1.2728914) q[3];
sx q[3];
rz(-0.26510173) q[3];
sx q[3];
rz(0.42241514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0203005) q[0];
sx q[0];
rz(-1.8676119) q[0];
sx q[0];
rz(-2.2913388) q[0];
rz(1.1212768) q[1];
sx q[1];
rz(-1.3669776) q[1];
sx q[1];
rz(2.3698295) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92895011) q[0];
sx q[0];
rz(-1.54689) q[0];
sx q[0];
rz(2.2646409) q[0];
rz(-pi) q[1];
rz(-0.1807801) q[2];
sx q[2];
rz(-1.639699) q[2];
sx q[2];
rz(0.034417424) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5215251) q[1];
sx q[1];
rz(-1.8142533) q[1];
sx q[1];
rz(0.056774898) q[1];
rz(-pi) q[2];
x q[2];
rz(0.87987514) q[3];
sx q[3];
rz(-2.7304683) q[3];
sx q[3];
rz(3.0539235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0261592) q[2];
sx q[2];
rz(-1.5966281) q[2];
sx q[2];
rz(-2.3010632) q[2];
rz(-0.080951512) q[3];
sx q[3];
rz(-2.5637124) q[3];
sx q[3];
rz(-0.24278434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7131272) q[0];
sx q[0];
rz(-2.9436538) q[0];
sx q[0];
rz(-1.9848829) q[0];
rz(-2.1139862) q[1];
sx q[1];
rz(-1.7637858) q[1];
sx q[1];
rz(0.81046945) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5549703) q[0];
sx q[0];
rz(-1.7867309) q[0];
sx q[0];
rz(-0.76540375) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2357622) q[2];
sx q[2];
rz(-1.050569) q[2];
sx q[2];
rz(-0.65641415) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9864784) q[1];
sx q[1];
rz(-2.7271075) q[1];
sx q[1];
rz(2.3398967) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4224206) q[3];
sx q[3];
rz(-1.0142438) q[3];
sx q[3];
rz(-0.30239964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.345574) q[2];
sx q[2];
rz(-1.2556262) q[2];
sx q[2];
rz(-1.6176809) q[2];
rz(-1.1003305) q[3];
sx q[3];
rz(-1.1805781) q[3];
sx q[3];
rz(-2.9617917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0236459) q[0];
sx q[0];
rz(-1.4143586) q[0];
sx q[0];
rz(2.216862) q[0];
rz(0.042451518) q[1];
sx q[1];
rz(-1.1141384) q[1];
sx q[1];
rz(1.3420807) q[1];
rz(1.4412075) q[2];
sx q[2];
rz(-0.40934206) q[2];
sx q[2];
rz(-0.48624292) q[2];
rz(-3.068559) q[3];
sx q[3];
rz(-2.6831153) q[3];
sx q[3];
rz(-1.4014184) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
