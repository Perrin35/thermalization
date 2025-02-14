OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.429739) q[0];
sx q[0];
rz(-1.9965594) q[0];
sx q[0];
rz(2.1249007) q[0];
rz(-3.9859803) q[1];
sx q[1];
rz(3.9581668) q[1];
sx q[1];
rz(8.858455) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53439683) q[0];
sx q[0];
rz(-1.6677987) q[0];
sx q[0];
rz(-2.1792381) q[0];
rz(1.9400305) q[2];
sx q[2];
rz(-2.9297631) q[2];
sx q[2];
rz(-0.43583696) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8493718) q[1];
sx q[1];
rz(-1.148889) q[1];
sx q[1];
rz(1.6590236) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.58523607) q[3];
sx q[3];
rz(-1.7437616) q[3];
sx q[3];
rz(-1.1542785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6301959) q[2];
sx q[2];
rz(-2.2679195) q[2];
sx q[2];
rz(-2.013618) q[2];
rz(-1.6389716) q[3];
sx q[3];
rz(-2.2962544) q[3];
sx q[3];
rz(-2.716841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2750435) q[0];
sx q[0];
rz(-2.6847222) q[0];
sx q[0];
rz(-0.04091111) q[0];
rz(-2.5919137) q[1];
sx q[1];
rz(-2.7626541) q[1];
sx q[1];
rz(-3.0612225) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.550404) q[0];
sx q[0];
rz(-2.4192717) q[0];
sx q[0];
rz(-0.0063615464) q[0];
rz(-pi) q[1];
rz(-0.88826142) q[2];
sx q[2];
rz(-2.3739034) q[2];
sx q[2];
rz(-1.3554128) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6192542) q[1];
sx q[1];
rz(-1.0864746) q[1];
sx q[1];
rz(-3.0794512) q[1];
rz(-2.8714058) q[3];
sx q[3];
rz(-0.98434292) q[3];
sx q[3];
rz(-1.0932066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4139159) q[2];
sx q[2];
rz(-2.2985986) q[2];
sx q[2];
rz(-2.8576039) q[2];
rz(-1.5208288) q[3];
sx q[3];
rz(-1.5512543) q[3];
sx q[3];
rz(0.76333299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82336998) q[0];
sx q[0];
rz(-2.9871873) q[0];
sx q[0];
rz(-2.7247317) q[0];
rz(0.93337026) q[1];
sx q[1];
rz(-2.2552172) q[1];
sx q[1];
rz(-2.0534168) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25167187) q[0];
sx q[0];
rz(-1.8150738) q[0];
sx q[0];
rz(1.0142465) q[0];
x q[1];
rz(-2.5809885) q[2];
sx q[2];
rz(-1.5365155) q[2];
sx q[2];
rz(-2.3115668) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8510906) q[1];
sx q[1];
rz(-0.41176418) q[1];
sx q[1];
rz(-0.74300933) q[1];
x q[2];
rz(2.8540594) q[3];
sx q[3];
rz(-1.5085847) q[3];
sx q[3];
rz(-2.2287318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0205959) q[2];
sx q[2];
rz(-0.72046295) q[2];
sx q[2];
rz(0.32568112) q[2];
rz(-2.7459512) q[3];
sx q[3];
rz(-0.49042693) q[3];
sx q[3];
rz(-0.71780786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91738236) q[0];
sx q[0];
rz(-2.2069187) q[0];
sx q[0];
rz(-0.14064661) q[0];
rz(2.7463101) q[1];
sx q[1];
rz(-1.544516) q[1];
sx q[1];
rz(-2.7093754) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89807876) q[0];
sx q[0];
rz(-1.9742516) q[0];
sx q[0];
rz(-1.0591169) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0555058) q[2];
sx q[2];
rz(-2.0755092) q[2];
sx q[2];
rz(1.1881669) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.85104698) q[1];
sx q[1];
rz(-1.6509027) q[1];
sx q[1];
rz(-2.8677031) q[1];
rz(2.3234576) q[3];
sx q[3];
rz(-2.3664118) q[3];
sx q[3];
rz(-1.1213746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.846659) q[2];
sx q[2];
rz(-1.1120956) q[2];
sx q[2];
rz(-0.46889949) q[2];
rz(-2.9851959) q[3];
sx q[3];
rz(-2.1726435) q[3];
sx q[3];
rz(1.8739353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0676607) q[0];
sx q[0];
rz(-1.1336552) q[0];
sx q[0];
rz(-2.3611948) q[0];
rz(0.91267768) q[1];
sx q[1];
rz(-0.78881216) q[1];
sx q[1];
rz(1.7524293) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75994223) q[0];
sx q[0];
rz(-1.3033322) q[0];
sx q[0];
rz(2.777369) q[0];
rz(-pi) q[1];
rz(1.4413799) q[2];
sx q[2];
rz(-2.1485818) q[2];
sx q[2];
rz(2.8313864) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1300182) q[1];
sx q[1];
rz(-2.5352074) q[1];
sx q[1];
rz(-0.19488402) q[1];
rz(-2.1209929) q[3];
sx q[3];
rz(-2.3781099) q[3];
sx q[3];
rz(2.965345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3461561) q[2];
sx q[2];
rz(-2.4861591) q[2];
sx q[2];
rz(-0.1795086) q[2];
rz(-1.6620212) q[3];
sx q[3];
rz(-2.447465) q[3];
sx q[3];
rz(-3.1365385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1805304) q[0];
sx q[0];
rz(-1.625165) q[0];
sx q[0];
rz(1.6777212) q[0];
rz(3.1278817) q[1];
sx q[1];
rz(-1.6746215) q[1];
sx q[1];
rz(-0.94508583) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22064117) q[0];
sx q[0];
rz(-2.2837229) q[0];
sx q[0];
rz(0.46895194) q[0];
x q[1];
rz(1.7795947) q[2];
sx q[2];
rz(-1.2660625) q[2];
sx q[2];
rz(-1.7727675) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.764115) q[1];
sx q[1];
rz(-2.4903653) q[1];
sx q[1];
rz(0.73173827) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.97894815) q[3];
sx q[3];
rz(-2.0267539) q[3];
sx q[3];
rz(0.78612529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1549687) q[2];
sx q[2];
rz(-1.0633609) q[2];
sx q[2];
rz(2.092579) q[2];
rz(-1.6074041) q[3];
sx q[3];
rz(-2.1395855) q[3];
sx q[3];
rz(1.775942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5904215) q[0];
sx q[0];
rz(-2.6541002) q[0];
sx q[0];
rz(2.3792939) q[0];
rz(0.47362348) q[1];
sx q[1];
rz(-1.381424) q[1];
sx q[1];
rz(2.2332938) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98971924) q[0];
sx q[0];
rz(-2.3494968) q[0];
sx q[0];
rz(0.0072974722) q[0];
x q[1];
rz(-1.3495096) q[2];
sx q[2];
rz(-1.7402116) q[2];
sx q[2];
rz(0.34299037) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.087203793) q[1];
sx q[1];
rz(-1.8076732) q[1];
sx q[1];
rz(0.014976784) q[1];
rz(1.9929816) q[3];
sx q[3];
rz(-3.1109538) q[3];
sx q[3];
rz(0.8493166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.121754) q[2];
sx q[2];
rz(-0.25889954) q[2];
sx q[2];
rz(0.961595) q[2];
rz(0.52729765) q[3];
sx q[3];
rz(-1.7617825) q[3];
sx q[3];
rz(0.56078792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2335662) q[0];
sx q[0];
rz(-0.59995025) q[0];
sx q[0];
rz(0.34616923) q[0];
rz(1.8442122) q[1];
sx q[1];
rz(-0.66645122) q[1];
sx q[1];
rz(0.60874879) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2045743) q[0];
sx q[0];
rz(-1.6547583) q[0];
sx q[0];
rz(-0.89421009) q[0];
rz(-pi) q[1];
rz(-2.8738451) q[2];
sx q[2];
rz(-1.7742947) q[2];
sx q[2];
rz(0.45058077) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.94683719) q[1];
sx q[1];
rz(-1.5018592) q[1];
sx q[1];
rz(-0.26534715) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5777588) q[3];
sx q[3];
rz(-2.2632416) q[3];
sx q[3];
rz(3.0094224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.79494563) q[2];
sx q[2];
rz(-0.55507675) q[2];
sx q[2];
rz(0.65183276) q[2];
rz(0.73976222) q[3];
sx q[3];
rz(-2.0756105) q[3];
sx q[3];
rz(0.8968001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5480492) q[0];
sx q[0];
rz(-0.47320047) q[0];
sx q[0];
rz(2.5439673) q[0];
rz(-1.8400486) q[1];
sx q[1];
rz(-1.5664682) q[1];
sx q[1];
rz(2.8155933) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2151011) q[0];
sx q[0];
rz(-2.282906) q[0];
sx q[0];
rz(2.7551921) q[0];
x q[1];
rz(-2.4818899) q[2];
sx q[2];
rz(-2.7355425) q[2];
sx q[2];
rz(2.7999807) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.9218775) q[1];
sx q[1];
rz(-1.829521) q[1];
sx q[1];
rz(2.2398276) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6889309) q[3];
sx q[3];
rz(-1.6836327) q[3];
sx q[3];
rz(-1.090534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.43805435) q[2];
sx q[2];
rz(-1.9376829) q[2];
sx q[2];
rz(0.34029141) q[2];
rz(1.5002286) q[3];
sx q[3];
rz(-0.54268018) q[3];
sx q[3];
rz(-2.5435508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.930645) q[0];
sx q[0];
rz(-2.0238545) q[0];
sx q[0];
rz(-0.043638226) q[0];
rz(-2.4234407) q[1];
sx q[1];
rz(-1.3190045) q[1];
sx q[1];
rz(1.4290379) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6788504) q[0];
sx q[0];
rz(-1.7587816) q[0];
sx q[0];
rz(1.3975271) q[0];
x q[1];
rz(-1.4347836) q[2];
sx q[2];
rz(-2.2439661) q[2];
sx q[2];
rz(1.7409131) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.43914) q[1];
sx q[1];
rz(-1.0826775) q[1];
sx q[1];
rz(1.6678651) q[1];
rz(-pi) q[2];
x q[2];
rz(0.36504443) q[3];
sx q[3];
rz(-1.2170346) q[3];
sx q[3];
rz(0.34788528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4622197) q[2];
sx q[2];
rz(-2.2023109) q[2];
sx q[2];
rz(-0.77793724) q[2];
rz(1.0434307) q[3];
sx q[3];
rz(-1.5305887) q[3];
sx q[3];
rz(1.9061609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0935681) q[0];
sx q[0];
rz(-1.0014191) q[0];
sx q[0];
rz(0.77145664) q[0];
rz(-1.7780766) q[1];
sx q[1];
rz(-1.9693146) q[1];
sx q[1];
rz(1.6065425) q[1];
rz(-1.0422937) q[2];
sx q[2];
rz(-0.92005554) q[2];
sx q[2];
rz(-1.0307606) q[2];
rz(2.1322973) q[3];
sx q[3];
rz(-2.3875368) q[3];
sx q[3];
rz(-0.33215678) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
