OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.45286173) q[0];
sx q[0];
rz(-0.17012574) q[0];
sx q[0];
rz(2.3556019) q[0];
rz(0.6056447) q[1];
sx q[1];
rz(3.7925386) q[1];
sx q[1];
rz(8.7906919) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4579826) q[0];
sx q[0];
rz(-1.542122) q[0];
sx q[0];
rz(-2.3258414) q[0];
x q[1];
rz(3.0787266) q[2];
sx q[2];
rz(-2.2241484) q[2];
sx q[2];
rz(-1.1970929) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4418728) q[1];
sx q[1];
rz(-0.77077121) q[1];
sx q[1];
rz(2.2295879) q[1];
rz(-pi) q[2];
rz(-0.4001873) q[3];
sx q[3];
rz(-1.549198) q[3];
sx q[3];
rz(2.9633629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2259851) q[2];
sx q[2];
rz(-2.6245124) q[2];
sx q[2];
rz(-1.263164) q[2];
rz(1.2849215) q[3];
sx q[3];
rz(-1.4611171) q[3];
sx q[3];
rz(0.29299709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15853515) q[0];
sx q[0];
rz(-1.2936658) q[0];
sx q[0];
rz(0.43757004) q[0];
rz(0.63105398) q[1];
sx q[1];
rz(-2.8129306) q[1];
sx q[1];
rz(-0.22110573) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76925812) q[0];
sx q[0];
rz(-2.6466469) q[0];
sx q[0];
rz(1.1713722) q[0];
rz(-pi) q[1];
rz(-1.8684623) q[2];
sx q[2];
rz(-2.7781099) q[2];
sx q[2];
rz(3.1306981) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4202538) q[1];
sx q[1];
rz(-1.6920648) q[1];
sx q[1];
rz(0.46055693) q[1];
rz(-0.50237327) q[3];
sx q[3];
rz(-1.6564192) q[3];
sx q[3];
rz(1.8300213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9235886) q[2];
sx q[2];
rz(-1.4322586) q[2];
sx q[2];
rz(-0.58829266) q[2];
rz(-0.44899392) q[3];
sx q[3];
rz(-0.42789999) q[3];
sx q[3];
rz(2.0935521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56851971) q[0];
sx q[0];
rz(-3.0439215) q[0];
sx q[0];
rz(-1.0004689) q[0];
rz(2.4160642) q[1];
sx q[1];
rz(-2.0670481) q[1];
sx q[1];
rz(-0.75769889) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0453148) q[0];
sx q[0];
rz(-1.0490388) q[0];
sx q[0];
rz(1.2009215) q[0];
rz(-pi) q[1];
rz(2.9341142) q[2];
sx q[2];
rz(-1.7277272) q[2];
sx q[2];
rz(-2.4776138) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0211027) q[1];
sx q[1];
rz(-2.0334525) q[1];
sx q[1];
rz(-0.91336577) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5030701) q[3];
sx q[3];
rz(-1.6604074) q[3];
sx q[3];
rz(1.5479969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2086601) q[2];
sx q[2];
rz(-2.9186086) q[2];
sx q[2];
rz(-2.3051252) q[2];
rz(1.4423192) q[3];
sx q[3];
rz(-1.6675555) q[3];
sx q[3];
rz(-2.9377655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.02012415) q[0];
sx q[0];
rz(-0.080105372) q[0];
sx q[0];
rz(0.51112038) q[0];
rz(0.077443667) q[1];
sx q[1];
rz(-0.42235342) q[1];
sx q[1];
rz(1.3285332) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15645813) q[0];
sx q[0];
rz(-1.4159604) q[0];
sx q[0];
rz(-1.3550718) q[0];
rz(0.51283299) q[2];
sx q[2];
rz(-2.3094258) q[2];
sx q[2];
rz(-1.8682478) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1991785) q[1];
sx q[1];
rz(-1.5702827) q[1];
sx q[1];
rz(-1.8375988) q[1];
x q[2];
rz(0.10947157) q[3];
sx q[3];
rz(-1.7913892) q[3];
sx q[3];
rz(2.4601439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0478583) q[2];
sx q[2];
rz(-1.9231984) q[2];
sx q[2];
rz(-2.034534) q[2];
rz(2.6691061) q[3];
sx q[3];
rz(-1.5269591) q[3];
sx q[3];
rz(2.2842177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.040314019) q[0];
sx q[0];
rz(-2.2450228) q[0];
sx q[0];
rz(-2.8033946) q[0];
rz(1.2942554) q[1];
sx q[1];
rz(-1.5599524) q[1];
sx q[1];
rz(-2.6370874) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.550068) q[0];
sx q[0];
rz(-1.803777) q[0];
sx q[0];
rz(-2.1165119) q[0];
rz(-1.9485103) q[2];
sx q[2];
rz(-1.1337122) q[2];
sx q[2];
rz(1.9405685) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.63197631) q[1];
sx q[1];
rz(-2.7512433) q[1];
sx q[1];
rz(0.901464) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.65977804) q[3];
sx q[3];
rz(-2.1187083) q[3];
sx q[3];
rz(-1.1783311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0304886) q[2];
sx q[2];
rz(-2.3031074) q[2];
sx q[2];
rz(1.4939235) q[2];
rz(-1.6882287) q[3];
sx q[3];
rz(-0.93795347) q[3];
sx q[3];
rz(-1.7593613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21022739) q[0];
sx q[0];
rz(-0.87451044) q[0];
sx q[0];
rz(-1.0908303) q[0];
rz(0.53972721) q[1];
sx q[1];
rz(-1.9878186) q[1];
sx q[1];
rz(-0.18879034) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.328619) q[0];
sx q[0];
rz(-2.6972174) q[0];
sx q[0];
rz(2.5286753) q[0];
x q[1];
rz(1.8421474) q[2];
sx q[2];
rz(-1.6711298) q[2];
sx q[2];
rz(-3.0999822) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.69953883) q[1];
sx q[1];
rz(-0.25146723) q[1];
sx q[1];
rz(-0.26775189) q[1];
x q[2];
rz(2.4981899) q[3];
sx q[3];
rz(-0.2989558) q[3];
sx q[3];
rz(1.9770196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9699036) q[2];
sx q[2];
rz(-2.1128555) q[2];
sx q[2];
rz(-1.3151273) q[2];
rz(-1.6361489) q[3];
sx q[3];
rz(-2.7841778) q[3];
sx q[3];
rz(-1.2833387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.090102) q[0];
sx q[0];
rz(-2.3110456) q[0];
sx q[0];
rz(0.32399696) q[0];
rz(1.3011159) q[1];
sx q[1];
rz(-2.3062861) q[1];
sx q[1];
rz(-1.3791929) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3430816) q[0];
sx q[0];
rz(-1.3268952) q[0];
sx q[0];
rz(-3.0573465) q[0];
x q[1];
rz(-0.46160134) q[2];
sx q[2];
rz(-1.201655) q[2];
sx q[2];
rz(0.75418562) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.52723253) q[1];
sx q[1];
rz(-2.7990606) q[1];
sx q[1];
rz(0.50369461) q[1];
rz(-pi) q[2];
rz(-2.2961388) q[3];
sx q[3];
rz(-0.77973706) q[3];
sx q[3];
rz(-2.9065135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1088915) q[2];
sx q[2];
rz(-1.4543433) q[2];
sx q[2];
rz(-0.94318715) q[2];
rz(-2.8105248) q[3];
sx q[3];
rz(-1.7739242) q[3];
sx q[3];
rz(0.29512063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.023495) q[0];
sx q[0];
rz(-2.5546615) q[0];
sx q[0];
rz(2.3773637) q[0];
rz(0.14097342) q[1];
sx q[1];
rz(-2.7455175) q[1];
sx q[1];
rz(-1.2932628) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9613567) q[0];
sx q[0];
rz(-2.0016252) q[0];
sx q[0];
rz(-0.21813099) q[0];
rz(-pi) q[1];
x q[1];
rz(0.13883491) q[2];
sx q[2];
rz(-2.5359557) q[2];
sx q[2];
rz(0.56251898) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.075899374) q[1];
sx q[1];
rz(-1.432044) q[1];
sx q[1];
rz(1.1636415) q[1];
rz(-pi) q[2];
rz(-1.771365) q[3];
sx q[3];
rz(-1.9631533) q[3];
sx q[3];
rz(-1.4810824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3747037) q[2];
sx q[2];
rz(-1.039144) q[2];
sx q[2];
rz(-2.5602706) q[2];
rz(2.2733722) q[3];
sx q[3];
rz(-2.7243082) q[3];
sx q[3];
rz(-1.8301615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.593489) q[0];
sx q[0];
rz(-2.5618401) q[0];
sx q[0];
rz(-1.8909489) q[0];
rz(2.1447694) q[1];
sx q[1];
rz(-2.2579028) q[1];
sx q[1];
rz(1.2876127) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9171137) q[0];
sx q[0];
rz(-1.0815485) q[0];
sx q[0];
rz(3.1130303) q[0];
rz(-2.655517) q[2];
sx q[2];
rz(-1.4176148) q[2];
sx q[2];
rz(-1.7024405) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1907562) q[1];
sx q[1];
rz(-1.3637929) q[1];
sx q[1];
rz(0.91916577) q[1];
x q[2];
rz(2.3171114) q[3];
sx q[3];
rz(-1.5319676) q[3];
sx q[3];
rz(-0.83394921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8210956) q[2];
sx q[2];
rz(-0.38791502) q[2];
sx q[2];
rz(-1.256475) q[2];
rz(-2.9750032) q[3];
sx q[3];
rz(-1.5649256) q[3];
sx q[3];
rz(1.0725718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2607516) q[0];
sx q[0];
rz(-0.80417019) q[0];
sx q[0];
rz(-0.32521954) q[0];
rz(-1.1351769) q[1];
sx q[1];
rz(-2.4644641) q[1];
sx q[1];
rz(2.7744055) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31736483) q[0];
sx q[0];
rz(-1.5396848) q[0];
sx q[0];
rz(3.0935862) q[0];
x q[1];
rz(2.514421) q[2];
sx q[2];
rz(-1.3417202) q[2];
sx q[2];
rz(-2.0723745) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.90875188) q[1];
sx q[1];
rz(-0.68896657) q[1];
sx q[1];
rz(2.8244551) q[1];
rz(-2.1807947) q[3];
sx q[3];
rz(-2.4028117) q[3];
sx q[3];
rz(-1.184954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8964768) q[2];
sx q[2];
rz(-0.8090691) q[2];
sx q[2];
rz(1.0661351) q[2];
rz(-0.079244763) q[3];
sx q[3];
rz(-0.60278046) q[3];
sx q[3];
rz(-1.4762896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29522482) q[0];
sx q[0];
rz(-1.5548779) q[0];
sx q[0];
rz(1.5665733) q[0];
rz(0.29905839) q[1];
sx q[1];
rz(-2.538264) q[1];
sx q[1];
rz(-2.6187142) q[1];
rz(1.6499741) q[2];
sx q[2];
rz(-2.3186602) q[2];
sx q[2];
rz(3.0849948) q[2];
rz(-0.15300898) q[3];
sx q[3];
rz(-2.2867793) q[3];
sx q[3];
rz(2.2147562) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
