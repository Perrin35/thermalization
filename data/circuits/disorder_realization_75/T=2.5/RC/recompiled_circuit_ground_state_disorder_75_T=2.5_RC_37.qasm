OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2158382) q[0];
sx q[0];
rz(-2.821142) q[0];
sx q[0];
rz(-3.0132063) q[0];
rz(-2.1865891) q[1];
sx q[1];
rz(-0.65989143) q[1];
sx q[1];
rz(2.554472) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4383343) q[0];
sx q[0];
rz(-1.4263784) q[0];
sx q[0];
rz(-2.6718706) q[0];
rz(-0.70694114) q[2];
sx q[2];
rz(-0.64919186) q[2];
sx q[2];
rz(-0.45237088) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.94226664) q[1];
sx q[1];
rz(-2.0055416) q[1];
sx q[1];
rz(-0.64202692) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1863094) q[3];
sx q[3];
rz(-1.5732854) q[3];
sx q[3];
rz(-1.3979504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.2653653) q[2];
sx q[2];
rz(-2.9076125) q[2];
sx q[2];
rz(-0.9210251) q[2];
rz(1.5246576) q[3];
sx q[3];
rz(-1.8504668) q[3];
sx q[3];
rz(0.18572148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8601473) q[0];
sx q[0];
rz(-2.7466819) q[0];
sx q[0];
rz(-1.7922147) q[0];
rz(0.34740627) q[1];
sx q[1];
rz(-1.1888622) q[1];
sx q[1];
rz(1.7385534) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13578519) q[0];
sx q[0];
rz(-2.5134235) q[0];
sx q[0];
rz(2.7676393) q[0];
x q[1];
rz(-2.5591056) q[2];
sx q[2];
rz(-0.91160027) q[2];
sx q[2];
rz(2.2381353) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.385437) q[1];
sx q[1];
rz(-1.7240579) q[1];
sx q[1];
rz(-2.4244244) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2469134) q[3];
sx q[3];
rz(-2.4270456) q[3];
sx q[3];
rz(2.0216048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3226402) q[2];
sx q[2];
rz(-1.3078657) q[2];
sx q[2];
rz(-2.9827706) q[2];
rz(-0.21765503) q[3];
sx q[3];
rz(-1.7526151) q[3];
sx q[3];
rz(1.3505107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36667103) q[0];
sx q[0];
rz(-1.3008302) q[0];
sx q[0];
rz(1.2694673) q[0];
rz(0.41995755) q[1];
sx q[1];
rz(-1.0008413) q[1];
sx q[1];
rz(-1.6023844) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29572648) q[0];
sx q[0];
rz(-0.55385607) q[0];
sx q[0];
rz(0.94828301) q[0];
rz(1.7544657) q[2];
sx q[2];
rz(-1.3546126) q[2];
sx q[2];
rz(0.64858299) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6091385) q[1];
sx q[1];
rz(-1.6241606) q[1];
sx q[1];
rz(-1.1751925) q[1];
rz(-pi) q[2];
x q[2];
rz(0.56681239) q[3];
sx q[3];
rz(-1.2897964) q[3];
sx q[3];
rz(2.5917201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.082108214) q[2];
sx q[2];
rz(-2.2993645) q[2];
sx q[2];
rz(-0.19169894) q[2];
rz(-0.55255237) q[3];
sx q[3];
rz(-1.5250165) q[3];
sx q[3];
rz(1.0244757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7315652) q[0];
sx q[0];
rz(-0.64842328) q[0];
sx q[0];
rz(-0.60436526) q[0];
rz(1.8961689) q[1];
sx q[1];
rz(-0.75029293) q[1];
sx q[1];
rz(0.85974685) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.031279292) q[0];
sx q[0];
rz(-2.1826751) q[0];
sx q[0];
rz(-0.41637929) q[0];
x q[1];
rz(0.014195125) q[2];
sx q[2];
rz(-1.6728587) q[2];
sx q[2];
rz(-1.0400187) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8670876) q[1];
sx q[1];
rz(-2.9531859) q[1];
sx q[1];
rz(2.8727358) q[1];
rz(2.4212461) q[3];
sx q[3];
rz(-1.3784153) q[3];
sx q[3];
rz(-1.9494353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8613646) q[2];
sx q[2];
rz(-2.1275529) q[2];
sx q[2];
rz(1.7029943) q[2];
rz(1.8493308) q[3];
sx q[3];
rz(-2.3137213) q[3];
sx q[3];
rz(2.0869702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0629405) q[0];
sx q[0];
rz(-2.8684454) q[0];
sx q[0];
rz(2.8237421) q[0];
rz(1.3392797) q[1];
sx q[1];
rz(-0.41792089) q[1];
sx q[1];
rz(0.97871614) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1237549) q[0];
sx q[0];
rz(-2.7394501) q[0];
sx q[0];
rz(0.077362424) q[0];
rz(2.9462325) q[2];
sx q[2];
rz(-1.6580515) q[2];
sx q[2];
rz(-1.7785787) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5086245) q[1];
sx q[1];
rz(-2.3787254) q[1];
sx q[1];
rz(-2.9588225) q[1];
x q[2];
rz(-1.7634912) q[3];
sx q[3];
rz(-1.3808625) q[3];
sx q[3];
rz(0.99270051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9432482) q[2];
sx q[2];
rz(-2.1965616) q[2];
sx q[2];
rz(1.1172969) q[2];
rz(2.9605588) q[3];
sx q[3];
rz(-1.2369316) q[3];
sx q[3];
rz(-1.9443289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4365874) q[0];
sx q[0];
rz(-0.26708189) q[0];
sx q[0];
rz(-1.8319112) q[0];
rz(-2.1197223) q[1];
sx q[1];
rz(-0.836687) q[1];
sx q[1];
rz(1.6530564) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2383136) q[0];
sx q[0];
rz(-2.2198703) q[0];
sx q[0];
rz(-0.42239503) q[0];
x q[1];
rz(2.8544442) q[2];
sx q[2];
rz(-2.4427772) q[2];
sx q[2];
rz(-0.1255807) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0206581) q[1];
sx q[1];
rz(-1.7132732) q[1];
sx q[1];
rz(-0.49748904) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11373489) q[3];
sx q[3];
rz(-1.0128504) q[3];
sx q[3];
rz(0.21765366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4788907) q[2];
sx q[2];
rz(-1.3871437) q[2];
sx q[2];
rz(2.877511) q[2];
rz(-1.4763907) q[3];
sx q[3];
rz(-2.4614406) q[3];
sx q[3];
rz(2.7217854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61178094) q[0];
sx q[0];
rz(-1.6407069) q[0];
sx q[0];
rz(-1.06426) q[0];
rz(3.0423959) q[1];
sx q[1];
rz(-1.3490889) q[1];
sx q[1];
rz(-1.4605716) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49584606) q[0];
sx q[0];
rz(-0.48862132) q[0];
sx q[0];
rz(0.13617985) q[0];
rz(-pi) q[1];
rz(1.869603) q[2];
sx q[2];
rz(-1.5332216) q[2];
sx q[2];
rz(2.7731098) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.32930347) q[1];
sx q[1];
rz(-1.9871622) q[1];
sx q[1];
rz(-1.8805481) q[1];
rz(-0.52525446) q[3];
sx q[3];
rz(-1.958985) q[3];
sx q[3];
rz(-1.3850016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3508241) q[2];
sx q[2];
rz(-2.0227573) q[2];
sx q[2];
rz(0.59949818) q[2];
rz(0.89764578) q[3];
sx q[3];
rz(-1.7220327) q[3];
sx q[3];
rz(0.038014855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-0.81383234) q[0];
sx q[0];
rz(-2.0812415) q[0];
sx q[0];
rz(1.4014442) q[0];
rz(0.071488149) q[1];
sx q[1];
rz(-1.46336) q[1];
sx q[1];
rz(2.1404526) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6342696) q[0];
sx q[0];
rz(-1.3183013) q[0];
sx q[0];
rz(-2.6802313) q[0];
x q[1];
rz(0.66364786) q[2];
sx q[2];
rz(-1.8816684) q[2];
sx q[2];
rz(0.03229095) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0467172) q[1];
sx q[1];
rz(-2.4946199) q[1];
sx q[1];
rz(2.9217974) q[1];
x q[2];
rz(-0.70434477) q[3];
sx q[3];
rz(-0.80432361) q[3];
sx q[3];
rz(-2.3856304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9643758) q[2];
sx q[2];
rz(-1.5528677) q[2];
sx q[2];
rz(2.4231518) q[2];
rz(-0.046772379) q[3];
sx q[3];
rz(-0.7074357) q[3];
sx q[3];
rz(-2.5717946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8319594) q[0];
sx q[0];
rz(-2.9871812) q[0];
sx q[0];
rz(-0.69044789) q[0];
rz(-2.9613103) q[1];
sx q[1];
rz(-2.0803662) q[1];
sx q[1];
rz(-0.32275018) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13234102) q[0];
sx q[0];
rz(-1.7242431) q[0];
sx q[0];
rz(1.4645534) q[0];
x q[1];
rz(-3.1351219) q[2];
sx q[2];
rz(-1.0186884) q[2];
sx q[2];
rz(3.0478549) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1230939) q[1];
sx q[1];
rz(-2.1040175) q[1];
sx q[1];
rz(-0.64792222) q[1];
x q[2];
rz(-0.68122562) q[3];
sx q[3];
rz(-2.109057) q[3];
sx q[3];
rz(-0.97506879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.47101578) q[2];
sx q[2];
rz(-1.7936423) q[2];
sx q[2];
rz(2.1237109) q[2];
rz(3.1380623) q[3];
sx q[3];
rz(-0.96960932) q[3];
sx q[3];
rz(-1.9297011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8135391) q[0];
sx q[0];
rz(-2.4040451) q[0];
sx q[0];
rz(2.2156583) q[0];
rz(-3.0042341) q[1];
sx q[1];
rz(-2.4691212) q[1];
sx q[1];
rz(-0.59991178) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5515251) q[0];
sx q[0];
rz(-0.54393629) q[0];
sx q[0];
rz(-0.005225709) q[0];
rz(-pi) q[1];
rz(0.7057759) q[2];
sx q[2];
rz(-0.76294357) q[2];
sx q[2];
rz(1.1102138) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9194738) q[1];
sx q[1];
rz(-1.6012234) q[1];
sx q[1];
rz(-0.98201507) q[1];
rz(-2.152485) q[3];
sx q[3];
rz(-0.55053655) q[3];
sx q[3];
rz(-1.4359695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7291193) q[2];
sx q[2];
rz(-1.0189265) q[2];
sx q[2];
rz(2.4998383) q[2];
rz(1.1563835) q[3];
sx q[3];
rz(-0.76120794) q[3];
sx q[3];
rz(-1.6183287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0625576) q[0];
sx q[0];
rz(-2.2618444) q[0];
sx q[0];
rz(-2.3656144) q[0];
rz(1.7017801) q[1];
sx q[1];
rz(-1.6701313) q[1];
sx q[1];
rz(-3.0727542) q[1];
rz(0.65563249) q[2];
sx q[2];
rz(-0.31972319) q[2];
sx q[2];
rz(2.2421851) q[2];
rz(-0.84550459) q[3];
sx q[3];
rz(-0.68544023) q[3];
sx q[3];
rz(-1.1123085) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
