OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.051209) q[0];
sx q[0];
rz(-0.78855711) q[0];
sx q[0];
rz(2.8330044) q[0];
rz(0.42449549) q[1];
sx q[1];
rz(4.2257809) q[1];
sx q[1];
rz(9.0492166) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5121927) q[0];
sx q[0];
rz(-1.6464656) q[0];
sx q[0];
rz(-1.2887495) q[0];
rz(-pi) q[1];
rz(-3.0388366) q[2];
sx q[2];
rz(-2.8381172) q[2];
sx q[2];
rz(-3.0927401) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.019770456) q[1];
sx q[1];
rz(-2.0149061) q[1];
sx q[1];
rz(-1.7239718) q[1];
rz(1.5442763) q[3];
sx q[3];
rz(-0.14355583) q[3];
sx q[3];
rz(-1.7100818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8033119) q[2];
sx q[2];
rz(-2.8901926) q[2];
sx q[2];
rz(-0.54331642) q[2];
rz(-1.4079037) q[3];
sx q[3];
rz(-1.7270154) q[3];
sx q[3];
rz(-0.89731115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5498085) q[0];
sx q[0];
rz(-1.9928638) q[0];
sx q[0];
rz(-0.59709221) q[0];
rz(-2.0789371) q[1];
sx q[1];
rz(-0.63908827) q[1];
sx q[1];
rz(-0.39892453) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0380533) q[0];
sx q[0];
rz(-2.1415496) q[0];
sx q[0];
rz(2.6817003) q[0];
x q[1];
rz(0.024256134) q[2];
sx q[2];
rz(-2.127671) q[2];
sx q[2];
rz(2.4766937) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0555206) q[1];
sx q[1];
rz(-0.87213736) q[1];
sx q[1];
rz(2.3585206) q[1];
rz(2.2116304) q[3];
sx q[3];
rz(-2.0130664) q[3];
sx q[3];
rz(1.363036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4661633) q[2];
sx q[2];
rz(-1.9025981) q[2];
sx q[2];
rz(-2.0002401) q[2];
rz(-0.91836786) q[3];
sx q[3];
rz(-0.71056241) q[3];
sx q[3];
rz(2.5520798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68701768) q[0];
sx q[0];
rz(-2.7641251) q[0];
sx q[0];
rz(-0.17534176) q[0];
rz(1.8270127) q[1];
sx q[1];
rz(-1.250123) q[1];
sx q[1];
rz(-2.4741727) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9475896) q[0];
sx q[0];
rz(-1.075248) q[0];
sx q[0];
rz(3.010514) q[0];
x q[1];
rz(2.5481648) q[2];
sx q[2];
rz(-1.9847514) q[2];
sx q[2];
rz(-2.7155795) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2167425) q[1];
sx q[1];
rz(-1.5713673) q[1];
sx q[1];
rz(-2.8257829) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5273058) q[3];
sx q[3];
rz(-1.7899872) q[3];
sx q[3];
rz(3.1295071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1158925) q[2];
sx q[2];
rz(-2.1860055) q[2];
sx q[2];
rz(2.058775) q[2];
rz(1.8358021) q[3];
sx q[3];
rz(-0.81597733) q[3];
sx q[3];
rz(0.79735565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1378655) q[0];
sx q[0];
rz(-0.19345134) q[0];
sx q[0];
rz(0.61099148) q[0];
rz(-2.5773279) q[1];
sx q[1];
rz(-1.0287501) q[1];
sx q[1];
rz(-2.4345523) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96933041) q[0];
sx q[0];
rz(-2.3670787) q[0];
sx q[0];
rz(-1.5277083) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66989278) q[2];
sx q[2];
rz(-0.50648738) q[2];
sx q[2];
rz(2.2746925) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0174071) q[1];
sx q[1];
rz(-1.1356135) q[1];
sx q[1];
rz(-2.6098677) q[1];
rz(-pi) q[2];
rz(1.8537997) q[3];
sx q[3];
rz(-1.0866829) q[3];
sx q[3];
rz(-3.1023539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6881312) q[2];
sx q[2];
rz(-0.42669272) q[2];
sx q[2];
rz(2.1235662) q[2];
rz(2.7536143) q[3];
sx q[3];
rz(-1.6951025) q[3];
sx q[3];
rz(-2.80262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46868789) q[0];
sx q[0];
rz(-2.0141116) q[0];
sx q[0];
rz(2.6197523) q[0];
rz(0.28779596) q[1];
sx q[1];
rz(-0.58240533) q[1];
sx q[1];
rz(-0.60307455) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1408242) q[0];
sx q[0];
rz(-0.10061564) q[0];
sx q[0];
rz(-2.6314497) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0126883) q[2];
sx q[2];
rz(-2.1753722) q[2];
sx q[2];
rz(2.9964436) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6213563) q[1];
sx q[1];
rz(-1.4870315) q[1];
sx q[1];
rz(0.4256953) q[1];
rz(-pi) q[2];
rz(2.2207845) q[3];
sx q[3];
rz(-0.42144708) q[3];
sx q[3];
rz(1.8325266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7120984) q[2];
sx q[2];
rz(-2.3488729) q[2];
sx q[2];
rz(-2.402044) q[2];
rz(-1.0150917) q[3];
sx q[3];
rz(-2.6138217) q[3];
sx q[3];
rz(0.019006193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58448142) q[0];
sx q[0];
rz(-2.8040573) q[0];
sx q[0];
rz(0.8648411) q[0];
rz(2.5014014) q[1];
sx q[1];
rz(-2.4914111) q[1];
sx q[1];
rz(2.6472299) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1304754) q[0];
sx q[0];
rz(-1.1251161) q[0];
sx q[0];
rz(2.9167239) q[0];
x q[1];
rz(0.041448822) q[2];
sx q[2];
rz(-0.96872708) q[2];
sx q[2];
rz(-0.96884251) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3013005) q[1];
sx q[1];
rz(-1.5363785) q[1];
sx q[1];
rz(0.057386668) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0406144) q[3];
sx q[3];
rz(-1.1270071) q[3];
sx q[3];
rz(0.63999301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2690027) q[2];
sx q[2];
rz(-0.63429093) q[2];
sx q[2];
rz(2.1799901) q[2];
rz(1.5087992) q[3];
sx q[3];
rz(-1.5140994) q[3];
sx q[3];
rz(0.076920286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85840571) q[0];
sx q[0];
rz(-2.7981693) q[0];
sx q[0];
rz(-1.0431694) q[0];
rz(0.6768325) q[1];
sx q[1];
rz(-2.4969641) q[1];
sx q[1];
rz(2.4136995) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60715577) q[0];
sx q[0];
rz(-0.42947665) q[0];
sx q[0];
rz(1.3017824) q[0];
x q[1];
rz(2.3995326) q[2];
sx q[2];
rz(-1.5734317) q[2];
sx q[2];
rz(-1.9085086) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.54503578) q[1];
sx q[1];
rz(-1.2378197) q[1];
sx q[1];
rz(0.92988981) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2292101) q[3];
sx q[3];
rz(-1.0083958) q[3];
sx q[3];
rz(-1.7096568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4595043) q[2];
sx q[2];
rz(-2.3900718) q[2];
sx q[2];
rz(-2.6085243) q[2];
rz(1.3047949) q[3];
sx q[3];
rz(-2.6663836) q[3];
sx q[3];
rz(-2.3615725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.042453658) q[0];
sx q[0];
rz(-3.070153) q[0];
sx q[0];
rz(3.1169917) q[0];
rz(-0.050361659) q[1];
sx q[1];
rz(-0.55614007) q[1];
sx q[1];
rz(0.76739001) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8883945) q[0];
sx q[0];
rz(-2.1108187) q[0];
sx q[0];
rz(1.5408976) q[0];
rz(2.6219756) q[2];
sx q[2];
rz(-2.765968) q[2];
sx q[2];
rz(-0.56359282) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.142006) q[1];
sx q[1];
rz(-2.1697727) q[1];
sx q[1];
rz(2.0525648) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9629836) q[3];
sx q[3];
rz(-1.0309891) q[3];
sx q[3];
rz(1.7843877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6699827) q[2];
sx q[2];
rz(-2.6675384) q[2];
sx q[2];
rz(-1.6101884) q[2];
rz(-2.1286185) q[3];
sx q[3];
rz(-1.6771202) q[3];
sx q[3];
rz(2.5392635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7811964) q[0];
sx q[0];
rz(-2.4735326) q[0];
sx q[0];
rz(0.49041954) q[0];
rz(-0.40052739) q[1];
sx q[1];
rz(-0.35919765) q[1];
sx q[1];
rz(-1.5708956) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2119031) q[0];
sx q[0];
rz(-2.4725998) q[0];
sx q[0];
rz(2.5446135) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2068858) q[2];
sx q[2];
rz(-1.2635974) q[2];
sx q[2];
rz(-0.72774923) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1410349) q[1];
sx q[1];
rz(-0.36492202) q[1];
sx q[1];
rz(1.8561897) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3491507) q[3];
sx q[3];
rz(-2.1413448) q[3];
sx q[3];
rz(1.5287409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5639497) q[2];
sx q[2];
rz(-0.22259139) q[2];
sx q[2];
rz(0.27601784) q[2];
rz(0.66463071) q[3];
sx q[3];
rz(-2.1641927) q[3];
sx q[3];
rz(0.19653921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4589602) q[0];
sx q[0];
rz(-0.17550547) q[0];
sx q[0];
rz(-1.578414) q[0];
rz(0.75421929) q[1];
sx q[1];
rz(-1.027532) q[1];
sx q[1];
rz(-3.077502) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10326021) q[0];
sx q[0];
rz(-2.0287559) q[0];
sx q[0];
rz(0.026491212) q[0];
rz(2.5245177) q[2];
sx q[2];
rz(-1.1960746) q[2];
sx q[2];
rz(2.9385757) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8386993) q[1];
sx q[1];
rz(-2.6461227) q[1];
sx q[1];
rz(-1.7458781) q[1];
rz(-pi) q[2];
rz(-2.1782297) q[3];
sx q[3];
rz(-0.86403908) q[3];
sx q[3];
rz(-0.53234311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2854707) q[2];
sx q[2];
rz(-2.8557114) q[2];
sx q[2];
rz(0.90606436) q[2];
rz(-1.9237349) q[3];
sx q[3];
rz(-1.2713615) q[3];
sx q[3];
rz(2.0876032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8080407) q[0];
sx q[0];
rz(-1.522949) q[0];
sx q[0];
rz(-0.20725313) q[0];
rz(2.8605657) q[1];
sx q[1];
rz(-1.3916176) q[1];
sx q[1];
rz(2.4236046) q[1];
rz(-1.4886552) q[2];
sx q[2];
rz(-1.1700656) q[2];
sx q[2];
rz(-1.1943371) q[2];
rz(-1.5850375) q[3];
sx q[3];
rz(-1.0938627) q[3];
sx q[3];
rz(0.59577019) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
