OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6501939) q[0];
sx q[0];
rz(-2.8770652) q[0];
sx q[0];
rz(-2.7471623) q[0];
rz(0.0061622942) q[1];
sx q[1];
rz(-0.34024629) q[1];
sx q[1];
rz(1.9415829) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7317176) q[0];
sx q[0];
rz(-3.0933676) q[0];
sx q[0];
rz(-1.4531141) q[0];
rz(-pi) q[1];
rz(-2.1985487) q[2];
sx q[2];
rz(-0.53484166) q[2];
sx q[2];
rz(0.56378555) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.66558054) q[1];
sx q[1];
rz(-1.1886547) q[1];
sx q[1];
rz(-1.0784472) q[1];
rz(-pi) q[2];
rz(1.5267738) q[3];
sx q[3];
rz(-1.9485954) q[3];
sx q[3];
rz(1.6954741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3068984) q[2];
sx q[2];
rz(-1.4935741) q[2];
sx q[2];
rz(-2.8519894) q[2];
rz(-2.2662207) q[3];
sx q[3];
rz(-1.0062199) q[3];
sx q[3];
rz(-3.0818821) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61525476) q[0];
sx q[0];
rz(-0.86741388) q[0];
sx q[0];
rz(-2.7976024) q[0];
rz(0.084331766) q[1];
sx q[1];
rz(-2.4721959) q[1];
sx q[1];
rz(1.7864236) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8665168) q[0];
sx q[0];
rz(-1.5059024) q[0];
sx q[0];
rz(0.95922031) q[0];
rz(-pi) q[1];
rz(-2.9789574) q[2];
sx q[2];
rz(-1.4504823) q[2];
sx q[2];
rz(-1.6441117) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5835186) q[1];
sx q[1];
rz(-1.1221702) q[1];
sx q[1];
rz(3.0706057) q[1];
x q[2];
rz(1.9676898) q[3];
sx q[3];
rz(-0.19860425) q[3];
sx q[3];
rz(1.7484776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8460059) q[2];
sx q[2];
rz(-2.3082374) q[2];
sx q[2];
rz(-0.53768349) q[2];
rz(-0.50283557) q[3];
sx q[3];
rz(-0.42481315) q[3];
sx q[3];
rz(-2.1311549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66353345) q[0];
sx q[0];
rz(-2.6038267) q[0];
sx q[0];
rz(0.64087254) q[0];
rz(0.74869853) q[1];
sx q[1];
rz(-2.1332032) q[1];
sx q[1];
rz(2.0764988) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0707866) q[0];
sx q[0];
rz(-0.35937989) q[0];
sx q[0];
rz(1.552836) q[0];
rz(-pi) q[1];
rz(-0.89945729) q[2];
sx q[2];
rz(-0.81842917) q[2];
sx q[2];
rz(-0.98086548) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6667337) q[1];
sx q[1];
rz(-1.6068646) q[1];
sx q[1];
rz(-0.95400793) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.25032708) q[3];
sx q[3];
rz(-0.86391376) q[3];
sx q[3];
rz(-1.7771429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.37725267) q[2];
sx q[2];
rz(-1.179402) q[2];
sx q[2];
rz(0.71737814) q[2];
rz(2.453089) q[3];
sx q[3];
rz(-0.62763667) q[3];
sx q[3];
rz(-2.8440516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82729572) q[0];
sx q[0];
rz(-1.2001487) q[0];
sx q[0];
rz(2.8919019) q[0];
rz(-2.1266134) q[1];
sx q[1];
rz(-2.8453638) q[1];
sx q[1];
rz(-3.1304741) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.714689) q[0];
sx q[0];
rz(-2.5137797) q[0];
sx q[0];
rz(-2.3054302) q[0];
rz(-2.501802) q[2];
sx q[2];
rz(-0.94546972) q[2];
sx q[2];
rz(-1.421979) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8324082) q[1];
sx q[1];
rz(-2.5384181) q[1];
sx q[1];
rz(0.96997728) q[1];
x q[2];
rz(2.6075881) q[3];
sx q[3];
rz(-0.36557331) q[3];
sx q[3];
rz(3.0766069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.49542385) q[2];
sx q[2];
rz(-1.8858706) q[2];
sx q[2];
rz(-0.28309506) q[2];
rz(0.66343534) q[3];
sx q[3];
rz(-2.5484271) q[3];
sx q[3];
rz(2.2535113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0304612) q[0];
sx q[0];
rz(-2.9009394) q[0];
sx q[0];
rz(-2.8097613) q[0];
rz(2.6470673) q[1];
sx q[1];
rz(-1.2978413) q[1];
sx q[1];
rz(-1.3269075) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3834284) q[0];
sx q[0];
rz(-0.82793068) q[0];
sx q[0];
rz(1.3616256) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3438498) q[2];
sx q[2];
rz(-2.0125055) q[2];
sx q[2];
rz(1.9517348) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.27969589) q[1];
sx q[1];
rz(-1.8021291) q[1];
sx q[1];
rz(-0.33613236) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9782449) q[3];
sx q[3];
rz(-1.3031928) q[3];
sx q[3];
rz(-2.101055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1422687) q[2];
sx q[2];
rz(-0.48196718) q[2];
sx q[2];
rz(-1.2456606) q[2];
rz(-1.2549531) q[3];
sx q[3];
rz(-0.84147036) q[3];
sx q[3];
rz(-0.83827034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44928837) q[0];
sx q[0];
rz(-3.1384387) q[0];
sx q[0];
rz(-0.6814878) q[0];
rz(0.20755126) q[1];
sx q[1];
rz(-2.6782268) q[1];
sx q[1];
rz(2.025827) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.241248) q[0];
sx q[0];
rz(-1.9755409) q[0];
sx q[0];
rz(0.83165283) q[0];
rz(-pi) q[1];
rz(1.0461147) q[2];
sx q[2];
rz(-0.74138734) q[2];
sx q[2];
rz(1.807883) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4286194) q[1];
sx q[1];
rz(-1.553112) q[1];
sx q[1];
rz(0.58847217) q[1];
rz(-1.7498155) q[3];
sx q[3];
rz(-2.2324623) q[3];
sx q[3];
rz(-2.9359407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9289124) q[2];
sx q[2];
rz(-1.8447515) q[2];
sx q[2];
rz(2.7872655) q[2];
rz(2.8220693) q[3];
sx q[3];
rz(-1.9890246) q[3];
sx q[3];
rz(2.7697146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1324683) q[0];
sx q[0];
rz(-0.23129825) q[0];
sx q[0];
rz(-2.4672467) q[0];
rz(2.0293503) q[1];
sx q[1];
rz(-2.4770885) q[1];
sx q[1];
rz(-0.56232125) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4903957) q[0];
sx q[0];
rz(-0.90264747) q[0];
sx q[0];
rz(0.018880318) q[0];
x q[1];
rz(1.6717031) q[2];
sx q[2];
rz(-1.9011874) q[2];
sx q[2];
rz(-0.23320564) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.454969) q[1];
sx q[1];
rz(-1.3775871) q[1];
sx q[1];
rz(-1.4229694) q[1];
rz(3.0543442) q[3];
sx q[3];
rz(-0.97594075) q[3];
sx q[3];
rz(1.113254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.101863) q[2];
sx q[2];
rz(-2.1477284) q[2];
sx q[2];
rz(-2.8015461) q[2];
rz(-0.21751054) q[3];
sx q[3];
rz(-0.9483996) q[3];
sx q[3];
rz(0.31869179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44889221) q[0];
sx q[0];
rz(-3.0954439) q[0];
sx q[0];
rz(-0.39644077) q[0];
rz(3.0026644) q[1];
sx q[1];
rz(-0.46008343) q[1];
sx q[1];
rz(-1.5213535) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5312708) q[0];
sx q[0];
rz(-1.5025508) q[0];
sx q[0];
rz(2.9444429) q[0];
x q[1];
rz(0.40201681) q[2];
sx q[2];
rz(-2.0094299) q[2];
sx q[2];
rz(-1.5356262) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.13008598) q[1];
sx q[1];
rz(-2.3475921) q[1];
sx q[1];
rz(1.4873051) q[1];
rz(-pi) q[2];
x q[2];
rz(1.429854) q[3];
sx q[3];
rz(-0.68416506) q[3];
sx q[3];
rz(-1.9977026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.56269318) q[2];
sx q[2];
rz(-1.0572628) q[2];
sx q[2];
rz(2.8472624) q[2];
rz(-2.0108022) q[3];
sx q[3];
rz(-1.3785988) q[3];
sx q[3];
rz(2.14595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49333736) q[0];
sx q[0];
rz(-2.2264037) q[0];
sx q[0];
rz(-2.7822568) q[0];
rz(-2.1954779) q[1];
sx q[1];
rz(-0.39603907) q[1];
sx q[1];
rz(-0.27063453) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.04979241) q[0];
sx q[0];
rz(-1.5349749) q[0];
sx q[0];
rz(1.6019751) q[0];
rz(-2.052202) q[2];
sx q[2];
rz(-2.319616) q[2];
sx q[2];
rz(-2.7871728) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2808025) q[1];
sx q[1];
rz(-1.5250912) q[1];
sx q[1];
rz(1.6153107) q[1];
x q[2];
rz(-0.98468303) q[3];
sx q[3];
rz(-0.95041785) q[3];
sx q[3];
rz(-1.8970722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3140807) q[2];
sx q[2];
rz(-0.22958799) q[2];
sx q[2];
rz(0.45544004) q[2];
rz(-0.81196249) q[3];
sx q[3];
rz(-2.2596695) q[3];
sx q[3];
rz(2.5922095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51046002) q[0];
sx q[0];
rz(-1.6332508) q[0];
sx q[0];
rz(0.73927885) q[0];
rz(2.9108858) q[1];
sx q[1];
rz(-2.6735327) q[1];
sx q[1];
rz(0.49490067) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98289821) q[0];
sx q[0];
rz(-1.6086846) q[0];
sx q[0];
rz(-2.057103) q[0];
rz(-pi) q[1];
rz(1.5636744) q[2];
sx q[2];
rz(-1.2429503) q[2];
sx q[2];
rz(2.7210826) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.478185) q[1];
sx q[1];
rz(-2.2067398) q[1];
sx q[1];
rz(3.0031167) q[1];
rz(-pi) q[2];
rz(0.87402113) q[3];
sx q[3];
rz(-1.4975944) q[3];
sx q[3];
rz(-0.30833581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.13359244) q[2];
sx q[2];
rz(-2.0501037) q[2];
sx q[2];
rz(-2.8137394) q[2];
rz(2.949529) q[3];
sx q[3];
rz(-0.24774194) q[3];
sx q[3];
rz(2.1081934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.0192169) q[0];
sx q[0];
rz(-1.6288971) q[0];
sx q[0];
rz(-1.6678641) q[0];
rz(-0.83256759) q[1];
sx q[1];
rz(-1.4214129) q[1];
sx q[1];
rz(1.3690154) q[1];
rz(-2.9359948) q[2];
sx q[2];
rz(-1.0369221) q[2];
sx q[2];
rz(2.592929) q[2];
rz(-0.26631793) q[3];
sx q[3];
rz(-1.3238293) q[3];
sx q[3];
rz(-2.7662591) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
