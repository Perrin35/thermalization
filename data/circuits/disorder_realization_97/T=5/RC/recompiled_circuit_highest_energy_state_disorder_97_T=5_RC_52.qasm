OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.834637) q[0];
sx q[0];
rz(-0.87943465) q[0];
sx q[0];
rz(2.6479794) q[0];
rz(1.8618795) q[1];
sx q[1];
rz(-0.6266098) q[1];
sx q[1];
rz(-1.6630747) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0042838) q[0];
sx q[0];
rz(-2.2861436) q[0];
sx q[0];
rz(2.2527931) q[0];
rz(-pi) q[1];
rz(2.2343687) q[2];
sx q[2];
rz(-1.3665082) q[2];
sx q[2];
rz(2.8641618) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.263674) q[1];
sx q[1];
rz(-1.3527737) q[1];
sx q[1];
rz(-2.430357) q[1];
x q[2];
rz(-2.0527203) q[3];
sx q[3];
rz(-0.69447172) q[3];
sx q[3];
rz(-2.0884909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0066068) q[2];
sx q[2];
rz(-1.1400433) q[2];
sx q[2];
rz(-2.0331649) q[2];
rz(-1.3125575) q[3];
sx q[3];
rz(-2.8969942) q[3];
sx q[3];
rz(0.31497064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7600064) q[0];
sx q[0];
rz(-2.3158323) q[0];
sx q[0];
rz(-0.015856892) q[0];
rz(0.99041692) q[1];
sx q[1];
rz(-0.76737338) q[1];
sx q[1];
rz(2.2784065) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47967692) q[0];
sx q[0];
rz(-2.4004395) q[0];
sx q[0];
rz(-2.9984498) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.241774) q[2];
sx q[2];
rz(-1.8890427) q[2];
sx q[2];
rz(1.1877738) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6241315) q[1];
sx q[1];
rz(-1.3966113) q[1];
sx q[1];
rz(2.6427173) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7402418) q[3];
sx q[3];
rz(-0.72899216) q[3];
sx q[3];
rz(0.21641009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.70011675) q[2];
sx q[2];
rz(-2.8642544) q[2];
sx q[2];
rz(0.46879834) q[2];
rz(0.68764728) q[3];
sx q[3];
rz(-1.5117437) q[3];
sx q[3];
rz(2.8414753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-1.0013292) q[0];
sx q[0];
rz(-2.2183473) q[0];
sx q[0];
rz(-0.73626751) q[0];
rz(-2.7983792) q[1];
sx q[1];
rz(-2.2540269) q[1];
sx q[1];
rz(0.18377486) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60584468) q[0];
sx q[0];
rz(-2.3130406) q[0];
sx q[0];
rz(1.6607619) q[0];
rz(-1.0509346) q[2];
sx q[2];
rz(-2.1987872) q[2];
sx q[2];
rz(1.7597212) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.50810516) q[1];
sx q[1];
rz(-0.87292665) q[1];
sx q[1];
rz(1.9644323) q[1];
rz(-0.36005693) q[3];
sx q[3];
rz(-1.0601794) q[3];
sx q[3];
rz(-2.1329481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.78202128) q[2];
sx q[2];
rz(-2.6247793) q[2];
sx q[2];
rz(2.6256631) q[2];
rz(-1.1082209) q[3];
sx q[3];
rz(-2.5907232) q[3];
sx q[3];
rz(1.9903323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16294031) q[0];
sx q[0];
rz(-1.2115703) q[0];
sx q[0];
rz(-0.51280713) q[0];
rz(-2.6817952) q[1];
sx q[1];
rz(-1.5922981) q[1];
sx q[1];
rz(2.3777681) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0054454) q[0];
sx q[0];
rz(-0.17373611) q[0];
sx q[0];
rz(2.4197015) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.725179) q[2];
sx q[2];
rz(-1.8685307) q[2];
sx q[2];
rz(2.1418051) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3556174) q[1];
sx q[1];
rz(-1.2316868) q[1];
sx q[1];
rz(0.29735844) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4277589) q[3];
sx q[3];
rz(-1.9120354) q[3];
sx q[3];
rz(2.9590497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0189455) q[2];
sx q[2];
rz(-0.15572369) q[2];
sx q[2];
rz(2.8141008) q[2];
rz(-0.28199768) q[3];
sx q[3];
rz(-1.3982541) q[3];
sx q[3];
rz(1.5544844) q[3];
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
rz(pi/2) q[3];
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
rz(0.17664385) q[0];
sx q[0];
rz(-2.7556941) q[0];
sx q[0];
rz(-2.1642245) q[0];
rz(0.36879677) q[1];
sx q[1];
rz(-0.86124033) q[1];
sx q[1];
rz(-1.4126973) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0127129) q[0];
sx q[0];
rz(-1.0815718) q[0];
sx q[0];
rz(-2.9262744) q[0];
rz(-pi) q[1];
rz(2.0875889) q[2];
sx q[2];
rz(-1.5828504) q[2];
sx q[2];
rz(-0.64753676) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.818644) q[1];
sx q[1];
rz(-1.6020157) q[1];
sx q[1];
rz(1.853456) q[1];
rz(-pi) q[2];
rz(1.2431954) q[3];
sx q[3];
rz(-2.1246111) q[3];
sx q[3];
rz(-2.6462951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8480924) q[2];
sx q[2];
rz(-0.81587452) q[2];
sx q[2];
rz(-0.14981848) q[2];
rz(-2.7430429) q[3];
sx q[3];
rz(-1.5236676) q[3];
sx q[3];
rz(-2.985305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77001101) q[0];
sx q[0];
rz(-3.0122029) q[0];
sx q[0];
rz(0.55321252) q[0];
rz(-0.54168701) q[1];
sx q[1];
rz(-1.0463511) q[1];
sx q[1];
rz(0.4479301) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49419379) q[0];
sx q[0];
rz(-1.9564391) q[0];
sx q[0];
rz(1.7269887) q[0];
rz(-pi) q[1];
x q[1];
rz(0.42945736) q[2];
sx q[2];
rz(-1.1633368) q[2];
sx q[2];
rz(2.7384659) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.74400157) q[1];
sx q[1];
rz(-0.41644704) q[1];
sx q[1];
rz(-1.2214425) q[1];
rz(-pi) q[2];
rz(0.33282307) q[3];
sx q[3];
rz(-1.5170005) q[3];
sx q[3];
rz(-1.4383565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.64122671) q[2];
sx q[2];
rz(-1.1643103) q[2];
sx q[2];
rz(0.79916239) q[2];
rz(0.3847807) q[3];
sx q[3];
rz(-0.32618263) q[3];
sx q[3];
rz(1.0001812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36050972) q[0];
sx q[0];
rz(-1.0436844) q[0];
sx q[0];
rz(0.59180301) q[0];
rz(0.38161713) q[1];
sx q[1];
rz(-2.7615669) q[1];
sx q[1];
rz(2.9891678) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7397933) q[0];
sx q[0];
rz(-1.2619979) q[0];
sx q[0];
rz(-1.0303506) q[0];
rz(-2.811829) q[2];
sx q[2];
rz(-0.50055365) q[2];
sx q[2];
rz(-0.19937521) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3445216) q[1];
sx q[1];
rz(-2.7943369) q[1];
sx q[1];
rz(-0.053649112) q[1];
rz(-pi) q[2];
rz(0.8922718) q[3];
sx q[3];
rz(-2.0990058) q[3];
sx q[3];
rz(-1.5367374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3064208) q[2];
sx q[2];
rz(-2.1495337) q[2];
sx q[2];
rz(-1.6888118) q[2];
rz(-0.49550223) q[3];
sx q[3];
rz(-0.82363868) q[3];
sx q[3];
rz(-2.5775094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
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
rz(2.751916) q[0];
sx q[0];
rz(-0.037411995) q[0];
sx q[0];
rz(0.16146846) q[0];
rz(-3.1096733) q[1];
sx q[1];
rz(-0.64161623) q[1];
sx q[1];
rz(-1.2772824) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.227018) q[0];
sx q[0];
rz(-0.76138568) q[0];
sx q[0];
rz(-2.7622591) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0525682) q[2];
sx q[2];
rz(-0.46899295) q[2];
sx q[2];
rz(-0.8587786) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5885791) q[1];
sx q[1];
rz(-1.2672851) q[1];
sx q[1];
rz(1.4576333) q[1];
rz(3.1031514) q[3];
sx q[3];
rz(-1.2769967) q[3];
sx q[3];
rz(-2.163559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6892467) q[2];
sx q[2];
rz(-0.9114868) q[2];
sx q[2];
rz(-0.39608836) q[2];
rz(-0.39997697) q[3];
sx q[3];
rz(-0.59498274) q[3];
sx q[3];
rz(2.2649435) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9006627) q[0];
sx q[0];
rz(-0.018095896) q[0];
sx q[0];
rz(2.990429) q[0];
rz(-2.1647272) q[1];
sx q[1];
rz(-2.7604389) q[1];
sx q[1];
rz(-2.3968598) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1158651) q[0];
sx q[0];
rz(-0.34894279) q[0];
sx q[0];
rz(-2.2141488) q[0];
rz(-pi) q[1];
rz(1.9836203) q[2];
sx q[2];
rz(-1.7016439) q[2];
sx q[2];
rz(-1.3589588) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3754282) q[1];
sx q[1];
rz(-0.64830983) q[1];
sx q[1];
rz(0.0023771087) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3296534) q[3];
sx q[3];
rz(-1.9134812) q[3];
sx q[3];
rz(-1.8543275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5186844) q[2];
sx q[2];
rz(-1.238995) q[2];
sx q[2];
rz(2.1935479) q[2];
rz(-1.0575804) q[3];
sx q[3];
rz(-1.6653929) q[3];
sx q[3];
rz(2.0675366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0483765) q[0];
sx q[0];
rz(-0.26234782) q[0];
sx q[0];
rz(2.9076305) q[0];
rz(1.9562862) q[1];
sx q[1];
rz(-0.92820853) q[1];
sx q[1];
rz(-1.8918461) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3435182) q[0];
sx q[0];
rz(-1.2815223) q[0];
sx q[0];
rz(-0.57855655) q[0];
x q[1];
rz(0.28567254) q[2];
sx q[2];
rz(-0.64027363) q[2];
sx q[2];
rz(-1.7264186) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1162303) q[1];
sx q[1];
rz(-0.90060189) q[1];
sx q[1];
rz(-0.7109364) q[1];
rz(-pi) q[2];
rz(2.9261111) q[3];
sx q[3];
rz(-0.69284791) q[3];
sx q[3];
rz(-2.1699435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.027111) q[2];
sx q[2];
rz(-1.9805084) q[2];
sx q[2];
rz(-1.7005881) q[2];
rz(3.0103185) q[3];
sx q[3];
rz(-0.461853) q[3];
sx q[3];
rz(-2.89768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.048653614) q[0];
sx q[0];
rz(-0.96554148) q[0];
sx q[0];
rz(1.1337793) q[0];
rz(-1.0406021) q[1];
sx q[1];
rz(-2.2941209) q[1];
sx q[1];
rz(2.6170731) q[1];
rz(1.9807627) q[2];
sx q[2];
rz(-0.66902918) q[2];
sx q[2];
rz(1.4794028) q[2];
rz(0.63325044) q[3];
sx q[3];
rz(-1.6926822) q[3];
sx q[3];
rz(2.5183986) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
