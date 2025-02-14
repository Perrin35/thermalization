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
rz(-0.96002785) q[0];
sx q[0];
rz(-2.5224944) q[0];
sx q[0];
rz(1.7382789) q[0];
rz(2.4461441) q[1];
sx q[1];
rz(-0.56013501) q[1];
sx q[1];
rz(0.6518031) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2346828) q[0];
sx q[0];
rz(-1.2359338) q[0];
sx q[0];
rz(0.85483179) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5648834) q[2];
sx q[2];
rz(-1.2893218) q[2];
sx q[2];
rz(-0.23879063) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5320325) q[1];
sx q[1];
rz(-1.802562) q[1];
sx q[1];
rz(1.547481) q[1];
x q[2];
rz(-1.8599717) q[3];
sx q[3];
rz(-0.77176982) q[3];
sx q[3];
rz(-1.4909596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.92242509) q[2];
sx q[2];
rz(-1.7520011) q[2];
sx q[2];
rz(0.31828848) q[2];
rz(-2.2394032) q[3];
sx q[3];
rz(-0.92034322) q[3];
sx q[3];
rz(2.2778146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1282463) q[0];
sx q[0];
rz(-2.8885169) q[0];
sx q[0];
rz(0.52023649) q[0];
rz(3.0546313) q[1];
sx q[1];
rz(-2.088701) q[1];
sx q[1];
rz(-2.3435977) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8142321) q[0];
sx q[0];
rz(-1.8810187) q[0];
sx q[0];
rz(-1.7807175) q[0];
rz(-pi) q[1];
rz(-0.17391674) q[2];
sx q[2];
rz(-0.7510646) q[2];
sx q[2];
rz(-0.30663315) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6317111) q[1];
sx q[1];
rz(-2.7260511) q[1];
sx q[1];
rz(-2.7870537) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6949072) q[3];
sx q[3];
rz(-1.2036714) q[3];
sx q[3];
rz(2.0172831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.308455) q[2];
sx q[2];
rz(-1.1546346) q[2];
sx q[2];
rz(2.7412565) q[2];
rz(2.9588302) q[3];
sx q[3];
rz(-0.57923135) q[3];
sx q[3];
rz(1.8496752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0188145) q[0];
sx q[0];
rz(-3.0570539) q[0];
sx q[0];
rz(1.9345181) q[0];
rz(1.5030376) q[1];
sx q[1];
rz(-1.259558) q[1];
sx q[1];
rz(-2.9685793) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6956222) q[0];
sx q[0];
rz(-2.7467331) q[0];
sx q[0];
rz(-1.7941712) q[0];
rz(-pi) q[1];
rz(0.11008693) q[2];
sx q[2];
rz(-0.99186691) q[2];
sx q[2];
rz(0.26012173) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.26721482) q[1];
sx q[1];
rz(-2.3372531) q[1];
sx q[1];
rz(-0.72093236) q[1];
rz(-pi) q[2];
rz(-2.9152619) q[3];
sx q[3];
rz(-0.76068288) q[3];
sx q[3];
rz(2.5193391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5485237) q[2];
sx q[2];
rz(-2.5544281) q[2];
sx q[2];
rz(-2.1068088) q[2];
rz(-0.89094025) q[3];
sx q[3];
rz(-0.8689298) q[3];
sx q[3];
rz(1.2263891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064875038) q[0];
sx q[0];
rz(-2.1893976) q[0];
sx q[0];
rz(0.12884831) q[0];
rz(-2.7128362) q[1];
sx q[1];
rz(-1.7196451) q[1];
sx q[1];
rz(-1.6645924) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77177202) q[0];
sx q[0];
rz(-1.5025398) q[0];
sx q[0];
rz(2.7959156) q[0];
rz(-1.1294133) q[2];
sx q[2];
rz(-1.980482) q[2];
sx q[2];
rz(0.022155174) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4739371) q[1];
sx q[1];
rz(-1.154711) q[1];
sx q[1];
rz(2.9101084) q[1];
x q[2];
rz(-1.5691532) q[3];
sx q[3];
rz(-0.59125671) q[3];
sx q[3];
rz(1.7776684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.19026549) q[2];
sx q[2];
rz(-1.4172047) q[2];
sx q[2];
rz(1.3235271) q[2];
rz(1.6773112) q[3];
sx q[3];
rz(-2.5530294) q[3];
sx q[3];
rz(-2.467449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4071963) q[0];
sx q[0];
rz(-0.36235991) q[0];
sx q[0];
rz(1.6998442) q[0];
rz(2.7184519) q[1];
sx q[1];
rz(-2.1839881) q[1];
sx q[1];
rz(-0.29022455) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5231756) q[0];
sx q[0];
rz(-1.2783861) q[0];
sx q[0];
rz(-0.66177701) q[0];
rz(-0.96137233) q[2];
sx q[2];
rz(-1.0842241) q[2];
sx q[2];
rz(-0.79604641) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4729974) q[1];
sx q[1];
rz(-1.5221704) q[1];
sx q[1];
rz(0.12811382) q[1];
rz(-pi) q[2];
rz(2.4423744) q[3];
sx q[3];
rz(-2.7479246) q[3];
sx q[3];
rz(0.75692219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.71771249) q[2];
sx q[2];
rz(-0.78705698) q[2];
sx q[2];
rz(0.94669739) q[2];
rz(2.8801019) q[3];
sx q[3];
rz(-1.3638834) q[3];
sx q[3];
rz(0.3869032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70447266) q[0];
sx q[0];
rz(-1.3334714) q[0];
sx q[0];
rz(-0.76960808) q[0];
rz(-0.17598027) q[1];
sx q[1];
rz(-1.8006005) q[1];
sx q[1];
rz(1.5716858) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.044445) q[0];
sx q[0];
rz(-2.189069) q[0];
sx q[0];
rz(0.78938214) q[0];
rz(-pi) q[1];
x q[1];
rz(0.069934037) q[2];
sx q[2];
rz(-0.85288793) q[2];
sx q[2];
rz(-0.96558918) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.6750727) q[1];
sx q[1];
rz(-0.41374712) q[1];
sx q[1];
rz(0.4725666) q[1];
rz(-0.50521781) q[3];
sx q[3];
rz(-1.1105624) q[3];
sx q[3];
rz(1.203618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0104597) q[2];
sx q[2];
rz(-0.95644462) q[2];
sx q[2];
rz(-1.825911) q[2];
rz(-1.2223505) q[3];
sx q[3];
rz(-1.1244011) q[3];
sx q[3];
rz(-2.840672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-2.9199801) q[0];
sx q[0];
rz(-1.4497117) q[0];
sx q[0];
rz(1.8020887) q[0];
rz(-1.9446531) q[1];
sx q[1];
rz(-2.0868128) q[1];
sx q[1];
rz(2.1766677) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5366936) q[0];
sx q[0];
rz(-1.6654547) q[0];
sx q[0];
rz(2.6306334) q[0];
x q[1];
rz(-2.8469725) q[2];
sx q[2];
rz(-1.9082532) q[2];
sx q[2];
rz(-1.1883512) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9274018) q[1];
sx q[1];
rz(-1.3510696) q[1];
sx q[1];
rz(-1.6564547) q[1];
x q[2];
rz(-2.6581702) q[3];
sx q[3];
rz(-1.7258378) q[3];
sx q[3];
rz(-2.8446236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7397466) q[2];
sx q[2];
rz(-1.9332989) q[2];
sx q[2];
rz(-1.9604663) q[2];
rz(0.014569672) q[3];
sx q[3];
rz(-2.7404692) q[3];
sx q[3];
rz(-1.9692263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3907617) q[0];
sx q[0];
rz(-1.9347235) q[0];
sx q[0];
rz(2.3857351) q[0];
rz(-2.4769056) q[1];
sx q[1];
rz(-1.1988147) q[1];
sx q[1];
rz(-1.2896779) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3493746) q[0];
sx q[0];
rz(-1.5070033) q[0];
sx q[0];
rz(2.347058) q[0];
rz(-2.6014156) q[2];
sx q[2];
rz(-1.8312757) q[2];
sx q[2];
rz(-1.8249701) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2665796) q[1];
sx q[1];
rz(-1.8266051) q[1];
sx q[1];
rz(-1.2775879) q[1];
rz(-pi) q[2];
rz(-2.8323402) q[3];
sx q[3];
rz(-1.4746631) q[3];
sx q[3];
rz(-2.1949196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1332625) q[2];
sx q[2];
rz(-1.8066758) q[2];
sx q[2];
rz(1.5745715) q[2];
rz(-1.3957006) q[3];
sx q[3];
rz(-1.402366) q[3];
sx q[3];
rz(2.375026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.72346) q[0];
sx q[0];
rz(-0.27579951) q[0];
sx q[0];
rz(-0.3904528) q[0];
rz(-2.0402724) q[1];
sx q[1];
rz(-1.7951671) q[1];
sx q[1];
rz(0.93213814) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9207805) q[0];
sx q[0];
rz(-1.5340065) q[0];
sx q[0];
rz(-3.1139741) q[0];
rz(0.53932346) q[2];
sx q[2];
rz(-1.612609) q[2];
sx q[2];
rz(1.5780748) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.336672) q[1];
sx q[1];
rz(-1.763736) q[1];
sx q[1];
rz(-2.8590917) q[1];
rz(1.325439) q[3];
sx q[3];
rz(-2.0319157) q[3];
sx q[3];
rz(2.2803015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.091888) q[2];
sx q[2];
rz(-0.69685093) q[2];
sx q[2];
rz(1.6799124) q[2];
rz(-1.0852496) q[3];
sx q[3];
rz(-2.0094252) q[3];
sx q[3];
rz(-2.5011044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.853448) q[0];
sx q[0];
rz(-0.65971056) q[0];
sx q[0];
rz(-2.0154542) q[0];
rz(2.0938734) q[1];
sx q[1];
rz(-0.62647096) q[1];
sx q[1];
rz(-1.610021) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1140014) q[0];
sx q[0];
rz(-0.023311255) q[0];
sx q[0];
rz(-3.0816881) q[0];
rz(-2.6833186) q[2];
sx q[2];
rz(-1.7802139) q[2];
sx q[2];
rz(-0.50814123) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.46343741) q[1];
sx q[1];
rz(-3.040294) q[1];
sx q[1];
rz(-0.8763635) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5793581) q[3];
sx q[3];
rz(-1.4959129) q[3];
sx q[3];
rz(0.39505388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5213726) q[2];
sx q[2];
rz(-2.3282101) q[2];
sx q[2];
rz(-0.20624557) q[2];
rz(0.5591875) q[3];
sx q[3];
rz(-1.5234448) q[3];
sx q[3];
rz(2.0172113) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6013721) q[0];
sx q[0];
rz(-1.843353) q[0];
sx q[0];
rz(-2.699615) q[0];
rz(1.1454918) q[1];
sx q[1];
rz(-1.3064697) q[1];
sx q[1];
rz(-1.3189955) q[1];
rz(-1.8505605) q[2];
sx q[2];
rz(-1.8342809) q[2];
sx q[2];
rz(-0.4763436) q[2];
rz(-0.72412422) q[3];
sx q[3];
rz(-1.4403314) q[3];
sx q[3];
rz(-2.452313) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
