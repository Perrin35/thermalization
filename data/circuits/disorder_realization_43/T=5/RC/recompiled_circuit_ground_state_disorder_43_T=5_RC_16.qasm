OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2744098) q[0];
sx q[0];
rz(2.6743968) q[0];
sx q[0];
rz(13.462486) q[0];
rz(2.3236302) q[1];
sx q[1];
rz(-1.5939413) q[1];
sx q[1];
rz(0.98734468) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2786699) q[0];
sx q[0];
rz(-0.76919829) q[0];
sx q[0];
rz(-2.132036) q[0];
x q[1];
rz(-0.57968037) q[2];
sx q[2];
rz(-2.2319921) q[2];
sx q[2];
rz(1.5722317) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3707177) q[1];
sx q[1];
rz(-2.8159119) q[1];
sx q[1];
rz(2.5925255) q[1];
rz(0.39405502) q[3];
sx q[3];
rz(-2.7480304) q[3];
sx q[3];
rz(1.8436028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2200615) q[2];
sx q[2];
rz(-1.0342197) q[2];
sx q[2];
rz(2.6424778) q[2];
rz(2.3257997) q[3];
sx q[3];
rz(-2.54839) q[3];
sx q[3];
rz(-0.35713404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65207425) q[0];
sx q[0];
rz(-2.7777785) q[0];
sx q[0];
rz(2.9079085) q[0];
rz(-3.0417327) q[1];
sx q[1];
rz(-1.4870817) q[1];
sx q[1];
rz(1.5335836) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8312992) q[0];
sx q[0];
rz(-1.6673894) q[0];
sx q[0];
rz(1.4486758) q[0];
x q[1];
rz(-1.3855782) q[2];
sx q[2];
rz(-1.1389009) q[2];
sx q[2];
rz(-0.76755953) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.74190226) q[1];
sx q[1];
rz(-2.5292225) q[1];
sx q[1];
rz(1.4111817) q[1];
rz(1.0879806) q[3];
sx q[3];
rz(-1.8466966) q[3];
sx q[3];
rz(1.538572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0443772) q[2];
sx q[2];
rz(-2.2558687) q[2];
sx q[2];
rz(0.054917939) q[2];
rz(-0.6461668) q[3];
sx q[3];
rz(-1.61444) q[3];
sx q[3];
rz(-3.0687148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9570479) q[0];
sx q[0];
rz(-0.2453198) q[0];
sx q[0];
rz(-2.3967337) q[0];
rz(-1.2767731) q[1];
sx q[1];
rz(-2.1121912) q[1];
sx q[1];
rz(-0.863711) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0833763) q[0];
sx q[0];
rz(-1.5124707) q[0];
sx q[0];
rz(2.7013679) q[0];
rz(-pi) q[1];
rz(2.6361739) q[2];
sx q[2];
rz(-0.8552455) q[2];
sx q[2];
rz(2.9723013) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.59434592) q[1];
sx q[1];
rz(-0.9556831) q[1];
sx q[1];
rz(-2.1662461) q[1];
rz(1.0705804) q[3];
sx q[3];
rz(-2.1693834) q[3];
sx q[3];
rz(-2.9974496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.471571) q[2];
sx q[2];
rz(-1.5315285) q[2];
sx q[2];
rz(-0.395533) q[2];
rz(2.1192571) q[3];
sx q[3];
rz(-2.6774355) q[3];
sx q[3];
rz(0.49183229) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1763497) q[0];
sx q[0];
rz(-1.9597766) q[0];
sx q[0];
rz(3.0856207) q[0];
rz(2.5296027) q[1];
sx q[1];
rz(-1.5490218) q[1];
sx q[1];
rz(-2.2956119) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2200945) q[0];
sx q[0];
rz(-1.3125696) q[0];
sx q[0];
rz(1.0357712) q[0];
rz(-pi) q[1];
rz(-0.34807713) q[2];
sx q[2];
rz(-2.161059) q[2];
sx q[2];
rz(-1.2952309) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3624907) q[1];
sx q[1];
rz(-1.7622708) q[1];
sx q[1];
rz(-2.8250241) q[1];
rz(-pi) q[2];
x q[2];
rz(0.37644551) q[3];
sx q[3];
rz(-2.0168224) q[3];
sx q[3];
rz(-2.6335187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1299639) q[2];
sx q[2];
rz(-1.446898) q[2];
sx q[2];
rz(2.4212627) q[2];
rz(0.35308853) q[3];
sx q[3];
rz(-1.1938286) q[3];
sx q[3];
rz(2.8928355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40097749) q[0];
sx q[0];
rz(-1.688513) q[0];
sx q[0];
rz(0.98689669) q[0];
rz(1.997021) q[1];
sx q[1];
rz(-2.3026376) q[1];
sx q[1];
rz(-2.1867337) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3026149) q[0];
sx q[0];
rz(-0.90388008) q[0];
sx q[0];
rz(-0.033292183) q[0];
x q[1];
rz(2.8138732) q[2];
sx q[2];
rz(-1.0340889) q[2];
sx q[2];
rz(-2.2959054) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8870626) q[1];
sx q[1];
rz(-1.945244) q[1];
sx q[1];
rz(2.6545054) q[1];
rz(-0.14820672) q[3];
sx q[3];
rz(-1.8216405) q[3];
sx q[3];
rz(-0.36484066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5130875) q[2];
sx q[2];
rz(-0.91431618) q[2];
sx q[2];
rz(2.0950192) q[2];
rz(2.4322677) q[3];
sx q[3];
rz(-1.9966634) q[3];
sx q[3];
rz(-2.5078702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3777622) q[0];
sx q[0];
rz(-1.4077633) q[0];
sx q[0];
rz(2.8705226) q[0];
rz(-0.079004869) q[1];
sx q[1];
rz(-2.7075691) q[1];
sx q[1];
rz(1.7105506) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1880497) q[0];
sx q[0];
rz(-2.8859647) q[0];
sx q[0];
rz(-2.5357657) q[0];
x q[1];
rz(-0.56366326) q[2];
sx q[2];
rz(-1.0978292) q[2];
sx q[2];
rz(1.2207292) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5137427) q[1];
sx q[1];
rz(-1.8166891) q[1];
sx q[1];
rz(-0.69468433) q[1];
rz(-pi) q[2];
rz(-1.2215516) q[3];
sx q[3];
rz(-0.442083) q[3];
sx q[3];
rz(-0.52607049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0926823) q[2];
sx q[2];
rz(-1.9581257) q[2];
sx q[2];
rz(3.1177706) q[2];
rz(-1.228099) q[3];
sx q[3];
rz(-2.7619599) q[3];
sx q[3];
rz(-2.9983799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.1517076) q[0];
sx q[0];
rz(-0.30549529) q[0];
sx q[0];
rz(0.09224961) q[0];
rz(2.8406738) q[1];
sx q[1];
rz(-1.9124799) q[1];
sx q[1];
rz(-1.0801962) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0940822) q[0];
sx q[0];
rz(-1.9157529) q[0];
sx q[0];
rz(-3.0408919) q[0];
rz(-pi) q[1];
rz(0.36209468) q[2];
sx q[2];
rz(-1.9685192) q[2];
sx q[2];
rz(2.5366572) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.791886) q[1];
sx q[1];
rz(-2.9936446) q[1];
sx q[1];
rz(-0.011758864) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8870513) q[3];
sx q[3];
rz(-2.2407534) q[3];
sx q[3];
rz(3.0681821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.13176189) q[2];
sx q[2];
rz(-2.7907382) q[2];
sx q[2];
rz(0.62823137) q[2];
rz(2.175711) q[3];
sx q[3];
rz(-1.5647669) q[3];
sx q[3];
rz(0.32111827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0308762) q[0];
sx q[0];
rz(-2.4346011) q[0];
sx q[0];
rz(1.734717) q[0];
rz(3.0137317) q[1];
sx q[1];
rz(-1.7117932) q[1];
sx q[1];
rz(-1.1258639) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4425621) q[0];
sx q[0];
rz(-1.3872212) q[0];
sx q[0];
rz(-0.72893127) q[0];
x q[1];
rz(-3.1229805) q[2];
sx q[2];
rz(-2.0224704) q[2];
sx q[2];
rz(-1.1455331) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.023214528) q[1];
sx q[1];
rz(-1.5074308) q[1];
sx q[1];
rz(-2.9411282) q[1];
rz(-0.61922686) q[3];
sx q[3];
rz(-1.9398749) q[3];
sx q[3];
rz(-0.31821966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.13059482) q[2];
sx q[2];
rz(-2.0086292) q[2];
sx q[2];
rz(0.093718378) q[2];
rz(-1.7081918) q[3];
sx q[3];
rz(-0.283537) q[3];
sx q[3];
rz(-1.7482429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.691064) q[0];
sx q[0];
rz(-2.0557623) q[0];
sx q[0];
rz(1.0040671) q[0];
rz(1.1035236) q[1];
sx q[1];
rz(-0.48201489) q[1];
sx q[1];
rz(1.7335266) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.133014) q[0];
sx q[0];
rz(-0.51515085) q[0];
sx q[0];
rz(-0.012284474) q[0];
rz(-2.8292921) q[2];
sx q[2];
rz(-1.4648629) q[2];
sx q[2];
rz(2.9779129) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4346659) q[1];
sx q[1];
rz(-1.3916755) q[1];
sx q[1];
rz(1.1302901) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.086661913) q[3];
sx q[3];
rz(-0.52059697) q[3];
sx q[3];
rz(-0.39256061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.39190009) q[2];
sx q[2];
rz(-1.5211279) q[2];
sx q[2];
rz(-0.54769546) q[2];
rz(0.31217602) q[3];
sx q[3];
rz(-1.7269937) q[3];
sx q[3];
rz(1.4148022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80426973) q[0];
sx q[0];
rz(-1.6707358) q[0];
sx q[0];
rz(0.7789337) q[0];
rz(-0.3745105) q[1];
sx q[1];
rz(-1.6284527) q[1];
sx q[1];
rz(0.83190727) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13369416) q[0];
sx q[0];
rz(-1.8323032) q[0];
sx q[0];
rz(0.71870872) q[0];
rz(-pi) q[1];
rz(0.7829297) q[2];
sx q[2];
rz(-2.1179869) q[2];
sx q[2];
rz(1.3732571) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.1268046) q[1];
sx q[1];
rz(-0.90378692) q[1];
sx q[1];
rz(-2.2875525) q[1];
rz(-pi) q[2];
rz(-3.1111701) q[3];
sx q[3];
rz(-0.796954) q[3];
sx q[3];
rz(1.4642293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1067918) q[2];
sx q[2];
rz(-0.35934862) q[2];
sx q[2];
rz(-0.0084776004) q[2];
rz(-2.7360385) q[3];
sx q[3];
rz(-1.6885992) q[3];
sx q[3];
rz(-1.4439553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99210284) q[0];
sx q[0];
rz(-1.6076417) q[0];
sx q[0];
rz(-2.3544307) q[0];
rz(-1.0943195) q[1];
sx q[1];
rz(-0.37847395) q[1];
sx q[1];
rz(2.7056221) q[1];
rz(-3.0158184) q[2];
sx q[2];
rz(-1.3604506) q[2];
sx q[2];
rz(2.4764555) q[2];
rz(-0.052362818) q[3];
sx q[3];
rz(-3.047245) q[3];
sx q[3];
rz(0.17518763) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
