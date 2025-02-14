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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40872657) q[0];
sx q[0];
rz(-0.48983296) q[0];
sx q[0];
rz(0.31087713) q[0];
rz(-pi) q[1];
rz(-0.52337661) q[2];
sx q[2];
rz(-1.1672772) q[2];
sx q[2];
rz(1.7159107) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.11553227) q[1];
sx q[1];
rz(-2.3838413) q[1];
sx q[1];
rz(2.4819786) q[1];
rz(-pi) q[2];
rz(3.1389075) q[3];
sx q[3];
rz(-1.1863106) q[3];
sx q[3];
rz(2.9677396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.2653653) q[2];
sx q[2];
rz(-0.23398016) q[2];
sx q[2];
rz(0.9210251) q[2];
rz(1.6169351) q[3];
sx q[3];
rz(-1.8504668) q[3];
sx q[3];
rz(2.9558712) q[3];
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
rz(2.8601473) q[0];
sx q[0];
rz(-2.7466819) q[0];
sx q[0];
rz(-1.349378) q[0];
rz(0.34740627) q[1];
sx q[1];
rz(-1.9527304) q[1];
sx q[1];
rz(1.4030392) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0058075) q[0];
sx q[0];
rz(-0.62816915) q[0];
sx q[0];
rz(0.37395333) q[0];
rz(-2.5591056) q[2];
sx q[2];
rz(-0.91160027) q[2];
sx q[2];
rz(-0.90345736) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.053005347) q[1];
sx q[1];
rz(-2.2777875) q[1];
sx q[1];
rz(-1.7729575) q[1];
x q[2];
rz(-2.8722309) q[3];
sx q[3];
rz(-0.90051631) q[3];
sx q[3];
rz(-1.6034077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3226402) q[2];
sx q[2];
rz(-1.3078657) q[2];
sx q[2];
rz(-0.15882203) q[2];
rz(-2.9239376) q[3];
sx q[3];
rz(-1.7526151) q[3];
sx q[3];
rz(1.7910819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7749216) q[0];
sx q[0];
rz(-1.8407624) q[0];
sx q[0];
rz(1.2694673) q[0];
rz(-0.41995755) q[1];
sx q[1];
rz(-2.1407514) q[1];
sx q[1];
rz(-1.6023844) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7363971) q[0];
sx q[0];
rz(-2.0123082) q[0];
sx q[0];
rz(-0.34607743) q[0];
rz(-pi) q[1];
rz(-2.4478618) q[2];
sx q[2];
rz(-2.8588534) q[2];
sx q[2];
rz(3.0764584) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0809787) q[1];
sx q[1];
rz(-1.1757869) q[1];
sx q[1];
rz(0.057821349) q[1];
rz(-pi) q[2];
rz(0.56681239) q[3];
sx q[3];
rz(-1.2897964) q[3];
sx q[3];
rz(2.5917201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0594844) q[2];
sx q[2];
rz(-2.2993645) q[2];
sx q[2];
rz(0.19169894) q[2];
rz(2.5890403) q[3];
sx q[3];
rz(-1.6165761) q[3];
sx q[3];
rz(-1.0244757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4100274) q[0];
sx q[0];
rz(-2.4931694) q[0];
sx q[0];
rz(2.5372274) q[0];
rz(1.2454237) q[1];
sx q[1];
rz(-2.3912997) q[1];
sx q[1];
rz(-2.2818458) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1103134) q[0];
sx q[0];
rz(-0.95891751) q[0];
sx q[0];
rz(-2.7252134) q[0];
x q[1];
rz(0.014195125) q[2];
sx q[2];
rz(-1.6728587) q[2];
sx q[2];
rz(-1.0400187) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.27450505) q[1];
sx q[1];
rz(-0.18840677) q[1];
sx q[1];
rz(0.26885689) q[1];
x q[2];
rz(-2.4212461) q[3];
sx q[3];
rz(-1.7631774) q[3];
sx q[3];
rz(-1.9494353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2802281) q[2];
sx q[2];
rz(-2.1275529) q[2];
sx q[2];
rz(-1.7029943) q[2];
rz(1.8493308) q[3];
sx q[3];
rz(-0.82787138) q[3];
sx q[3];
rz(1.0546225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0786521) q[0];
sx q[0];
rz(-2.8684454) q[0];
sx q[0];
rz(-0.31785059) q[0];
rz(1.3392797) q[1];
sx q[1];
rz(-0.41792089) q[1];
sx q[1];
rz(-2.1628765) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.933799) q[0];
sx q[0];
rz(-1.1699256) q[0];
sx q[0];
rz(1.6036556) q[0];
rz(0.42338223) q[2];
sx q[2];
rz(-0.21373323) q[2];
sx q[2];
rz(0.20704421) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3826344) q[1];
sx q[1];
rz(-2.3178604) q[1];
sx q[1];
rz(-1.7428223) q[1];
rz(0.19342761) q[3];
sx q[3];
rz(-1.3816091) q[3];
sx q[3];
rz(2.6003169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.19834441) q[2];
sx q[2];
rz(-0.94503108) q[2];
sx q[2];
rz(-2.0242958) q[2];
rz(2.9605588) q[3];
sx q[3];
rz(-1.9046611) q[3];
sx q[3];
rz(1.9443289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.70500526) q[0];
sx q[0];
rz(-0.26708189) q[0];
sx q[0];
rz(1.8319112) q[0];
rz(1.0218703) q[1];
sx q[1];
rz(-2.3049057) q[1];
sx q[1];
rz(-1.6530564) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7393417) q[0];
sx q[0];
rz(-1.9034804) q[0];
sx q[0];
rz(-0.87693033) q[0];
rz(-pi) q[1];
rz(1.8044293) q[2];
sx q[2];
rz(-2.2356845) q[2];
sx q[2];
rz(2.6478752) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.19410876) q[1];
sx q[1];
rz(-2.6257537) q[1];
sx q[1];
rz(2.8495925) q[1];
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
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.66270193) q[2];
sx q[2];
rz(-1.7544489) q[2];
sx q[2];
rz(-2.877511) q[2];
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
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5298117) q[0];
sx q[0];
rz(-1.5008858) q[0];
sx q[0];
rz(1.06426) q[0];
rz(-3.0423959) q[1];
sx q[1];
rz(-1.7925037) q[1];
sx q[1];
rz(-1.4605716) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.195358) q[0];
sx q[0];
rz(-1.6345662) q[0];
sx q[0];
rz(-2.6568165) q[0];
x q[1];
rz(1.6978092) q[2];
sx q[2];
rz(-0.30108967) q[2];
sx q[2];
rz(1.3236698) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3702195) q[1];
sx q[1];
rz(-1.2882731) q[1];
sx q[1];
rz(-2.706896) q[1];
rz(-pi) q[2];
rz(-2.0123294) q[3];
sx q[3];
rz(-2.0534424) q[3];
sx q[3];
rz(-3.1114374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3508241) q[2];
sx q[2];
rz(-1.1188353) q[2];
sx q[2];
rz(-2.5420945) q[2];
rz(2.2439469) q[3];
sx q[3];
rz(-1.41956) q[3];
sx q[3];
rz(-3.1035778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3277603) q[0];
sx q[0];
rz(-2.0812415) q[0];
sx q[0];
rz(-1.4014442) q[0];
rz(-3.0701045) q[1];
sx q[1];
rz(-1.46336) q[1];
sx q[1];
rz(2.1404526) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4021804) q[0];
sx q[0];
rz(-0.52153679) q[0];
sx q[0];
rz(2.6163382) q[0];
rz(-pi) q[1];
rz(0.48076081) q[2];
sx q[2];
rz(-0.72274739) q[2];
sx q[2];
rz(-1.2302081) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.65232507) q[1];
sx q[1];
rz(-1.4389924) q[1];
sx q[1];
rz(-2.5062923) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1628405) q[3];
sx q[3];
rz(-2.1518937) q[3];
sx q[3];
rz(1.6421902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9643758) q[2];
sx q[2];
rz(-1.588725) q[2];
sx q[2];
rz(-0.71844086) q[2];
rz(-0.046772379) q[3];
sx q[3];
rz(-2.434157) q[3];
sx q[3];
rz(2.5717946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8319594) q[0];
sx q[0];
rz(-2.9871812) q[0];
sx q[0];
rz(2.4511448) q[0];
rz(-0.18028232) q[1];
sx q[1];
rz(-2.0803662) q[1];
sx q[1];
rz(0.32275018) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6868387) q[0];
sx q[0];
rz(-1.4658064) q[0];
sx q[0];
rz(-2.9872894) q[0];
rz(-pi) q[1];
rz(1.5602925) q[2];
sx q[2];
rz(-0.5521419) q[2];
sx q[2];
rz(-0.10607468) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9195936) q[1];
sx q[1];
rz(-2.1173491) q[1];
sx q[1];
rz(2.2080457) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.75876816) q[3];
sx q[3];
rz(-2.3009217) q[3];
sx q[3];
rz(-1.981995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.47101578) q[2];
sx q[2];
rz(-1.7936423) q[2];
sx q[2];
rz(-1.0178817) q[2];
rz(-0.0035303591) q[3];
sx q[3];
rz(-0.96960932) q[3];
sx q[3];
rz(-1.9297011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.8135391) q[0];
sx q[0];
rz(-2.4040451) q[0];
sx q[0];
rz(-0.92593431) q[0];
rz(-0.13735859) q[1];
sx q[1];
rz(-2.4691212) q[1];
sx q[1];
rz(0.59991178) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5900676) q[0];
sx q[0];
rz(-2.5976564) q[0];
sx q[0];
rz(3.1363669) q[0];
rz(-2.5125347) q[2];
sx q[2];
rz(-1.1060083) q[2];
sx q[2];
rz(-1.0126707) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.39420262) q[1];
sx q[1];
rz(-0.58947403) q[1];
sx q[1];
rz(-1.516045) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8163054) q[3];
sx q[3];
rz(-1.1184177) q[3];
sx q[3];
rz(0.77879209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7291193) q[2];
sx q[2];
rz(-1.0189265) q[2];
sx q[2];
rz(2.4998383) q[2];
rz(1.9852091) q[3];
sx q[3];
rz(-0.76120794) q[3];
sx q[3];
rz(1.6183287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.079035096) q[0];
sx q[0];
rz(-2.2618444) q[0];
sx q[0];
rz(-2.3656144) q[0];
rz(-1.7017801) q[1];
sx q[1];
rz(-1.4714614) q[1];
sx q[1];
rz(0.068838483) q[1];
rz(-1.769968) q[2];
sx q[2];
rz(-1.3190075) q[2];
sx q[2];
rz(-1.5803303) q[2];
rz(-1.0216807) q[3];
sx q[3];
rz(-1.1374478) q[3];
sx q[3];
rz(2.9986387) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
