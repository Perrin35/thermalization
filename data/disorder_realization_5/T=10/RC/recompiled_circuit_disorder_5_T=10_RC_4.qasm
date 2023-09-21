OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5390227) q[0];
sx q[0];
rz(-2.5780926) q[0];
sx q[0];
rz(-0.45698693) q[0];
rz(2.5198088) q[1];
sx q[1];
rz(-2.4609202) q[1];
sx q[1];
rz(-1.2759804) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7514483) q[0];
sx q[0];
rz(-1.7062618) q[0];
sx q[0];
rz(-0.22060237) q[0];
rz(3.0511191) q[2];
sx q[2];
rz(-1.4960559) q[2];
sx q[2];
rz(-0.59404101) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4133271) q[1];
sx q[1];
rz(-1.773355) q[1];
sx q[1];
rz(-3.0800746) q[1];
rz(2.7492417) q[3];
sx q[3];
rz(-0.96071834) q[3];
sx q[3];
rz(0.073079212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0322545) q[2];
sx q[2];
rz(-0.87301746) q[2];
sx q[2];
rz(-2.9764552) q[2];
rz(-1.6312381) q[3];
sx q[3];
rz(-0.32838467) q[3];
sx q[3];
rz(2.3993649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6363268) q[0];
sx q[0];
rz(-0.16997448) q[0];
sx q[0];
rz(0.045036137) q[0];
rz(-0.33915195) q[1];
sx q[1];
rz(-1.76666) q[1];
sx q[1];
rz(0.80274686) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.734056) q[0];
sx q[0];
rz(-1.3450087) q[0];
sx q[0];
rz(1.9682103) q[0];
x q[1];
rz(-2.0741834) q[2];
sx q[2];
rz(-1.479584) q[2];
sx q[2];
rz(2.3561321) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.75623679) q[1];
sx q[1];
rz(-1.6223575) q[1];
sx q[1];
rz(-2.1897584) q[1];
rz(1.734415) q[3];
sx q[3];
rz(-0.67577261) q[3];
sx q[3];
rz(-2.8676652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7654483) q[2];
sx q[2];
rz(-0.62952289) q[2];
sx q[2];
rz(1.7155898) q[2];
rz(-2.5358893) q[3];
sx q[3];
rz(-0.96174812) q[3];
sx q[3];
rz(3.0595996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1894492) q[0];
sx q[0];
rz(-1.209963) q[0];
sx q[0];
rz(-2.8564575) q[0];
rz(-0.51672283) q[1];
sx q[1];
rz(-1.6500094) q[1];
sx q[1];
rz(1.3396938) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8782608) q[0];
sx q[0];
rz(-2.2790331) q[0];
sx q[0];
rz(-1.0017298) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2013024) q[2];
sx q[2];
rz(-0.72417799) q[2];
sx q[2];
rz(-2.5232814) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7741751) q[1];
sx q[1];
rz(-1.7304825) q[1];
sx q[1];
rz(1.7596485) q[1];
rz(-2.4524868) q[3];
sx q[3];
rz(-1.8854685) q[3];
sx q[3];
rz(-0.93022197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4703935) q[2];
sx q[2];
rz(-0.931804) q[2];
sx q[2];
rz(-1.7335256) q[2];
rz(-0.18313289) q[3];
sx q[3];
rz(-2.1268842) q[3];
sx q[3];
rz(0.071921913) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6614439) q[0];
sx q[0];
rz(-2.792795) q[0];
sx q[0];
rz(-0.21425042) q[0];
rz(0.42901531) q[1];
sx q[1];
rz(-2.0206101) q[1];
sx q[1];
rz(2.4750211) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4144856) q[0];
sx q[0];
rz(-1.7761199) q[0];
sx q[0];
rz(1.3798825) q[0];
x q[1];
rz(1.8604943) q[2];
sx q[2];
rz(-2.2640267) q[2];
sx q[2];
rz(1.9300269) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9208593) q[1];
sx q[1];
rz(-1.024106) q[1];
sx q[1];
rz(2.2842555) q[1];
rz(-pi) q[2];
rz(-0.781226) q[3];
sx q[3];
rz(-2.9328212) q[3];
sx q[3];
rz(-0.77199948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.2864705) q[2];
sx q[2];
rz(-2.9292967) q[2];
sx q[2];
rz(-0.66037035) q[2];
rz(-1.1874229) q[3];
sx q[3];
rz(-1.238845) q[3];
sx q[3];
rz(-0.56046265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51263556) q[0];
sx q[0];
rz(-1.9218788) q[0];
sx q[0];
rz(0.83706013) q[0];
rz(2.191026) q[1];
sx q[1];
rz(-0.66851139) q[1];
sx q[1];
rz(1.9821092) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54407185) q[0];
sx q[0];
rz(-2.263875) q[0];
sx q[0];
rz(-1.9835299) q[0];
rz(-pi) q[1];
rz(-1.7972838) q[2];
sx q[2];
rz(-0.61033568) q[2];
sx q[2];
rz(-3.0938873) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0817464) q[1];
sx q[1];
rz(-1.0330327) q[1];
sx q[1];
rz(0.50077011) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5043143) q[3];
sx q[3];
rz(-1.1253469) q[3];
sx q[3];
rz(-0.74419903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.435047) q[2];
sx q[2];
rz(-2.0303625) q[2];
sx q[2];
rz(-1.2716028) q[2];
rz(1.050625) q[3];
sx q[3];
rz(-1.8836861) q[3];
sx q[3];
rz(-0.064855382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5669252) q[0];
sx q[0];
rz(-2.8820679) q[0];
sx q[0];
rz(-2.0264453) q[0];
rz(-0.043958157) q[1];
sx q[1];
rz(-1.3479439) q[1];
sx q[1];
rz(0.10087068) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2642919) q[0];
sx q[0];
rz(-1.4192686) q[0];
sx q[0];
rz(0.63440462) q[0];
rz(1.2783373) q[2];
sx q[2];
rz(-1.1012226) q[2];
sx q[2];
rz(-1.8243607) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.74360352) q[1];
sx q[1];
rz(-1.3685703) q[1];
sx q[1];
rz(-2.9823751) q[1];
rz(1.8039861) q[3];
sx q[3];
rz(-1.9662074) q[3];
sx q[3];
rz(0.49664859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4042525) q[2];
sx q[2];
rz(-1.6254144) q[2];
sx q[2];
rz(2.2163088) q[2];
rz(0.18151367) q[3];
sx q[3];
rz(-2.4003024) q[3];
sx q[3];
rz(1.8484176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7589384) q[0];
sx q[0];
rz(-2.721334) q[0];
sx q[0];
rz(0.96310258) q[0];
rz(-1.5513647) q[1];
sx q[1];
rz(-1.0179049) q[1];
sx q[1];
rz(0.68774736) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4884268) q[0];
sx q[0];
rz(-1.352172) q[0];
sx q[0];
rz(1.9125008) q[0];
rz(-2.9155832) q[2];
sx q[2];
rz(-1.2603716) q[2];
sx q[2];
rz(-1.4357391) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3878165) q[1];
sx q[1];
rz(-3.0101335) q[1];
sx q[1];
rz(0.17863518) q[1];
rz(-pi) q[2];
rz(0.01597605) q[3];
sx q[3];
rz(-2.139058) q[3];
sx q[3];
rz(2.361134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8035651) q[2];
sx q[2];
rz(-1.7269208) q[2];
sx q[2];
rz(2.2040099) q[2];
rz(-2.1607416) q[3];
sx q[3];
rz(-2.8068481) q[3];
sx q[3];
rz(-0.78398314) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6136318) q[0];
sx q[0];
rz(-2.1870446) q[0];
sx q[0];
rz(3.0623073) q[0];
rz(-2.0203363) q[1];
sx q[1];
rz(-2.2032578) q[1];
sx q[1];
rz(-1.942873) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6262561) q[0];
sx q[0];
rz(-1.5021828) q[0];
sx q[0];
rz(-2.7068044) q[0];
x q[1];
rz(1.2077246) q[2];
sx q[2];
rz(-1.6512401) q[2];
sx q[2];
rz(2.9094839) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3938013) q[1];
sx q[1];
rz(-2.0732905) q[1];
sx q[1];
rz(-1.5138813) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6783887) q[3];
sx q[3];
rz(-1.3363046) q[3];
sx q[3];
rz(-0.47298688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2856059) q[2];
sx q[2];
rz(-0.14582835) q[2];
sx q[2];
rz(0.78906995) q[2];
rz(0.35813913) q[3];
sx q[3];
rz(-1.5877682) q[3];
sx q[3];
rz(-1.8607128) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57551861) q[0];
sx q[0];
rz(-1.7552567) q[0];
sx q[0];
rz(-0.58832204) q[0];
rz(0.026780216) q[1];
sx q[1];
rz(-2.2356922) q[1];
sx q[1];
rz(-0.9334329) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3073472) q[0];
sx q[0];
rz(-1.5114307) q[0];
sx q[0];
rz(3.0349915) q[0];
rz(-pi) q[1];
rz(-0.8546631) q[2];
sx q[2];
rz(-0.48644201) q[2];
sx q[2];
rz(-1.9180627) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.48298207) q[1];
sx q[1];
rz(-2.254199) q[1];
sx q[1];
rz(1.618209) q[1];
x q[2];
rz(-0.45456072) q[3];
sx q[3];
rz(-0.70485669) q[3];
sx q[3];
rz(-2.9398033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0788706) q[2];
sx q[2];
rz(-0.82010078) q[2];
sx q[2];
rz(2.771647) q[2];
rz(-2.2423819) q[3];
sx q[3];
rz(-1.4827385) q[3];
sx q[3];
rz(-2.1883011) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9545492) q[0];
sx q[0];
rz(-1.8148913) q[0];
sx q[0];
rz(-2.488234) q[0];
rz(1.6409142) q[1];
sx q[1];
rz(-1.4705642) q[1];
sx q[1];
rz(0.18383372) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.375658) q[0];
sx q[0];
rz(-1.4544393) q[0];
sx q[0];
rz(-3.0202306) q[0];
rz(0.83689883) q[2];
sx q[2];
rz(-0.95015929) q[2];
sx q[2];
rz(-1.9377973) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2059584) q[1];
sx q[1];
rz(-1.4050254) q[1];
sx q[1];
rz(0.88840719) q[1];
rz(-1.5552181) q[3];
sx q[3];
rz(-0.9874978) q[3];
sx q[3];
rz(0.52535666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.9188149) q[2];
sx q[2];
rz(-0.99514014) q[2];
sx q[2];
rz(2.9392021) q[2];
rz(1.0732132) q[3];
sx q[3];
rz(-1.9066633) q[3];
sx q[3];
rz(1.7741268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4595173) q[0];
sx q[0];
rz(-0.85953241) q[0];
sx q[0];
rz(-0.55534242) q[0];
rz(-0.36322414) q[1];
sx q[1];
rz(-1.8050615) q[1];
sx q[1];
rz(-0.25711679) q[1];
rz(2.6828962) q[2];
sx q[2];
rz(-1.9445322) q[2];
sx q[2];
rz(-1.4945488) q[2];
rz(-0.17584569) q[3];
sx q[3];
rz(-0.98777117) q[3];
sx q[3];
rz(-2.5143757) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];