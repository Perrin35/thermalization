OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.00081113022) q[0];
sx q[0];
rz(-0.82511628) q[0];
sx q[0];
rz(0.42186475) q[0];
rz(-0.74692625) q[1];
sx q[1];
rz(3.7318228) q[1];
sx q[1];
rz(11.745315) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3481188) q[0];
sx q[0];
rz(-0.095190053) q[0];
sx q[0];
rz(-0.4918672) q[0];
rz(2.9157588) q[2];
sx q[2];
rz(-1.3703114) q[2];
sx q[2];
rz(3.0124913) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6118879) q[1];
sx q[1];
rz(-2.6571353) q[1];
sx q[1];
rz(2.2147708) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9299292) q[3];
sx q[3];
rz(-2.7756423) q[3];
sx q[3];
rz(-2.8173878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.39516285) q[2];
sx q[2];
rz(-2.2300827) q[2];
sx q[2];
rz(-0.54473031) q[2];
rz(1.644545) q[3];
sx q[3];
rz(-1.112273) q[3];
sx q[3];
rz(1.8018319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4561975) q[0];
sx q[0];
rz(-2.39769) q[0];
sx q[0];
rz(0.5734545) q[0];
rz(-3.0990797) q[1];
sx q[1];
rz(-1.7492234) q[1];
sx q[1];
rz(-1.858985) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79270482) q[0];
sx q[0];
rz(-0.12785873) q[0];
sx q[0];
rz(-0.16221817) q[0];
rz(-pi) q[1];
rz(0.2942652) q[2];
sx q[2];
rz(-1.8163101) q[2];
sx q[2];
rz(0.15799274) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6635884) q[1];
sx q[1];
rz(-1.4299031) q[1];
sx q[1];
rz(0.53014596) q[1];
rz(-pi) q[2];
x q[2];
rz(0.76000472) q[3];
sx q[3];
rz(-2.3541321) q[3];
sx q[3];
rz(-1.852688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2635228) q[2];
sx q[2];
rz(-1.1493378) q[2];
sx q[2];
rz(-2.6980706) q[2];
rz(1.589132) q[3];
sx q[3];
rz(-2.7510721) q[3];
sx q[3];
rz(0.21903567) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.012101128) q[0];
sx q[0];
rz(-2.2624367) q[0];
sx q[0];
rz(2.7398859) q[0];
rz(-1.7237639) q[1];
sx q[1];
rz(-0.44770733) q[1];
sx q[1];
rz(-1.1161425) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7340461) q[0];
sx q[0];
rz(-0.73483682) q[0];
sx q[0];
rz(-0.46291344) q[0];
x q[1];
rz(2.2625917) q[2];
sx q[2];
rz(-1.689581) q[2];
sx q[2];
rz(2.9853161) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1765343) q[1];
sx q[1];
rz(-1.985834) q[1];
sx q[1];
rz(-2.1188291) q[1];
x q[2];
rz(0.12401993) q[3];
sx q[3];
rz(-2.5895666) q[3];
sx q[3];
rz(0.13904143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9125354) q[2];
sx q[2];
rz(-2.8450862) q[2];
sx q[2];
rz(1.0079481) q[2];
rz(-1.5063162) q[3];
sx q[3];
rz(-1.9624036) q[3];
sx q[3];
rz(2.262871) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2849543) q[0];
sx q[0];
rz(-0.027218787) q[0];
sx q[0];
rz(-0.83754367) q[0];
rz(-1.8800927) q[1];
sx q[1];
rz(-1.7439758) q[1];
sx q[1];
rz(-1.2417485) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.094114) q[0];
sx q[0];
rz(-0.54007184) q[0];
sx q[0];
rz(-2.5940226) q[0];
x q[1];
rz(1.2707082) q[2];
sx q[2];
rz(-2.5113912) q[2];
sx q[2];
rz(0.84033191) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1832402) q[1];
sx q[1];
rz(-1.1179525) q[1];
sx q[1];
rz(1.1151821) q[1];
rz(2.1513274) q[3];
sx q[3];
rz(-1.255569) q[3];
sx q[3];
rz(-2.6516556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9813098) q[2];
sx q[2];
rz(-2.1872988) q[2];
sx q[2];
rz(1.0120288) q[2];
rz(-2.1435598) q[3];
sx q[3];
rz(-2.4055552) q[3];
sx q[3];
rz(-1.6691104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0618458) q[0];
sx q[0];
rz(-2.340402) q[0];
sx q[0];
rz(-0.27859846) q[0];
rz(-2.8246236) q[1];
sx q[1];
rz(-1.7520889) q[1];
sx q[1];
rz(0.46151361) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7636895) q[0];
sx q[0];
rz(-0.17556854) q[0];
sx q[0];
rz(1.8493091) q[0];
rz(-pi) q[1];
rz(-0.95038173) q[2];
sx q[2];
rz(-1.912552) q[2];
sx q[2];
rz(-2.1555962) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7263067) q[1];
sx q[1];
rz(-1.4771645) q[1];
sx q[1];
rz(-1.3312396) q[1];
rz(-pi) q[2];
rz(1.4123437) q[3];
sx q[3];
rz(-2.5907443) q[3];
sx q[3];
rz(0.85516632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3658112) q[2];
sx q[2];
rz(-1.3373988) q[2];
sx q[2];
rz(2.865045) q[2];
rz(0.60025275) q[3];
sx q[3];
rz(-2.4252031) q[3];
sx q[3];
rz(-2.6071809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14715956) q[0];
sx q[0];
rz(-2.5614547) q[0];
sx q[0];
rz(-0.68122) q[0];
rz(0.45411202) q[1];
sx q[1];
rz(-1.157016) q[1];
sx q[1];
rz(-0.62201321) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2300782) q[0];
sx q[0];
rz(-2.921836) q[0];
sx q[0];
rz(2.6917372) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63779442) q[2];
sx q[2];
rz(-2.0908815) q[2];
sx q[2];
rz(2.3931723) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4677538) q[1];
sx q[1];
rz(-2.2430393) q[1];
sx q[1];
rz(-2.6121187) q[1];
x q[2];
rz(1.9196904) q[3];
sx q[3];
rz(-1.8233646) q[3];
sx q[3];
rz(3.0608197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5155718) q[2];
sx q[2];
rz(-2.5280648) q[2];
sx q[2];
rz(0.75508368) q[2];
rz(-0.10609047) q[3];
sx q[3];
rz(-1.875149) q[3];
sx q[3];
rz(2.6809926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.047091529) q[0];
sx q[0];
rz(-0.61151183) q[0];
sx q[0];
rz(0.066135429) q[0];
rz(-0.051008929) q[1];
sx q[1];
rz(-0.74791932) q[1];
sx q[1];
rz(1.8775108) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28461489) q[0];
sx q[0];
rz(-0.98941411) q[0];
sx q[0];
rz(0.435243) q[0];
x q[1];
rz(-2.4263229) q[2];
sx q[2];
rz(-1.8085305) q[2];
sx q[2];
rz(-1.0843074) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3466929) q[1];
sx q[1];
rz(-0.82603077) q[1];
sx q[1];
rz(-0.54462437) q[1];
rz(-pi) q[2];
rz(2.7080688) q[3];
sx q[3];
rz(-1.2085087) q[3];
sx q[3];
rz(1.9725245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2833726) q[2];
sx q[2];
rz(-1.1274575) q[2];
sx q[2];
rz(-2.9952725) q[2];
rz(0.60394168) q[3];
sx q[3];
rz(-2.5950409) q[3];
sx q[3];
rz(1.4124136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42501763) q[0];
sx q[0];
rz(-1.9442433) q[0];
sx q[0];
rz(-0.091751598) q[0];
rz(-3.0692406) q[1];
sx q[1];
rz(-1.3183343) q[1];
sx q[1];
rz(-0.66973698) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3689174) q[0];
sx q[0];
rz(-1.2784504) q[0];
sx q[0];
rz(1.4293074) q[0];
rz(-1.7013427) q[2];
sx q[2];
rz(-1.1645082) q[2];
sx q[2];
rz(-2.8216336) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3460541) q[1];
sx q[1];
rz(-2.8020952) q[1];
sx q[1];
rz(0.44793753) q[1];
rz(-1.2708951) q[3];
sx q[3];
rz(-2.6811973) q[3];
sx q[3];
rz(3.0391191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3503795) q[2];
sx q[2];
rz(-0.71724856) q[2];
sx q[2];
rz(-0.72163248) q[2];
rz(2.3676938) q[3];
sx q[3];
rz(-2.1589203) q[3];
sx q[3];
rz(-2.6017046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72188193) q[0];
sx q[0];
rz(-1.9419436) q[0];
sx q[0];
rz(-0.53959674) q[0];
rz(-2.3609912) q[1];
sx q[1];
rz(-1.8949948) q[1];
sx q[1];
rz(0.4221198) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84655428) q[0];
sx q[0];
rz(-2.1384412) q[0];
sx q[0];
rz(0.66935434) q[0];
rz(-pi) q[1];
rz(-1.5109748) q[2];
sx q[2];
rz(-1.079353) q[2];
sx q[2];
rz(-2.7697542) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6680702) q[1];
sx q[1];
rz(-2.0501185) q[1];
sx q[1];
rz(-0.042906656) q[1];
rz(-pi) q[2];
rz(1.6160746) q[3];
sx q[3];
rz(-2.1589734) q[3];
sx q[3];
rz(0.28376337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0427986) q[2];
sx q[2];
rz(-1.021794) q[2];
sx q[2];
rz(-0.54879028) q[2];
rz(0.15268606) q[3];
sx q[3];
rz(-2.1137674) q[3];
sx q[3];
rz(-2.4299183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3720836) q[0];
sx q[0];
rz(-0.39396572) q[0];
sx q[0];
rz(-2.5373996) q[0];
rz(-1.4671885) q[1];
sx q[1];
rz(-1.3507651) q[1];
sx q[1];
rz(-1.1473354) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7580099) q[0];
sx q[0];
rz(-2.1910163) q[0];
sx q[0];
rz(2.0877286) q[0];
rz(-pi) q[1];
rz(-0.34530039) q[2];
sx q[2];
rz(-2.3757114) q[2];
sx q[2];
rz(1.6102546) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5903476) q[1];
sx q[1];
rz(-0.4830803) q[1];
sx q[1];
rz(0.03309588) q[1];
rz(-pi) q[2];
rz(2.2828034) q[3];
sx q[3];
rz(-2.6325573) q[3];
sx q[3];
rz(1.3648175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4597822) q[2];
sx q[2];
rz(-1.9777538) q[2];
sx q[2];
rz(2.5414844) q[2];
rz(-2.4188304) q[3];
sx q[3];
rz(-0.99812752) q[3];
sx q[3];
rz(0.37989894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73846524) q[0];
sx q[0];
rz(-1.6275788) q[0];
sx q[0];
rz(1.7049261) q[0];
rz(-1.5858831) q[1];
sx q[1];
rz(-1.7441505) q[1];
sx q[1];
rz(2.2081262) q[1];
rz(-2.2384833) q[2];
sx q[2];
rz(-1.9676859) q[2];
sx q[2];
rz(-3.0953593) q[2];
rz(-1.7865576) q[3];
sx q[3];
rz(-1.8117306) q[3];
sx q[3];
rz(2.4227052) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
