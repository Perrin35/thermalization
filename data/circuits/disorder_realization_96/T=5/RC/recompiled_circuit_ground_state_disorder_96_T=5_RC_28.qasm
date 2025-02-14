OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.95028967) q[0];
sx q[0];
rz(6.0130881) q[0];
sx q[0];
rz(10.313378) q[0];
rz(1.8285881) q[1];
sx q[1];
rz(-1.5421901) q[1];
sx q[1];
rz(1.3808274) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61533538) q[0];
sx q[0];
rz(-1.3852296) q[0];
sx q[0];
rz(-2.944817) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6216952) q[2];
sx q[2];
rz(-0.87018229) q[2];
sx q[2];
rz(-2.0127279) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.12871582) q[1];
sx q[1];
rz(-0.33564645) q[1];
sx q[1];
rz(0.0291834) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2065776) q[3];
sx q[3];
rz(-1.7813588) q[3];
sx q[3];
rz(-0.59277804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8218653) q[2];
sx q[2];
rz(-1.3250019) q[2];
sx q[2];
rz(-0.74903178) q[2];
rz(0.15394112) q[3];
sx q[3];
rz(-1.0356244) q[3];
sx q[3];
rz(3.0416987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.446949) q[0];
sx q[0];
rz(-1.6634989) q[0];
sx q[0];
rz(-1.9248167) q[0];
rz(2.079839) q[1];
sx q[1];
rz(-0.80445015) q[1];
sx q[1];
rz(0.42713508) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9574643) q[0];
sx q[0];
rz(-1.2250568) q[0];
sx q[0];
rz(2.5347802) q[0];
rz(1.8893625) q[2];
sx q[2];
rz(-0.67020352) q[2];
sx q[2];
rz(0.88662749) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1698902) q[1];
sx q[1];
rz(-2.8198543) q[1];
sx q[1];
rz(1.1119026) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3668381) q[3];
sx q[3];
rz(-1.5195623) q[3];
sx q[3];
rz(2.8683087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.079166807) q[2];
sx q[2];
rz(-1.4126974) q[2];
sx q[2];
rz(-1.742935) q[2];
rz(0.22377293) q[3];
sx q[3];
rz(-0.90884915) q[3];
sx q[3];
rz(-1.5980501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-3.1200714) q[0];
sx q[0];
rz(-1.7912309) q[0];
sx q[0];
rz(-3.0960826) q[0];
rz(-1.9728164) q[1];
sx q[1];
rz(-1.2380506) q[1];
sx q[1];
rz(-2.5659836) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5110014) q[0];
sx q[0];
rz(-1.5819893) q[0];
sx q[0];
rz(2.5024947) q[0];
x q[1];
rz(-1.5928629) q[2];
sx q[2];
rz(-0.66078545) q[2];
sx q[2];
rz(1.2118076) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.61316669) q[1];
sx q[1];
rz(-1.3537145) q[1];
sx q[1];
rz(-1.0124799) q[1];
rz(-0.91797249) q[3];
sx q[3];
rz(-1.5906207) q[3];
sx q[3];
rz(1.6719847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1316954) q[2];
sx q[2];
rz(-2.3991149) q[2];
sx q[2];
rz(0.95727813) q[2];
rz(2.5668528) q[3];
sx q[3];
rz(-1.5301907) q[3];
sx q[3];
rz(1.3055698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9855373) q[0];
sx q[0];
rz(-0.95273459) q[0];
sx q[0];
rz(1.2611457) q[0];
rz(-0.082322923) q[1];
sx q[1];
rz(-1.0876834) q[1];
sx q[1];
rz(-1.3538768) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.026431) q[0];
sx q[0];
rz(-1.3360093) q[0];
sx q[0];
rz(-0.26682667) q[0];
rz(-pi) q[1];
x q[1];
rz(0.48575966) q[2];
sx q[2];
rz(-1.4076715) q[2];
sx q[2];
rz(-0.32147929) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2980068) q[1];
sx q[1];
rz(-1.6248967) q[1];
sx q[1];
rz(-2.6463616) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1534377) q[3];
sx q[3];
rz(-0.76863063) q[3];
sx q[3];
rz(0.60291327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.42133078) q[2];
sx q[2];
rz(-2.223184) q[2];
sx q[2];
rz(1.8850231) q[2];
rz(-0.71150696) q[3];
sx q[3];
rz(-1.3308176) q[3];
sx q[3];
rz(-0.28212696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8307777) q[0];
sx q[0];
rz(-0.84596914) q[0];
sx q[0];
rz(-0.34314439) q[0];
rz(-3.0798196) q[1];
sx q[1];
rz(-2.1715178) q[1];
sx q[1];
rz(-1.4168581) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2767216) q[0];
sx q[0];
rz(-1.7407932) q[0];
sx q[0];
rz(1.8525339) q[0];
x q[1];
rz(1.2447276) q[2];
sx q[2];
rz(-0.89674258) q[2];
sx q[2];
rz(1.774396) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0793905) q[1];
sx q[1];
rz(-0.50143999) q[1];
sx q[1];
rz(-1.2859243) q[1];
rz(-2.9067578) q[3];
sx q[3];
rz(-1.1560625) q[3];
sx q[3];
rz(2.107634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5938277) q[2];
sx q[2];
rz(-1.4676899) q[2];
sx q[2];
rz(1.0859547) q[2];
rz(-2.2222399) q[3];
sx q[3];
rz(-1.7707526) q[3];
sx q[3];
rz(-0.61387387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3980961) q[0];
sx q[0];
rz(-1.2092104) q[0];
sx q[0];
rz(-0.79291517) q[0];
rz(2.4010557) q[1];
sx q[1];
rz(-2.1426327) q[1];
sx q[1];
rz(0.83121306) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75126264) q[0];
sx q[0];
rz(-0.4340651) q[0];
sx q[0];
rz(2.084226) q[0];
rz(-pi) q[1];
rz(3.1341427) q[2];
sx q[2];
rz(-1.2937355) q[2];
sx q[2];
rz(-1.3598833) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2789687) q[1];
sx q[1];
rz(-1.2662856) q[1];
sx q[1];
rz(0.83359756) q[1];
rz(0.020757787) q[3];
sx q[3];
rz(-1.100173) q[3];
sx q[3];
rz(-0.94604674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.071216019) q[2];
sx q[2];
rz(-1.8698317) q[2];
sx q[2];
rz(-2.0224723) q[2];
rz(1.8708771) q[3];
sx q[3];
rz(-1.1071353) q[3];
sx q[3];
rz(-1.3814111) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1374461) q[0];
sx q[0];
rz(-0.16682145) q[0];
sx q[0];
rz(1.5420089) q[0];
rz(-2.1265325) q[1];
sx q[1];
rz(-1.5856182) q[1];
sx q[1];
rz(0.17280811) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8625582) q[0];
sx q[0];
rz(-2.1365215) q[0];
sx q[0];
rz(-2.5507257) q[0];
rz(-pi) q[1];
rz(-0.44266959) q[2];
sx q[2];
rz(-2.4076862) q[2];
sx q[2];
rz(1.9706322) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4932996) q[1];
sx q[1];
rz(-1.9956335) q[1];
sx q[1];
rz(0.27314911) q[1];
x q[2];
rz(2.1234346) q[3];
sx q[3];
rz(-0.88200906) q[3];
sx q[3];
rz(-1.8398726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.479594) q[2];
sx q[2];
rz(-0.84344232) q[2];
sx q[2];
rz(2.1638828) q[2];
rz(0.57502037) q[3];
sx q[3];
rz(-1.9366555) q[3];
sx q[3];
rz(-1.9112816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8136895) q[0];
sx q[0];
rz(-0.29569018) q[0];
sx q[0];
rz(3.1258702) q[0];
rz(-1.9533336) q[1];
sx q[1];
rz(-0.63242042) q[1];
sx q[1];
rz(-1.1994919) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6576783) q[0];
sx q[0];
rz(-1.8059314) q[0];
sx q[0];
rz(1.1170618) q[0];
x q[1];
rz(1.079915) q[2];
sx q[2];
rz(-1.5245617) q[2];
sx q[2];
rz(2.2904928) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3856628) q[1];
sx q[1];
rz(-2.0634868) q[1];
sx q[1];
rz(1.0316959) q[1];
rz(-2.8335081) q[3];
sx q[3];
rz(-1.020806) q[3];
sx q[3];
rz(0.96984449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1477995) q[2];
sx q[2];
rz(-1.2387929) q[2];
sx q[2];
rz(-1.3851059) q[2];
rz(-2.5177453) q[3];
sx q[3];
rz(-1.8892989) q[3];
sx q[3];
rz(-1.0998211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95595908) q[0];
sx q[0];
rz(-2.3210242) q[0];
sx q[0];
rz(2.4483335) q[0];
rz(2.8765826) q[1];
sx q[1];
rz(-0.82844228) q[1];
sx q[1];
rz(-1.6185435) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5361621) q[0];
sx q[0];
rz(-0.79728617) q[0];
sx q[0];
rz(1.5553586) q[0];
x q[1];
rz(-2.1586645) q[2];
sx q[2];
rz(-1.0318021) q[2];
sx q[2];
rz(1.3608152) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5793188) q[1];
sx q[1];
rz(-0.34046158) q[1];
sx q[1];
rz(-0.071694362) q[1];
x q[2];
rz(-0.070793666) q[3];
sx q[3];
rz(-2.1507743) q[3];
sx q[3];
rz(-1.1759315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.30274621) q[2];
sx q[2];
rz(-1.6281717) q[2];
sx q[2];
rz(1.4576853) q[2];
rz(-2.1863106) q[3];
sx q[3];
rz(-0.49946076) q[3];
sx q[3];
rz(-2.3063851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7569698) q[0];
sx q[0];
rz(-2.5769825) q[0];
sx q[0];
rz(-0.48267522) q[0];
rz(-2.2118498) q[1];
sx q[1];
rz(-2.1371806) q[1];
sx q[1];
rz(2.7453056) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9549462) q[0];
sx q[0];
rz(-0.56342331) q[0];
sx q[0];
rz(-0.71840973) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.70102923) q[2];
sx q[2];
rz(-2.9091479) q[2];
sx q[2];
rz(1.503141) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0490108) q[1];
sx q[1];
rz(-1.3900458) q[1];
sx q[1];
rz(0.84047079) q[1];
rz(-pi) q[2];
rz(2.4195005) q[3];
sx q[3];
rz(-0.51755899) q[3];
sx q[3];
rz(-2.232377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3489909) q[2];
sx q[2];
rz(-1.7020117) q[2];
sx q[2];
rz(-2.8144042) q[2];
rz(0.78091019) q[3];
sx q[3];
rz(-1.9802997) q[3];
sx q[3];
rz(1.4769295) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0384211) q[0];
sx q[0];
rz(-2.1506943) q[0];
sx q[0];
rz(-2.8839169) q[0];
rz(2.7720263) q[1];
sx q[1];
rz(-0.8820487) q[1];
sx q[1];
rz(-0.65912156) q[1];
rz(1.6160028) q[2];
sx q[2];
rz(-0.87799413) q[2];
sx q[2];
rz(-2.9064657) q[2];
rz(-2.7575708) q[3];
sx q[3];
rz(-0.99746084) q[3];
sx q[3];
rz(1.8538047) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
