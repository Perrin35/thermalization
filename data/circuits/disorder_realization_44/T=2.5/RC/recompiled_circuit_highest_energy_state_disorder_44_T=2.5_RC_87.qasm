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
rz(-0.62701464) q[0];
sx q[0];
rz(3.7778683) q[0];
sx q[0];
rz(10.502622) q[0];
rz(4.2181349) q[1];
sx q[1];
rz(3.7706673) q[1];
sx q[1];
rz(7.315048) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6451749) q[0];
sx q[0];
rz(-1.5749567) q[0];
sx q[0];
rz(1.5711938) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5450942) q[2];
sx q[2];
rz(-2.7257127) q[2];
sx q[2];
rz(0.88230995) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.87776041) q[1];
sx q[1];
rz(-1.9139557) q[1];
sx q[1];
rz(-0.43107098) q[1];
rz(-1.4937819) q[3];
sx q[3];
rz(-1.2359914) q[3];
sx q[3];
rz(2.176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2304113) q[2];
sx q[2];
rz(-1.1831256) q[2];
sx q[2];
rz(-0.46768701) q[2];
rz(0.45105252) q[3];
sx q[3];
rz(-2.7428198) q[3];
sx q[3];
rz(-2.8635645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(0.69538799) q[0];
sx q[0];
rz(-0.22286335) q[0];
sx q[0];
rz(-2.9872802) q[0];
rz(0.34126869) q[1];
sx q[1];
rz(-0.35366615) q[1];
sx q[1];
rz(-0.42010677) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.882928) q[0];
sx q[0];
rz(-2.1999584) q[0];
sx q[0];
rz(-1.8101519) q[0];
rz(-pi) q[1];
x q[1];
rz(0.86067274) q[2];
sx q[2];
rz(-1.5352852) q[2];
sx q[2];
rz(-1.6954317) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1359033) q[1];
sx q[1];
rz(-1.5721675) q[1];
sx q[1];
rz(0.47079177) q[1];
rz(-pi) q[2];
rz(2.9734334) q[3];
sx q[3];
rz(-2.3871867) q[3];
sx q[3];
rz(-1.2446294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.61098617) q[2];
sx q[2];
rz(-2.2425118) q[2];
sx q[2];
rz(0.77303028) q[2];
rz(2.2981339) q[3];
sx q[3];
rz(-1.5777028) q[3];
sx q[3];
rz(0.57190603) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.968349) q[0];
sx q[0];
rz(-1.0272212) q[0];
sx q[0];
rz(-2.9395043) q[0];
rz(-2.659722) q[1];
sx q[1];
rz(-0.20129573) q[1];
sx q[1];
rz(-0.906382) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4545091) q[0];
sx q[0];
rz(-0.89620279) q[0];
sx q[0];
rz(2.1620511) q[0];
x q[1];
rz(-2.0252805) q[2];
sx q[2];
rz(-1.8826199) q[2];
sx q[2];
rz(-2.9224967) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9962768) q[1];
sx q[1];
rz(-1.4106803) q[1];
sx q[1];
rz(-0.54460633) q[1];
x q[2];
rz(2.0790948) q[3];
sx q[3];
rz(-0.33514272) q[3];
sx q[3];
rz(2.3329874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.071094461) q[2];
sx q[2];
rz(-1.1134032) q[2];
sx q[2];
rz(-2.5615198) q[2];
rz(-1.5602559) q[3];
sx q[3];
rz(-1.3908849) q[3];
sx q[3];
rz(-1.7177379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73109126) q[0];
sx q[0];
rz(-3.0803362) q[0];
sx q[0];
rz(0.047792338) q[0];
rz(1.2740678) q[1];
sx q[1];
rz(-1.85227) q[1];
sx q[1];
rz(-3.0963669) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84789373) q[0];
sx q[0];
rz(-1.7195104) q[0];
sx q[0];
rz(2.2151164) q[0];
rz(-pi) q[1];
rz(2.0038744) q[2];
sx q[2];
rz(-0.68223729) q[2];
sx q[2];
rz(-3.0909958) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3395244) q[1];
sx q[1];
rz(-1.5763723) q[1];
sx q[1];
rz(-3.1385723) q[1];
rz(-0.43680059) q[3];
sx q[3];
rz(-1.068401) q[3];
sx q[3];
rz(-0.36239195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4101326) q[2];
sx q[2];
rz(-1.2612017) q[2];
sx q[2];
rz(2.251808) q[2];
rz(-0.54640031) q[3];
sx q[3];
rz(-1.0011287) q[3];
sx q[3];
rz(-2.8304097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(-2.6379717) q[0];
sx q[0];
rz(-1.6872971) q[0];
sx q[0];
rz(-1.7244435) q[0];
rz(1.1196989) q[1];
sx q[1];
rz(-1.7599301) q[1];
sx q[1];
rz(-1.959257) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1954297) q[0];
sx q[0];
rz(-2.7481005) q[0];
sx q[0];
rz(2.959743) q[0];
rz(-pi) q[1];
x q[1];
rz(0.77967092) q[2];
sx q[2];
rz(-2.5354249) q[2];
sx q[2];
rz(-0.066628284) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.87403832) q[1];
sx q[1];
rz(-0.8505162) q[1];
sx q[1];
rz(-1.505386) q[1];
rz(-pi) q[2];
rz(2.3919258) q[3];
sx q[3];
rz(-1.9056068) q[3];
sx q[3];
rz(1.0286758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.22415367) q[2];
sx q[2];
rz(-2.3771693) q[2];
sx q[2];
rz(-2.7177641) q[2];
rz(1.1118927) q[3];
sx q[3];
rz(-2.6424776) q[3];
sx q[3];
rz(-0.65139884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79913419) q[0];
sx q[0];
rz(-3.0693711) q[0];
sx q[0];
rz(-0.95241958) q[0];
rz(-0.77596387) q[1];
sx q[1];
rz(-2.2064078) q[1];
sx q[1];
rz(-1.5752569) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1492406) q[0];
sx q[0];
rz(-1.5910205) q[0];
sx q[0];
rz(3.1035191) q[0];
rz(-pi) q[1];
rz(0.96409728) q[2];
sx q[2];
rz(-1.5390054) q[2];
sx q[2];
rz(-2.3111985) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1426865) q[1];
sx q[1];
rz(-0.64780462) q[1];
sx q[1];
rz(-2.7728989) q[1];
x q[2];
rz(2.9569472) q[3];
sx q[3];
rz(-2.054541) q[3];
sx q[3];
rz(0.32102206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7427407) q[2];
sx q[2];
rz(-0.80790085) q[2];
sx q[2];
rz(0.74637949) q[2];
rz(1.0910723) q[3];
sx q[3];
rz(-1.4336136) q[3];
sx q[3];
rz(-2.5892042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93254507) q[0];
sx q[0];
rz(-0.046367558) q[0];
sx q[0];
rz(-2.5982502) q[0];
rz(2.6047193) q[1];
sx q[1];
rz(-2.8165292) q[1];
sx q[1];
rz(-3.0184025) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7808852) q[0];
sx q[0];
rz(-2.505143) q[0];
sx q[0];
rz(1.6223915) q[0];
rz(-pi) q[1];
rz(-2.6307879) q[2];
sx q[2];
rz(-1.9419365) q[2];
sx q[2];
rz(-0.57194158) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.035288485) q[1];
sx q[1];
rz(-2.114944) q[1];
sx q[1];
rz(-0.72334163) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3770272) q[3];
sx q[3];
rz(-0.95029059) q[3];
sx q[3];
rz(2.4303049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8366375) q[2];
sx q[2];
rz(-1.7540437) q[2];
sx q[2];
rz(0.48509625) q[2];
rz(-1.1525611) q[3];
sx q[3];
rz(-1.3644812) q[3];
sx q[3];
rz(0.047957234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7570067) q[0];
sx q[0];
rz(-0.93757367) q[0];
sx q[0];
rz(2.6458929) q[0];
rz(3.0777625) q[1];
sx q[1];
rz(-0.7494691) q[1];
sx q[1];
rz(0.071050342) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4249504) q[0];
sx q[0];
rz(-1.1122238) q[0];
sx q[0];
rz(1.3032662) q[0];
rz(-pi) q[1];
rz(-1.7456648) q[2];
sx q[2];
rz(-0.9728469) q[2];
sx q[2];
rz(-1.1451461) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.20375511) q[1];
sx q[1];
rz(-1.8120288) q[1];
sx q[1];
rz(-2.9822442) q[1];
rz(-pi) q[2];
rz(-0.17545731) q[3];
sx q[3];
rz(-0.4532632) q[3];
sx q[3];
rz(2.3153265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.12585982) q[2];
sx q[2];
rz(-1.4733529) q[2];
sx q[2];
rz(-0.86137548) q[2];
rz(1.3043978) q[3];
sx q[3];
rz(-2.5671037) q[3];
sx q[3];
rz(1.5931574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17643377) q[0];
sx q[0];
rz(-2.6506944) q[0];
sx q[0];
rz(1.7023671) q[0];
rz(-3.0610415) q[1];
sx q[1];
rz(-1.4713902) q[1];
sx q[1];
rz(-1.7505987) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8497365) q[0];
sx q[0];
rz(-2.448659) q[0];
sx q[0];
rz(0.87397184) q[0];
rz(-pi) q[1];
rz(-0.020419196) q[2];
sx q[2];
rz(-1.2660742) q[2];
sx q[2];
rz(-2.8430276) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7841107) q[1];
sx q[1];
rz(-2.3598089) q[1];
sx q[1];
rz(-0.89872054) q[1];
x q[2];
rz(-0.37064044) q[3];
sx q[3];
rz(-1.9270867) q[3];
sx q[3];
rz(-0.72539893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3535658) q[2];
sx q[2];
rz(-1.2949508) q[2];
sx q[2];
rz(-3.0143152) q[2];
rz(0.92653972) q[3];
sx q[3];
rz(-0.24181952) q[3];
sx q[3];
rz(-1.658879) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7552898) q[0];
sx q[0];
rz(-2.4918064) q[0];
sx q[0];
rz(-1.6408828) q[0];
rz(-0.81037784) q[1];
sx q[1];
rz(-2.2893298) q[1];
sx q[1];
rz(-0.77218974) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0836853) q[0];
sx q[0];
rz(-0.9389239) q[0];
sx q[0];
rz(-2.4993012) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5224546) q[2];
sx q[2];
rz(-2.7840814) q[2];
sx q[2];
rz(-2.641397) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4160847) q[1];
sx q[1];
rz(-2.5920752) q[1];
sx q[1];
rz(-1.548442) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2507319) q[3];
sx q[3];
rz(-1.5308342) q[3];
sx q[3];
rz(3.0639632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7605674) q[2];
sx q[2];
rz(-2.3367391) q[2];
sx q[2];
rz(2.5123361) q[2];
rz(-2.4106846) q[3];
sx q[3];
rz(-0.84354246) q[3];
sx q[3];
rz(-1.6985016) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6998941) q[0];
sx q[0];
rz(-0.48342539) q[0];
sx q[0];
rz(-0.40406686) q[0];
rz(0.40456698) q[1];
sx q[1];
rz(-1.4063944) q[1];
sx q[1];
rz(-0.688596) q[1];
rz(-2.8597833) q[2];
sx q[2];
rz(-1.6360313) q[2];
sx q[2];
rz(-2.4665063) q[2];
rz(2.6062905) q[3];
sx q[3];
rz(-1.3131159) q[3];
sx q[3];
rz(2.5406607) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
