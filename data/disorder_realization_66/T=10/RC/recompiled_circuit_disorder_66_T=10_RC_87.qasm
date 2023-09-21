OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.80823094) q[0];
sx q[0];
rz(4.1132676) q[0];
sx q[0];
rz(10.905807) q[0];
rz(-0.12227585) q[1];
sx q[1];
rz(-0.08637698) q[1];
sx q[1];
rz(-3.1187305) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45863736) q[0];
sx q[0];
rz(-1.6874755) q[0];
sx q[0];
rz(1.4215683) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0055662) q[2];
sx q[2];
rz(-2.044894) q[2];
sx q[2];
rz(2.9413162) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.37576807) q[1];
sx q[1];
rz(-1.4210912) q[1];
sx q[1];
rz(1.2234115) q[1];
rz(0.25306563) q[3];
sx q[3];
rz(-0.96732891) q[3];
sx q[3];
rz(-2.8367708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.4188529) q[2];
sx q[2];
rz(-0.8272233) q[2];
sx q[2];
rz(0.76618761) q[2];
rz(-2.9700759) q[3];
sx q[3];
rz(-0.73232108) q[3];
sx q[3];
rz(2.5550301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3294753) q[0];
sx q[0];
rz(-2.7024039) q[0];
sx q[0];
rz(-0.086659327) q[0];
rz(0.86241972) q[1];
sx q[1];
rz(-2.5984867) q[1];
sx q[1];
rz(0.08509732) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25227308) q[0];
sx q[0];
rz(-0.62250455) q[0];
sx q[0];
rz(0.028153367) q[0];
rz(-pi) q[1];
rz(-2.8367963) q[2];
sx q[2];
rz(-2.1509503) q[2];
sx q[2];
rz(0.09588974) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.18499204) q[1];
sx q[1];
rz(-1.4586095) q[1];
sx q[1];
rz(0.086602028) q[1];
x q[2];
rz(-2.8220903) q[3];
sx q[3];
rz(-1.1211683) q[3];
sx q[3];
rz(0.055469661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.26248419) q[2];
sx q[2];
rz(-1.5510677) q[2];
sx q[2];
rz(-1.4651728) q[2];
rz(1.9654467) q[3];
sx q[3];
rz(-2.2361103) q[3];
sx q[3];
rz(2.8951077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6572606) q[0];
sx q[0];
rz(-3.0069139) q[0];
sx q[0];
rz(2.235967) q[0];
rz(-2.8088645) q[1];
sx q[1];
rz(-0.85488027) q[1];
sx q[1];
rz(1.7571626) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2565838) q[0];
sx q[0];
rz(-2.6589767) q[0];
sx q[0];
rz(0.17871876) q[0];
rz(2.076782) q[2];
sx q[2];
rz(-1.9188606) q[2];
sx q[2];
rz(-2.6315174) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3123734) q[1];
sx q[1];
rz(-2.3260498) q[1];
sx q[1];
rz(-1.6014598) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1372041) q[3];
sx q[3];
rz(-1.1096138) q[3];
sx q[3];
rz(-3.0571292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.24421346) q[2];
sx q[2];
rz(-0.34003568) q[2];
sx q[2];
rz(1.8133694) q[2];
rz(1.7437079) q[3];
sx q[3];
rz(-2.6013247) q[3];
sx q[3];
rz(-0.11169294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.06552799) q[0];
sx q[0];
rz(-2.9632443) q[0];
sx q[0];
rz(-0.67101014) q[0];
rz(-1.3440075) q[1];
sx q[1];
rz(-1.9799415) q[1];
sx q[1];
rz(-0.20733325) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3610791) q[0];
sx q[0];
rz(-1.5797537) q[0];
sx q[0];
rz(1.6121959) q[0];
rz(-1.7447128) q[2];
sx q[2];
rz(-1.5408278) q[2];
sx q[2];
rz(-2.7001691) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9878848) q[1];
sx q[1];
rz(-2.8578836) q[1];
sx q[1];
rz(2.1464286) q[1];
rz(-0.025709318) q[3];
sx q[3];
rz(-2.6298012) q[3];
sx q[3];
rz(0.6230841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.31071445) q[2];
sx q[2];
rz(-1.0881311) q[2];
sx q[2];
rz(2.3332398) q[2];
rz(-2.3338142) q[3];
sx q[3];
rz(-2.419796) q[3];
sx q[3];
rz(-0.19367735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5508674) q[0];
sx q[0];
rz(-0.71722513) q[0];
sx q[0];
rz(2.2461058) q[0];
rz(-0.26516178) q[1];
sx q[1];
rz(-0.83509713) q[1];
sx q[1];
rz(-2.6079544) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9981209) q[0];
sx q[0];
rz(-2.7406685) q[0];
sx q[0];
rz(-2.3401767) q[0];
rz(-pi) q[1];
rz(2.3737714) q[2];
sx q[2];
rz(-1.7349093) q[2];
sx q[2];
rz(-2.5196911) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4644949) q[1];
sx q[1];
rz(-0.2650731) q[1];
sx q[1];
rz(0.58880834) q[1];
rz(-pi) q[2];
x q[2];
rz(2.93612) q[3];
sx q[3];
rz(-0.95147248) q[3];
sx q[3];
rz(-1.6831786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2395893) q[2];
sx q[2];
rz(-1.5782372) q[2];
sx q[2];
rz(-2.5732102) q[2];
rz(1.9249632) q[3];
sx q[3];
rz(-0.47074461) q[3];
sx q[3];
rz(2.7980428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50826532) q[0];
sx q[0];
rz(-0.87843043) q[0];
sx q[0];
rz(-1.9653962) q[0];
rz(-1.5559224) q[1];
sx q[1];
rz(-2.1891179) q[1];
sx q[1];
rz(1.0046545) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83197901) q[0];
sx q[0];
rz(-1.2011315) q[0];
sx q[0];
rz(1.2020822) q[0];
x q[1];
rz(2.7447694) q[2];
sx q[2];
rz(-1.6453711) q[2];
sx q[2];
rz(-0.65149388) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.55014729) q[1];
sx q[1];
rz(-2.1131556) q[1];
sx q[1];
rz(-1.28246) q[1];
rz(1.5435013) q[3];
sx q[3];
rz(-2.7276162) q[3];
sx q[3];
rz(-2.9881145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4040318) q[2];
sx q[2];
rz(-0.99016756) q[2];
sx q[2];
rz(0.44815865) q[2];
rz(0.29799497) q[3];
sx q[3];
rz(-2.4500676) q[3];
sx q[3];
rz(-2.7929849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77666831) q[0];
sx q[0];
rz(-1.8402599) q[0];
sx q[0];
rz(0.65814322) q[0];
rz(1.8572042) q[1];
sx q[1];
rz(-2.6874266) q[1];
sx q[1];
rz(-0.67214322) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9854239) q[0];
sx q[0];
rz(-0.16616136) q[0];
sx q[0];
rz(1.6155433) q[0];
rz(-pi) q[1];
rz(0.80588801) q[2];
sx q[2];
rz(-0.58008367) q[2];
sx q[2];
rz(2.9758331) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0828515) q[1];
sx q[1];
rz(-2.5141659) q[1];
sx q[1];
rz(2.4861366) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4403205) q[3];
sx q[3];
rz(-1.6810732) q[3];
sx q[3];
rz(2.9993204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.60823524) q[2];
sx q[2];
rz(-2.50864) q[2];
sx q[2];
rz(-0.39880025) q[2];
rz(2.3170025) q[3];
sx q[3];
rz(-1.3922858) q[3];
sx q[3];
rz(-2.5430172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7010715) q[0];
sx q[0];
rz(-2.7109787) q[0];
sx q[0];
rz(-2.9833802) q[0];
rz(2.054706) q[1];
sx q[1];
rz(-0.67651665) q[1];
sx q[1];
rz(-0.80668443) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0643443) q[0];
sx q[0];
rz(-1.6883388) q[0];
sx q[0];
rz(-2.006152) q[0];
rz(-pi) q[1];
rz(-2.9358747) q[2];
sx q[2];
rz(-1.6763902) q[2];
sx q[2];
rz(-0.92668698) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.403277) q[1];
sx q[1];
rz(-0.74790955) q[1];
sx q[1];
rz(0.43055375) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1095803) q[3];
sx q[3];
rz(-1.8714928) q[3];
sx q[3];
rz(-0.34919448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5640101) q[2];
sx q[2];
rz(-1.492604) q[2];
sx q[2];
rz(-2.1419443) q[2];
rz(-0.10351652) q[3];
sx q[3];
rz(-3.0045356) q[3];
sx q[3];
rz(2.0172393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12061159) q[0];
sx q[0];
rz(-2.3636901) q[0];
sx q[0];
rz(-2.4840684) q[0];
rz(2.855037) q[1];
sx q[1];
rz(-0.92266881) q[1];
sx q[1];
rz(0.34067571) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9690902) q[0];
sx q[0];
rz(-1.2524676) q[0];
sx q[0];
rz(0.83842917) q[0];
rz(-pi) q[1];
rz(1.0057783) q[2];
sx q[2];
rz(-1.3823595) q[2];
sx q[2];
rz(-0.68067683) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.77010158) q[1];
sx q[1];
rz(-1.0189459) q[1];
sx q[1];
rz(0.21367195) q[1];
rz(1.7394003) q[3];
sx q[3];
rz(-1.7662893) q[3];
sx q[3];
rz(0.80381264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.75491536) q[2];
sx q[2];
rz(-1.200054) q[2];
sx q[2];
rz(0.49794751) q[2];
rz(2.8399816) q[3];
sx q[3];
rz(-2.7409654) q[3];
sx q[3];
rz(2.6224459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4122445) q[0];
sx q[0];
rz(-3.0817741) q[0];
sx q[0];
rz(-0.27858946) q[0];
rz(2.5623698) q[1];
sx q[1];
rz(-0.93943739) q[1];
sx q[1];
rz(0.07671193) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.840832) q[0];
sx q[0];
rz(-1.6941841) q[0];
sx q[0];
rz(1.8873909) q[0];
rz(-0.99761988) q[2];
sx q[2];
rz(-1.4263037) q[2];
sx q[2];
rz(-1.0695374) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0127276) q[1];
sx q[1];
rz(-2.7755133) q[1];
sx q[1];
rz(-0.058458316) q[1];
rz(-pi) q[2];
x q[2];
rz(0.17681392) q[3];
sx q[3];
rz(-1.483695) q[3];
sx q[3];
rz(0.84482312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.77438337) q[2];
sx q[2];
rz(-2.4343906) q[2];
sx q[2];
rz(0.45551604) q[2];
rz(-2.7101743) q[3];
sx q[3];
rz(-0.12532561) q[3];
sx q[3];
rz(3.07807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1800304) q[0];
sx q[0];
rz(-1.4008235) q[0];
sx q[0];
rz(0.93749198) q[0];
rz(0.58615276) q[1];
sx q[1];
rz(-1.502232) q[1];
sx q[1];
rz(-1.4486817) q[1];
rz(-2.0417487) q[2];
sx q[2];
rz(-1.3146613) q[2];
sx q[2];
rz(1.4801499) q[2];
rz(-1.9361817) q[3];
sx q[3];
rz(-1.972651) q[3];
sx q[3];
rz(-2.7046711) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];