OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.20733362) q[0];
sx q[0];
rz(3.7319558) q[0];
sx q[0];
rz(9.0537602) q[0];
rz(-0.38129216) q[1];
sx q[1];
rz(2.5420904) q[1];
sx q[1];
rz(11.190344) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3289514) q[0];
sx q[0];
rz(-2.0094123) q[0];
sx q[0];
rz(-0.89226512) q[0];
rz(2.5901428) q[2];
sx q[2];
rz(-0.79556888) q[2];
sx q[2];
rz(-1.431682) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8810597) q[1];
sx q[1];
rz(-0.36306371) q[1];
sx q[1];
rz(2.0102324) q[1];
x q[2];
rz(-1.514643) q[3];
sx q[3];
rz(-0.75222844) q[3];
sx q[3];
rz(0.55263954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.084289) q[2];
sx q[2];
rz(-2.7412582) q[2];
sx q[2];
rz(-2.1526745) q[2];
rz(-2.3890498) q[3];
sx q[3];
rz(-1.1458594) q[3];
sx q[3];
rz(-0.83077103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9343524) q[0];
sx q[0];
rz(-3.0293284) q[0];
sx q[0];
rz(1.1799312) q[0];
rz(-0.99769366) q[1];
sx q[1];
rz(-1.8583863) q[1];
sx q[1];
rz(0.72431272) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9162826) q[0];
sx q[0];
rz(-2.3819469) q[0];
sx q[0];
rz(-1.1067953) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.486859) q[2];
sx q[2];
rz(-1.7695904) q[2];
sx q[2];
rz(2.4059911) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.95485479) q[1];
sx q[1];
rz(-1.9222944) q[1];
sx q[1];
rz(2.9224706) q[1];
rz(-pi) q[2];
rz(-0.23726666) q[3];
sx q[3];
rz(-1.7627343) q[3];
sx q[3];
rz(-2.3290079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2362242) q[2];
sx q[2];
rz(-1.3304109) q[2];
sx q[2];
rz(2.779707) q[2];
rz(0.13606717) q[3];
sx q[3];
rz(-2.5858904) q[3];
sx q[3];
rz(-3.0959685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1419462) q[0];
sx q[0];
rz(-2.0479585) q[0];
sx q[0];
rz(1.746159) q[0];
rz(2.6793001) q[1];
sx q[1];
rz(-2.7170083) q[1];
sx q[1];
rz(-1.2190855) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65788668) q[0];
sx q[0];
rz(-1.833796) q[0];
sx q[0];
rz(-1.2067243) q[0];
x q[1];
rz(-0.8692603) q[2];
sx q[2];
rz(-1.7995036) q[2];
sx q[2];
rz(1.9321835) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0137274) q[1];
sx q[1];
rz(-1.8121108) q[1];
sx q[1];
rz(-1.024854) q[1];
x q[2];
rz(-1.6270646) q[3];
sx q[3];
rz(-2.8731822) q[3];
sx q[3];
rz(-2.066582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7200155) q[2];
sx q[2];
rz(-2.4013459) q[2];
sx q[2];
rz(-1.6960309) q[2];
rz(-2.5727663) q[3];
sx q[3];
rz(-0.84972644) q[3];
sx q[3];
rz(-0.11051699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22359426) q[0];
sx q[0];
rz(-0.44819865) q[0];
sx q[0];
rz(-2.6233327) q[0];
rz(2.4261684) q[1];
sx q[1];
rz(-1.1162858) q[1];
sx q[1];
rz(0.82675654) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2430192) q[0];
sx q[0];
rz(-1.9802226) q[0];
sx q[0];
rz(2.5893674) q[0];
x q[1];
rz(3.0024198) q[2];
sx q[2];
rz(-2.3690802) q[2];
sx q[2];
rz(3.006209) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3863694) q[1];
sx q[1];
rz(-1.1896903) q[1];
sx q[1];
rz(-1.2632881) q[1];
rz(2.7923146) q[3];
sx q[3];
rz(-0.75540245) q[3];
sx q[3];
rz(-2.3625771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.15239079) q[2];
sx q[2];
rz(-0.19897904) q[2];
sx q[2];
rz(1.3789122) q[2];
rz(-0.072323024) q[3];
sx q[3];
rz(-0.81243378) q[3];
sx q[3];
rz(1.6453843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.3146661) q[0];
sx q[0];
rz(-2.5214654) q[0];
sx q[0];
rz(1.1258874) q[0];
rz(2.2391438) q[1];
sx q[1];
rz(-2.1676962) q[1];
sx q[1];
rz(2.856423) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0641159) q[0];
sx q[0];
rz(-0.47668326) q[0];
sx q[0];
rz(-1.478273) q[0];
rz(-0.6638078) q[2];
sx q[2];
rz(-1.8823349) q[2];
sx q[2];
rz(0.89154348) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3045029) q[1];
sx q[1];
rz(-1.2577406) q[1];
sx q[1];
rz(-3.1411527) q[1];
rz(-pi) q[2];
rz(0.40163715) q[3];
sx q[3];
rz(-1.4816227) q[3];
sx q[3];
rz(-3.0107486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8482762) q[2];
sx q[2];
rz(-2.6360376) q[2];
sx q[2];
rz(-2.6021393) q[2];
rz(-2.8347677) q[3];
sx q[3];
rz(-2.2570733) q[3];
sx q[3];
rz(-0.45421281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8834615) q[0];
sx q[0];
rz(-2.3911609) q[0];
sx q[0];
rz(-0.15701292) q[0];
rz(2.4482588) q[1];
sx q[1];
rz(-0.88070977) q[1];
sx q[1];
rz(1.7745811) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1098423) q[0];
sx q[0];
rz(-1.6037918) q[0];
sx q[0];
rz(-1.6389636) q[0];
rz(2.5622257) q[2];
sx q[2];
rz(-2.2420792) q[2];
sx q[2];
rz(0.52057779) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4821266) q[1];
sx q[1];
rz(-0.30059338) q[1];
sx q[1];
rz(-1.8468922) q[1];
rz(-pi) q[2];
rz(-1.4671765) q[3];
sx q[3];
rz(-1.2195671) q[3];
sx q[3];
rz(-2.2675089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8453025) q[2];
sx q[2];
rz(-0.85001105) q[2];
sx q[2];
rz(2.7381251) q[2];
rz(0.48163313) q[3];
sx q[3];
rz(-1.0721595) q[3];
sx q[3];
rz(2.6223555) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8994609) q[0];
sx q[0];
rz(-2.2583028) q[0];
sx q[0];
rz(-2.2677299) q[0];
rz(-2.6938687) q[1];
sx q[1];
rz(-2.402585) q[1];
sx q[1];
rz(-1.9708995) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0543594) q[0];
sx q[0];
rz(-1.5938252) q[0];
sx q[0];
rz(-1.5769632) q[0];
x q[1];
rz(-2.7295693) q[2];
sx q[2];
rz(-3.014866) q[2];
sx q[2];
rz(2.3840981) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2504187) q[1];
sx q[1];
rz(-2.182057) q[1];
sx q[1];
rz(-0.82959081) q[1];
rz(-2.6550754) q[3];
sx q[3];
rz(-2.5518637) q[3];
sx q[3];
rz(0.041134838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0447023) q[2];
sx q[2];
rz(-2.5723852) q[2];
sx q[2];
rz(0.61075413) q[2];
rz(2.6664873) q[3];
sx q[3];
rz(-2.0510309) q[3];
sx q[3];
rz(-0.92774123) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24191813) q[0];
sx q[0];
rz(-3.0245259) q[0];
sx q[0];
rz(-2.8444667) q[0];
rz(-1.7469453) q[1];
sx q[1];
rz(-1.1509832) q[1];
sx q[1];
rz(-2.4954605) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1387678) q[0];
sx q[0];
rz(-1.3647623) q[0];
sx q[0];
rz(3.0384008) q[0];
x q[1];
rz(-0.32142873) q[2];
sx q[2];
rz(-2.6476963) q[2];
sx q[2];
rz(-2.0733881) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7568126) q[1];
sx q[1];
rz(-1.9128748) q[1];
sx q[1];
rz(0.53201075) q[1];
x q[2];
rz(0.93692245) q[3];
sx q[3];
rz(-2.2735032) q[3];
sx q[3];
rz(-3.1261409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0769161) q[2];
sx q[2];
rz(-0.95305324) q[2];
sx q[2];
rz(0.45483744) q[2];
rz(0.70139766) q[3];
sx q[3];
rz(-1.0304136) q[3];
sx q[3];
rz(1.1340244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4421473) q[0];
sx q[0];
rz(-3*pi/16) q[0];
sx q[0];
rz(-0.79750693) q[0];
rz(-2.6240255) q[1];
sx q[1];
rz(-0.8126173) q[1];
sx q[1];
rz(3.033175) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78128821) q[0];
sx q[0];
rz(-0.66626781) q[0];
sx q[0];
rz(-1.6643307) q[0];
rz(0.046594521) q[2];
sx q[2];
rz(-1.3923044) q[2];
sx q[2];
rz(-0.5878693) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0260967) q[1];
sx q[1];
rz(-0.27448389) q[1];
sx q[1];
rz(0.91067578) q[1];
rz(-pi) q[2];
rz(-1.7897723) q[3];
sx q[3];
rz(-0.7162381) q[3];
sx q[3];
rz(2.4718474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1291528) q[2];
sx q[2];
rz(-1.7947349) q[2];
sx q[2];
rz(0.33995315) q[2];
rz(2.7231976) q[3];
sx q[3];
rz(-2.5451626) q[3];
sx q[3];
rz(2.4160014) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6091992) q[0];
sx q[0];
rz(-2.7476855) q[0];
sx q[0];
rz(-0.6788196) q[0];
rz(0.36418307) q[1];
sx q[1];
rz(-1.697425) q[1];
sx q[1];
rz(3.0864339) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77627742) q[0];
sx q[0];
rz(-0.3826097) q[0];
sx q[0];
rz(0.17392735) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8137384) q[2];
sx q[2];
rz(-0.58460669) q[2];
sx q[2];
rz(0.84529982) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2436279) q[1];
sx q[1];
rz(-0.98476714) q[1];
sx q[1];
rz(0.90378739) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.03667128) q[3];
sx q[3];
rz(-1.1150556) q[3];
sx q[3];
rz(-0.40487056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.98383343) q[2];
sx q[2];
rz(-1.021421) q[2];
sx q[2];
rz(2.6514163) q[2];
rz(-0.13752078) q[3];
sx q[3];
rz(-1.0995882) q[3];
sx q[3];
rz(0.93808758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72538439) q[0];
sx q[0];
rz(-1.2513456) q[0];
sx q[0];
rz(-0.67847897) q[0];
rz(-0.2086808) q[1];
sx q[1];
rz(-1.122767) q[1];
sx q[1];
rz(-1.541419) q[1];
rz(0.29253929) q[2];
sx q[2];
rz(-0.32513041) q[2];
sx q[2];
rz(0.61376094) q[2];
rz(2.901554) q[3];
sx q[3];
rz(-1.6143027) q[3];
sx q[3];
rz(-1.3054813) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
