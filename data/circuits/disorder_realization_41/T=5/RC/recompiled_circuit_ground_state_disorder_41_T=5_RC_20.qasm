OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5557033) q[0];
sx q[0];
rz(4.3309431) q[0];
sx q[0];
rz(10.588607) q[0];
rz(-3.8883348) q[1];
sx q[1];
rz(0.26489869) q[1];
sx q[1];
rz(12.737552) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2585325) q[0];
sx q[0];
rz(-2.8730321) q[0];
sx q[0];
rz(0.95897755) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2581882) q[2];
sx q[2];
rz(-2.2210178) q[2];
sx q[2];
rz(-1.5387077) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.491558) q[1];
sx q[1];
rz(-1.2753873) q[1];
sx q[1];
rz(3.1409688) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.016841932) q[3];
sx q[3];
rz(-1.5653658) q[3];
sx q[3];
rz(3.0433486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.050194) q[2];
sx q[2];
rz(-1.564743) q[2];
sx q[2];
rz(-2.0719299) q[2];
rz(1.7012677) q[3];
sx q[3];
rz(-1.0245208) q[3];
sx q[3];
rz(1.9694156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1759724) q[0];
sx q[0];
rz(-1.4485285) q[0];
sx q[0];
rz(0.58797055) q[0];
rz(-1.2929471) q[1];
sx q[1];
rz(-2.1473532) q[1];
sx q[1];
rz(2.9867244) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.800134) q[0];
sx q[0];
rz(-1.4636369) q[0];
sx q[0];
rz(0.5002621) q[0];
rz(0.54041086) q[2];
sx q[2];
rz(-0.27920461) q[2];
sx q[2];
rz(2.0824733) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.57666393) q[1];
sx q[1];
rz(-1.1047078) q[1];
sx q[1];
rz(-0.376815) q[1];
rz(2.5741379) q[3];
sx q[3];
rz(-2.8914872) q[3];
sx q[3];
rz(-2.9222744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.64314848) q[2];
sx q[2];
rz(-0.73068205) q[2];
sx q[2];
rz(-0.11239642) q[2];
rz(0.82320881) q[3];
sx q[3];
rz(-1.221849) q[3];
sx q[3];
rz(0.52197758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4560029) q[0];
sx q[0];
rz(-2.0528448) q[0];
sx q[0];
rz(2.1734557) q[0];
rz(-2.0227506) q[1];
sx q[1];
rz(-1.5348996) q[1];
sx q[1];
rz(1.3333837) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3385425) q[0];
sx q[0];
rz(-2.1678939) q[0];
sx q[0];
rz(-0.19495585) q[0];
x q[1];
rz(2.5081464) q[2];
sx q[2];
rz(-0.70745969) q[2];
sx q[2];
rz(-0.58282575) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1044429) q[1];
sx q[1];
rz(-0.57174505) q[1];
sx q[1];
rz(1.3859831) q[1];
rz(-3.0564752) q[3];
sx q[3];
rz(-1.9276351) q[3];
sx q[3];
rz(2.3841355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.94990388) q[2];
sx q[2];
rz(-1.0368985) q[2];
sx q[2];
rz(2.8766768) q[2];
rz(2.8252937) q[3];
sx q[3];
rz(-0.33360544) q[3];
sx q[3];
rz(-1.3914289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.815627) q[0];
sx q[0];
rz(-2.9209324) q[0];
sx q[0];
rz(-1.4937481) q[0];
rz(-1.4119459) q[1];
sx q[1];
rz(-2.1144046) q[1];
sx q[1];
rz(0.6573917) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43577584) q[0];
sx q[0];
rz(-3.0141493) q[0];
sx q[0];
rz(-2.5246922) q[0];
x q[1];
rz(-0.78423402) q[2];
sx q[2];
rz(-2.0080697) q[2];
sx q[2];
rz(-2.0511049) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.79815976) q[1];
sx q[1];
rz(-1.533759) q[1];
sx q[1];
rz(-1.7658893) q[1];
x q[2];
rz(-2.3393624) q[3];
sx q[3];
rz(-2.1539237) q[3];
sx q[3];
rz(2.1183573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.21939453) q[2];
sx q[2];
rz(-1.2224835) q[2];
sx q[2];
rz(0.40327367) q[2];
rz(-1.9379617) q[3];
sx q[3];
rz(-1.5091242) q[3];
sx q[3];
rz(0.50886124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68172425) q[0];
sx q[0];
rz(-2.7175856) q[0];
sx q[0];
rz(0.62393171) q[0];
rz(0.722018) q[1];
sx q[1];
rz(-1.7941509) q[1];
sx q[1];
rz(-3.0375979) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34891221) q[0];
sx q[0];
rz(-2.4390872) q[0];
sx q[0];
rz(0.2200495) q[0];
x q[1];
rz(-0.76461069) q[2];
sx q[2];
rz(-0.48073623) q[2];
sx q[2];
rz(1.3388456) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.38253575) q[1];
sx q[1];
rz(-1.0604572) q[1];
sx q[1];
rz(-2.3063763) q[1];
x q[2];
rz(-2.4071619) q[3];
sx q[3];
rz(-1.0356257) q[3];
sx q[3];
rz(-1.9246235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7907052) q[2];
sx q[2];
rz(-0.71983379) q[2];
sx q[2];
rz(0.78901115) q[2];
rz(-1.2645432) q[3];
sx q[3];
rz(-1.0435373) q[3];
sx q[3];
rz(2.1544382) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8378545) q[0];
sx q[0];
rz(-0.721295) q[0];
sx q[0];
rz(-0.20027941) q[0];
rz(-1.6750977) q[1];
sx q[1];
rz(-1.3267582) q[1];
sx q[1];
rz(1.6237578) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12833327) q[0];
sx q[0];
rz(-0.64494123) q[0];
sx q[0];
rz(-0.63887424) q[0];
x q[1];
rz(-1.3369722) q[2];
sx q[2];
rz(-0.5153729) q[2];
sx q[2];
rz(-0.17763868) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1788097) q[1];
sx q[1];
rz(-0.97985044) q[1];
sx q[1];
rz(-1.5750118) q[1];
rz(-2.9041457) q[3];
sx q[3];
rz(-2.3288764) q[3];
sx q[3];
rz(2.0456912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.62310654) q[2];
sx q[2];
rz(-0.97272626) q[2];
sx q[2];
rz(-2.8002807) q[2];
rz(0.85706472) q[3];
sx q[3];
rz(-1.2271481) q[3];
sx q[3];
rz(-2.6763693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(1.051006) q[0];
sx q[0];
rz(-1.0336646) q[0];
sx q[0];
rz(1.2475913) q[0];
rz(3.0233851) q[1];
sx q[1];
rz(-1.8659614) q[1];
sx q[1];
rz(0.035621312) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44862952) q[0];
sx q[0];
rz(-1.5933196) q[0];
sx q[0];
rz(-1.453214) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9526677) q[2];
sx q[2];
rz(-1.1729203) q[2];
sx q[2];
rz(-2.973345) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5069711) q[1];
sx q[1];
rz(-1.8160607) q[1];
sx q[1];
rz(-2.5140106) q[1];
rz(3.0183718) q[3];
sx q[3];
rz(-0.78500063) q[3];
sx q[3];
rz(0.30320019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.15141618) q[2];
sx q[2];
rz(-1.6665062) q[2];
sx q[2];
rz(2.0659633) q[2];
rz(-0.52078024) q[3];
sx q[3];
rz(-2.5213089) q[3];
sx q[3];
rz(-1.5889997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-2.2369279) q[0];
sx q[0];
rz(-2.1187145) q[0];
sx q[0];
rz(-2.2383595) q[0];
rz(-1.6385993) q[1];
sx q[1];
rz(-1.6939949) q[1];
sx q[1];
rz(1.8797871) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1539858) q[0];
sx q[0];
rz(-1.7075451) q[0];
sx q[0];
rz(1.8022949) q[0];
x q[1];
rz(-0.52623279) q[2];
sx q[2];
rz(-1.1723435) q[2];
sx q[2];
rz(0.72771074) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.256989) q[1];
sx q[1];
rz(-1.6620364) q[1];
sx q[1];
rz(-1.6216519) q[1];
x q[2];
rz(0.66916211) q[3];
sx q[3];
rz(-0.72297308) q[3];
sx q[3];
rz(-2.1694136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5624076) q[2];
sx q[2];
rz(-2.4895442) q[2];
sx q[2];
rz(3.0273666) q[2];
rz(-1.7646029) q[3];
sx q[3];
rz(-2.0012794) q[3];
sx q[3];
rz(2.637114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1219516) q[0];
sx q[0];
rz(-2.1421102) q[0];
sx q[0];
rz(-1.9695388) q[0];
rz(1.7578341) q[1];
sx q[1];
rz(-1.9969321) q[1];
sx q[1];
rz(0.11464548) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5929554) q[0];
sx q[0];
rz(-1.0186188) q[0];
sx q[0];
rz(2.1594285) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9681712) q[2];
sx q[2];
rz(-0.49931128) q[2];
sx q[2];
rz(2.2713646) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0681359) q[1];
sx q[1];
rz(-2.8836859) q[1];
sx q[1];
rz(-1.8515481) q[1];
rz(1.4632567) q[3];
sx q[3];
rz(-0.77725526) q[3];
sx q[3];
rz(1.4449256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5961479) q[2];
sx q[2];
rz(-2.9344276) q[2];
sx q[2];
rz(-0.88085112) q[2];
rz(-2.7018069) q[3];
sx q[3];
rz(-1.1652911) q[3];
sx q[3];
rz(-1.1254651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83774829) q[0];
sx q[0];
rz(-1.4864018) q[0];
sx q[0];
rz(-3.0434171) q[0];
rz(-0.57442609) q[1];
sx q[1];
rz(-1.6822633) q[1];
sx q[1];
rz(-2.548545) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7710073) q[0];
sx q[0];
rz(-0.86261504) q[0];
sx q[0];
rz(2.7561491) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1109652) q[2];
sx q[2];
rz(-1.9367276) q[2];
sx q[2];
rz(2.1618787) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.39506158) q[1];
sx q[1];
rz(-1.5575002) q[1];
sx q[1];
rz(-0.011237267) q[1];
rz(-1.2393059) q[3];
sx q[3];
rz(-0.96486366) q[3];
sx q[3];
rz(0.87406033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9888931) q[2];
sx q[2];
rz(-1.9579192) q[2];
sx q[2];
rz(0.93842554) q[2];
rz(-1.7484131) q[3];
sx q[3];
rz(-1.8311071) q[3];
sx q[3];
rz(0.2383298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68168454) q[0];
sx q[0];
rz(-2.5354698) q[0];
sx q[0];
rz(1.348362) q[0];
rz(1.0472736) q[1];
sx q[1];
rz(-0.92887639) q[1];
sx q[1];
rz(2.1705719) q[1];
rz(0.96228308) q[2];
sx q[2];
rz(-3.0627927) q[2];
sx q[2];
rz(0.40712955) q[2];
rz(-0.88225928) q[3];
sx q[3];
rz(-2.8591446) q[3];
sx q[3];
rz(2.1459116) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
