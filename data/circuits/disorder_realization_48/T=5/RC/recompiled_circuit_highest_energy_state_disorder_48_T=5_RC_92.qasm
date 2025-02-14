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
rz(0.91520619) q[0];
sx q[0];
rz(2.6829166) q[0];
sx q[0];
rz(12.248326) q[0];
rz(1.3660499) q[1];
sx q[1];
rz(-2.4429758) q[1];
sx q[1];
rz(-3.0787992) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2554277) q[0];
sx q[0];
rz(-2.0768099) q[0];
sx q[0];
rz(-0.54765986) q[0];
rz(1.0452966) q[2];
sx q[2];
rz(-0.7514779) q[2];
sx q[2];
rz(2.1906302) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.54431984) q[1];
sx q[1];
rz(-1.8493127) q[1];
sx q[1];
rz(-1.0611978) q[1];
x q[2];
rz(-2.1413689) q[3];
sx q[3];
rz(-0.65815845) q[3];
sx q[3];
rz(-0.47129813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.723145) q[2];
sx q[2];
rz(-1.8258839) q[2];
sx q[2];
rz(-2.1212228) q[2];
rz(-0.72888914) q[3];
sx q[3];
rz(-0.43292361) q[3];
sx q[3];
rz(1.5322022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.51434022) q[0];
sx q[0];
rz(-1.9753375) q[0];
sx q[0];
rz(-0.038272055) q[0];
rz(3.017784) q[1];
sx q[1];
rz(-2.6688711) q[1];
sx q[1];
rz(1.570787) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8077652) q[0];
sx q[0];
rz(-0.021919576) q[0];
sx q[0];
rz(1.0967361) q[0];
x q[1];
rz(1.6699617) q[2];
sx q[2];
rz(-2.3395774) q[2];
sx q[2];
rz(-0.85814171) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.60197406) q[1];
sx q[1];
rz(-1.5718196) q[1];
sx q[1];
rz(0.00042331528) q[1];
x q[2];
rz(-1.6917211) q[3];
sx q[3];
rz(-2.0837415) q[3];
sx q[3];
rz(0.67739048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6191285) q[2];
sx q[2];
rz(-1.8121441) q[2];
sx q[2];
rz(0.5788571) q[2];
rz(2.0959496) q[3];
sx q[3];
rz(-2.3561616) q[3];
sx q[3];
rz(0.75700179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4724562) q[0];
sx q[0];
rz(-2.0643015) q[0];
sx q[0];
rz(-0.45397595) q[0];
rz(2.9767735) q[1];
sx q[1];
rz(-2.3655128) q[1];
sx q[1];
rz(-1.4705315) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8225518) q[0];
sx q[0];
rz(-1.7421573) q[0];
sx q[0];
rz(-1.1771855) q[0];
x q[1];
rz(2.8118089) q[2];
sx q[2];
rz(-2.0726786) q[2];
sx q[2];
rz(1.9374073) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0297894) q[1];
sx q[1];
rz(-1.8743427) q[1];
sx q[1];
rz(-0.70036594) q[1];
x q[2];
rz(2.4824247) q[3];
sx q[3];
rz(-1.9999095) q[3];
sx q[3];
rz(1.3339213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7620324) q[2];
sx q[2];
rz(-1.6635868) q[2];
sx q[2];
rz(1.7851768) q[2];
rz(2.6421269) q[3];
sx q[3];
rz(-1.5288589) q[3];
sx q[3];
rz(2.0904026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22719638) q[0];
sx q[0];
rz(-2.2375186) q[0];
sx q[0];
rz(2.0830578) q[0];
rz(-1.9805485) q[1];
sx q[1];
rz(-2.5807022) q[1];
sx q[1];
rz(-3.0243691) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49699052) q[0];
sx q[0];
rz(-1.0287026) q[0];
sx q[0];
rz(1.5250456) q[0];
rz(-0.92553301) q[2];
sx q[2];
rz(-1.4129854) q[2];
sx q[2];
rz(0.35343808) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1438751) q[1];
sx q[1];
rz(-1.2923601) q[1];
sx q[1];
rz(-0.77634676) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53026188) q[3];
sx q[3];
rz(-2.432309) q[3];
sx q[3];
rz(-1.5269296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.78493541) q[2];
sx q[2];
rz(-1.6931809) q[2];
sx q[2];
rz(1.3680722) q[2];
rz(-1.8562227) q[3];
sx q[3];
rz(-1.1077489) q[3];
sx q[3];
rz(0.72803298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0990937) q[0];
sx q[0];
rz(-2.2137401) q[0];
sx q[0];
rz(1.7294783) q[0];
rz(-2.1389351) q[1];
sx q[1];
rz(-0.24191562) q[1];
sx q[1];
rz(-1.5589421) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6928685) q[0];
sx q[0];
rz(-2.4131615) q[0];
sx q[0];
rz(2.4850575) q[0];
rz(-pi) q[1];
rz(0.095076128) q[2];
sx q[2];
rz(-0.72059599) q[2];
sx q[2];
rz(2.528185) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6465368) q[1];
sx q[1];
rz(-0.55083129) q[1];
sx q[1];
rz(-1.6158094) q[1];
x q[2];
rz(-1.2491426) q[3];
sx q[3];
rz(-1.6603004) q[3];
sx q[3];
rz(0.22150515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.58854181) q[2];
sx q[2];
rz(-1.7091227) q[2];
sx q[2];
rz(-2.7480965) q[2];
rz(-2.1816175) q[3];
sx q[3];
rz(-2.4656651) q[3];
sx q[3];
rz(2.1737607) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8680962) q[0];
sx q[0];
rz(-0.12729004) q[0];
sx q[0];
rz(-0.18336503) q[0];
rz(2.6496437) q[1];
sx q[1];
rz(-2.349647) q[1];
sx q[1];
rz(-1.9221745) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47226071) q[0];
sx q[0];
rz(-0.12024256) q[0];
sx q[0];
rz(-0.73445102) q[0];
rz(-pi) q[1];
rz(1.9398795) q[2];
sx q[2];
rz(-3.0574692) q[2];
sx q[2];
rz(0.28757986) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.26039429) q[1];
sx q[1];
rz(-1.9250084) q[1];
sx q[1];
rz(-2.7071035) q[1];
rz(1.0991092) q[3];
sx q[3];
rz(-1.6460643) q[3];
sx q[3];
rz(-1.4510297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0101953) q[2];
sx q[2];
rz(-1.4364028) q[2];
sx q[2];
rz(0.08237002) q[2];
rz(2.2149337) q[3];
sx q[3];
rz(-0.72946531) q[3];
sx q[3];
rz(-1.6389219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53967404) q[0];
sx q[0];
rz(-2.813756) q[0];
sx q[0];
rz(2.4251921) q[0];
rz(-0.40388233) q[1];
sx q[1];
rz(-1.5733746) q[1];
sx q[1];
rz(1.7049888) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5157455) q[0];
sx q[0];
rz(-1.4270009) q[0];
sx q[0];
rz(1.189144) q[0];
rz(-1.4384957) q[2];
sx q[2];
rz(-1.0463593) q[2];
sx q[2];
rz(-0.99614267) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.77715141) q[1];
sx q[1];
rz(-1.6702515) q[1];
sx q[1];
rz(-0.39945368) q[1];
rz(-0.8248803) q[3];
sx q[3];
rz(-0.91478225) q[3];
sx q[3];
rz(-3.0716621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7776103) q[2];
sx q[2];
rz(-2.1668375) q[2];
sx q[2];
rz(-3.0762365) q[2];
rz(0.86723793) q[3];
sx q[3];
rz(-1.5012274) q[3];
sx q[3];
rz(-1.8170099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6548178) q[0];
sx q[0];
rz(-1.7205394) q[0];
sx q[0];
rz(-2.4105893) q[0];
rz(0.77516088) q[1];
sx q[1];
rz(-2.8072) q[1];
sx q[1];
rz(0.41466546) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1486007) q[0];
sx q[0];
rz(-0.96382695) q[0];
sx q[0];
rz(-0.66120045) q[0];
rz(0.35279775) q[2];
sx q[2];
rz(-2.5830306) q[2];
sx q[2];
rz(1.1272023) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6952299) q[1];
sx q[1];
rz(-1.2712423) q[1];
sx q[1];
rz(-1.3079554) q[1];
rz(-pi) q[2];
rz(-0.13402588) q[3];
sx q[3];
rz(-1.7334337) q[3];
sx q[3];
rz(0.37510474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.23184648) q[2];
sx q[2];
rz(-0.48603386) q[2];
sx q[2];
rz(3.0492142) q[2];
rz(-2.3653024) q[3];
sx q[3];
rz(-1.4624566) q[3];
sx q[3];
rz(0.20364729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48731503) q[0];
sx q[0];
rz(-1.4947816) q[0];
sx q[0];
rz(0.017729433) q[0];
rz(-2.0954633) q[1];
sx q[1];
rz(-1.7283864) q[1];
sx q[1];
rz(1.0606934) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26955206) q[0];
sx q[0];
rz(-1.5020796) q[0];
sx q[0];
rz(-1.5848573) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1589375) q[2];
sx q[2];
rz(-1.7870296) q[2];
sx q[2];
rz(0.71046358) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2613147) q[1];
sx q[1];
rz(-2.3572817) q[1];
sx q[1];
rz(-1.8159291) q[1];
rz(-pi) q[2];
rz(1.5117731) q[3];
sx q[3];
rz(-2.0938329) q[3];
sx q[3];
rz(-0.28429261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1549418) q[2];
sx q[2];
rz(-0.81601802) q[2];
sx q[2];
rz(-2.4972534) q[2];
rz(0.84195697) q[3];
sx q[3];
rz(-1.569845) q[3];
sx q[3];
rz(-2.2346066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3802721) q[0];
sx q[0];
rz(-0.43571061) q[0];
sx q[0];
rz(2.6897588) q[0];
rz(-2.1322346) q[1];
sx q[1];
rz(-1.511907) q[1];
sx q[1];
rz(1.5712646) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1399978) q[0];
sx q[0];
rz(-2.9712935) q[0];
sx q[0];
rz(-3.1055121) q[0];
rz(1.0958448) q[2];
sx q[2];
rz(-1.6745076) q[2];
sx q[2];
rz(-0.46670676) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.81262368) q[1];
sx q[1];
rz(-1.6996034) q[1];
sx q[1];
rz(1.9616496) q[1];
rz(-0.7223981) q[3];
sx q[3];
rz(-1.5649475) q[3];
sx q[3];
rz(-1.4129782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5181804) q[2];
sx q[2];
rz(-2.0698915) q[2];
sx q[2];
rz(-1.5706459) q[2];
rz(3.0359641) q[3];
sx q[3];
rz(-0.94625866) q[3];
sx q[3];
rz(0.51641974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.327772) q[0];
sx q[0];
rz(-1.0186503) q[0];
sx q[0];
rz(-3.0042197) q[0];
rz(2.9227921) q[1];
sx q[1];
rz(-1.7411502) q[1];
sx q[1];
rz(2.2314744) q[1];
rz(-0.21135052) q[2];
sx q[2];
rz(-2.0049145) q[2];
sx q[2];
rz(2.66578) q[2];
rz(2.4841251) q[3];
sx q[3];
rz(-1.2764662) q[3];
sx q[3];
rz(2.488967) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
