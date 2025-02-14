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
rz(0.49993604) q[0];
sx q[0];
rz(1.9498107) q[0];
sx q[0];
rz(10.796588) q[0];
rz(1.6504047) q[1];
sx q[1];
rz(-0.86725441) q[1];
sx q[1];
rz(0.39630085) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9797452) q[0];
sx q[0];
rz(-1.7801741) q[0];
sx q[0];
rz(-1.4896979) q[0];
rz(-pi) q[1];
rz(-1.9849586) q[2];
sx q[2];
rz(-2.7348619) q[2];
sx q[2];
rz(1.0309645) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.635879) q[1];
sx q[1];
rz(-1.6649655) q[1];
sx q[1];
rz(-1.0416743) q[1];
rz(-pi) q[2];
rz(2.5175573) q[3];
sx q[3];
rz(-2.6665832) q[3];
sx q[3];
rz(1.6866682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0013915) q[2];
sx q[2];
rz(-0.80696693) q[2];
sx q[2];
rz(1.8001455) q[2];
rz(-1.0691102) q[3];
sx q[3];
rz(-0.1461229) q[3];
sx q[3];
rz(1.7599546) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57217252) q[0];
sx q[0];
rz(-2.2241346) q[0];
sx q[0];
rz(-0.65895748) q[0];
rz(0.072703687) q[1];
sx q[1];
rz(-2.0854988) q[1];
sx q[1];
rz(1.2603849) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2598686) q[0];
sx q[0];
rz(-1.4894852) q[0];
sx q[0];
rz(-2.0242244) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0837696) q[2];
sx q[2];
rz(-1.1758846) q[2];
sx q[2];
rz(2.9213443) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1998493) q[1];
sx q[1];
rz(-1.3056271) q[1];
sx q[1];
rz(-1.6380235) q[1];
rz(-2.1815592) q[3];
sx q[3];
rz(-1.472578) q[3];
sx q[3];
rz(-1.4032422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.8947123) q[2];
sx q[2];
rz(-2.3983045) q[2];
sx q[2];
rz(1.7903719) q[2];
rz(2.5249935) q[3];
sx q[3];
rz(-1.0298046) q[3];
sx q[3];
rz(2.8504347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5093444) q[0];
sx q[0];
rz(-0.46724874) q[0];
sx q[0];
rz(0.21009357) q[0];
rz(-3.0585152) q[1];
sx q[1];
rz(-1.2385085) q[1];
sx q[1];
rz(-0.067213623) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34850304) q[0];
sx q[0];
rz(-1.5572892) q[0];
sx q[0];
rz(-1.5828787) q[0];
rz(-2.1234197) q[2];
sx q[2];
rz(-0.45028307) q[2];
sx q[2];
rz(2.1076815) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1695425) q[1];
sx q[1];
rz(-1.2360555) q[1];
sx q[1];
rz(-2.2856345) q[1];
rz(-2.4359853) q[3];
sx q[3];
rz(-1.5021245) q[3];
sx q[3];
rz(-2.9182383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.82789603) q[2];
sx q[2];
rz(-1.4915497) q[2];
sx q[2];
rz(-1.776604) q[2];
rz(2.7374173) q[3];
sx q[3];
rz(-1.3173236) q[3];
sx q[3];
rz(1.9425758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2760524) q[0];
sx q[0];
rz(-1.4006389) q[0];
sx q[0];
rz(-0.47819594) q[0];
rz(-2.1887691) q[1];
sx q[1];
rz(-0.15004221) q[1];
sx q[1];
rz(2.731954) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7573625) q[0];
sx q[0];
rz(-2.2262947) q[0];
sx q[0];
rz(1.537498) q[0];
rz(-pi) q[1];
x q[1];
rz(0.45646052) q[2];
sx q[2];
rz(-1.4829516) q[2];
sx q[2];
rz(-2.8473471) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8781902) q[1];
sx q[1];
rz(-2.1904084) q[1];
sx q[1];
rz(-2.3248169) q[1];
rz(-pi) q[2];
rz(-1.8919935) q[3];
sx q[3];
rz(-1.7068591) q[3];
sx q[3];
rz(0.67067671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.59565583) q[2];
sx q[2];
rz(-1.7031534) q[2];
sx q[2];
rz(-1.1379918) q[2];
rz(0.99500895) q[3];
sx q[3];
rz(-0.55750877) q[3];
sx q[3];
rz(3.0549684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74395162) q[0];
sx q[0];
rz(-1.8331563) q[0];
sx q[0];
rz(-2.8635136) q[0];
rz(1.2644348) q[1];
sx q[1];
rz(-1.2828628) q[1];
sx q[1];
rz(1.5778731) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0456384) q[0];
sx q[0];
rz(-1.6976893) q[0];
sx q[0];
rz(-1.9357827) q[0];
x q[1];
rz(0.35985882) q[2];
sx q[2];
rz(-0.82170502) q[2];
sx q[2];
rz(2.6287959) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.63999004) q[1];
sx q[1];
rz(-1.4516561) q[1];
sx q[1];
rz(0.57699012) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4117354) q[3];
sx q[3];
rz(-1.4791227) q[3];
sx q[3];
rz(-1.5605469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9331253) q[2];
sx q[2];
rz(-1.5310023) q[2];
sx q[2];
rz(2.7739286) q[2];
rz(2.0465046) q[3];
sx q[3];
rz(-2.118066) q[3];
sx q[3];
rz(2.9663405) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6943738) q[0];
sx q[0];
rz(-0.29964888) q[0];
sx q[0];
rz(-1.3391986) q[0];
rz(0.12204349) q[1];
sx q[1];
rz(-1.3587147) q[1];
sx q[1];
rz(2.7809714) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0568389) q[0];
sx q[0];
rz(-0.90843102) q[0];
sx q[0];
rz(-1.0785036) q[0];
rz(0.039741411) q[2];
sx q[2];
rz(-1.1510599) q[2];
sx q[2];
rz(2.6750203) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.569446) q[1];
sx q[1];
rz(-2.0098369) q[1];
sx q[1];
rz(2.8296058) q[1];
rz(0.020645647) q[3];
sx q[3];
rz(-0.25293487) q[3];
sx q[3];
rz(1.3731352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.7755166) q[2];
sx q[2];
rz(-1.4630432) q[2];
sx q[2];
rz(2.0029081) q[2];
rz(0.25389296) q[3];
sx q[3];
rz(-2.4963899) q[3];
sx q[3];
rz(-2.3937288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.89895407) q[0];
sx q[0];
rz(-2.2518318) q[0];
sx q[0];
rz(-2.1658072) q[0];
rz(2.0121393) q[1];
sx q[1];
rz(-0.44810805) q[1];
sx q[1];
rz(-1.7879558) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0746552) q[0];
sx q[0];
rz(-0.88464175) q[0];
sx q[0];
rz(-1.1891436) q[0];
rz(-pi) q[1];
rz(3.0616272) q[2];
sx q[2];
rz(-1.1478394) q[2];
sx q[2];
rz(-2.1048196) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2372386) q[1];
sx q[1];
rz(-1.2356346) q[1];
sx q[1];
rz(-1.4314913) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4134992) q[3];
sx q[3];
rz(-0.59641664) q[3];
sx q[3];
rz(0.28486681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5740616) q[2];
sx q[2];
rz(-0.83063829) q[2];
sx q[2];
rz(0.91721025) q[2];
rz(-1.5926825) q[3];
sx q[3];
rz(-0.33532381) q[3];
sx q[3];
rz(2.3622021) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9354644) q[0];
sx q[0];
rz(-1.520282) q[0];
sx q[0];
rz(-0.024854831) q[0];
rz(1.286233) q[1];
sx q[1];
rz(-1.7846466) q[1];
sx q[1];
rz(1.6110427) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0162189) q[0];
sx q[0];
rz(-1.1825996) q[0];
sx q[0];
rz(2.6632705) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31168657) q[2];
sx q[2];
rz(-2.4309845) q[2];
sx q[2];
rz(-2.54541) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3250585) q[1];
sx q[1];
rz(-1.5360502) q[1];
sx q[1];
rz(-2.860022) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1389109) q[3];
sx q[3];
rz(-1.362097) q[3];
sx q[3];
rz(-0.025996093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.020869104) q[2];
sx q[2];
rz(-2.1376231) q[2];
sx q[2];
rz(-1.2845767) q[2];
rz(-2.8203216) q[3];
sx q[3];
rz(-0.64906859) q[3];
sx q[3];
rz(-2.9200714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.016841737) q[0];
sx q[0];
rz(-1.1469954) q[0];
sx q[0];
rz(2.2220213) q[0];
rz(0.59182566) q[1];
sx q[1];
rz(-2.070919) q[1];
sx q[1];
rz(0.33669534) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3580456) q[0];
sx q[0];
rz(-1.3033195) q[0];
sx q[0];
rz(0.13937431) q[0];
rz(-pi) q[1];
rz(-0.36328237) q[2];
sx q[2];
rz(-1.6840625) q[2];
sx q[2];
rz(-0.47081468) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0398061) q[1];
sx q[1];
rz(-2.617604) q[1];
sx q[1];
rz(-2.8807959) q[1];
x q[2];
rz(2.982364) q[3];
sx q[3];
rz(-1.8223636) q[3];
sx q[3];
rz(0.23345527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0261953) q[2];
sx q[2];
rz(-0.43033174) q[2];
sx q[2];
rz(2.3587312) q[2];
rz(3.03249) q[3];
sx q[3];
rz(-1.0166054) q[3];
sx q[3];
rz(1.2469863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1054909) q[0];
sx q[0];
rz(-1.6721268) q[0];
sx q[0];
rz(-2.2176657) q[0];
rz(0.52664122) q[1];
sx q[1];
rz(-1.8662607) q[1];
sx q[1];
rz(0.99871666) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.963388) q[0];
sx q[0];
rz(-2.9792157) q[0];
sx q[0];
rz(2.8078662) q[0];
rz(1.8615103) q[2];
sx q[2];
rz(-1.4913342) q[2];
sx q[2];
rz(3.0018011) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0622039) q[1];
sx q[1];
rz(-2.1583589) q[1];
sx q[1];
rz(3.1000137) q[1];
rz(-pi) q[2];
rz(2.4630354) q[3];
sx q[3];
rz(-1.863409) q[3];
sx q[3];
rz(-1.7998526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0052789) q[2];
sx q[2];
rz(-0.88999358) q[2];
sx q[2];
rz(2.8078553) q[2];
rz(-0.8606832) q[3];
sx q[3];
rz(-1.656683) q[3];
sx q[3];
rz(0.0066283289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.6776047) q[0];
sx q[0];
rz(-1.1900359) q[0];
sx q[0];
rz(0.66587454) q[0];
rz(2.5869276) q[1];
sx q[1];
rz(-1.278109) q[1];
sx q[1];
rz(0.14229933) q[1];
rz(-1.3210422) q[2];
sx q[2];
rz(-1.142923) q[2];
sx q[2];
rz(-0.10071071) q[2];
rz(2.4246115) q[3];
sx q[3];
rz(-0.76631279) q[3];
sx q[3];
rz(-2.8240639) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
