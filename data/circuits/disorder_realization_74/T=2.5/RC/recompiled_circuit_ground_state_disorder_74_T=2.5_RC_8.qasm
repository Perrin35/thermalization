OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.35535204) q[0];
sx q[0];
rz(2.8347637) q[0];
sx q[0];
rz(9.1260202) q[0];
rz(-0.6660676) q[1];
sx q[1];
rz(2.5112285) q[1];
sx q[1];
rz(11.093333) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10160343) q[0];
sx q[0];
rz(-2.371006) q[0];
sx q[0];
rz(1.5250852) q[0];
x q[1];
rz(-1.2700419) q[2];
sx q[2];
rz(-1.5879121) q[2];
sx q[2];
rz(2.764876) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4358878) q[1];
sx q[1];
rz(-2.0667564) q[1];
sx q[1];
rz(2.3889747) q[1];
rz(1.6410337) q[3];
sx q[3];
rz(-1.8459709) q[3];
sx q[3];
rz(0.17632139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1851958) q[2];
sx q[2];
rz(-2.7608725) q[2];
sx q[2];
rz(1.8308651) q[2];
rz(3.0500566) q[3];
sx q[3];
rz(-0.59713489) q[3];
sx q[3];
rz(-2.6105647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.071534261) q[0];
sx q[0];
rz(-2.0508843) q[0];
sx q[0];
rz(0.48010215) q[0];
rz(1.1473038) q[1];
sx q[1];
rz(-2.5105748) q[1];
sx q[1];
rz(0.70153418) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40532058) q[0];
sx q[0];
rz(-0.87001409) q[0];
sx q[0];
rz(2.2474656) q[0];
rz(-pi) q[1];
rz(0.84757324) q[2];
sx q[2];
rz(-0.59047359) q[2];
sx q[2];
rz(2.5440885) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.638826) q[1];
sx q[1];
rz(-2.354489) q[1];
sx q[1];
rz(0.9197286) q[1];
x q[2];
rz(1.2626889) q[3];
sx q[3];
rz(-2.1854221) q[3];
sx q[3];
rz(-2.7648787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.39701617) q[2];
sx q[2];
rz(-1.1822367) q[2];
sx q[2];
rz(-1.5783295) q[2];
rz(-2.8619316) q[3];
sx q[3];
rz(-0.71664387) q[3];
sx q[3];
rz(-0.40951148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.79769832) q[0];
sx q[0];
rz(-0.41223031) q[0];
sx q[0];
rz(1.718234) q[0];
rz(-3.0176945) q[1];
sx q[1];
rz(-2.2975477) q[1];
sx q[1];
rz(-1.557225) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9950969) q[0];
sx q[0];
rz(-2.2687758) q[0];
sx q[0];
rz(2.7510038) q[0];
rz(-2.7742375) q[2];
sx q[2];
rz(-2.2599788) q[2];
sx q[2];
rz(-2.286498) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8516006) q[1];
sx q[1];
rz(-0.25892913) q[1];
sx q[1];
rz(-1.2570791) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.91482678) q[3];
sx q[3];
rz(-0.97712356) q[3];
sx q[3];
rz(2.3542924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.065217) q[2];
sx q[2];
rz(-0.6535483) q[2];
sx q[2];
rz(-2.2067113) q[2];
rz(0.30385083) q[3];
sx q[3];
rz(-0.6128208) q[3];
sx q[3];
rz(-0.8479619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(2.3728751) q[0];
sx q[0];
rz(-0.77943742) q[0];
sx q[0];
rz(0.26707643) q[0];
rz(0.13485394) q[1];
sx q[1];
rz(-0.57661533) q[1];
sx q[1];
rz(-0.83736247) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42138153) q[0];
sx q[0];
rz(-2.1108339) q[0];
sx q[0];
rz(-2.8296986) q[0];
rz(2.3100987) q[2];
sx q[2];
rz(-1.1165501) q[2];
sx q[2];
rz(1.9902094) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4761047) q[1];
sx q[1];
rz(-1.7649516) q[1];
sx q[1];
rz(0.23745115) q[1];
x q[2];
rz(0.92045201) q[3];
sx q[3];
rz(-0.077215791) q[3];
sx q[3];
rz(1.6473351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.59552938) q[2];
sx q[2];
rz(-0.17543051) q[2];
sx q[2];
rz(-2.9626633) q[2];
rz(1.9723802) q[3];
sx q[3];
rz(-1.365265) q[3];
sx q[3];
rz(-0.21807142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1201852) q[0];
sx q[0];
rz(-2.6299801) q[0];
sx q[0];
rz(-2.2688493) q[0];
rz(-1.9841638) q[1];
sx q[1];
rz(-2.2667784) q[1];
sx q[1];
rz(-3.1246368) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32497787) q[0];
sx q[0];
rz(-2.3053753) q[0];
sx q[0];
rz(1.5543544) q[0];
x q[1];
rz(1.5808231) q[2];
sx q[2];
rz(-2.1705049) q[2];
sx q[2];
rz(-1.6184652) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.53345799) q[1];
sx q[1];
rz(-0.82913405) q[1];
sx q[1];
rz(-1.9749179) q[1];
rz(-pi) q[2];
rz(2.3153002) q[3];
sx q[3];
rz(-0.38681627) q[3];
sx q[3];
rz(1.4329965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.056328) q[2];
sx q[2];
rz(-1.7076098) q[2];
sx q[2];
rz(-2.8223574) q[2];
rz(-0.6790092) q[3];
sx q[3];
rz(-1.7947936) q[3];
sx q[3];
rz(2.3398248) q[3];
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
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3646669) q[0];
sx q[0];
rz(-0.32061446) q[0];
sx q[0];
rz(2.4989682) q[0];
rz(-1.7721666) q[1];
sx q[1];
rz(-0.6232999) q[1];
sx q[1];
rz(-0.97698897) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60212712) q[0];
sx q[0];
rz(-1.4475679) q[0];
sx q[0];
rz(-1.4719523) q[0];
rz(0.80987038) q[2];
sx q[2];
rz(-2.4634482) q[2];
sx q[2];
rz(-2.7026351) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9639324) q[1];
sx q[1];
rz(-2.2626082) q[1];
sx q[1];
rz(3.1236137) q[1];
rz(-pi) q[2];
rz(-2.1565735) q[3];
sx q[3];
rz(-1.5259229) q[3];
sx q[3];
rz(2.6834727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6151108) q[2];
sx q[2];
rz(-1.4075764) q[2];
sx q[2];
rz(-1.2283481) q[2];
rz(2.2743716) q[3];
sx q[3];
rz(-0.4129748) q[3];
sx q[3];
rz(0.76798463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41384554) q[0];
sx q[0];
rz(-0.21738805) q[0];
sx q[0];
rz(-0.17284285) q[0];
rz(0.8051644) q[1];
sx q[1];
rz(-0.54904896) q[1];
sx q[1];
rz(-0.13596143) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14912389) q[0];
sx q[0];
rz(-0.71399105) q[0];
sx q[0];
rz(-1.7538422) q[0];
rz(-0.82201652) q[2];
sx q[2];
rz(-1.4944585) q[2];
sx q[2];
rz(-2.1427296) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0522642) q[1];
sx q[1];
rz(-1.58984) q[1];
sx q[1];
rz(0.11472265) q[1];
rz(-pi) q[2];
rz(-0.32973955) q[3];
sx q[3];
rz(-1.4434333) q[3];
sx q[3];
rz(1.4500666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.92676306) q[2];
sx q[2];
rz(-1.5462993) q[2];
sx q[2];
rz(2.9570441) q[2];
rz(-0.18937011) q[3];
sx q[3];
rz(-0.53818494) q[3];
sx q[3];
rz(2.2034933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.9823343) q[0];
sx q[0];
rz(-2.8026447) q[0];
sx q[0];
rz(0.92639297) q[0];
rz(0.039208086) q[1];
sx q[1];
rz(-0.67752939) q[1];
sx q[1];
rz(2.1047986) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77039546) q[0];
sx q[0];
rz(-1.3950893) q[0];
sx q[0];
rz(1.8699339) q[0];
x q[1];
rz(-2.7095058) q[2];
sx q[2];
rz(-1.7361904) q[2];
sx q[2];
rz(-1.8036606) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1649303) q[1];
sx q[1];
rz(-1.9769985) q[1];
sx q[1];
rz(1.3066533) q[1];
rz(-pi) q[2];
rz(2.7091712) q[3];
sx q[3];
rz(-0.42539551) q[3];
sx q[3];
rz(-2.2063125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2490273) q[2];
sx q[2];
rz(-2.0165063) q[2];
sx q[2];
rz(0.79891515) q[2];
rz(-0.51698452) q[3];
sx q[3];
rz(-0.40701443) q[3];
sx q[3];
rz(0.5504722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7901881) q[0];
sx q[0];
rz(-0.0049954448) q[0];
sx q[0];
rz(-1.6222401) q[0];
rz(-1.5478569) q[1];
sx q[1];
rz(-2.3212815) q[1];
sx q[1];
rz(-2.4511852) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6192157) q[0];
sx q[0];
rz(-1.3478312) q[0];
sx q[0];
rz(-2.7821543) q[0];
rz(-pi) q[1];
rz(0.099748513) q[2];
sx q[2];
rz(-2.3398551) q[2];
sx q[2];
rz(-0.24101098) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.71296299) q[1];
sx q[1];
rz(-0.74364122) q[1];
sx q[1];
rz(0.78674591) q[1];
rz(2.1097095) q[3];
sx q[3];
rz(-1.6953985) q[3];
sx q[3];
rz(2.3923158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7030299) q[2];
sx q[2];
rz(-2.5327693) q[2];
sx q[2];
rz(0.969886) q[2];
rz(-0.25740933) q[3];
sx q[3];
rz(-0.22203797) q[3];
sx q[3];
rz(-1.7820057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6064706) q[0];
sx q[0];
rz(-2.4115925) q[0];
sx q[0];
rz(0.085513376) q[0];
rz(-2.8657148) q[1];
sx q[1];
rz(-2.52067) q[1];
sx q[1];
rz(1.1891018) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5082423) q[0];
sx q[0];
rz(-1.2062643) q[0];
sx q[0];
rz(-2.7056498) q[0];
x q[1];
rz(-3.1136572) q[2];
sx q[2];
rz(-2.8314674) q[2];
sx q[2];
rz(-1.5512229) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.72451776) q[1];
sx q[1];
rz(-2.7226884) q[1];
sx q[1];
rz(-1.275283) q[1];
rz(2.7932634) q[3];
sx q[3];
rz(-1.4831721) q[3];
sx q[3];
rz(-1.5548116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6659866) q[2];
sx q[2];
rz(-2.2624272) q[2];
sx q[2];
rz(-2.6574668) q[2];
rz(-3.0020946) q[3];
sx q[3];
rz(-0.77553427) q[3];
sx q[3];
rz(2.9805396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.465268) q[0];
sx q[0];
rz(-1.9857255) q[0];
sx q[0];
rz(2.3406512) q[0];
rz(1.4576661) q[1];
sx q[1];
rz(-1.6438345) q[1];
sx q[1];
rz(1.8297292) q[1];
rz(-1.3670078) q[2];
sx q[2];
rz(-0.93175722) q[2];
sx q[2];
rz(2.5358806) q[2];
rz(-2.3845354) q[3];
sx q[3];
rz(-1.5499877) q[3];
sx q[3];
rz(-3.1201759) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
