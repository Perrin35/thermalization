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
rz(-0.30682895) q[0];
sx q[0];
rz(0.29875779) q[0];
rz(-0.6660676) q[1];
sx q[1];
rz(2.5112285) q[1];
sx q[1];
rz(11.093333) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.436384) q[0];
sx q[0];
rz(-1.602631) q[0];
sx q[0];
rz(-0.80073204) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5130784) q[2];
sx q[2];
rz(-2.8403663) q[2];
sx q[2];
rz(-2.0026375) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4358878) q[1];
sx q[1];
rz(-1.0748362) q[1];
sx q[1];
rz(-2.3889747) q[1];
x q[2];
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
rz(0.95639688) q[2];
sx q[2];
rz(-2.7608725) q[2];
sx q[2];
rz(-1.8308651) q[2];
rz(-0.091536097) q[3];
sx q[3];
rz(-2.5444578) q[3];
sx q[3];
rz(-0.53102791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0700584) q[0];
sx q[0];
rz(-1.0907084) q[0];
sx q[0];
rz(2.6614905) q[0];
rz(1.1473038) q[1];
sx q[1];
rz(-0.6310178) q[1];
sx q[1];
rz(2.4400585) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6525878) q[0];
sx q[0];
rz(-2.2092144) q[0];
sx q[0];
rz(0.63853635) q[0];
rz(-pi) q[1];
rz(-2.7240997) q[2];
sx q[2];
rz(-1.1402297) q[2];
sx q[2];
rz(-2.9232962) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5614723) q[1];
sx q[1];
rz(-2.0144723) q[1];
sx q[1];
rz(2.2444112) q[1];
rz(-pi) q[2];
rz(-1.2626889) q[3];
sx q[3];
rz(-0.95617056) q[3];
sx q[3];
rz(0.37671396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.39701617) q[2];
sx q[2];
rz(-1.959356) q[2];
sx q[2];
rz(1.5632632) q[2];
rz(-2.8619316) q[3];
sx q[3];
rz(-0.71664387) q[3];
sx q[3];
rz(2.7320812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79769832) q[0];
sx q[0];
rz(-0.41223031) q[0];
sx q[0];
rz(1.4233587) q[0];
rz(-0.1238981) q[1];
sx q[1];
rz(-2.2975477) q[1];
sx q[1];
rz(-1.5843676) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7182962) q[0];
sx q[0];
rz(-0.78351057) q[0];
sx q[0];
rz(1.1447203) q[0];
x q[1];
rz(2.7742375) q[2];
sx q[2];
rz(-0.88161385) q[2];
sx q[2];
rz(-2.286498) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6137984) q[1];
sx q[1];
rz(-1.3247733) q[1];
sx q[1];
rz(3.0600344) q[1];
x q[2];
rz(0.70544542) q[3];
sx q[3];
rz(-2.100832) q[3];
sx q[3];
rz(-2.7646661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.065217) q[2];
sx q[2];
rz(-0.6535483) q[2];
sx q[2];
rz(-2.2067113) q[2];
rz(0.30385083) q[3];
sx q[3];
rz(-0.6128208) q[3];
sx q[3];
rz(2.2936308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3728751) q[0];
sx q[0];
rz(-0.77943742) q[0];
sx q[0];
rz(-0.26707643) q[0];
rz(0.13485394) q[1];
sx q[1];
rz(-0.57661533) q[1];
sx q[1];
rz(2.3042302) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7202111) q[0];
sx q[0];
rz(-1.0307587) q[0];
sx q[0];
rz(-2.8296986) q[0];
rz(-pi) q[1];
rz(-2.5576511) q[2];
sx q[2];
rz(-2.2211233) q[2];
sx q[2];
rz(0.038829858) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.76772753) q[1];
sx q[1];
rz(-0.30555913) q[1];
sx q[1];
rz(2.4453246) q[1];
rz(0.046810026) q[3];
sx q[3];
rz(-1.6322281) q[3];
sx q[3];
rz(0.84247473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.59552938) q[2];
sx q[2];
rz(-0.17543051) q[2];
sx q[2];
rz(-2.9626633) q[2];
rz(1.1692125) q[3];
sx q[3];
rz(-1.365265) q[3];
sx q[3];
rz(0.21807142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1201852) q[0];
sx q[0];
rz(-2.6299801) q[0];
sx q[0];
rz(-2.2688493) q[0];
rz(1.9841638) q[1];
sx q[1];
rz(-0.87481421) q[1];
sx q[1];
rz(-3.1246368) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2568396) q[0];
sx q[0];
rz(-1.5829979) q[0];
sx q[0];
rz(0.73464616) q[0];
rz(-pi) q[1];
rz(-3.1269286) q[2];
sx q[2];
rz(-2.5418106) q[2];
sx q[2];
rz(-1.5408915) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3185721) q[1];
sx q[1];
rz(-1.2766395) q[1];
sx q[1];
rz(-2.3579954) q[1];
x q[2];
rz(-1.8618552) q[3];
sx q[3];
rz(-1.3123056) q[3];
sx q[3];
rz(0.84421009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.056328) q[2];
sx q[2];
rz(-1.4339829) q[2];
sx q[2];
rz(-2.8223574) q[2];
rz(-2.4625835) q[3];
sx q[3];
rz(-1.346799) q[3];
sx q[3];
rz(-0.80176789) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7769258) q[0];
sx q[0];
rz(-0.32061446) q[0];
sx q[0];
rz(0.6426245) q[0];
rz(-1.7721666) q[1];
sx q[1];
rz(-0.6232999) q[1];
sx q[1];
rz(-0.97698897) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8606023) q[0];
sx q[0];
rz(-2.9837768) q[0];
sx q[0];
rz(-0.67276038) q[0];
rz(-0.80987038) q[2];
sx q[2];
rz(-2.4634482) q[2];
sx q[2];
rz(-0.43895753) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.20583892) q[1];
sx q[1];
rz(-0.69200695) q[1];
sx q[1];
rz(1.5924953) q[1];
x q[2];
rz(0.98501916) q[3];
sx q[3];
rz(-1.6156697) q[3];
sx q[3];
rz(-2.6834727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6151108) q[2];
sx q[2];
rz(-1.4075764) q[2];
sx q[2];
rz(1.2283481) q[2];
rz(-2.2743716) q[3];
sx q[3];
rz(-0.4129748) q[3];
sx q[3];
rz(2.373608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41384554) q[0];
sx q[0];
rz(-0.21738805) q[0];
sx q[0];
rz(2.9687498) q[0];
rz(-2.3364283) q[1];
sx q[1];
rz(-2.5925437) q[1];
sx q[1];
rz(0.13596143) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8589218) q[0];
sx q[0];
rz(-1.4513119) q[0];
sx q[0];
rz(2.2764599) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3195761) q[2];
sx q[2];
rz(-1.4944585) q[2];
sx q[2];
rz(2.1427296) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0522642) q[1];
sx q[1];
rz(-1.5517527) q[1];
sx q[1];
rz(0.11472265) q[1];
x q[2];
rz(2.7649859) q[3];
sx q[3];
rz(-2.7889502) q[3];
sx q[3];
rz(-2.6654748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.92676306) q[2];
sx q[2];
rz(-1.5952933) q[2];
sx q[2];
rz(2.9570441) q[2];
rz(0.18937011) q[3];
sx q[3];
rz(-0.53818494) q[3];
sx q[3];
rz(0.93809938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1592584) q[0];
sx q[0];
rz(-0.33894798) q[0];
sx q[0];
rz(2.2151997) q[0];
rz(0.039208086) q[1];
sx q[1];
rz(-0.67752939) q[1];
sx q[1];
rz(2.1047986) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3711972) q[0];
sx q[0];
rz(-1.7465034) q[0];
sx q[0];
rz(1.8699339) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.43208684) q[2];
sx q[2];
rz(-1.4054023) q[2];
sx q[2];
rz(-1.8036606) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8421887) q[1];
sx q[1];
rz(-1.8129947) q[1];
sx q[1];
rz(-2.7224225) q[1];
rz(-pi) q[2];
x q[2];
rz(1.758427) q[3];
sx q[3];
rz(-1.1867282) q[3];
sx q[3];
rz(1.404286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.89256531) q[2];
sx q[2];
rz(-1.1250863) q[2];
sx q[2];
rz(-2.3426775) q[2];
rz(2.6246081) q[3];
sx q[3];
rz(-0.40701443) q[3];
sx q[3];
rz(0.5504722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35140458) q[0];
sx q[0];
rz(-3.1365972) q[0];
sx q[0];
rz(-1.5193526) q[0];
rz(1.5937357) q[1];
sx q[1];
rz(-0.82031119) q[1];
sx q[1];
rz(2.4511852) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.522377) q[0];
sx q[0];
rz(-1.3478312) q[0];
sx q[0];
rz(-2.7821543) q[0];
rz(-pi) q[1];
rz(1.6733273) q[2];
sx q[2];
rz(-2.3674115) q[2];
sx q[2];
rz(0.098086327) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.71296299) q[1];
sx q[1];
rz(-2.3979514) q[1];
sx q[1];
rz(-0.78674591) q[1];
x q[2];
rz(-0.14491187) q[3];
sx q[3];
rz(-1.0365126) q[3];
sx q[3];
rz(-0.89569672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7030299) q[2];
sx q[2];
rz(-2.5327693) q[2];
sx q[2];
rz(-0.969886) q[2];
rz(-2.8841833) q[3];
sx q[3];
rz(-0.22203797) q[3];
sx q[3];
rz(-1.359587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6064706) q[0];
sx q[0];
rz(-0.73000014) q[0];
sx q[0];
rz(-3.0560793) q[0];
rz(-2.8657148) q[1];
sx q[1];
rz(-0.62092263) q[1];
sx q[1];
rz(1.9524908) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10201926) q[0];
sx q[0];
rz(-1.1652526) q[0];
sx q[0];
rz(-1.1723551) q[0];
rz(2.8315808) q[2];
sx q[2];
rz(-1.5793206) q[2];
sx q[2];
rz(-0.046176813) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.738742) q[1];
sx q[1];
rz(-1.9704809) q[1];
sx q[1];
rz(-3.012639) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6639878) q[3];
sx q[3];
rz(-1.9177327) q[3];
sx q[3];
rz(-0.047752927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6659866) q[2];
sx q[2];
rz(-2.2624272) q[2];
sx q[2];
rz(-0.48412588) q[2];
rz(-0.13949805) q[3];
sx q[3];
rz(-2.3660584) q[3];
sx q[3];
rz(2.9805396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6763247) q[0];
sx q[0];
rz(-1.9857255) q[0];
sx q[0];
rz(2.3406512) q[0];
rz(1.6839266) q[1];
sx q[1];
rz(-1.4977581) q[1];
sx q[1];
rz(-1.3118634) q[1];
rz(-1.3670078) q[2];
sx q[2];
rz(-0.93175722) q[2];
sx q[2];
rz(2.5358806) q[2];
rz(2.3845354) q[3];
sx q[3];
rz(-1.591605) q[3];
sx q[3];
rz(0.021416728) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
