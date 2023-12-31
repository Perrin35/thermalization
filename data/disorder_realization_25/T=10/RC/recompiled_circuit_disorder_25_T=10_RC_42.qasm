OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0796354) q[0];
sx q[0];
rz(-2.893653) q[0];
sx q[0];
rz(-2.3214582) q[0];
rz(0.040854383) q[1];
sx q[1];
rz(-2.3526469) q[1];
sx q[1];
rz(0.087892428) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1577507) q[0];
sx q[0];
rz(-2.0294667) q[0];
sx q[0];
rz(1.0755324) q[0];
x q[1];
rz(2.7152039) q[2];
sx q[2];
rz(-1.083056) q[2];
sx q[2];
rz(0.79977712) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.81401285) q[1];
sx q[1];
rz(-2.6373133) q[1];
sx q[1];
rz(-3.1370647) q[1];
x q[2];
rz(-1.5486693) q[3];
sx q[3];
rz(-0.25875729) q[3];
sx q[3];
rz(0.72507897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0007881) q[2];
sx q[2];
rz(-1.452383) q[2];
sx q[2];
rz(0.56935707) q[2];
rz(-1.5287483) q[3];
sx q[3];
rz(-0.6362392) q[3];
sx q[3];
rz(1.7830085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.59812087) q[0];
sx q[0];
rz(-2.9765029) q[0];
sx q[0];
rz(-0.55364451) q[0];
rz(-1.2373295) q[1];
sx q[1];
rz(-1.3668704) q[1];
sx q[1];
rz(-1.9083317) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90399088) q[0];
sx q[0];
rz(-2.1008137) q[0];
sx q[0];
rz(-2.3474098) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4670062) q[2];
sx q[2];
rz(-1.470675) q[2];
sx q[2];
rz(1.9763725) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1439647) q[1];
sx q[1];
rz(-1.7768754) q[1];
sx q[1];
rz(-0.016184316) q[1];
rz(-pi) q[2];
x q[2];
rz(0.66744653) q[3];
sx q[3];
rz(-1.7597886) q[3];
sx q[3];
rz(-0.28039704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.0040940293) q[2];
sx q[2];
rz(-1.6823021) q[2];
sx q[2];
rz(-0.0022350524) q[2];
rz(0.8301174) q[3];
sx q[3];
rz(-2.3486962) q[3];
sx q[3];
rz(1.2494276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8925791) q[0];
sx q[0];
rz(-2.5254288) q[0];
sx q[0];
rz(2.2900443) q[0];
rz(0.7710723) q[1];
sx q[1];
rz(-1.4665736) q[1];
sx q[1];
rz(-3.070389) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6547346) q[0];
sx q[0];
rz(-1.808597) q[0];
sx q[0];
rz(0.20240692) q[0];
rz(1.0686915) q[2];
sx q[2];
rz(-0.50902589) q[2];
sx q[2];
rz(0.44391649) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5350266) q[1];
sx q[1];
rz(-0.92381322) q[1];
sx q[1];
rz(1.4024629) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86665385) q[3];
sx q[3];
rz(-0.94368499) q[3];
sx q[3];
rz(-1.8303527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3123793) q[2];
sx q[2];
rz(-2.2177314) q[2];
sx q[2];
rz(-0.52345792) q[2];
rz(0.33189014) q[3];
sx q[3];
rz(-3.0811946) q[3];
sx q[3];
rz(-2.6383242) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7317384) q[0];
sx q[0];
rz(-0.24895746) q[0];
sx q[0];
rz(2.3160034) q[0];
rz(-1.9748953) q[1];
sx q[1];
rz(-0.86667934) q[1];
sx q[1];
rz(-1.694214) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3281876) q[0];
sx q[0];
rz(-1.4077328) q[0];
sx q[0];
rz(1.8036611) q[0];
rz(-1.6214192) q[2];
sx q[2];
rz(-2.9630337) q[2];
sx q[2];
rz(0.52390097) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.96781603) q[1];
sx q[1];
rz(-2.2602091) q[1];
sx q[1];
rz(0.20676989) q[1];
x q[2];
rz(-1.8131687) q[3];
sx q[3];
rz(-1.216785) q[3];
sx q[3];
rz(-2.6995475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6772785) q[2];
sx q[2];
rz(-1.7769287) q[2];
sx q[2];
rz(0.14492598) q[2];
rz(1.011301) q[3];
sx q[3];
rz(-0.62711182) q[3];
sx q[3];
rz(-3.1269126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99575627) q[0];
sx q[0];
rz(-0.053481426) q[0];
sx q[0];
rz(-0.68403912) q[0];
rz(1.1902635) q[1];
sx q[1];
rz(-2.4763156) q[1];
sx q[1];
rz(1.9285944) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0048464674) q[0];
sx q[0];
rz(-2.6803225) q[0];
sx q[0];
rz(2.524551) q[0];
x q[1];
rz(0.62972516) q[2];
sx q[2];
rz(-1.0169528) q[2];
sx q[2];
rz(0.50309203) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4124406) q[1];
sx q[1];
rz(-1.8262595) q[1];
sx q[1];
rz(1.4490119) q[1];
rz(-pi) q[2];
rz(0.17739399) q[3];
sx q[3];
rz(-1.3131485) q[3];
sx q[3];
rz(2.859476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.50903901) q[2];
sx q[2];
rz(-1.2712487) q[2];
sx q[2];
rz(2.5115013) q[2];
rz(-0.57224327) q[3];
sx q[3];
rz(-2.4961491) q[3];
sx q[3];
rz(2.2563289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0260139) q[0];
sx q[0];
rz(-0.84340874) q[0];
sx q[0];
rz(2.8552326) q[0];
rz(0.59965602) q[1];
sx q[1];
rz(-1.1279794) q[1];
sx q[1];
rz(2.5568331) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50305784) q[0];
sx q[0];
rz(-2.1612744) q[0];
sx q[0];
rz(-1.9496586) q[0];
rz(-pi) q[1];
rz(-2.5271687) q[2];
sx q[2];
rz(-1.5493869) q[2];
sx q[2];
rz(-1.1504088) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6905745) q[1];
sx q[1];
rz(-1.283053) q[1];
sx q[1];
rz(-2.2698127) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9752475) q[3];
sx q[3];
rz(-1.5208941) q[3];
sx q[3];
rz(2.9970616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.16453234) q[2];
sx q[2];
rz(-0.88447276) q[2];
sx q[2];
rz(1.5131081) q[2];
rz(0.96308723) q[3];
sx q[3];
rz(-1.77308) q[3];
sx q[3];
rz(-2.517038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.088783711) q[0];
sx q[0];
rz(-1.4316906) q[0];
sx q[0];
rz(-2.5928296) q[0];
rz(-0.051963003) q[1];
sx q[1];
rz(-2.6497662) q[1];
sx q[1];
rz(-2.5351977) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.711395) q[0];
sx q[0];
rz(-2.093962) q[0];
sx q[0];
rz(2.8682312) q[0];
rz(-pi) q[1];
x q[1];
rz(0.37574558) q[2];
sx q[2];
rz(-0.94259113) q[2];
sx q[2];
rz(1.9919765) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4013306) q[1];
sx q[1];
rz(-1.4583424) q[1];
sx q[1];
rz(2.3349891) q[1];
rz(-pi) q[2];
rz(-0.86319478) q[3];
sx q[3];
rz(-0.51897012) q[3];
sx q[3];
rz(0.88245813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2018044) q[2];
sx q[2];
rz(-0.99210343) q[2];
sx q[2];
rz(-0.50841224) q[2];
rz(-2.0243747) q[3];
sx q[3];
rz(-0.17337392) q[3];
sx q[3];
rz(-1.4846876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8023119) q[0];
sx q[0];
rz(-0.34582129) q[0];
sx q[0];
rz(-2.3876277) q[0];
rz(2.1915961) q[1];
sx q[1];
rz(-1.4818622) q[1];
sx q[1];
rz(2.9139013) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4906625) q[0];
sx q[0];
rz(-0.41914808) q[0];
sx q[0];
rz(-1.7463669) q[0];
x q[1];
rz(-1.20485) q[2];
sx q[2];
rz(-1.9908675) q[2];
sx q[2];
rz(2.7811188) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.55351174) q[1];
sx q[1];
rz(-1.3462102) q[1];
sx q[1];
rz(2.5976546) q[1];
rz(0.065683059) q[3];
sx q[3];
rz(-0.55995299) q[3];
sx q[3];
rz(0.63510676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.078538744) q[2];
sx q[2];
rz(-0.20919122) q[2];
sx q[2];
rz(0.83855808) q[2];
rz(-2.7770384) q[3];
sx q[3];
rz(-1.7611327) q[3];
sx q[3];
rz(2.3013039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39695981) q[0];
sx q[0];
rz(-1.7663706) q[0];
sx q[0];
rz(-2.3378085) q[0];
rz(-2.0869758) q[1];
sx q[1];
rz(-0.63916731) q[1];
sx q[1];
rz(-1.9931591) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52590695) q[0];
sx q[0];
rz(-1.6427759) q[0];
sx q[0];
rz(0.023948897) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7788413) q[2];
sx q[2];
rz(-1.0748378) q[2];
sx q[2];
rz(-0.19778684) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8879537) q[1];
sx q[1];
rz(-1.1432853) q[1];
sx q[1];
rz(1.0484139) q[1];
x q[2];
rz(-2.7896342) q[3];
sx q[3];
rz(-1.5378012) q[3];
sx q[3];
rz(-1.9381423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0917197) q[2];
sx q[2];
rz(-2.1240978) q[2];
sx q[2];
rz(-0.77825528) q[2];
rz(0.19566472) q[3];
sx q[3];
rz(-2.1019432) q[3];
sx q[3];
rz(-1.0218609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2739928) q[0];
sx q[0];
rz(-2.1430528) q[0];
sx q[0];
rz(2.8883949) q[0];
rz(-2.4018438) q[1];
sx q[1];
rz(-0.7363798) q[1];
sx q[1];
rz(-2.443312) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2343276) q[0];
sx q[0];
rz(-2.8303879) q[0];
sx q[0];
rz(-0.23647232) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.024784879) q[2];
sx q[2];
rz(-1.6449271) q[2];
sx q[2];
rz(0.63400808) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0776978) q[1];
sx q[1];
rz(-2.6297036) q[1];
sx q[1];
rz(2.0236532) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9000557) q[3];
sx q[3];
rz(-1.0165443) q[3];
sx q[3];
rz(2.9709771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0239821) q[2];
sx q[2];
rz(-2.2587946) q[2];
sx q[2];
rz(2.3013766) q[2];
rz(-2.501287) q[3];
sx q[3];
rz(-0.87602031) q[3];
sx q[3];
rz(-1.1673814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.8828076) q[0];
sx q[0];
rz(-2.2866645) q[0];
sx q[0];
rz(-1.4186161) q[0];
rz(-0.029126833) q[1];
sx q[1];
rz(-3.0283785) q[1];
sx q[1];
rz(1.821847) q[1];
rz(2.6559033) q[2];
sx q[2];
rz(-2.4808241) q[2];
sx q[2];
rz(-0.11447699) q[2];
rz(-2.2610353) q[3];
sx q[3];
rz(-1.6869443) q[3];
sx q[3];
rz(-1.501367) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
