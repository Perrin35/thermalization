OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6157827) q[0];
sx q[0];
rz(-1.4178185) q[0];
sx q[0];
rz(-2.5807227) q[0];
rz(1.1129192) q[1];
sx q[1];
rz(-1.7634044) q[1];
sx q[1];
rz(1.2150432) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62984798) q[0];
sx q[0];
rz(-1.6748322) q[0];
sx q[0];
rz(-1.3826136) q[0];
x q[1];
rz(1.2826074) q[2];
sx q[2];
rz(-2.211314) q[2];
sx q[2];
rz(0.033601947) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3455968) q[1];
sx q[1];
rz(-1.0618292) q[1];
sx q[1];
rz(0.3791581) q[1];
x q[2];
rz(2.9925572) q[3];
sx q[3];
rz(-1.8334853) q[3];
sx q[3];
rz(1.5460154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3540196) q[2];
sx q[2];
rz(-0.95280567) q[2];
sx q[2];
rz(-0.18307486) q[2];
rz(2.7637774) q[3];
sx q[3];
rz(-2.0928045) q[3];
sx q[3];
rz(-0.29418501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-0.29782444) q[0];
sx q[0];
rz(-2.4968708) q[0];
sx q[0];
rz(0.077117292) q[0];
rz(-0.33879694) q[1];
sx q[1];
rz(-2.0270551) q[1];
sx q[1];
rz(1.5391301) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2167643) q[0];
sx q[0];
rz(-2.1568858) q[0];
sx q[0];
rz(3.0490962) q[0];
x q[1];
rz(0.23761959) q[2];
sx q[2];
rz(-2.3887206) q[2];
sx q[2];
rz(-1.2653637) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.96670818) q[1];
sx q[1];
rz(-2.0497353) q[1];
sx q[1];
rz(2.944988) q[1];
rz(-2.7178571) q[3];
sx q[3];
rz(-1.3556004) q[3];
sx q[3];
rz(-2.7620897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.845528) q[2];
sx q[2];
rz(-1.8639996) q[2];
sx q[2];
rz(-2.4831333) q[2];
rz(-2.9902839) q[3];
sx q[3];
rz(-2.1189809) q[3];
sx q[3];
rz(-2.4466799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028458683) q[0];
sx q[0];
rz(-0.78335339) q[0];
sx q[0];
rz(-2.7084896) q[0];
rz(-1.1921047) q[1];
sx q[1];
rz(-1.9299709) q[1];
sx q[1];
rz(2.5862397) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1486737) q[0];
sx q[0];
rz(-2.2745471) q[0];
sx q[0];
rz(-0.91627319) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0981512) q[2];
sx q[2];
rz(-0.61908365) q[2];
sx q[2];
rz(1.667779) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1173646) q[1];
sx q[1];
rz(-1.7324565) q[1];
sx q[1];
rz(2.2092186) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.20136307) q[3];
sx q[3];
rz(-1.0346197) q[3];
sx q[3];
rz(-0.45504967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.73734036) q[2];
sx q[2];
rz(-2.3534687) q[2];
sx q[2];
rz(1.2505442) q[2];
rz(-2.897443) q[3];
sx q[3];
rz(-1.282225) q[3];
sx q[3];
rz(1.4499433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26043949) q[0];
sx q[0];
rz(-2.6840211) q[0];
sx q[0];
rz(0.81480169) q[0];
rz(1.762215) q[1];
sx q[1];
rz(-0.35019362) q[1];
sx q[1];
rz(-2.8864158) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4677306) q[0];
sx q[0];
rz(-0.47463372) q[0];
sx q[0];
rz(0.69068308) q[0];
x q[1];
rz(-1.2404664) q[2];
sx q[2];
rz(-2.1836046) q[2];
sx q[2];
rz(1.4532879) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7123588) q[1];
sx q[1];
rz(-0.36839596) q[1];
sx q[1];
rz(0.952094) q[1];
rz(-pi) q[2];
rz(1.4820443) q[3];
sx q[3];
rz(-2.6188861) q[3];
sx q[3];
rz(-2.6356217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2531551) q[2];
sx q[2];
rz(-1.5506813) q[2];
sx q[2];
rz(-0.17318428) q[2];
rz(0.52982461) q[3];
sx q[3];
rz(-0.14557043) q[3];
sx q[3];
rz(-0.1023275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.2816876) q[0];
sx q[0];
rz(-1.6485933) q[0];
sx q[0];
rz(1.7657071) q[0];
rz(-1.2777404) q[1];
sx q[1];
rz(-2.3294096) q[1];
sx q[1];
rz(3.0854991) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1306886) q[0];
sx q[0];
rz(-2.4405257) q[0];
sx q[0];
rz(0.58347337) q[0];
rz(-pi) q[1];
rz(2.7943139) q[2];
sx q[2];
rz(-1.1355073) q[2];
sx q[2];
rz(-1.4594644) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1069378) q[1];
sx q[1];
rz(-1.5853303) q[1];
sx q[1];
rz(-0.30277877) q[1];
rz(-1.4291184) q[3];
sx q[3];
rz(-2.5767527) q[3];
sx q[3];
rz(2.7431938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.00099480199) q[2];
sx q[2];
rz(-2.2237015) q[2];
sx q[2];
rz(2.7094005) q[2];
rz(2.2473992) q[3];
sx q[3];
rz(-2.0420572) q[3];
sx q[3];
rz(-1.6754707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11319259) q[0];
sx q[0];
rz(-0.88554651) q[0];
sx q[0];
rz(-2.4940441) q[0];
rz(1.8796857) q[1];
sx q[1];
rz(-1.4636661) q[1];
sx q[1];
rz(-0.9544968) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0237085) q[0];
sx q[0];
rz(-0.55340289) q[0];
sx q[0];
rz(-1.8694359) q[0];
rz(-pi) q[1];
rz(1.5231832) q[2];
sx q[2];
rz(-0.63112586) q[2];
sx q[2];
rz(0.39436755) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8050025) q[1];
sx q[1];
rz(-1.0953566) q[1];
sx q[1];
rz(0.57979433) q[1];
rz(-pi) q[2];
x q[2];
rz(0.16320634) q[3];
sx q[3];
rz(-1.2833793) q[3];
sx q[3];
rz(0.67374574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.59297562) q[2];
sx q[2];
rz(-1.2334712) q[2];
sx q[2];
rz(-1.0423638) q[2];
rz(-0.43867612) q[3];
sx q[3];
rz(-1.0500267) q[3];
sx q[3];
rz(1.8235122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8577268) q[0];
sx q[0];
rz(-0.23290578) q[0];
sx q[0];
rz(-2.3983811) q[0];
rz(-1.6339533) q[1];
sx q[1];
rz(-2.4217024) q[1];
sx q[1];
rz(-0.61002237) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2110721) q[0];
sx q[0];
rz(-1.0809582) q[0];
sx q[0];
rz(-1.2852438) q[0];
rz(-pi) q[1];
rz(-0.1158175) q[2];
sx q[2];
rz(-0.91997416) q[2];
sx q[2];
rz(-1.4027558) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.66079084) q[1];
sx q[1];
rz(-1.5958438) q[1];
sx q[1];
rz(-2.2381496) q[1];
x q[2];
rz(-2.9355572) q[3];
sx q[3];
rz(-1.9427951) q[3];
sx q[3];
rz(-2.4601065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.33621776) q[2];
sx q[2];
rz(-1.4423794) q[2];
sx q[2];
rz(2.2231893) q[2];
rz(1.5504799) q[3];
sx q[3];
rz(-2.1912626) q[3];
sx q[3];
rz(-2.7526855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
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
rz(-0.36088762) q[0];
sx q[0];
rz(-0.66910678) q[0];
sx q[0];
rz(-1.5135182) q[0];
rz(-0.52945119) q[1];
sx q[1];
rz(-2.0748731) q[1];
sx q[1];
rz(-0.73658529) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.026222762) q[0];
sx q[0];
rz(-1.6008458) q[0];
sx q[0];
rz(0.010859246) q[0];
rz(-pi) q[1];
rz(-1.2543711) q[2];
sx q[2];
rz(-2.8036615) q[2];
sx q[2];
rz(-0.28883176) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.089162) q[1];
sx q[1];
rz(-2.6658635) q[1];
sx q[1];
rz(2.9176941) q[1];
x q[2];
rz(-0.32324507) q[3];
sx q[3];
rz(-1.9976227) q[3];
sx q[3];
rz(2.1904898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1071876) q[2];
sx q[2];
rz(-1.9341058) q[2];
sx q[2];
rz(-0.68391189) q[2];
rz(-1.2290139) q[3];
sx q[3];
rz(-1.7714272) q[3];
sx q[3];
rz(-1.3945403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9443611) q[0];
sx q[0];
rz(-1.5988388) q[0];
sx q[0];
rz(-2.9558682) q[0];
rz(-0.99705237) q[1];
sx q[1];
rz(-1.8763708) q[1];
sx q[1];
rz(2.396778) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1574402) q[0];
sx q[0];
rz(-2.3072349) q[0];
sx q[0];
rz(2.0122583) q[0];
rz(-pi) q[1];
rz(1.9717734) q[2];
sx q[2];
rz(-2.3224761) q[2];
sx q[2];
rz(-2.4405406) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3639431) q[1];
sx q[1];
rz(-1.3476106) q[1];
sx q[1];
rz(1.9294444) q[1];
x q[2];
rz(2.7471077) q[3];
sx q[3];
rz(-2.3186765) q[3];
sx q[3];
rz(2.8701973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0361438) q[2];
sx q[2];
rz(-0.75573409) q[2];
sx q[2];
rz(-1.194681) q[2];
rz(2.1448994) q[3];
sx q[3];
rz(-1.2160622) q[3];
sx q[3];
rz(-0.99635807) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3982518) q[0];
sx q[0];
rz(-0.78813362) q[0];
sx q[0];
rz(2.7375896) q[0];
rz(-3.1104654) q[1];
sx q[1];
rz(-1.4844091) q[1];
sx q[1];
rz(-1.9706479) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084328018) q[0];
sx q[0];
rz(-1.9976166) q[0];
sx q[0];
rz(-1.5469993) q[0];
x q[1];
rz(-0.77565907) q[2];
sx q[2];
rz(-1.7663029) q[2];
sx q[2];
rz(-2.069371) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3552637) q[1];
sx q[1];
rz(-1.5794282) q[1];
sx q[1];
rz(1.5481871) q[1];
rz(-pi) q[2];
rz(2.5578299) q[3];
sx q[3];
rz(-2.6423892) q[3];
sx q[3];
rz(-1.1399869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6955473) q[2];
sx q[2];
rz(-1.7798767) q[2];
sx q[2];
rz(-0.5919624) q[2];
rz(2.5752318) q[3];
sx q[3];
rz(-2.9768894) q[3];
sx q[3];
rz(-1.6177572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3175209) q[0];
sx q[0];
rz(-2.1614647) q[0];
sx q[0];
rz(1.9807057) q[0];
rz(3.042165) q[1];
sx q[1];
rz(-1.2482523) q[1];
sx q[1];
rz(-2.0773239) q[1];
rz(0.912491) q[2];
sx q[2];
rz(-0.9463263) q[2];
sx q[2];
rz(1.8778388) q[2];
rz(-0.016146544) q[3];
sx q[3];
rz(-1.2373677) q[3];
sx q[3];
rz(2.2131372) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];