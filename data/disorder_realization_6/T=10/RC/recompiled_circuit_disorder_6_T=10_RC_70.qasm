OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.52580994) q[0];
sx q[0];
rz(-1.7237741) q[0];
sx q[0];
rz(-0.56086993) q[0];
rz(-2.0286735) q[1];
sx q[1];
rz(-1.3781883) q[1];
sx q[1];
rz(-1.2150432) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96072223) q[0];
sx q[0];
rz(-1.3836432) q[0];
sx q[0];
rz(-3.0357009) q[0];
rz(2.4807793) q[2];
sx q[2];
rz(-1.3408957) q[2];
sx q[2];
rz(-1.7125318) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.58303761) q[1];
sx q[1];
rz(-1.8999294) q[1];
sx q[1];
rz(2.1117044) q[1];
rz(-pi) q[2];
rz(-1.3052985) q[3];
sx q[3];
rz(-1.4269097) q[3];
sx q[3];
rz(-3.0778411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3540196) q[2];
sx q[2];
rz(-2.188787) q[2];
sx q[2];
rz(-2.9585178) q[2];
rz(0.37781528) q[3];
sx q[3];
rz(-2.0928045) q[3];
sx q[3];
rz(-2.8474076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29782444) q[0];
sx q[0];
rz(-2.4968708) q[0];
sx q[0];
rz(3.0644754) q[0];
rz(2.8027957) q[1];
sx q[1];
rz(-2.0270551) q[1];
sx q[1];
rz(1.5391301) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3829271) q[0];
sx q[0];
rz(-2.5490952) q[0];
sx q[0];
rz(1.4325607) q[0];
rz(0.23761959) q[2];
sx q[2];
rz(-2.3887206) q[2];
sx q[2];
rz(1.876229) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3746678) q[1];
sx q[1];
rz(-2.6267849) q[1];
sx q[1];
rz(-1.9306081) q[1];
x q[2];
rz(-1.8061403) q[3];
sx q[3];
rz(-1.1574405) q[3];
sx q[3];
rz(-1.2873161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2960647) q[2];
sx q[2];
rz(-1.277593) q[2];
sx q[2];
rz(-2.4831333) q[2];
rz(0.15130875) q[3];
sx q[3];
rz(-1.0226117) q[3];
sx q[3];
rz(-0.69491274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028458683) q[0];
sx q[0];
rz(-0.78335339) q[0];
sx q[0];
rz(0.43310305) q[0];
rz(1.9494879) q[1];
sx q[1];
rz(-1.2116218) q[1];
sx q[1];
rz(-2.5862397) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1028324) q[0];
sx q[0];
rz(-2.0534678) q[0];
sx q[0];
rz(-2.32248) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31366445) q[2];
sx q[2];
rz(-2.1137538) q[2];
sx q[2];
rz(2.0344337) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8086116) q[1];
sx q[1];
rz(-2.4858027) q[1];
sx q[1];
rz(1.8379184) q[1];
rz(-pi) q[2];
x q[2];
rz(1.025612) q[3];
sx q[3];
rz(-1.7435929) q[3];
sx q[3];
rz(-1.9219414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.73734036) q[2];
sx q[2];
rz(-2.3534687) q[2];
sx q[2];
rz(-1.2505442) q[2];
rz(-2.897443) q[3];
sx q[3];
rz(-1.8593676) q[3];
sx q[3];
rz(-1.4499433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8811532) q[0];
sx q[0];
rz(-0.45757159) q[0];
sx q[0];
rz(0.81480169) q[0];
rz(1.762215) q[1];
sx q[1];
rz(-0.35019362) q[1];
sx q[1];
rz(0.25517685) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92514738) q[0];
sx q[0];
rz(-1.9307923) q[0];
sx q[0];
rz(-1.2544592) q[0];
rz(-pi) q[1];
rz(-0.63919477) q[2];
sx q[2];
rz(-1.8393469) q[2];
sx q[2];
rz(-3.064379) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5863122) q[1];
sx q[1];
rz(-1.781207) q[1];
sx q[1];
rz(1.2661238) q[1];
rz(3.0905686) q[3];
sx q[3];
rz(-1.050356) q[3];
sx q[3];
rz(-2.737962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.8884376) q[2];
sx q[2];
rz(-1.5909114) q[2];
sx q[2];
rz(0.17318428) q[2];
rz(-2.611768) q[3];
sx q[3];
rz(-2.9960222) q[3];
sx q[3];
rz(-3.0392652) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.859905) q[0];
sx q[0];
rz(-1.6485933) q[0];
sx q[0];
rz(-1.3758855) q[0];
rz(-1.8638523) q[1];
sx q[1];
rz(-0.81218305) q[1];
sx q[1];
rz(3.0854991) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2984021) q[0];
sx q[0];
rz(-1.0023596) q[0];
sx q[0];
rz(2.0060904) q[0];
rz(1.1115083) q[2];
sx q[2];
rz(-1.2570724) q[2];
sx q[2];
rz(2.8788061) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1069378) q[1];
sx q[1];
rz(-1.5853303) q[1];
sx q[1];
rz(-2.8388139) q[1];
rz(-pi) q[2];
rz(-1.7124743) q[3];
sx q[3];
rz(-0.56483993) q[3];
sx q[3];
rz(-0.39839881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1405979) q[2];
sx q[2];
rz(-2.2237015) q[2];
sx q[2];
rz(0.43219217) q[2];
rz(2.2473992) q[3];
sx q[3];
rz(-2.0420572) q[3];
sx q[3];
rz(-1.6754707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11319259) q[0];
sx q[0];
rz(-2.2560461) q[0];
sx q[0];
rz(-2.4940441) q[0];
rz(1.2619069) q[1];
sx q[1];
rz(-1.4636661) q[1];
sx q[1];
rz(0.9544968) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0237085) q[0];
sx q[0];
rz(-0.55340289) q[0];
sx q[0];
rz(-1.2721567) q[0];
rz(-pi) q[1];
rz(-0.94021057) q[2];
sx q[2];
rz(-1.5988837) q[2];
sx q[2];
rz(-1.1379776) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8050025) q[1];
sx q[1];
rz(-2.046236) q[1];
sx q[1];
rz(0.57979433) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2797221) q[3];
sx q[3];
rz(-1.4143412) q[3];
sx q[3];
rz(-0.94369704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.59297562) q[2];
sx q[2];
rz(-1.2334712) q[2];
sx q[2];
rz(2.0992289) q[2];
rz(-0.43867612) q[3];
sx q[3];
rz(-1.0500267) q[3];
sx q[3];
rz(1.8235122) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2838659) q[0];
sx q[0];
rz(-0.23290578) q[0];
sx q[0];
rz(2.3983811) q[0];
rz(-1.5076393) q[1];
sx q[1];
rz(-2.4217024) q[1];
sx q[1];
rz(-2.5315703) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2110721) q[0];
sx q[0];
rz(-2.0606344) q[0];
sx q[0];
rz(1.2852438) q[0];
x q[1];
rz(1.4201944) q[2];
sx q[2];
rz(-2.4820231) q[2];
sx q[2];
rz(1.9285551) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.66079084) q[1];
sx q[1];
rz(-1.5958438) q[1];
sx q[1];
rz(-2.2381496) q[1];
x q[2];
rz(2.9355572) q[3];
sx q[3];
rz(-1.1987975) q[3];
sx q[3];
rz(-2.4601065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.33621776) q[2];
sx q[2];
rz(-1.6992133) q[2];
sx q[2];
rz(0.91840333) q[2];
rz(-1.5911128) q[3];
sx q[3];
rz(-0.95033002) q[3];
sx q[3];
rz(-0.38890719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.780705) q[0];
sx q[0];
rz(-2.4724859) q[0];
sx q[0];
rz(1.5135182) q[0];
rz(-0.52945119) q[1];
sx q[1];
rz(-1.0667195) q[1];
sx q[1];
rz(0.73658529) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7685331) q[0];
sx q[0];
rz(-3.1096418) q[0];
sx q[0];
rz(1.2241227) q[0];
rz(-3.0326764) q[2];
sx q[2];
rz(-1.250259) q[2];
sx q[2];
rz(0.045217302) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.089162) q[1];
sx q[1];
rz(-2.6658635) q[1];
sx q[1];
rz(-2.9176941) q[1];
rz(-pi) q[2];
rz(-2.1804817) q[3];
sx q[3];
rz(-2.612252) q[3];
sx q[3];
rz(1.6314268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0344051) q[2];
sx q[2];
rz(-1.2074869) q[2];
sx q[2];
rz(0.68391189) q[2];
rz(-1.9125787) q[3];
sx q[3];
rz(-1.3701655) q[3];
sx q[3];
rz(-1.3945403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9443611) q[0];
sx q[0];
rz(-1.5988388) q[0];
sx q[0];
rz(-0.18572447) q[0];
rz(-0.99705237) q[1];
sx q[1];
rz(-1.2652218) q[1];
sx q[1];
rz(-2.396778) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72070044) q[0];
sx q[0];
rz(-1.8928796) q[0];
sx q[0];
rz(2.3548404) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.79297519) q[2];
sx q[2];
rz(-1.2816396) q[2];
sx q[2];
rz(0.58794978) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4312268) q[1];
sx q[1];
rz(-1.2214298) q[1];
sx q[1];
rz(0.23780312) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7471077) q[3];
sx q[3];
rz(-2.3186765) q[3];
sx q[3];
rz(2.8701973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0361438) q[2];
sx q[2];
rz(-0.75573409) q[2];
sx q[2];
rz(1.194681) q[2];
rz(-0.99669325) q[3];
sx q[3];
rz(-1.2160622) q[3];
sx q[3];
rz(-0.99635807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3982518) q[0];
sx q[0];
rz(-0.78813362) q[0];
sx q[0];
rz(-2.7375896) q[0];
rz(-0.031127302) q[1];
sx q[1];
rz(-1.4844091) q[1];
sx q[1];
rz(-1.1709447) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0572646) q[0];
sx q[0];
rz(-1.9976166) q[0];
sx q[0];
rz(-1.5469993) q[0];
x q[1];
rz(-0.77565907) q[2];
sx q[2];
rz(-1.7663029) q[2];
sx q[2];
rz(1.0722216) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3569301) q[1];
sx q[1];
rz(-1.5934048) q[1];
sx q[1];
rz(3.1329586) q[1];
rz(-1.8627432) q[3];
sx q[3];
rz(-1.1598831) q[3];
sx q[3];
rz(1.3565855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6955473) q[2];
sx q[2];
rz(-1.3617159) q[2];
sx q[2];
rz(0.5919624) q[2];
rz(-0.56636089) q[3];
sx q[3];
rz(-0.16470328) q[3];
sx q[3];
rz(1.6177572) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82407172) q[0];
sx q[0];
rz(-0.98012797) q[0];
sx q[0];
rz(-1.160887) q[0];
rz(-0.099427632) q[1];
sx q[1];
rz(-1.2482523) q[1];
sx q[1];
rz(-2.0773239) q[1];
rz(2.4377433) q[2];
sx q[2];
rz(-2.2675632) q[2];
sx q[2];
rz(2.8013196) q[2];
rz(1.2373274) q[3];
sx q[3];
rz(-1.5860535) q[3];
sx q[3];
rz(0.64762583) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
