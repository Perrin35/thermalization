OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7140952) q[0];
sx q[0];
rz(3.7035898) q[0];
sx q[0];
rz(9.193767) q[0];
rz(0.24569874) q[1];
sx q[1];
rz(2.6872771) q[1];
sx q[1];
rz(8.1375577) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19350138) q[0];
sx q[0];
rz(-2.5768902) q[0];
sx q[0];
rz(-2.094784) q[0];
rz(-2.8962249) q[2];
sx q[2];
rz(-1.7984219) q[2];
sx q[2];
rz(0.99492225) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.27632984) q[1];
sx q[1];
rz(-1.4722605) q[1];
sx q[1];
rz(0.16698412) q[1];
rz(1.7460398) q[3];
sx q[3];
rz(-1.6449494) q[3];
sx q[3];
rz(-1.4369304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7370558) q[2];
sx q[2];
rz(-0.70713592) q[2];
sx q[2];
rz(-2.7810968) q[2];
rz(-2.0170085) q[3];
sx q[3];
rz(-2.093061) q[3];
sx q[3];
rz(1.8726965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-2.8812113) q[0];
sx q[0];
rz(-2.9660048) q[0];
sx q[0];
rz(2.4847109) q[0];
rz(2.8086713) q[1];
sx q[1];
rz(-2.0491144) q[1];
sx q[1];
rz(2.2064256) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1293949) q[0];
sx q[0];
rz(-0.10888616) q[0];
sx q[0];
rz(2.8588141) q[0];
rz(-pi) q[1];
rz(0.25721154) q[2];
sx q[2];
rz(-1.4818483) q[2];
sx q[2];
rz(2.6572029) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3352901) q[1];
sx q[1];
rz(-2.9893901) q[1];
sx q[1];
rz(-2.8782513) q[1];
rz(1.7788497) q[3];
sx q[3];
rz(-1.2223772) q[3];
sx q[3];
rz(1.6121685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8344581) q[2];
sx q[2];
rz(-1.5156526) q[2];
sx q[2];
rz(1.2163986) q[2];
rz(2.6136716) q[3];
sx q[3];
rz(-1.0390176) q[3];
sx q[3];
rz(1.886604) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5288178) q[0];
sx q[0];
rz(-0.82614326) q[0];
sx q[0];
rz(2.5566027) q[0];
rz(-0.12475573) q[1];
sx q[1];
rz(-0.59586066) q[1];
sx q[1];
rz(1.4488719) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9609279) q[0];
sx q[0];
rz(-1.573473) q[0];
sx q[0];
rz(-1.5737783) q[0];
rz(-pi) q[1];
rz(-1.5572433) q[2];
sx q[2];
rz(-1.6429735) q[2];
sx q[2];
rz(-0.37301317) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8068778) q[1];
sx q[1];
rz(-0.92744614) q[1];
sx q[1];
rz(-2.592464) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8336867) q[3];
sx q[3];
rz(-1.8975199) q[3];
sx q[3];
rz(-2.9430091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2584194) q[2];
sx q[2];
rz(-1.0018188) q[2];
sx q[2];
rz(-1.6955356) q[2];
rz(-2.3840733) q[3];
sx q[3];
rz(-1.5816553) q[3];
sx q[3];
rz(-2.3022046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3676753) q[0];
sx q[0];
rz(-2.2125419) q[0];
sx q[0];
rz(-2.0539334) q[0];
rz(-2.7663973) q[1];
sx q[1];
rz(-2.0162069) q[1];
sx q[1];
rz(-0.96022022) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86726928) q[0];
sx q[0];
rz(-1.835863) q[0];
sx q[0];
rz(-0.041170711) q[0];
rz(-pi) q[1];
rz(-1.2795936) q[2];
sx q[2];
rz(-3.0377977) q[2];
sx q[2];
rz(1.4196996) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8254777) q[1];
sx q[1];
rz(-1.1255472) q[1];
sx q[1];
rz(0.47834217) q[1];
x q[2];
rz(0.26953617) q[3];
sx q[3];
rz(-1.4785267) q[3];
sx q[3];
rz(0.59449457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9356161) q[2];
sx q[2];
rz(-1.2674067) q[2];
sx q[2];
rz(-1.6440294) q[2];
rz(-2.2802672) q[3];
sx q[3];
rz(-0.70278168) q[3];
sx q[3];
rz(-2.6845045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1645666) q[0];
sx q[0];
rz(-2.470546) q[0];
sx q[0];
rz(-2.5373051) q[0];
rz(-1.1445649) q[1];
sx q[1];
rz(-1.6555758) q[1];
sx q[1];
rz(2.7191275) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9086994) q[0];
sx q[0];
rz(-2.9541203) q[0];
sx q[0];
rz(1.1195681) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.56296641) q[2];
sx q[2];
rz(-2.0694695) q[2];
sx q[2];
rz(0.50567217) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8026655) q[1];
sx q[1];
rz(-2.4484854) q[1];
sx q[1];
rz(3.0285804) q[1];
rz(0.38099184) q[3];
sx q[3];
rz(-1.89011) q[3];
sx q[3];
rz(1.8007985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8233238) q[2];
sx q[2];
rz(-0.6508998) q[2];
sx q[2];
rz(0.65625119) q[2];
rz(-1.6107791) q[3];
sx q[3];
rz(-1.2698413) q[3];
sx q[3];
rz(-2.5608565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73606473) q[0];
sx q[0];
rz(-1.0117714) q[0];
sx q[0];
rz(2.5166125) q[0];
rz(2.3896353) q[1];
sx q[1];
rz(-1.0686921) q[1];
sx q[1];
rz(-1.5350852) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67095516) q[0];
sx q[0];
rz(-0.29330253) q[0];
sx q[0];
rz(0.17099149) q[0];
x q[1];
rz(2.4322832) q[2];
sx q[2];
rz(-2.0481718) q[2];
sx q[2];
rz(-2.2717182) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.62404406) q[1];
sx q[1];
rz(-0.66429783) q[1];
sx q[1];
rz(0.37599998) q[1];
x q[2];
rz(2.296359) q[3];
sx q[3];
rz(-0.74016011) q[3];
sx q[3];
rz(-0.46421591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9525166) q[2];
sx q[2];
rz(-1.7837046) q[2];
sx q[2];
rz(-0.9291741) q[2];
rz(0.82516986) q[3];
sx q[3];
rz(-1.5114096) q[3];
sx q[3];
rz(1.1024124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55649844) q[0];
sx q[0];
rz(-0.80699054) q[0];
sx q[0];
rz(1.7671385) q[0];
rz(0.51042405) q[1];
sx q[1];
rz(-2.3511032) q[1];
sx q[1];
rz(0.62320954) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5730815) q[0];
sx q[0];
rz(-2.7740736) q[0];
sx q[0];
rz(2.0042242) q[0];
rz(-pi) q[1];
rz(-2.6323967) q[2];
sx q[2];
rz(-1.7422581) q[2];
sx q[2];
rz(0.96656884) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.53155073) q[1];
sx q[1];
rz(-1.7167713) q[1];
sx q[1];
rz(0.11256217) q[1];
x q[2];
rz(-2.9233452) q[3];
sx q[3];
rz(-0.75957662) q[3];
sx q[3];
rz(1.4324566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0509384) q[2];
sx q[2];
rz(-0.81523681) q[2];
sx q[2];
rz(1.9742924) q[2];
rz(-3.1189611) q[3];
sx q[3];
rz(-0.52112094) q[3];
sx q[3];
rz(-1.1289271) q[3];
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
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8959494) q[0];
sx q[0];
rz(-1.1351981) q[0];
sx q[0];
rz(1.0173215) q[0];
rz(1.3941049) q[1];
sx q[1];
rz(-0.21109763) q[1];
sx q[1];
rz(2.5140433) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8857972) q[0];
sx q[0];
rz(-1.1748519) q[0];
sx q[0];
rz(-2.2084479) q[0];
rz(1.6553788) q[2];
sx q[2];
rz(-2.248317) q[2];
sx q[2];
rz(-1.6135474) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.095762091) q[1];
sx q[1];
rz(-1.374829) q[1];
sx q[1];
rz(-2.8992871) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5379337) q[3];
sx q[3];
rz(-0.55760477) q[3];
sx q[3];
rz(1.4035567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5158186) q[2];
sx q[2];
rz(-0.63483441) q[2];
sx q[2];
rz(-0.41075692) q[2];
rz(0.050203236) q[3];
sx q[3];
rz(-1.2529195) q[3];
sx q[3];
rz(-3.0661809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5589767) q[0];
sx q[0];
rz(-0.89685431) q[0];
sx q[0];
rz(0.26891747) q[0];
rz(2.8315663) q[1];
sx q[1];
rz(-2.0545484) q[1];
sx q[1];
rz(0.1300098) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5101679) q[0];
sx q[0];
rz(-1.893259) q[0];
sx q[0];
rz(2.9084951) q[0];
rz(-pi) q[1];
x q[1];
rz(1.214005) q[2];
sx q[2];
rz(-1.2399106) q[2];
sx q[2];
rz(-0.39168229) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0669032) q[1];
sx q[1];
rz(-0.58614158) q[1];
sx q[1];
rz(2.8476004) q[1];
x q[2];
rz(-2.1653529) q[3];
sx q[3];
rz(-0.77139445) q[3];
sx q[3];
rz(-2.6736128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2076063) q[2];
sx q[2];
rz(-0.76675582) q[2];
sx q[2];
rz(3.0450191) q[2];
rz(1.6377595) q[3];
sx q[3];
rz(-1.8042754) q[3];
sx q[3];
rz(0.79969978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12829517) q[0];
sx q[0];
rz(-0.18823637) q[0];
sx q[0];
rz(0.53303322) q[0];
rz(-3.10532) q[1];
sx q[1];
rz(-0.78250042) q[1];
sx q[1];
rz(-1.689555) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5851327) q[0];
sx q[0];
rz(-1.0074573) q[0];
sx q[0];
rz(2.6908532) q[0];
x q[1];
rz(-1.4617986) q[2];
sx q[2];
rz(-1.8882635) q[2];
sx q[2];
rz(-0.017680971) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9033135) q[1];
sx q[1];
rz(-0.36064816) q[1];
sx q[1];
rz(-0.084899501) q[1];
rz(-pi) q[2];
rz(-0.17590268) q[3];
sx q[3];
rz(-0.47187343) q[3];
sx q[3];
rz(2.7681818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4286246) q[2];
sx q[2];
rz(-0.73178256) q[2];
sx q[2];
rz(-0.0326322) q[2];
rz(0.84515682) q[3];
sx q[3];
rz(-1.9149575) q[3];
sx q[3];
rz(2.0613861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7964771) q[0];
sx q[0];
rz(-0.65237541) q[0];
sx q[0];
rz(2.8271578) q[0];
rz(0.79700094) q[1];
sx q[1];
rz(-0.92756699) q[1];
sx q[1];
rz(-3.0753593) q[1];
rz(-1.2751515) q[2];
sx q[2];
rz(-2.0870123) q[2];
sx q[2];
rz(-0.43475702) q[2];
rz(1.6975523) q[3];
sx q[3];
rz(-0.76382617) q[3];
sx q[3];
rz(2.7504117) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
