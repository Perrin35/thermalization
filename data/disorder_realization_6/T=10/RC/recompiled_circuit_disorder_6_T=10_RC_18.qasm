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
rz(4.5594112) q[0];
sx q[0];
rz(8.863908) q[0];
rz(4.2545118) q[1];
sx q[1];
rz(1.7634044) q[1];
sx q[1];
rz(7.4982285) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44170609) q[0];
sx q[0];
rz(-2.9268648) q[0];
sx q[0];
rz(1.0617274) q[0];
rz(-pi) q[1];
rz(-1.8589852) q[2];
sx q[2];
rz(-0.93027861) q[2];
sx q[2];
rz(-0.033601947) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.79599586) q[1];
sx q[1];
rz(-2.0797634) q[1];
sx q[1];
rz(2.7624346) q[1];
x q[2];
rz(1.3052985) q[3];
sx q[3];
rz(-1.4269097) q[3];
sx q[3];
rz(3.0778411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3540196) q[2];
sx q[2];
rz(-2.188787) q[2];
sx q[2];
rz(-0.18307486) q[2];
rz(2.7637774) q[3];
sx q[3];
rz(-2.0928045) q[3];
sx q[3];
rz(2.8474076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29782444) q[0];
sx q[0];
rz(-0.64472187) q[0];
sx q[0];
rz(0.077117292) q[0];
rz(-0.33879694) q[1];
sx q[1];
rz(-2.0270551) q[1];
sx q[1];
rz(1.5391301) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30277006) q[0];
sx q[0];
rz(-1.6478224) q[0];
sx q[0];
rz(2.1588615) q[0];
rz(-pi) q[1];
x q[1];
rz(2.402926) q[2];
sx q[2];
rz(-1.4091485) q[2];
sx q[2];
rz(-2.6612298) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1748845) q[1];
sx q[1];
rz(-1.0918573) q[1];
sx q[1];
rz(-2.944988) q[1];
x q[2];
rz(-2.7178571) q[3];
sx q[3];
rz(-1.3556004) q[3];
sx q[3];
rz(0.37950294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.845528) q[2];
sx q[2];
rz(-1.277593) q[2];
sx q[2];
rz(0.65845931) q[2];
rz(-2.9902839) q[3];
sx q[3];
rz(-1.0226117) q[3];
sx q[3];
rz(2.4466799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.028458683) q[0];
sx q[0];
rz(-2.3582393) q[0];
sx q[0];
rz(0.43310305) q[0];
rz(-1.9494879) q[1];
sx q[1];
rz(-1.2116218) q[1];
sx q[1];
rz(-0.55535299) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0387602) q[0];
sx q[0];
rz(-1.0881249) q[0];
sx q[0];
rz(-2.32248) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0981512) q[2];
sx q[2];
rz(-2.522509) q[2];
sx q[2];
rz(1.4738136) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.66545031) q[1];
sx q[1];
rz(-2.1995771) q[1];
sx q[1];
rz(0.20035845) q[1];
rz(-pi) q[2];
x q[2];
rz(0.20136307) q[3];
sx q[3];
rz(-2.106973) q[3];
sx q[3];
rz(2.686543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4042523) q[2];
sx q[2];
rz(-0.78812391) q[2];
sx q[2];
rz(1.2505442) q[2];
rz(0.2441497) q[3];
sx q[3];
rz(-1.282225) q[3];
sx q[3];
rz(-1.6916493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(-0.26043949) q[0];
sx q[0];
rz(-0.45757159) q[0];
sx q[0];
rz(2.326791) q[0];
rz(-1.3793777) q[1];
sx q[1];
rz(-2.791399) q[1];
sx q[1];
rz(-0.25517685) q[1];
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
rz(-2.5023979) q[2];
sx q[2];
rz(-1.3022458) q[2];
sx q[2];
rz(0.07721363) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5863122) q[1];
sx q[1];
rz(-1.3603856) q[1];
sx q[1];
rz(-1.8754688) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4820443) q[3];
sx q[3];
rz(-0.52270652) q[3];
sx q[3];
rz(0.50597092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2531551) q[2];
sx q[2];
rz(-1.5506813) q[2];
sx q[2];
rz(-2.9684084) q[2];
rz(0.52982461) q[3];
sx q[3];
rz(-2.9960222) q[3];
sx q[3];
rz(-3.0392652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(1.859905) q[0];
sx q[0];
rz(-1.6485933) q[0];
sx q[0];
rz(-1.7657071) q[0];
rz(1.2777404) q[1];
sx q[1];
rz(-0.81218305) q[1];
sx q[1];
rz(-0.056093562) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2984021) q[0];
sx q[0];
rz(-2.139233) q[0];
sx q[0];
rz(2.0060904) q[0];
x q[1];
rz(0.9390097) q[2];
sx q[2];
rz(-0.54982215) q[2];
sx q[2];
rz(0.75013559) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6519421) q[1];
sx q[1];
rz(-0.30311668) q[1];
sx q[1];
rz(0.048708212) q[1];
rz(-pi) q[2];
rz(1.7124743) q[3];
sx q[3];
rz(-0.56483993) q[3];
sx q[3];
rz(-2.7431938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.00099480199) q[2];
sx q[2];
rz(-0.91789118) q[2];
sx q[2];
rz(-2.7094005) q[2];
rz(0.8941935) q[3];
sx q[3];
rz(-1.0995355) q[3];
sx q[3];
rz(-1.6754707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11319259) q[0];
sx q[0];
rz(-2.2560461) q[0];
sx q[0];
rz(-0.64754852) q[0];
rz(-1.2619069) q[1];
sx q[1];
rz(-1.6779265) q[1];
sx q[1];
rz(0.9544968) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.944825) q[0];
sx q[0];
rz(-1.4155354) q[0];
sx q[0];
rz(1.0374271) q[0];
x q[1];
rz(-3.1068222) q[2];
sx q[2];
rz(-0.94049847) q[2];
sx q[2];
rz(-2.6882753) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8442321) q[1];
sx q[1];
rz(-2.4095222) q[1];
sx q[1];
rz(2.3872603) q[1];
rz(-pi) q[2];
rz(2.0733842) q[3];
sx q[3];
rz(-2.8121901) q[3];
sx q[3];
rz(0.1474895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.59297562) q[2];
sx q[2];
rz(-1.9081215) q[2];
sx q[2];
rz(-1.0423638) q[2];
rz(2.7029165) q[3];
sx q[3];
rz(-2.091566) q[3];
sx q[3];
rz(-1.8235122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8577268) q[0];
sx q[0];
rz(-2.9086869) q[0];
sx q[0];
rz(2.3983811) q[0];
rz(-1.6339533) q[1];
sx q[1];
rz(-0.71989027) q[1];
sx q[1];
rz(0.61002237) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9191223) q[0];
sx q[0];
rz(-1.3195992) q[0];
sx q[0];
rz(-0.50719502) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7213983) q[2];
sx q[2];
rz(-2.4820231) q[2];
sx q[2];
rz(1.9285551) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2513189) q[1];
sx q[1];
rz(-2.2379025) q[1];
sx q[1];
rz(-3.1097079) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.20603541) q[3];
sx q[3];
rz(-1.1987975) q[3];
sx q[3];
rz(0.68148617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.780705) q[0];
sx q[0];
rz(-0.66910678) q[0];
sx q[0];
rz(-1.6280744) q[0];
rz(0.52945119) q[1];
sx q[1];
rz(-2.0748731) q[1];
sx q[1];
rz(-2.4050074) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7685331) q[0];
sx q[0];
rz(-3.1096418) q[0];
sx q[0];
rz(1.2241227) q[0];
rz(-pi) q[1];
rz(0.10891624) q[2];
sx q[2];
rz(-1.250259) q[2];
sx q[2];
rz(0.045217302) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0524307) q[1];
sx q[1];
rz(-2.6658635) q[1];
sx q[1];
rz(0.22389852) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96111091) q[3];
sx q[3];
rz(-2.612252) q[3];
sx q[3];
rz(-1.6314268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0344051) q[2];
sx q[2];
rz(-1.9341058) q[2];
sx q[2];
rz(2.4576808) q[2];
rz(-1.2290139) q[3];
sx q[3];
rz(-1.3701655) q[3];
sx q[3];
rz(1.3945403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9443611) q[0];
sx q[0];
rz(-1.5427538) q[0];
sx q[0];
rz(-0.18572447) q[0];
rz(-0.99705237) q[1];
sx q[1];
rz(-1.2652218) q[1];
sx q[1];
rz(-2.396778) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54430994) q[0];
sx q[0];
rz(-2.3047857) q[0];
sx q[0];
rz(-0.44041667) q[0];
rz(-0.79297519) q[2];
sx q[2];
rz(-1.2816396) q[2];
sx q[2];
rz(0.58794978) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4312268) q[1];
sx q[1];
rz(-1.9201628) q[1];
sx q[1];
rz(0.23780312) q[1];
rz(-pi) q[2];
rz(-0.39448491) q[3];
sx q[3];
rz(-2.3186765) q[3];
sx q[3];
rz(-0.2713954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0361438) q[2];
sx q[2];
rz(-2.3858586) q[2];
sx q[2];
rz(1.194681) q[2];
rz(2.1448994) q[3];
sx q[3];
rz(-1.2160622) q[3];
sx q[3];
rz(-0.99635807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74334082) q[0];
sx q[0];
rz(-0.78813362) q[0];
sx q[0];
rz(-0.40400305) q[0];
rz(-0.031127302) q[1];
sx q[1];
rz(-1.6571836) q[1];
sx q[1];
rz(1.1709447) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.02689657) q[0];
sx q[0];
rz(-2.7141502) q[0];
sx q[0];
rz(-0.052274152) q[0];
rz(-pi) q[1];
rz(-2.3659336) q[2];
sx q[2];
rz(-1.7663029) q[2];
sx q[2];
rz(2.069371) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3552637) q[1];
sx q[1];
rz(-1.5794282) q[1];
sx q[1];
rz(-1.5934056) q[1];
x q[2];
rz(-2.714614) q[3];
sx q[3];
rz(-1.3037762) q[3];
sx q[3];
rz(3.0468575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6955473) q[2];
sx q[2];
rz(-1.7798767) q[2];
sx q[2];
rz(2.5496303) q[2];
rz(2.5752318) q[3];
sx q[3];
rz(-0.16470328) q[3];
sx q[3];
rz(-1.5238354) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3175209) q[0];
sx q[0];
rz(-0.98012797) q[0];
sx q[0];
rz(-1.160887) q[0];
rz(3.042165) q[1];
sx q[1];
rz(-1.2482523) q[1];
sx q[1];
rz(-2.0773239) q[1];
rz(0.70384937) q[2];
sx q[2];
rz(-0.8740295) q[2];
sx q[2];
rz(-0.340273) q[2];
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
