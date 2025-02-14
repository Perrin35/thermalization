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
rz(0.089224815) q[0];
sx q[0];
rz(2.3951946) q[0];
sx q[0];
rz(8.7328773) q[0];
rz(-2.6785985) q[1];
sx q[1];
rz(-1.5249335) q[1];
sx q[1];
rz(-2.1114299) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97447831) q[0];
sx q[0];
rz(-1.8016707) q[0];
sx q[0];
rz(2.5690325) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9857731) q[2];
sx q[2];
rz(-0.43104592) q[2];
sx q[2];
rz(0.31722122) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9889415) q[1];
sx q[1];
rz(-1.2458774) q[1];
sx q[1];
rz(-0.12293651) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.452677) q[3];
sx q[3];
rz(-1.6647881) q[3];
sx q[3];
rz(1.9436364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2618711) q[2];
sx q[2];
rz(-2.3679831) q[2];
sx q[2];
rz(0.36641463) q[2];
rz(2.053818) q[3];
sx q[3];
rz(-0.78947869) q[3];
sx q[3];
rz(-2.4521949) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3159897) q[0];
sx q[0];
rz(-0.082254224) q[0];
sx q[0];
rz(0.91973037) q[0];
rz(-2.3986744) q[1];
sx q[1];
rz(-1.2580322) q[1];
sx q[1];
rz(1.119335) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9539328) q[0];
sx q[0];
rz(-1.7155678) q[0];
sx q[0];
rz(-0.06186084) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3886003) q[2];
sx q[2];
rz(-2.2815653) q[2];
sx q[2];
rz(1.3533392) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.38779681) q[1];
sx q[1];
rz(-0.8555853) q[1];
sx q[1];
rz(-1.254735) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1503937) q[3];
sx q[3];
rz(-2.0127735) q[3];
sx q[3];
rz(0.22721618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4982345) q[2];
sx q[2];
rz(-1.1946119) q[2];
sx q[2];
rz(-0.052113459) q[2];
rz(-2.0335967) q[3];
sx q[3];
rz(-0.85927695) q[3];
sx q[3];
rz(0.064854709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68010083) q[0];
sx q[0];
rz(-1.9100186) q[0];
sx q[0];
rz(1.8517866) q[0];
rz(-0.75311226) q[1];
sx q[1];
rz(-1.6878004) q[1];
sx q[1];
rz(0.49711102) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38579861) q[0];
sx q[0];
rz(-2.442474) q[0];
sx q[0];
rz(-2.505245) q[0];
rz(-2.2154664) q[2];
sx q[2];
rz(-2.1417649) q[2];
sx q[2];
rz(2.1059011) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.063035251) q[1];
sx q[1];
rz(-2.7655256) q[1];
sx q[1];
rz(0.1034115) q[1];
rz(-1.4744954) q[3];
sx q[3];
rz(-0.9203921) q[3];
sx q[3];
rz(1.4092829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6734753) q[2];
sx q[2];
rz(-0.27286369) q[2];
sx q[2];
rz(0.72431481) q[2];
rz(2.916548) q[3];
sx q[3];
rz(-1.2757755) q[3];
sx q[3];
rz(-0.79506522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1555136) q[0];
sx q[0];
rz(-0.88548311) q[0];
sx q[0];
rz(-0.28050637) q[0];
rz(1.4836503) q[1];
sx q[1];
rz(-1.2018485) q[1];
sx q[1];
rz(-0.4712421) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5007809) q[0];
sx q[0];
rz(-2.7735595) q[0];
sx q[0];
rz(-1.5076007) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6837101) q[2];
sx q[2];
rz(-1.5051923) q[2];
sx q[2];
rz(-2.0173617) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.60741456) q[1];
sx q[1];
rz(-1.8414268) q[1];
sx q[1];
rz(3.1128008) q[1];
rz(2.332451) q[3];
sx q[3];
rz(-1.1365206) q[3];
sx q[3];
rz(-0.44444042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9435141) q[2];
sx q[2];
rz(-2.5835865) q[2];
sx q[2];
rz(-2.1264709) q[2];
rz(0.73457581) q[3];
sx q[3];
rz(-1.7769122) q[3];
sx q[3];
rz(-0.53810292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9839142) q[0];
sx q[0];
rz(-1.9240802) q[0];
sx q[0];
rz(-0.96018803) q[0];
rz(3.0958815) q[1];
sx q[1];
rz(-0.8784596) q[1];
sx q[1];
rz(0.033230573) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3094745) q[0];
sx q[0];
rz(-0.72369472) q[0];
sx q[0];
rz(-2.3242818) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8383547) q[2];
sx q[2];
rz(-1.5853527) q[2];
sx q[2];
rz(2.6637258) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7575562) q[1];
sx q[1];
rz(-1.23856) q[1];
sx q[1];
rz(0.13793034) q[1];
rz(-pi) q[2];
rz(-0.73992336) q[3];
sx q[3];
rz(-0.56906521) q[3];
sx q[3];
rz(2.8875433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9783322) q[2];
sx q[2];
rz(-2.1425118) q[2];
sx q[2];
rz(-2.8283289) q[2];
rz(-0.35798171) q[3];
sx q[3];
rz(-2.101818) q[3];
sx q[3];
rz(0.09859214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.98704308) q[0];
sx q[0];
rz(-1.6313666) q[0];
sx q[0];
rz(-2.6556067) q[0];
rz(1.9588574) q[1];
sx q[1];
rz(-0.69570884) q[1];
sx q[1];
rz(2.8114496) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3257623) q[0];
sx q[0];
rz(-1.6005524) q[0];
sx q[0];
rz(-1.8682006) q[0];
rz(-pi) q[1];
rz(-2.3358104) q[2];
sx q[2];
rz(-1.0892765) q[2];
sx q[2];
rz(1.0891506) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.79336004) q[1];
sx q[1];
rz(-2.4168682) q[1];
sx q[1];
rz(-0.7591922) q[1];
rz(0.56056916) q[3];
sx q[3];
rz(-2.7383826) q[3];
sx q[3];
rz(-1.8016658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6351472) q[2];
sx q[2];
rz(-2.5250285) q[2];
sx q[2];
rz(0.35663024) q[2];
rz(0.59456524) q[3];
sx q[3];
rz(-1.6584572) q[3];
sx q[3];
rz(-0.73896343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67721382) q[0];
sx q[0];
rz(-2.3742299) q[0];
sx q[0];
rz(-2.5546524) q[0];
rz(-2.7691973) q[1];
sx q[1];
rz(-1.7814691) q[1];
sx q[1];
rz(0.24758235) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42668396) q[0];
sx q[0];
rz(-1.539402) q[0];
sx q[0];
rz(2.035729) q[0];
rz(1.3273986) q[2];
sx q[2];
rz(-0.95275022) q[2];
sx q[2];
rz(-0.40153654) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.54625) q[1];
sx q[1];
rz(-0.74186013) q[1];
sx q[1];
rz(-2.0764024) q[1];
rz(-0.34981899) q[3];
sx q[3];
rz(-1.3365776) q[3];
sx q[3];
rz(-1.8893591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.26329142) q[2];
sx q[2];
rz(-1.2211439) q[2];
sx q[2];
rz(0.2090052) q[2];
rz(0.59669295) q[3];
sx q[3];
rz(-2.7946819) q[3];
sx q[3];
rz(1.5359115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40808943) q[0];
sx q[0];
rz(-0.96036378) q[0];
sx q[0];
rz(0.39304131) q[0];
rz(2.4709002) q[1];
sx q[1];
rz(-1.9978943) q[1];
sx q[1];
rz(-1.447698) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2358952) q[0];
sx q[0];
rz(-1.6950771) q[0];
sx q[0];
rz(2.0688562) q[0];
x q[1];
rz(-1.1426439) q[2];
sx q[2];
rz(-0.65946666) q[2];
sx q[2];
rz(3.1127549) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5504584) q[1];
sx q[1];
rz(-2.0789642) q[1];
sx q[1];
rz(-0.9891467) q[1];
rz(-2.0441229) q[3];
sx q[3];
rz(-2.3103263) q[3];
sx q[3];
rz(2.748718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4286917) q[2];
sx q[2];
rz(-2.28076) q[2];
sx q[2];
rz(-3.0248896) q[2];
rz(1.0137001) q[3];
sx q[3];
rz(-1.2319535) q[3];
sx q[3];
rz(-1.3056171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.3882465) q[0];
sx q[0];
rz(-1.6755063) q[0];
sx q[0];
rz(1.4353132) q[0];
rz(-2.1460136) q[1];
sx q[1];
rz(-1.3187871) q[1];
sx q[1];
rz(-0.54042712) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89622241) q[0];
sx q[0];
rz(-2.5858552) q[0];
sx q[0];
rz(-2.627719) q[0];
x q[1];
rz(-0.90370205) q[2];
sx q[2];
rz(-2.0533178) q[2];
sx q[2];
rz(2.7806892) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.82481532) q[1];
sx q[1];
rz(-1.9476003) q[1];
sx q[1];
rz(2.3945859) q[1];
rz(-pi) q[2];
rz(1.1630658) q[3];
sx q[3];
rz(-2.1559048) q[3];
sx q[3];
rz(-1.9118904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.871375) q[2];
sx q[2];
rz(-1.8367218) q[2];
sx q[2];
rz(-0.77667856) q[2];
rz(-0.03260472) q[3];
sx q[3];
rz(-1.6083345) q[3];
sx q[3];
rz(1.5672654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8107574) q[0];
sx q[0];
rz(-2.7726655) q[0];
sx q[0];
rz(2.708013) q[0];
rz(2.4724204) q[1];
sx q[1];
rz(-1.9586261) q[1];
sx q[1];
rz(0.42632595) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2992914) q[0];
sx q[0];
rz(-0.29336818) q[0];
sx q[0];
rz(-0.82404739) q[0];
x q[1];
rz(1.2157006) q[2];
sx q[2];
rz(-0.85596687) q[2];
sx q[2];
rz(-2.1260043) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9712228) q[1];
sx q[1];
rz(-2.055238) q[1];
sx q[1];
rz(-1.1722527) q[1];
rz(-pi) q[2];
x q[2];
rz(0.58482184) q[3];
sx q[3];
rz(-1.1556781) q[3];
sx q[3];
rz(2.7266172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7508037) q[2];
sx q[2];
rz(-2.8533253) q[2];
sx q[2];
rz(-2.0354347) q[2];
rz(-0.11354167) q[3];
sx q[3];
rz(-2.2083486) q[3];
sx q[3];
rz(-3.09789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15585598) q[0];
sx q[0];
rz(-2.2423797) q[0];
sx q[0];
rz(2.5806497) q[0];
rz(0.38324311) q[1];
sx q[1];
rz(-1.1226729) q[1];
sx q[1];
rz(2.9164006) q[1];
rz(1.2651934) q[2];
sx q[2];
rz(-1.0355179) q[2];
sx q[2];
rz(-2.1260362) q[2];
rz(0.44456496) q[3];
sx q[3];
rz(-0.86935432) q[3];
sx q[3];
rz(-1.0039734) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
