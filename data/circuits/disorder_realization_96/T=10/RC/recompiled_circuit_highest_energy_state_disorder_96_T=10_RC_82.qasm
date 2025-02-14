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
rz(0.32956707) q[0];
sx q[0];
rz(-2.1962533) q[0];
sx q[0];
rz(-3.1402631) q[0];
rz(3.3450491) q[1];
sx q[1];
rz(3.3766881) q[1];
sx q[1];
rz(8.8535218) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.093319915) q[0];
sx q[0];
rz(-1.1128582) q[0];
sx q[0];
rz(0.67955525) q[0];
rz(-pi) q[1];
x q[1];
rz(0.18988256) q[2];
sx q[2];
rz(-1.7203091) q[2];
sx q[2];
rz(-0.25611311) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.23001901) q[1];
sx q[1];
rz(-2.0610448) q[1];
sx q[1];
rz(-1.5515627) q[1];
rz(-pi) q[2];
rz(2.4563229) q[3];
sx q[3];
rz(-2.5766386) q[3];
sx q[3];
rz(2.9034815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.098238952) q[2];
sx q[2];
rz(-1.2231772) q[2];
sx q[2];
rz(-2.5377048) q[2];
rz(-2.1308925) q[3];
sx q[3];
rz(-2.5778008) q[3];
sx q[3];
rz(0.92196661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0826223) q[0];
sx q[0];
rz(-2.048546) q[0];
sx q[0];
rz(-0.0086722886) q[0];
rz(1.9827093) q[1];
sx q[1];
rz(-0.68717879) q[1];
sx q[1];
rz(-3.029356) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6644635) q[0];
sx q[0];
rz(-0.57687974) q[0];
sx q[0];
rz(2.5974899) q[0];
x q[1];
rz(1.6084163) q[2];
sx q[2];
rz(-0.92447399) q[2];
sx q[2];
rz(2.1802156) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1232154) q[1];
sx q[1];
rz(-1.0800414) q[1];
sx q[1];
rz(1.823455) q[1];
rz(-pi) q[2];
rz(1.2848507) q[3];
sx q[3];
rz(-1.3166898) q[3];
sx q[3];
rz(-1.7320964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.23942854) q[2];
sx q[2];
rz(-1.3210693) q[2];
sx q[2];
rz(-0.93977896) q[2];
rz(-0.93722614) q[3];
sx q[3];
rz(-1.7669433) q[3];
sx q[3];
rz(-0.37041131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20060191) q[0];
sx q[0];
rz(-2.2414099) q[0];
sx q[0];
rz(-2.5584333) q[0];
rz(1.0519625) q[1];
sx q[1];
rz(-2.5556892) q[1];
sx q[1];
rz(0.016187035) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4154041) q[0];
sx q[0];
rz(-0.73999087) q[0];
sx q[0];
rz(-1.9053024) q[0];
rz(2.8444958) q[2];
sx q[2];
rz(-1.8645745) q[2];
sx q[2];
rz(2.370594) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.52798277) q[1];
sx q[1];
rz(-1.6787091) q[1];
sx q[1];
rz(-1.4069665) q[1];
rz(-0.82456581) q[3];
sx q[3];
rz(-2.0494132) q[3];
sx q[3];
rz(2.5156817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5200603) q[2];
sx q[2];
rz(-1.3935057) q[2];
sx q[2];
rz(0.072546093) q[2];
rz(-2.1324615) q[3];
sx q[3];
rz(-2.3435209) q[3];
sx q[3];
rz(0.23536853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16044727) q[0];
sx q[0];
rz(-0.35780847) q[0];
sx q[0];
rz(2.2571046) q[0];
rz(1.5011939) q[1];
sx q[1];
rz(-1.6694992) q[1];
sx q[1];
rz(0.59008682) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6654331) q[0];
sx q[0];
rz(-2.4976625) q[0];
sx q[0];
rz(2.2225999) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2954166) q[2];
sx q[2];
rz(-1.3510873) q[2];
sx q[2];
rz(-1.7412468) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6650538) q[1];
sx q[1];
rz(-1.0805942) q[1];
sx q[1];
rz(-0.3785822) q[1];
rz(-pi) q[2];
x q[2];
rz(1.45118) q[3];
sx q[3];
rz(-2.0053468) q[3];
sx q[3];
rz(1.0512607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0044452) q[2];
sx q[2];
rz(-1.4819757) q[2];
sx q[2];
rz(1.8972634) q[2];
rz(1.7826049) q[3];
sx q[3];
rz(-0.46670306) q[3];
sx q[3];
rz(-2.9984503) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2537848) q[0];
sx q[0];
rz(-0.65685993) q[0];
sx q[0];
rz(-2.7746871) q[0];
rz(-1.7985571) q[1];
sx q[1];
rz(-2.8470706) q[1];
sx q[1];
rz(-1.8273182) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99710354) q[0];
sx q[0];
rz(-0.95074749) q[0];
sx q[0];
rz(-1.0245566) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2814683) q[2];
sx q[2];
rz(-2.0380033) q[2];
sx q[2];
rz(-1.510646) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6180743) q[1];
sx q[1];
rz(-0.95361751) q[1];
sx q[1];
rz(-2.7568222) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.81277992) q[3];
sx q[3];
rz(-0.4970937) q[3];
sx q[3];
rz(1.7258096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.542995) q[2];
sx q[2];
rz(-0.44671217) q[2];
sx q[2];
rz(-0.57999769) q[2];
rz(-2.9567806) q[3];
sx q[3];
rz(-1.5671268) q[3];
sx q[3];
rz(-2.4549761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.065303236) q[0];
sx q[0];
rz(-2.6798601) q[0];
sx q[0];
rz(2.8522458) q[0];
rz(0.099960001) q[1];
sx q[1];
rz(-0.83671612) q[1];
sx q[1];
rz(-2.3438556) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.046043175) q[0];
sx q[0];
rz(-0.64402295) q[0];
sx q[0];
rz(1.3451856) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6835454) q[2];
sx q[2];
rz(-1.3703773) q[2];
sx q[2];
rz(-1.7419614) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.575036) q[1];
sx q[1];
rz(-2.8417491) q[1];
sx q[1];
rz(-0.17438247) q[1];
rz(-pi) q[2];
rz(1.6838668) q[3];
sx q[3];
rz(-0.73323876) q[3];
sx q[3];
rz(-1.3850141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.37819698) q[2];
sx q[2];
rz(-1.244647) q[2];
sx q[2];
rz(-0.81983105) q[2];
rz(1.4929006) q[3];
sx q[3];
rz(-2.7870352) q[3];
sx q[3];
rz(0.41560391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.019526871) q[0];
sx q[0];
rz(-0.30240914) q[0];
sx q[0];
rz(-2.8269738) q[0];
rz(-0.64612359) q[1];
sx q[1];
rz(-1.4433292) q[1];
sx q[1];
rz(-3.019928) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1878649) q[0];
sx q[0];
rz(-0.47891579) q[0];
sx q[0];
rz(-2.2035416) q[0];
x q[1];
rz(0.7971042) q[2];
sx q[2];
rz(-2.0893163) q[2];
sx q[2];
rz(2.4064877) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3683284) q[1];
sx q[1];
rz(-0.46302893) q[1];
sx q[1];
rz(0.29843389) q[1];
rz(-pi) q[2];
rz(-2.9470351) q[3];
sx q[3];
rz(-0.97134198) q[3];
sx q[3];
rz(-0.50473467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6803711) q[2];
sx q[2];
rz(-1.9663726) q[2];
sx q[2];
rz(-2.1628974) q[2];
rz(2.3693502) q[3];
sx q[3];
rz(-0.52670908) q[3];
sx q[3];
rz(-2.1286807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.0270281) q[0];
sx q[0];
rz(-1.1440729) q[0];
sx q[0];
rz(-1.8373328) q[0];
rz(0.14106855) q[1];
sx q[1];
rz(-2.8619659) q[1];
sx q[1];
rz(-1.28349) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9314037) q[0];
sx q[0];
rz(-0.21912665) q[0];
sx q[0];
rz(-1.7396084) q[0];
rz(1.8709932) q[2];
sx q[2];
rz(-1.4772479) q[2];
sx q[2];
rz(-1.5771505) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5775958) q[1];
sx q[1];
rz(-1.3118) q[1];
sx q[1];
rz(1.6902655) q[1];
rz(-1.5235376) q[3];
sx q[3];
rz(-0.79191557) q[3];
sx q[3];
rz(3.126865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2804602) q[2];
sx q[2];
rz(-0.91117636) q[2];
sx q[2];
rz(2.148441) q[2];
rz(-1.447621) q[3];
sx q[3];
rz(-1.2074892) q[3];
sx q[3];
rz(1.8705468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3250378) q[0];
sx q[0];
rz(-3.0386381) q[0];
sx q[0];
rz(-2.2651941) q[0];
rz(1.5986298) q[1];
sx q[1];
rz(-1.8537268) q[1];
sx q[1];
rz(1.1184982) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9063909) q[0];
sx q[0];
rz(-1.9287325) q[0];
sx q[0];
rz(-0.64235781) q[0];
x q[1];
rz(2.411941) q[2];
sx q[2];
rz(-2.0579946) q[2];
sx q[2];
rz(1.9097569) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4913018) q[1];
sx q[1];
rz(-0.87556616) q[1];
sx q[1];
rz(-0.14632605) q[1];
x q[2];
rz(-2.4966024) q[3];
sx q[3];
rz(-1.3324311) q[3];
sx q[3];
rz(1.1568943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2529651) q[2];
sx q[2];
rz(-2.1972563) q[2];
sx q[2];
rz(-1.2523119) q[2];
rz(2.0400932) q[3];
sx q[3];
rz(-1.3809729) q[3];
sx q[3];
rz(-1.8496877) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0445223) q[0];
sx q[0];
rz(-0.047053311) q[0];
sx q[0];
rz(-0.72730056) q[0];
rz(-0.76668382) q[1];
sx q[1];
rz(-0.84696451) q[1];
sx q[1];
rz(-1.9511706) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2461288) q[0];
sx q[0];
rz(-1.4448691) q[0];
sx q[0];
rz(-0.7271304) q[0];
rz(-pi) q[1];
x q[1];
rz(0.095049871) q[2];
sx q[2];
rz(-2.2288786) q[2];
sx q[2];
rz(0.50957818) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3695581) q[1];
sx q[1];
rz(-1.8716964) q[1];
sx q[1];
rz(1.8298261) q[1];
rz(2.2735209) q[3];
sx q[3];
rz(-2.4360354) q[3];
sx q[3];
rz(2.9265273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.35655725) q[2];
sx q[2];
rz(-1.0871474) q[2];
sx q[2];
rz(-1.9800775) q[2];
rz(-1.131743) q[3];
sx q[3];
rz(-0.83860207) q[3];
sx q[3];
rz(-2.13818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3840735) q[0];
sx q[0];
rz(-2.2390371) q[0];
sx q[0];
rz(-0.83691103) q[0];
rz(1.8867672) q[1];
sx q[1];
rz(-1.8950987) q[1];
sx q[1];
rz(-1.7991039) q[1];
rz(-2.135545) q[2];
sx q[2];
rz(-1.7972094) q[2];
sx q[2];
rz(-2.1564855) q[2];
rz(-1.6419717) q[3];
sx q[3];
rz(-2.945032) q[3];
sx q[3];
rz(-1.3155027) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
