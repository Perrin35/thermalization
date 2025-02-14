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
rz(0.0013295833) q[0];
rz(3.3450491) q[1];
sx q[1];
rz(3.3766881) q[1];
sx q[1];
rz(8.8535218) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1648096) q[0];
sx q[0];
rz(-2.342988) q[0];
sx q[0];
rz(2.476506) q[0];
rz(2.9517101) q[2];
sx q[2];
rz(-1.7203091) q[2];
sx q[2];
rz(-2.8854795) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3498342) q[1];
sx q[1];
rz(-1.5877643) q[1];
sx q[1];
rz(2.6512673) q[1];
rz(-pi) q[2];
rz(1.1892849) q[3];
sx q[3];
rz(-1.9982059) q[3];
sx q[3];
rz(-2.1344466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.098238952) q[2];
sx q[2];
rz(-1.2231772) q[2];
sx q[2];
rz(0.60388786) q[2];
rz(2.1308925) q[3];
sx q[3];
rz(-0.5637919) q[3];
sx q[3];
rz(0.92196661) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0826223) q[0];
sx q[0];
rz(-2.048546) q[0];
sx q[0];
rz(3.1329204) q[0];
rz(1.9827093) q[1];
sx q[1];
rz(-0.68717879) q[1];
sx q[1];
rz(-3.029356) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56298577) q[0];
sx q[0];
rz(-1.2845717) q[0];
sx q[0];
rz(2.6335825) q[0];
rz(-pi) q[1];
rz(-0.64666266) q[2];
sx q[2];
rz(-1.5407667) q[2];
sx q[2];
rz(-2.554837) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.018377233) q[1];
sx q[1];
rz(-1.0800414) q[1];
sx q[1];
rz(1.3181376) q[1];
x q[2];
rz(-0.82661622) q[3];
sx q[3];
rz(-0.38020779) q[3];
sx q[3];
rz(0.8687063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.23942854) q[2];
sx q[2];
rz(-1.3210693) q[2];
sx q[2];
rz(-0.93977896) q[2];
rz(0.93722614) q[3];
sx q[3];
rz(-1.3746494) q[3];
sx q[3];
rz(-0.37041131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9409907) q[0];
sx q[0];
rz(-0.90018278) q[0];
sx q[0];
rz(-0.58315939) q[0];
rz(1.0519625) q[1];
sx q[1];
rz(-0.58590341) q[1];
sx q[1];
rz(3.1254056) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0457387) q[0];
sx q[0];
rz(-1.7940137) q[0];
sx q[0];
rz(2.2824817) q[0];
rz(-1.2643686) q[2];
sx q[2];
rz(-1.8548022) q[2];
sx q[2];
rz(-2.4302182) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0250108) q[1];
sx q[1];
rz(-1.7336646) q[1];
sx q[1];
rz(-3.0322269) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3170268) q[3];
sx q[3];
rz(-1.0921794) q[3];
sx q[3];
rz(-0.625911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6215324) q[2];
sx q[2];
rz(-1.3935057) q[2];
sx q[2];
rz(3.0690466) q[2];
rz(1.0091311) q[3];
sx q[3];
rz(-0.79807177) q[3];
sx q[3];
rz(-0.23536853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9811454) q[0];
sx q[0];
rz(-2.7837842) q[0];
sx q[0];
rz(2.2571046) q[0];
rz(-1.6403987) q[1];
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
rz(-0.64393015) q[0];
sx q[0];
rz(0.91899271) q[0];
rz(-1.8461761) q[2];
sx q[2];
rz(-1.7905053) q[2];
sx q[2];
rz(1.7412468) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1780488) q[1];
sx q[1];
rz(-0.60979382) q[1];
sx q[1];
rz(-0.96505537) q[1];
x q[2];
rz(-0.25165148) q[3];
sx q[3];
rz(-0.44969895) q[3];
sx q[3];
rz(-2.3684199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.13714743) q[2];
sx q[2];
rz(-1.659617) q[2];
sx q[2];
rz(1.8972634) q[2];
rz(-1.3589877) q[3];
sx q[3];
rz(-2.6748896) q[3];
sx q[3];
rz(-0.14314237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2537848) q[0];
sx q[0];
rz(-0.65685993) q[0];
sx q[0];
rz(0.3669056) q[0];
rz(1.7985571) q[1];
sx q[1];
rz(-0.29452205) q[1];
sx q[1];
rz(1.3142745) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1444891) q[0];
sx q[0];
rz(-0.95074749) q[0];
sx q[0];
rz(2.117036) q[0];
x q[1];
rz(0.48449202) q[2];
sx q[2];
rz(-1.3132261) q[2];
sx q[2];
rz(0.07312873) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6180743) q[1];
sx q[1];
rz(-2.1879751) q[1];
sx q[1];
rz(-2.7568222) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9461104) q[3];
sx q[3];
rz(-1.2367782) q[3];
sx q[3];
rz(-0.53900146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.542995) q[2];
sx q[2];
rz(-0.44671217) q[2];
sx q[2];
rz(-2.561595) q[2];
rz(0.18481208) q[3];
sx q[3];
rz(-1.5671268) q[3];
sx q[3];
rz(-2.4549761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0762894) q[0];
sx q[0];
rz(-2.6798601) q[0];
sx q[0];
rz(-2.8522458) q[0];
rz(3.0416327) q[1];
sx q[1];
rz(-0.83671612) q[1];
sx q[1];
rz(-0.79773703) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3432309) q[0];
sx q[0];
rz(-1.4360748) q[0];
sx q[0];
rz(0.93905296) q[0];
rz(-pi) q[1];
rz(-1.4580472) q[2];
sx q[2];
rz(-1.7712153) q[2];
sx q[2];
rz(1.3996313) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.575036) q[1];
sx q[1];
rz(-0.29984353) q[1];
sx q[1];
rz(0.17438247) q[1];
rz(-2.3008505) q[3];
sx q[3];
rz(-1.4952097) q[3];
sx q[3];
rz(-0.10160916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.37819698) q[2];
sx q[2];
rz(-1.8969456) q[2];
sx q[2];
rz(2.3217616) q[2];
rz(-1.4929006) q[3];
sx q[3];
rz(-0.35455743) q[3];
sx q[3];
rz(0.41560391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1220658) q[0];
sx q[0];
rz(-2.8391835) q[0];
sx q[0];
rz(0.31461883) q[0];
rz(2.4954691) q[1];
sx q[1];
rz(-1.6982634) q[1];
sx q[1];
rz(3.019928) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80597654) q[0];
sx q[0];
rz(-1.2947963) q[0];
sx q[0];
rz(-1.1742623) q[0];
x q[1];
rz(-2.468277) q[2];
sx q[2];
rz(-2.2229869) q[2];
sx q[2];
rz(-2.7567692) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.104427) q[1];
sx q[1];
rz(-2.011877) q[1];
sx q[1];
rz(1.4250526) q[1];
rz(-0.96243919) q[3];
sx q[3];
rz(-1.4104868) q[3];
sx q[3];
rz(2.1862505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6803711) q[2];
sx q[2];
rz(-1.17522) q[2];
sx q[2];
rz(-2.1628974) q[2];
rz(0.77224246) q[3];
sx q[3];
rz(-2.6148836) q[3];
sx q[3];
rz(-2.1286807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0270281) q[0];
sx q[0];
rz(-1.9975198) q[0];
sx q[0];
rz(1.3042599) q[0];
rz(-0.14106855) q[1];
sx q[1];
rz(-0.27962676) q[1];
sx q[1];
rz(1.8581026) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0373223) q[0];
sx q[0];
rz(-1.3548343) q[0];
sx q[0];
rz(-3.1041935) q[0];
rz(0.097900585) q[2];
sx q[2];
rz(-1.86964) q[2];
sx q[2];
rz(-0.035261957) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.56399689) q[1];
sx q[1];
rz(-1.8297927) q[1];
sx q[1];
rz(1.4513272) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3621534) q[3];
sx q[3];
rz(-1.6044242) q[3];
sx q[3];
rz(-1.5228576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.86113247) q[2];
sx q[2];
rz(-0.91117636) q[2];
sx q[2];
rz(-0.99315161) q[2];
rz(1.447621) q[3];
sx q[3];
rz(-1.2074892) q[3];
sx q[3];
rz(-1.8705468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3250378) q[0];
sx q[0];
rz(-0.10295454) q[0];
sx q[0];
rz(2.2651941) q[0];
rz(-1.5429629) q[1];
sx q[1];
rz(-1.2878659) q[1];
sx q[1];
rz(-1.1184982) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7735104) q[0];
sx q[0];
rz(-0.72282571) q[0];
sx q[0];
rz(0.55814349) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.411941) q[2];
sx q[2];
rz(-2.0579946) q[2];
sx q[2];
rz(1.2318357) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4913018) q[1];
sx q[1];
rz(-2.2660265) q[1];
sx q[1];
rz(-2.9952666) q[1];
x q[2];
rz(-0.38409813) q[3];
sx q[3];
rz(-0.68168655) q[3];
sx q[3];
rz(0.1097928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8886275) q[2];
sx q[2];
rz(-2.1972563) q[2];
sx q[2];
rz(1.8892807) q[2];
rz(1.1014994) q[3];
sx q[3];
rz(-1.7606198) q[3];
sx q[3];
rz(1.291905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.0970704) q[0];
sx q[0];
rz(-0.047053311) q[0];
sx q[0];
rz(0.72730056) q[0];
rz(2.3749088) q[1];
sx q[1];
rz(-0.84696451) q[1];
sx q[1];
rz(1.1904221) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81556126) q[0];
sx q[0];
rz(-2.4056068) q[0];
sx q[0];
rz(-0.18819564) q[0];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.1228635) q[1];
sx q[1];
rz(-1.3236538) q[1];
sx q[1];
rz(0.31063947) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1471858) q[3];
sx q[3];
rz(-2.0032479) q[3];
sx q[3];
rz(-1.9284039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7850354) q[2];
sx q[2];
rz(-1.0871474) q[2];
sx q[2];
rz(-1.9800775) q[2];
rz(1.131743) q[3];
sx q[3];
rz(-0.83860207) q[3];
sx q[3];
rz(2.13818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7575191) q[0];
sx q[0];
rz(-0.90255559) q[0];
sx q[0];
rz(2.3046816) q[0];
rz(1.8867672) q[1];
sx q[1];
rz(-1.8950987) q[1];
sx q[1];
rz(-1.7991039) q[1];
rz(1.1643428) q[2];
sx q[2];
rz(-2.5377574) q[2];
sx q[2];
rz(-0.24518235) q[2];
rz(1.3747207) q[3];
sx q[3];
rz(-1.5846854) q[3];
sx q[3];
rz(-2.9561083) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
