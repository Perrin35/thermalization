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
rz(0.94533935) q[0];
sx q[0];
rz(9.4234484) q[0];
rz(3.3450491) q[1];
sx q[1];
rz(3.3766881) q[1];
sx q[1];
rz(8.8535218) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1648096) q[0];
sx q[0];
rz(-0.79860461) q[0];
sx q[0];
rz(2.476506) q[0];
rz(1.4185888) q[2];
sx q[2];
rz(-1.3830575) q[2];
sx q[2];
rz(-1.3433045) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9524051) q[1];
sx q[1];
rz(-2.6509977) q[1];
sx q[1];
rz(3.1055727) q[1];
rz(-pi) q[2];
x q[2];
rz(0.68526973) q[3];
sx q[3];
rz(-2.5766386) q[3];
sx q[3];
rz(-2.9034815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0433537) q[2];
sx q[2];
rz(-1.2231772) q[2];
sx q[2];
rz(-0.60388786) q[2];
rz(-1.0107001) q[3];
sx q[3];
rz(-0.5637919) q[3];
sx q[3];
rz(0.92196661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0826223) q[0];
sx q[0];
rz(-1.0930467) q[0];
sx q[0];
rz(-3.1329204) q[0];
rz(1.1588833) q[1];
sx q[1];
rz(-0.68717879) q[1];
sx q[1];
rz(3.029356) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56298577) q[0];
sx q[0];
rz(-1.2845717) q[0];
sx q[0];
rz(-0.5080101) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5331763) q[2];
sx q[2];
rz(-0.92447399) q[2];
sx q[2];
rz(0.96137709) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4680926) q[1];
sx q[1];
rz(-1.7930987) q[1];
sx q[1];
rz(-2.6373074) q[1];
rz(2.3149764) q[3];
sx q[3];
rz(-0.38020779) q[3];
sx q[3];
rz(-2.2728863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9021641) q[2];
sx q[2];
rz(-1.8205234) q[2];
sx q[2];
rz(0.93977896) q[2];
rz(-2.2043665) q[3];
sx q[3];
rz(-1.3746494) q[3];
sx q[3];
rz(2.7711813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20060191) q[0];
sx q[0];
rz(-0.90018278) q[0];
sx q[0];
rz(-0.58315939) q[0];
rz(2.0896301) q[1];
sx q[1];
rz(-2.5556892) q[1];
sx q[1];
rz(-0.016187035) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0457387) q[0];
sx q[0];
rz(-1.7940137) q[0];
sx q[0];
rz(-0.8591109) q[0];
rz(-0.80183713) q[2];
sx q[2];
rz(-2.7268598) q[2];
sx q[2];
rz(1.557359) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5213373) q[1];
sx q[1];
rz(-2.9456821) q[1];
sx q[1];
rz(-2.1571062) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2233769) q[3];
sx q[3];
rz(-0.86091061) q[3];
sx q[3];
rz(-0.48274279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6215324) q[2];
sx q[2];
rz(-1.748087) q[2];
sx q[2];
rz(3.0690466) q[2];
rz(-1.0091311) q[3];
sx q[3];
rz(-2.3435209) q[3];
sx q[3];
rz(-0.23536853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9811454) q[0];
sx q[0];
rz(-0.35780847) q[0];
sx q[0];
rz(0.88448802) q[0];
rz(1.6403987) q[1];
sx q[1];
rz(-1.6694992) q[1];
sx q[1];
rz(2.5515058) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5948759) q[0];
sx q[0];
rz(-1.1980431) q[0];
sx q[0];
rz(1.0327505) q[0];
rz(-1.8461761) q[2];
sx q[2];
rz(-1.7905053) q[2];
sx q[2];
rz(1.7412468) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1780488) q[1];
sx q[1];
rz(-0.60979382) q[1];
sx q[1];
rz(0.96505537) q[1];
rz(-2.7042974) q[3];
sx q[3];
rz(-1.6792494) q[3];
sx q[3];
rz(2.5714998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0044452) q[2];
sx q[2];
rz(-1.4819757) q[2];
sx q[2];
rz(1.8972634) q[2];
rz(1.7826049) q[3];
sx q[3];
rz(-2.6748896) q[3];
sx q[3];
rz(2.9984503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2537848) q[0];
sx q[0];
rz(-2.4847327) q[0];
sx q[0];
rz(2.7746871) q[0];
rz(-1.7985571) q[1];
sx q[1];
rz(-0.29452205) q[1];
sx q[1];
rz(1.8273182) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9132694) q[0];
sx q[0];
rz(-1.1342888) q[0];
sx q[0];
rz(-0.69605791) q[0];
rz(0.48449202) q[2];
sx q[2];
rz(-1.3132261) q[2];
sx q[2];
rz(-3.0684639) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.81707047) q[1];
sx q[1];
rz(-1.2597023) q[1];
sx q[1];
rz(2.2242311) q[1];
rz(-2.3288127) q[3];
sx q[3];
rz(-0.4970937) q[3];
sx q[3];
rz(-1.7258096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.542995) q[2];
sx q[2];
rz(-2.6948805) q[2];
sx q[2];
rz(2.561595) q[2];
rz(-0.18481208) q[3];
sx q[3];
rz(-1.5744659) q[3];
sx q[3];
rz(-2.4549761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0762894) q[0];
sx q[0];
rz(-0.46173254) q[0];
sx q[0];
rz(2.8522458) q[0];
rz(3.0416327) q[1];
sx q[1];
rz(-0.83671612) q[1];
sx q[1];
rz(-0.79773703) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3432309) q[0];
sx q[0];
rz(-1.7055178) q[0];
sx q[0];
rz(-2.2025397) q[0];
x q[1];
rz(-0.20166534) q[2];
sx q[2];
rz(-1.4603134) q[2];
sx q[2];
rz(0.1486272) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9790837) q[1];
sx q[1];
rz(-1.519527) q[1];
sx q[1];
rz(-2.8460345) q[1];
rz(3.0403071) q[3];
sx q[3];
rz(-2.2982979) q[3];
sx q[3];
rz(-1.6049214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7633957) q[2];
sx q[2];
rz(-1.8969456) q[2];
sx q[2];
rz(-2.3217616) q[2];
rz(1.648692) q[3];
sx q[3];
rz(-2.7870352) q[3];
sx q[3];
rz(-0.41560391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1220658) q[0];
sx q[0];
rz(-0.30240914) q[0];
sx q[0];
rz(0.31461883) q[0];
rz(-2.4954691) q[1];
sx q[1];
rz(-1.4433292) q[1];
sx q[1];
rz(3.019928) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1878649) q[0];
sx q[0];
rz(-0.47891579) q[0];
sx q[0];
rz(2.2035416) q[0];
rz(2.3444885) q[2];
sx q[2];
rz(-1.0522763) q[2];
sx q[2];
rz(2.4064877) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.77326425) q[1];
sx q[1];
rz(-2.6785637) q[1];
sx q[1];
rz(-0.29843389) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2950778) q[3];
sx q[3];
rz(-0.62653095) q[3];
sx q[3];
rz(2.3008418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6803711) q[2];
sx q[2];
rz(-1.9663726) q[2];
sx q[2];
rz(2.1628974) q[2];
rz(2.3693502) q[3];
sx q[3];
rz(-2.6148836) q[3];
sx q[3];
rz(2.1286807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.0270281) q[0];
sx q[0];
rz(-1.1440729) q[0];
sx q[0];
rz(1.3042599) q[0];
rz(-0.14106855) q[1];
sx q[1];
rz(-2.8619659) q[1];
sx q[1];
rz(-1.8581026) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.210189) q[0];
sx q[0];
rz(-2.922466) q[0];
sx q[0];
rz(-1.7396084) q[0];
x q[1];
rz(-1.8709932) q[2];
sx q[2];
rz(-1.4772479) q[2];
sx q[2];
rz(1.5771505) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5775958) q[1];
sx q[1];
rz(-1.3118) q[1];
sx q[1];
rz(1.4513272) q[1];
x q[2];
rz(0.77943927) q[3];
sx q[3];
rz(-1.6044242) q[3];
sx q[3];
rz(-1.618735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2804602) q[2];
sx q[2];
rz(-0.91117636) q[2];
sx q[2];
rz(-0.99315161) q[2];
rz(1.6939717) q[3];
sx q[3];
rz(-1.9341035) q[3];
sx q[3];
rz(-1.8705468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(1.3250378) q[0];
sx q[0];
rz(-0.10295454) q[0];
sx q[0];
rz(-0.8763985) q[0];
rz(-1.5429629) q[1];
sx q[1];
rz(-1.2878659) q[1];
sx q[1];
rz(-1.1184982) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0792313) q[0];
sx q[0];
rz(-0.97506279) q[0];
sx q[0];
rz(1.1337639) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1886979) q[2];
sx q[2];
rz(-0.94091533) q[2];
sx q[2];
rz(-0.057531051) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.42414819) q[1];
sx q[1];
rz(-2.433648) q[1];
sx q[1];
rz(-1.397754) q[1];
rz(-2.4966024) q[3];
sx q[3];
rz(-1.8091615) q[3];
sx q[3];
rz(-1.1568943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8886275) q[2];
sx q[2];
rz(-2.1972563) q[2];
sx q[2];
rz(1.8892807) q[2];
rz(-2.0400932) q[3];
sx q[3];
rz(-1.3809729) q[3];
sx q[3];
rz(1.8496877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0445223) q[0];
sx q[0];
rz(-0.047053311) q[0];
sx q[0];
rz(2.4142921) q[0];
rz(-2.3749088) q[1];
sx q[1];
rz(-2.2946281) q[1];
sx q[1];
rz(1.1904221) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81556126) q[0];
sx q[0];
rz(-0.73598586) q[0];
sx q[0];
rz(-0.18819564) q[0];
rz(3.0465428) q[2];
sx q[2];
rz(-2.2288786) q[2];
sx q[2];
rz(2.6320145) q[2];
x q[3];
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
rz(-pi) q[2];
x q[2];
rz(0.86807172) q[3];
sx q[3];
rz(-0.70555726) q[3];
sx q[3];
rz(-0.21506532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.35655725) q[2];
sx q[2];
rz(-2.0544453) q[2];
sx q[2];
rz(1.9800775) q[2];
rz(2.0098497) q[3];
sx q[3];
rz(-0.83860207) q[3];
sx q[3];
rz(1.0034126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7575191) q[0];
sx q[0];
rz(-2.2390371) q[0];
sx q[0];
rz(-0.83691103) q[0];
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
rz(-1.7668719) q[3];
sx q[3];
rz(-1.5846854) q[3];
sx q[3];
rz(-2.9561083) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
