OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.43006858) q[0];
sx q[0];
rz(-3.0741337) q[0];
sx q[0];
rz(2.467632) q[0];
rz(-0.31710467) q[1];
sx q[1];
rz(-1.6333406) q[1];
sx q[1];
rz(2.34692) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3622409) q[0];
sx q[0];
rz(-1.2138214) q[0];
sx q[0];
rz(-2.7428521) q[0];
rz(-pi) q[1];
rz(2.652466) q[2];
sx q[2];
rz(-1.7009652) q[2];
sx q[2];
rz(-0.39718539) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0920081) q[1];
sx q[1];
rz(-1.9132179) q[1];
sx q[1];
rz(-0.53978668) q[1];
rz(-pi) q[2];
rz(-1.2875597) q[3];
sx q[3];
rz(-2.6784416) q[3];
sx q[3];
rz(-1.8739665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.47544605) q[2];
sx q[2];
rz(-2.011516) q[2];
sx q[2];
rz(1.3249935) q[2];
rz(0.31630668) q[3];
sx q[3];
rz(-0.91747147) q[3];
sx q[3];
rz(-0.90707183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.9544202) q[0];
sx q[0];
rz(-0.81048727) q[0];
sx q[0];
rz(2.1818838) q[0];
rz(-1.6540487) q[1];
sx q[1];
rz(-0.68351775) q[1];
sx q[1];
rz(2.5141292) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7213631) q[0];
sx q[0];
rz(-3.1313167) q[0];
sx q[0];
rz(-0.72364877) q[0];
x q[1];
rz(0.10208315) q[2];
sx q[2];
rz(-0.21271579) q[2];
sx q[2];
rz(1.724147) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.97034772) q[1];
sx q[1];
rz(-1.3309609) q[1];
sx q[1];
rz(-0.14648267) q[1];
x q[2];
rz(-2.0262358) q[3];
sx q[3];
rz(-3.1133828) q[3];
sx q[3];
rz(1.3970323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.54005694) q[2];
sx q[2];
rz(-1.4011708) q[2];
sx q[2];
rz(-2.1616948) q[2];
rz(1.4510138) q[3];
sx q[3];
rz(-2.4308949) q[3];
sx q[3];
rz(-1.6987945) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5805805) q[0];
sx q[0];
rz(-3.09364) q[0];
sx q[0];
rz(1.7618435) q[0];
rz(-2.9648932) q[1];
sx q[1];
rz(-1.9243762) q[1];
sx q[1];
rz(-1.3476936) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4189258) q[0];
sx q[0];
rz(-1.8113266) q[0];
sx q[0];
rz(-0.47970432) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.548544) q[2];
sx q[2];
rz(-1.3705727) q[2];
sx q[2];
rz(-0.74860886) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.99464097) q[1];
sx q[1];
rz(-2.9839582) q[1];
sx q[1];
rz(-0.68871246) q[1];
x q[2];
rz(1.755581) q[3];
sx q[3];
rz(-1.4818958) q[3];
sx q[3];
rz(-0.88528663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1220864) q[2];
sx q[2];
rz(-2.2821189) q[2];
sx q[2];
rz(1.9879831) q[2];
rz(2.7518318) q[3];
sx q[3];
rz(-0.4699769) q[3];
sx q[3];
rz(-2.3526758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
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
rz(-2.6692114) q[0];
sx q[0];
rz(-2.5166469) q[0];
sx q[0];
rz(2.5770082) q[0];
rz(3.0584884) q[1];
sx q[1];
rz(-1.233498) q[1];
sx q[1];
rz(0.21534236) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4190061) q[0];
sx q[0];
rz(-0.79918396) q[0];
sx q[0];
rz(1.0976085) q[0];
rz(-1.3232857) q[2];
sx q[2];
rz(-0.60192054) q[2];
sx q[2];
rz(-1.9844696) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.035186471) q[1];
sx q[1];
rz(-1.3760929) q[1];
sx q[1];
rz(0.10722843) q[1];
rz(-2.5921949) q[3];
sx q[3];
rz(-1.975172) q[3];
sx q[3];
rz(0.79493633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5061491) q[2];
sx q[2];
rz(-2.052867) q[2];
sx q[2];
rz(-2.5396458) q[2];
rz(-0.69058949) q[3];
sx q[3];
rz(-0.75406051) q[3];
sx q[3];
rz(2.5155892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85201207) q[0];
sx q[0];
rz(-2.1471725) q[0];
sx q[0];
rz(-2.1767298) q[0];
rz(3.124974) q[1];
sx q[1];
rz(-2.7162572) q[1];
sx q[1];
rz(-0.66666493) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4283838) q[0];
sx q[0];
rz(-2.7292626) q[0];
sx q[0];
rz(1.8250188) q[0];
rz(-1.0275074) q[2];
sx q[2];
rz(-0.61858656) q[2];
sx q[2];
rz(-2.0999883) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.95496817) q[1];
sx q[1];
rz(-0.72997) q[1];
sx q[1];
rz(-0.43220046) q[1];
rz(-pi) q[2];
rz(-0.22294754) q[3];
sx q[3];
rz(-1.2811536) q[3];
sx q[3];
rz(1.851561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.91606402) q[2];
sx q[2];
rz(-1.7614044) q[2];
sx q[2];
rz(1.2529681) q[2];
rz(-1.3145674) q[3];
sx q[3];
rz(-1.1809177) q[3];
sx q[3];
rz(1.130828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1154293) q[0];
sx q[0];
rz(-2.207343) q[0];
sx q[0];
rz(1.3622267) q[0];
rz(-1.399614) q[1];
sx q[1];
rz(-1.5796966) q[1];
sx q[1];
rz(-1.5302352) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6513072) q[0];
sx q[0];
rz(-1.4978652) q[0];
sx q[0];
rz(-0.42071995) q[0];
x q[1];
rz(2.7026664) q[2];
sx q[2];
rz(-1.8669087) q[2];
sx q[2];
rz(0.93026464) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6758319) q[1];
sx q[1];
rz(-2.6869046) q[1];
sx q[1];
rz(-2.0783706) q[1];
x q[2];
rz(1.8239162) q[3];
sx q[3];
rz(-1.5217921) q[3];
sx q[3];
rz(-0.24310902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7851012) q[2];
sx q[2];
rz(-2.1337528) q[2];
sx q[2];
rz(-3.1121758) q[2];
rz(-1.4908837) q[3];
sx q[3];
rz(-0.95508948) q[3];
sx q[3];
rz(-1.9194549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-2.4239663) q[0];
sx q[0];
rz(-2.3762149) q[0];
sx q[0];
rz(0.18280612) q[0];
rz(0.43552661) q[1];
sx q[1];
rz(-2.2874449) q[1];
sx q[1];
rz(-0.51087728) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8683838) q[0];
sx q[0];
rz(-1.8330935) q[0];
sx q[0];
rz(1.457085) q[0];
rz(-pi) q[1];
rz(-2.4099318) q[2];
sx q[2];
rz(-2.1795142) q[2];
sx q[2];
rz(-0.70639709) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6888914) q[1];
sx q[1];
rz(-2.211262) q[1];
sx q[1];
rz(-2.074261) q[1];
rz(-pi) q[2];
rz(1.0484344) q[3];
sx q[3];
rz(-1.2921385) q[3];
sx q[3];
rz(3.0184961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2304948) q[2];
sx q[2];
rz(-0.7926597) q[2];
sx q[2];
rz(0.62057173) q[2];
rz(-0.42256045) q[3];
sx q[3];
rz(-1.7812604) q[3];
sx q[3];
rz(-1.871199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34185103) q[0];
sx q[0];
rz(-2.7809679) q[0];
sx q[0];
rz(-2.0525232) q[0];
rz(-2.0607121) q[1];
sx q[1];
rz(-1.2842478) q[1];
sx q[1];
rz(-0.63527766) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95406336) q[0];
sx q[0];
rz(-0.39217338) q[0];
sx q[0];
rz(-1.9847045) q[0];
rz(-2.9674171) q[2];
sx q[2];
rz(-1.0798228) q[2];
sx q[2];
rz(-0.27506405) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.63102555) q[1];
sx q[1];
rz(-2.4706565) q[1];
sx q[1];
rz(-0.4354233) q[1];
x q[2];
rz(2.1463296) q[3];
sx q[3];
rz(-1.4008153) q[3];
sx q[3];
rz(3.1135524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.13715956) q[2];
sx q[2];
rz(-1.6798423) q[2];
sx q[2];
rz(-3.1414202) q[2];
rz(-0.6862644) q[3];
sx q[3];
rz(-2.4133447) q[3];
sx q[3];
rz(-2.2837158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8453318) q[0];
sx q[0];
rz(-2.2289472) q[0];
sx q[0];
rz(-2.4849179) q[0];
rz(-0.36390057) q[1];
sx q[1];
rz(-2.6952422) q[1];
sx q[1];
rz(-0.90973967) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1014935) q[0];
sx q[0];
rz(-1.3001406) q[0];
sx q[0];
rz(2.8069242) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7934974) q[2];
sx q[2];
rz(-1.3653792) q[2];
sx q[2];
rz(1.5945438) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7418356) q[1];
sx q[1];
rz(-1.6941119) q[1];
sx q[1];
rz(-2.9305305) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2551366) q[3];
sx q[3];
rz(-1.3254032) q[3];
sx q[3];
rz(2.0592225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3108814) q[2];
sx q[2];
rz(-2.723366) q[2];
sx q[2];
rz(2.6549784) q[2];
rz(3.0330372) q[3];
sx q[3];
rz(-2.4311549) q[3];
sx q[3];
rz(1.9177115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.20031032) q[0];
sx q[0];
rz(-1.0815942) q[0];
sx q[0];
rz(0.51668984) q[0];
rz(-1.5962881) q[1];
sx q[1];
rz(-2.089112) q[1];
sx q[1];
rz(-0.89458481) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54297011) q[0];
sx q[0];
rz(-0.93790903) q[0];
sx q[0];
rz(-3.0892239) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.20608979) q[2];
sx q[2];
rz(-2.0317674) q[2];
sx q[2];
rz(3.0989625) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.91720944) q[1];
sx q[1];
rz(-1.50748) q[1];
sx q[1];
rz(1.9320095) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.92071269) q[3];
sx q[3];
rz(-1.4462785) q[3];
sx q[3];
rz(-1.6684106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9937925) q[2];
sx q[2];
rz(-2.1788969) q[2];
sx q[2];
rz(-1.0220035) q[2];
rz(1.3379124) q[3];
sx q[3];
rz(-1.6513848) q[3];
sx q[3];
rz(0.62906229) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.548303) q[0];
sx q[0];
rz(-1.6175445) q[0];
sx q[0];
rz(2.2486726) q[0];
rz(1.0409566) q[1];
sx q[1];
rz(-0.092408471) q[1];
sx q[1];
rz(1.668781) q[1];
rz(-2.2330166) q[2];
sx q[2];
rz(-0.77975811) q[2];
sx q[2];
rz(-2.1644885) q[2];
rz(-0.048687497) q[3];
sx q[3];
rz(-1.8713453) q[3];
sx q[3];
rz(0.71181675) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
