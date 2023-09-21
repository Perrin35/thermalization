OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.5716612) q[0];
sx q[0];
rz(-0.067458955) q[0];
sx q[0];
rz(10.098739) q[0];
rz(-0.31710467) q[1];
sx q[1];
rz(-1.6333406) q[1];
sx q[1];
rz(2.34692) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7793517) q[0];
sx q[0];
rz(-1.2138214) q[0];
sx q[0];
rz(-0.39874052) q[0];
rz(-pi) q[1];
x q[1];
rz(0.27172471) q[2];
sx q[2];
rz(-0.50479111) q[2];
sx q[2];
rz(1.4128078) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.049584576) q[1];
sx q[1];
rz(-1.9132179) q[1];
sx q[1];
rz(0.53978668) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2875597) q[3];
sx q[3];
rz(-2.6784416) q[3];
sx q[3];
rz(-1.2676261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.47544605) q[2];
sx q[2];
rz(-1.1300766) q[2];
sx q[2];
rz(-1.3249935) q[2];
rz(-0.31630668) q[3];
sx q[3];
rz(-0.91747147) q[3];
sx q[3];
rz(0.90707183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1871724) q[0];
sx q[0];
rz(-0.81048727) q[0];
sx q[0];
rz(2.1818838) q[0];
rz(-1.487544) q[1];
sx q[1];
rz(-2.4580749) q[1];
sx q[1];
rz(-0.62746343) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42022959) q[0];
sx q[0];
rz(-0.010275928) q[0];
sx q[0];
rz(-2.4179439) q[0];
rz(0.10208315) q[2];
sx q[2];
rz(-0.21271579) q[2];
sx q[2];
rz(-1.4174457) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5761766) q[1];
sx q[1];
rz(-1.7130573) q[1];
sx q[1];
rz(1.8131282) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1291817) q[3];
sx q[3];
rz(-1.5454626) q[3];
sx q[3];
rz(2.2001571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6015357) q[2];
sx q[2];
rz(-1.4011708) q[2];
sx q[2];
rz(2.1616948) q[2];
rz(1.6905789) q[3];
sx q[3];
rz(-2.4308949) q[3];
sx q[3];
rz(1.6987945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56101218) q[0];
sx q[0];
rz(-0.0479527) q[0];
sx q[0];
rz(-1.3797492) q[0];
rz(0.17669949) q[1];
sx q[1];
rz(-1.9243762) q[1];
sx q[1];
rz(-1.3476936) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.028570024) q[0];
sx q[0];
rz(-1.1060113) q[0];
sx q[0];
rz(-1.8405429) q[0];
x q[1];
rz(-1.5930487) q[2];
sx q[2];
rz(-1.3705727) q[2];
sx q[2];
rz(0.74860886) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8417931) q[1];
sx q[1];
rz(-1.6922957) q[1];
sx q[1];
rz(-1.6714765) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3860116) q[3];
sx q[3];
rz(-1.6596969) q[3];
sx q[3];
rz(-2.256306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0195062) q[2];
sx q[2];
rz(-2.2821189) q[2];
sx q[2];
rz(1.9879831) q[2];
rz(-2.7518318) q[3];
sx q[3];
rz(-2.6716158) q[3];
sx q[3];
rz(0.78891689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47238123) q[0];
sx q[0];
rz(-0.62494576) q[0];
sx q[0];
rz(-2.5770082) q[0];
rz(3.0584884) q[1];
sx q[1];
rz(-1.233498) q[1];
sx q[1];
rz(0.21534236) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1911083) q[0];
sx q[0];
rz(-1.2380301) q[0];
sx q[0];
rz(-0.82975181) q[0];
rz(1.3232857) q[2];
sx q[2];
rz(-0.60192054) q[2];
sx q[2];
rz(1.9844696) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1064062) q[1];
sx q[1];
rz(-1.3760929) q[1];
sx q[1];
rz(3.0343642) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5921949) q[3];
sx q[3];
rz(-1.975172) q[3];
sx q[3];
rz(-0.79493633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5061491) q[2];
sx q[2];
rz(-1.0887257) q[2];
sx q[2];
rz(-2.5396458) q[2];
rz(0.69058949) q[3];
sx q[3];
rz(-0.75406051) q[3];
sx q[3];
rz(0.62600342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2895806) q[0];
sx q[0];
rz(-2.1471725) q[0];
sx q[0];
rz(-0.96486282) q[0];
rz(0.016618641) q[1];
sx q[1];
rz(-0.42533541) q[1];
sx q[1];
rz(-0.66666493) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43685164) q[0];
sx q[0];
rz(-1.1724823) q[0];
sx q[0];
rz(-0.10956357) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0275074) q[2];
sx q[2];
rz(-0.61858656) q[2];
sx q[2];
rz(-1.0416043) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.9469229) q[1];
sx q[1];
rz(-1.2877081) q[1];
sx q[1];
rz(-0.6823632) q[1];
rz(-pi) q[2];
rz(-0.93249647) q[3];
sx q[3];
rz(-0.36358788) q[3];
sx q[3];
rz(2.5225085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.91606402) q[2];
sx q[2];
rz(-1.3801882) q[2];
sx q[2];
rz(1.2529681) q[2];
rz(1.8270252) q[3];
sx q[3];
rz(-1.1809177) q[3];
sx q[3];
rz(-2.0107646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1154293) q[0];
sx q[0];
rz(-0.9342497) q[0];
sx q[0];
rz(-1.3622267) q[0];
rz(1.399614) q[1];
sx q[1];
rz(-1.5796966) q[1];
sx q[1];
rz(1.5302352) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0936733) q[0];
sx q[0];
rz(-1.9903272) q[0];
sx q[0];
rz(-1.6506667) q[0];
rz(-pi) q[1];
rz(-0.43892626) q[2];
sx q[2];
rz(-1.8669087) q[2];
sx q[2];
rz(0.93026464) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6758319) q[1];
sx q[1];
rz(-2.6869046) q[1];
sx q[1];
rz(-2.0783706) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3774032) q[3];
sx q[3];
rz(-2.8838727) q[3];
sx q[3];
rz(1.6267488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7851012) q[2];
sx q[2];
rz(-2.1337528) q[2];
sx q[2];
rz(3.1121758) q[2];
rz(1.4908837) q[3];
sx q[3];
rz(-2.1865032) q[3];
sx q[3];
rz(1.2221378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4239663) q[0];
sx q[0];
rz(-0.76537776) q[0];
sx q[0];
rz(-0.18280612) q[0];
rz(0.43552661) q[1];
sx q[1];
rz(-0.85414779) q[1];
sx q[1];
rz(-2.6307154) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32719192) q[0];
sx q[0];
rz(-1.4609903) q[0];
sx q[0];
rz(-0.26392428) q[0];
rz(-2.4099318) q[2];
sx q[2];
rz(-2.1795142) q[2];
sx q[2];
rz(2.4351956) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6888914) q[1];
sx q[1];
rz(-0.93033067) q[1];
sx q[1];
rz(1.0673317) q[1];
rz(-pi) q[2];
rz(1.0501409) q[3];
sx q[3];
rz(-0.58590349) q[3];
sx q[3];
rz(-2.1396162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9110979) q[2];
sx q[2];
rz(-0.7926597) q[2];
sx q[2];
rz(-2.5210209) q[2];
rz(2.7190322) q[3];
sx q[3];
rz(-1.3603323) q[3];
sx q[3];
rz(-1.2703936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34185103) q[0];
sx q[0];
rz(-2.7809679) q[0];
sx q[0];
rz(-1.0890695) q[0];
rz(-2.0607121) q[1];
sx q[1];
rz(-1.8573449) q[1];
sx q[1];
rz(0.63527766) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95406336) q[0];
sx q[0];
rz(-0.39217338) q[0];
sx q[0];
rz(1.1568882) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8842472) q[2];
sx q[2];
rz(-2.6230276) q[2];
sx q[2];
rz(-2.5093362) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2891846) q[1];
sx q[1];
rz(-1.8361366) q[1];
sx q[1];
rz(-2.5177588) q[1];
x q[2];
rz(-0.20181228) q[3];
sx q[3];
rz(-1.0045856) q[3];
sx q[3];
rz(-1.7081529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.13715956) q[2];
sx q[2];
rz(-1.6798423) q[2];
sx q[2];
rz(-0.00017246406) q[2];
rz(2.4553283) q[3];
sx q[3];
rz(-0.72824794) q[3];
sx q[3];
rz(2.2837158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8453318) q[0];
sx q[0];
rz(-2.2289472) q[0];
sx q[0];
rz(-0.6566748) q[0];
rz(2.7776921) q[1];
sx q[1];
rz(-0.44635043) q[1];
sx q[1];
rz(0.90973967) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0400992) q[0];
sx q[0];
rz(-1.8414521) q[0];
sx q[0];
rz(0.33466848) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7934974) q[2];
sx q[2];
rz(-1.7762134) q[2];
sx q[2];
rz(1.5945438) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.64987953) q[1];
sx q[1];
rz(-0.24398206) q[1];
sx q[1];
rz(-2.6073543) q[1];
x q[2];
rz(-1.886456) q[3];
sx q[3];
rz(-1.3254032) q[3];
sx q[3];
rz(-1.0823702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.83071128) q[2];
sx q[2];
rz(-0.41822663) q[2];
sx q[2];
rz(-2.6549784) q[2];
rz(3.0330372) q[3];
sx q[3];
rz(-0.71043772) q[3];
sx q[3];
rz(1.2238812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9412823) q[0];
sx q[0];
rz(-2.0599984) q[0];
sx q[0];
rz(0.51668984) q[0];
rz(1.5453045) q[1];
sx q[1];
rz(-1.0524806) q[1];
sx q[1];
rz(-2.2470078) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6870118) q[0];
sx q[0];
rz(-0.63475383) q[0];
sx q[0];
rz(-1.4995585) q[0];
rz(-1.9616227) q[2];
sx q[2];
rz(-2.639689) q[2];
sx q[2];
rz(0.48197907) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.91720944) q[1];
sx q[1];
rz(-1.50748) q[1];
sx q[1];
rz(1.2095832) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9856332) q[3];
sx q[3];
rz(-0.92658639) q[3];
sx q[3];
rz(0.19176602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.14780012) q[2];
sx q[2];
rz(-0.96269572) q[2];
sx q[2];
rz(-2.1195892) q[2];
rz(1.8036802) q[3];
sx q[3];
rz(-1.6513848) q[3];
sx q[3];
rz(2.5125304) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.548303) q[0];
sx q[0];
rz(-1.6175445) q[0];
sx q[0];
rz(2.2486726) q[0];
rz(2.100636) q[1];
sx q[1];
rz(-3.0491842) q[1];
sx q[1];
rz(-1.4728117) q[1];
rz(2.2330877) q[2];
sx q[2];
rz(-2.0178595) q[2];
sx q[2];
rz(3.0541228) q[2];
rz(-1.4150423) q[3];
sx q[3];
rz(-2.8372436) q[3];
sx q[3];
rz(-2.5929034) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
