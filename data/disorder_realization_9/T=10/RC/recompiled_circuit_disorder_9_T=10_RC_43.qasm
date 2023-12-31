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
rz(-0.79467264) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.062382467) q[0];
sx q[0];
rz(-1.198472) q[0];
sx q[0];
rz(1.9553493) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.48912666) q[2];
sx q[2];
rz(-1.7009652) q[2];
sx q[2];
rz(2.7444073) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0104048) q[1];
sx q[1];
rz(-2.511575) q[1];
sx q[1];
rz(0.6063993) q[1];
x q[2];
rz(-1.854033) q[3];
sx q[3];
rz(-2.6784416) q[3];
sx q[3];
rz(1.8739665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6661466) q[2];
sx q[2];
rz(-2.011516) q[2];
sx q[2];
rz(-1.8165992) q[2];
rz(-0.31630668) q[3];
sx q[3];
rz(-2.2241212) q[3];
sx q[3];
rz(-0.90707183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9544202) q[0];
sx q[0];
rz(-2.3311054) q[0];
sx q[0];
rz(2.1818838) q[0];
rz(-1.487544) q[1];
sx q[1];
rz(-2.4580749) q[1];
sx q[1];
rz(2.5141292) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42022959) q[0];
sx q[0];
rz(-3.1313167) q[0];
sx q[0];
rz(2.4179439) q[0];
rz(-1.5487899) q[2];
sx q[2];
rz(-1.3592048) q[2];
sx q[2];
rz(-1.6197268) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5261425) q[1];
sx q[1];
rz(-2.8612988) q[1];
sx q[1];
rz(-2.1089094) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1153568) q[3];
sx q[3];
rz(-3.1133828) q[3];
sx q[3];
rz(1.3970323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.54005694) q[2];
sx q[2];
rz(-1.4011708) q[2];
sx q[2];
rz(-0.97989782) q[2];
rz(-1.6905789) q[3];
sx q[3];
rz(-2.4308949) q[3];
sx q[3];
rz(1.4427982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56101218) q[0];
sx q[0];
rz(-3.09364) q[0];
sx q[0];
rz(-1.7618435) q[0];
rz(0.17669949) q[1];
sx q[1];
rz(-1.9243762) q[1];
sx q[1];
rz(-1.3476936) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5603148) q[0];
sx q[0];
rz(-0.53239765) q[0];
sx q[0];
rz(-0.48849948) q[0];
rz(1.5930487) q[2];
sx q[2];
rz(-1.77102) q[2];
sx q[2];
rz(-2.3929838) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.29979953) q[1];
sx q[1];
rz(-1.449297) q[1];
sx q[1];
rz(1.4701162) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0224781) q[3];
sx q[3];
rz(-0.20483769) q[3];
sx q[3];
rz(-2.0126437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1220864) q[2];
sx q[2];
rz(-0.85947376) q[2];
sx q[2];
rz(-1.1536095) q[2];
rz(-0.38976088) q[3];
sx q[3];
rz(-0.4699769) q[3];
sx q[3];
rz(-2.3526758) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6692114) q[0];
sx q[0];
rz(-0.62494576) q[0];
sx q[0];
rz(-2.5770082) q[0];
rz(-3.0584884) q[1];
sx q[1];
rz(-1.233498) q[1];
sx q[1];
rz(-0.21534236) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089245361) q[0];
sx q[0];
rz(-0.87886341) q[0];
sx q[0];
rz(-2.7034876) q[0];
rz(-2.9748561) q[2];
sx q[2];
rz(-0.98966375) q[2];
sx q[2];
rz(-0.8596479) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.035186471) q[1];
sx q[1];
rz(-1.7654997) q[1];
sx q[1];
rz(3.0343642) q[1];
rz(-pi) q[2];
rz(-0.54939778) q[3];
sx q[3];
rz(-1.1664207) q[3];
sx q[3];
rz(0.79493633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5061491) q[2];
sx q[2];
rz(-2.052867) q[2];
sx q[2];
rz(2.5396458) q[2];
rz(-0.69058949) q[3];
sx q[3];
rz(-0.75406051) q[3];
sx q[3];
rz(-0.62600342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2895806) q[0];
sx q[0];
rz(-0.99442011) q[0];
sx q[0];
rz(0.96486282) q[0];
rz(-0.016618641) q[1];
sx q[1];
rz(-0.42533541) q[1];
sx q[1];
rz(-2.4749277) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0913038) q[0];
sx q[0];
rz(-1.6717523) q[0];
sx q[0];
rz(-1.9712649) q[0];
rz(-pi) q[1];
rz(1.0275074) q[2];
sx q[2];
rz(-0.61858656) q[2];
sx q[2];
rz(2.0999883) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.95496817) q[1];
sx q[1];
rz(-0.72997) q[1];
sx q[1];
rz(-0.43220046) q[1];
x q[2];
rz(-1.8673709) q[3];
sx q[3];
rz(-1.7843102) q[3];
sx q[3];
rz(0.34542686) q[3];
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
rz(1.8270252) q[3];
sx q[3];
rz(-1.9606749) q[3];
sx q[3];
rz(-1.130828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0261633) q[0];
sx q[0];
rz(-0.9342497) q[0];
sx q[0];
rz(-1.3622267) q[0];
rz(1.7419787) q[1];
sx q[1];
rz(-1.5796966) q[1];
sx q[1];
rz(-1.5302352) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8996457) q[0];
sx q[0];
rz(-2.7149704) q[0];
sx q[0];
rz(2.9645779) q[0];
rz(-pi) q[1];
rz(-1.2457232) q[2];
sx q[2];
rz(-1.9893861) q[2];
sx q[2];
rz(-0.77667728) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6758319) q[1];
sx q[1];
rz(-2.6869046) q[1];
sx q[1];
rz(-1.0632221) q[1];
rz(1.7641894) q[3];
sx q[3];
rz(-0.25771991) q[3];
sx q[3];
rz(1.5148439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7851012) q[2];
sx q[2];
rz(-2.1337528) q[2];
sx q[2];
rz(3.1121758) q[2];
rz(1.6507089) q[3];
sx q[3];
rz(-0.95508948) q[3];
sx q[3];
rz(1.2221378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4239663) q[0];
sx q[0];
rz(-2.3762149) q[0];
sx q[0];
rz(-0.18280612) q[0];
rz(0.43552661) q[1];
sx q[1];
rz(-2.2874449) q[1];
sx q[1];
rz(2.6307154) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32719192) q[0];
sx q[0];
rz(-1.4609903) q[0];
sx q[0];
rz(-2.8776684) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3235544) q[2];
sx q[2];
rz(-0.99070264) q[2];
sx q[2];
rz(1.8028508) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4527013) q[1];
sx q[1];
rz(-0.93033067) q[1];
sx q[1];
rz(2.074261) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8227311) q[3];
sx q[3];
rz(-1.0705035) q[3];
sx q[3];
rz(1.5368411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9110979) q[2];
sx q[2];
rz(-2.348933) q[2];
sx q[2];
rz(-2.5210209) q[2];
rz(-2.7190322) q[3];
sx q[3];
rz(-1.7812604) q[3];
sx q[3];
rz(-1.2703936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34185103) q[0];
sx q[0];
rz(-2.7809679) q[0];
sx q[0];
rz(1.0890695) q[0];
rz(2.0607121) q[1];
sx q[1];
rz(-1.8573449) q[1];
sx q[1];
rz(2.506315) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9104722) q[0];
sx q[0];
rz(-1.7251245) q[0];
sx q[0];
rz(1.9327823) q[0];
rz(-pi) q[1];
rz(-2.0681357) q[2];
sx q[2];
rz(-1.724223) q[2];
sx q[2];
rz(1.9286326) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0465225) q[1];
sx q[1];
rz(-0.97192837) q[1];
sx q[1];
rz(1.8938766) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.99526309) q[3];
sx q[3];
rz(-1.7407773) q[3];
sx q[3];
rz(0.02804027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0044331) q[2];
sx q[2];
rz(-1.6798423) q[2];
sx q[2];
rz(-3.1414202) q[2];
rz(0.6862644) q[3];
sx q[3];
rz(-0.72824794) q[3];
sx q[3];
rz(0.85787684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8453318) q[0];
sx q[0];
rz(-0.91264549) q[0];
sx q[0];
rz(0.6566748) q[0];
rz(-2.7776921) q[1];
sx q[1];
rz(-2.6952422) q[1];
sx q[1];
rz(-2.231853) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0400992) q[0];
sx q[0];
rz(-1.3001406) q[0];
sx q[0];
rz(-2.8069242) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7934974) q[2];
sx q[2];
rz(-1.3653792) q[2];
sx q[2];
rz(-1.5470488) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.39975702) q[1];
sx q[1];
rz(-1.6941119) q[1];
sx q[1];
rz(-0.21106212) q[1];
rz(-pi) q[2];
rz(0.25760381) q[3];
sx q[3];
rz(-1.2649049) q[3];
sx q[3];
rz(0.56759873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.83071128) q[2];
sx q[2];
rz(-0.41822663) q[2];
sx q[2];
rz(2.6549784) q[2];
rz(-0.1085554) q[3];
sx q[3];
rz(-0.71043772) q[3];
sx q[3];
rz(-1.9177115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9412823) q[0];
sx q[0];
rz(-1.0815942) q[0];
sx q[0];
rz(0.51668984) q[0];
rz(-1.5962881) q[1];
sx q[1];
rz(-1.0524806) q[1];
sx q[1];
rz(-2.2470078) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45458083) q[0];
sx q[0];
rz(-0.63475383) q[0];
sx q[0];
rz(-1.4995585) q[0];
rz(1.9616227) q[2];
sx q[2];
rz(-0.50190364) q[2];
sx q[2];
rz(0.48197907) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.48764187) q[1];
sx q[1];
rz(-0.36648053) q[1];
sx q[1];
rz(-1.7483064) q[1];
rz(2.9856332) q[3];
sx q[3];
rz(-2.2150063) q[3];
sx q[3];
rz(0.19176602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9937925) q[2];
sx q[2];
rz(-2.1788969) q[2];
sx q[2];
rz(-2.1195892) q[2];
rz(1.8036802) q[3];
sx q[3];
rz(-1.6513848) q[3];
sx q[3];
rz(-0.62906229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5932896) q[0];
sx q[0];
rz(-1.5240482) q[0];
sx q[0];
rz(-0.89292009) q[0];
rz(-1.0409566) q[1];
sx q[1];
rz(-3.0491842) q[1];
sx q[1];
rz(-1.4728117) q[1];
rz(-2.2330166) q[2];
sx q[2];
rz(-0.77975811) q[2];
sx q[2];
rz(-2.1644885) q[2];
rz(-1.2699119) q[3];
sx q[3];
rz(-1.6172998) q[3];
sx q[3];
rz(-0.87340364) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
