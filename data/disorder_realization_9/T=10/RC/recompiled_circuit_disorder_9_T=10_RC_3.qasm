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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7793517) q[0];
sx q[0];
rz(-1.2138214) q[0];
sx q[0];
rz(2.7428521) q[0];
rz(-2.8698679) q[2];
sx q[2];
rz(-2.6368015) q[2];
sx q[2];
rz(-1.4128078) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.049584576) q[1];
sx q[1];
rz(-1.9132179) q[1];
sx q[1];
rz(-2.601806) q[1];
rz(-0.13866339) q[3];
sx q[3];
rz(-1.1274459) q[3];
sx q[3];
rz(1.5594304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.47544605) q[2];
sx q[2];
rz(-2.011516) q[2];
sx q[2];
rz(1.3249935) q[2];
rz(2.825286) q[3];
sx q[3];
rz(-2.2241212) q[3];
sx q[3];
rz(2.2345208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9544202) q[0];
sx q[0];
rz(-2.3311054) q[0];
sx q[0];
rz(-0.95970884) q[0];
rz(1.487544) q[1];
sx q[1];
rz(-2.4580749) q[1];
sx q[1];
rz(0.62746343) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42694416) q[0];
sx q[0];
rz(-1.5639925) q[0];
sx q[0];
rz(3.1338918) q[0];
rz(-3.0395095) q[2];
sx q[2];
rz(-2.9288769) q[2];
sx q[2];
rz(1.4174457) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6154502) q[1];
sx q[1];
rz(-0.28029385) q[1];
sx q[1];
rz(1.0326833) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.596132) q[3];
sx q[3];
rz(-1.5832033) q[3];
sx q[3];
rz(-0.62904639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6015357) q[2];
sx q[2];
rz(-1.7404218) q[2];
sx q[2];
rz(0.97989782) q[2];
rz(-1.6905789) q[3];
sx q[3];
rz(-2.4308949) q[3];
sx q[3];
rz(1.4427982) q[3];
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
rz(2.5805805) q[0];
sx q[0];
rz(-0.0479527) q[0];
sx q[0];
rz(-1.3797492) q[0];
rz(-0.17669949) q[1];
sx q[1];
rz(-1.9243762) q[1];
sx q[1];
rz(1.3476936) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5603148) q[0];
sx q[0];
rz(-2.609195) q[0];
sx q[0];
rz(0.48849948) q[0];
x q[1];
rz(0.10920306) q[2];
sx q[2];
rz(-2.9401527) q[2];
sx q[2];
rz(-0.86004721) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1469517) q[1];
sx q[1];
rz(-2.9839582) q[1];
sx q[1];
rz(-2.4528802) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1191145) q[3];
sx q[3];
rz(-2.936755) q[3];
sx q[3];
rz(1.1289489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1220864) q[2];
sx q[2];
rz(-0.85947376) q[2];
sx q[2];
rz(-1.9879831) q[2];
rz(-0.38976088) q[3];
sx q[3];
rz(-2.6716158) q[3];
sx q[3];
rz(-0.78891689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6692114) q[0];
sx q[0];
rz(-2.5166469) q[0];
sx q[0];
rz(-2.5770082) q[0];
rz(3.0584884) q[1];
sx q[1];
rz(-1.233498) q[1];
sx q[1];
rz(-2.9262503) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1911083) q[0];
sx q[0];
rz(-1.9035625) q[0];
sx q[0];
rz(-2.3118408) q[0];
rz(2.1583546) q[2];
sx q[2];
rz(-1.431627) q[2];
sx q[2];
rz(-0.61901865) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.035186471) q[1];
sx q[1];
rz(-1.3760929) q[1];
sx q[1];
rz(-3.0343642) q[1];
rz(-2.5921949) q[3];
sx q[3];
rz(-1.1664207) q[3];
sx q[3];
rz(2.3466563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.6354436) q[2];
sx q[2];
rz(-1.0887257) q[2];
sx q[2];
rz(-2.5396458) q[2];
rz(-2.4510032) q[3];
sx q[3];
rz(-0.75406051) q[3];
sx q[3];
rz(-2.5155892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2895806) q[0];
sx q[0];
rz(-0.99442011) q[0];
sx q[0];
rz(-0.96486282) q[0];
rz(-0.016618641) q[1];
sx q[1];
rz(-0.42533541) q[1];
sx q[1];
rz(-2.4749277) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71320888) q[0];
sx q[0];
rz(-2.7292626) q[0];
sx q[0];
rz(1.8250188) q[0];
x q[1];
rz(-2.1140852) q[2];
sx q[2];
rz(-0.61858656) q[2];
sx q[2];
rz(2.0999883) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1946698) q[1];
sx q[1];
rz(-1.2877081) q[1];
sx q[1];
rz(-2.4592295) q[1];
x q[2];
rz(-0.93249647) q[3];
sx q[3];
rz(-0.36358788) q[3];
sx q[3];
rz(-0.61908412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.91606402) q[2];
sx q[2];
rz(-1.7614044) q[2];
sx q[2];
rz(-1.8886245) q[2];
rz(1.8270252) q[3];
sx q[3];
rz(-1.9606749) q[3];
sx q[3];
rz(2.0107646) q[3];
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
rz(-pi) q[3];
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
rz(1.1154293) q[0];
sx q[0];
rz(-2.207343) q[0];
sx q[0];
rz(1.3622267) q[0];
rz(-1.399614) q[1];
sx q[1];
rz(-1.5796966) q[1];
sx q[1];
rz(1.6113575) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0936733) q[0];
sx q[0];
rz(-1.9903272) q[0];
sx q[0];
rz(1.490926) q[0];
rz(2.5189581) q[2];
sx q[2];
rz(-0.52402516) q[2];
sx q[2];
rz(-3.0571836) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.02009) q[1];
sx q[1];
rz(-1.1768747) q[1];
sx q[1];
rz(-0.23328383) q[1];
rz(-pi) q[2];
rz(0.050614428) q[3];
sx q[3];
rz(-1.8236056) q[3];
sx q[3];
rz(-1.8265754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.35649148) q[2];
sx q[2];
rz(-1.0078398) q[2];
sx q[2];
rz(-0.029416857) q[2];
rz(1.6507089) q[3];
sx q[3];
rz(-0.95508948) q[3];
sx q[3];
rz(1.2221378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71762639) q[0];
sx q[0];
rz(-2.3762149) q[0];
sx q[0];
rz(2.9587865) q[0];
rz(-2.706066) q[1];
sx q[1];
rz(-0.85414779) q[1];
sx q[1];
rz(0.51087728) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32719192) q[0];
sx q[0];
rz(-1.4609903) q[0];
sx q[0];
rz(2.8776684) q[0];
rz(-pi) q[1];
x q[1];
rz(0.81803825) q[2];
sx q[2];
rz(-2.15089) q[2];
sx q[2];
rz(1.3387418) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4527013) q[1];
sx q[1];
rz(-0.93033067) q[1];
sx q[1];
rz(2.074261) q[1];
rz(-1.0484344) q[3];
sx q[3];
rz(-1.8494542) q[3];
sx q[3];
rz(-0.1230965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2304948) q[2];
sx q[2];
rz(-2.348933) q[2];
sx q[2];
rz(2.5210209) q[2];
rz(2.7190322) q[3];
sx q[3];
rz(-1.7812604) q[3];
sx q[3];
rz(1.2703936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(2.7997416) q[0];
sx q[0];
rz(-2.7809679) q[0];
sx q[0];
rz(1.0890695) q[0];
rz(2.0607121) q[1];
sx q[1];
rz(-1.8573449) q[1];
sx q[1];
rz(-0.63527766) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7437744) q[0];
sx q[0];
rz(-1.9282856) q[0];
sx q[0];
rz(-0.1648358) q[0];
x q[1];
rz(0.17417553) q[2];
sx q[2];
rz(-2.0617699) q[2];
sx q[2];
rz(0.27506405) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5105671) q[1];
sx q[1];
rz(-0.67093611) q[1];
sx q[1];
rz(0.4354233) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.99526309) q[3];
sx q[3];
rz(-1.4008153) q[3];
sx q[3];
rz(3.1135524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.13715956) q[2];
sx q[2];
rz(-1.6798423) q[2];
sx q[2];
rz(3.1414202) q[2];
rz(-0.6862644) q[3];
sx q[3];
rz(-0.72824794) q[3];
sx q[3];
rz(2.2837158) q[3];
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
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2962608) q[0];
sx q[0];
rz(-0.91264549) q[0];
sx q[0];
rz(2.4849179) q[0];
rz(-0.36390057) q[1];
sx q[1];
rz(-2.6952422) q[1];
sx q[1];
rz(2.231853) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0400992) q[0];
sx q[0];
rz(-1.3001406) q[0];
sx q[0];
rz(0.33466848) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7889195) q[2];
sx q[2];
rz(-1.2303196) q[2];
sx q[2];
rz(-3.043963) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7418356) q[1];
sx q[1];
rz(-1.4474807) q[1];
sx q[1];
rz(-0.21106212) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2496201) q[3];
sx q[3];
rz(-0.39728764) q[3];
sx q[3];
rz(0.15115034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3108814) q[2];
sx q[2];
rz(-2.723366) q[2];
sx q[2];
rz(-2.6549784) q[2];
rz(0.1085554) q[3];
sx q[3];
rz(-2.4311549) q[3];
sx q[3];
rz(-1.9177115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20031032) q[0];
sx q[0];
rz(-2.0599984) q[0];
sx q[0];
rz(-2.6249028) q[0];
rz(-1.5962881) q[1];
sx q[1];
rz(-1.0524806) q[1];
sx q[1];
rz(0.89458481) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54297011) q[0];
sx q[0];
rz(-2.2036836) q[0];
sx q[0];
rz(-3.0892239) q[0];
rz(1.1799699) q[2];
sx q[2];
rz(-0.50190364) q[2];
sx q[2];
rz(-0.48197907) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.91720944) q[1];
sx q[1];
rz(-1.6341126) q[1];
sx q[1];
rz(-1.9320095) q[1];
rz(-pi) q[2];
rz(0.15595943) q[3];
sx q[3];
rz(-2.2150063) q[3];
sx q[3];
rz(2.9498266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9937925) q[2];
sx q[2];
rz(-0.96269572) q[2];
sx q[2];
rz(-1.0220035) q[2];
rz(-1.8036802) q[3];
sx q[3];
rz(-1.6513848) q[3];
sx q[3];
rz(-2.5125304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5932896) q[0];
sx q[0];
rz(-1.5240482) q[0];
sx q[0];
rz(-0.89292009) q[0];
rz(1.0409566) q[1];
sx q[1];
rz(-0.092408471) q[1];
sx q[1];
rz(1.668781) q[1];
rz(-0.54626089) q[2];
sx q[2];
rz(-2.1585474) q[2];
sx q[2];
rz(-1.3331158) q[2];
rz(1.4150423) q[3];
sx q[3];
rz(-0.30434904) q[3];
sx q[3];
rz(0.54868922) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
