OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.8118892) q[0];
sx q[0];
rz(5.9726297) q[0];
sx q[0];
rz(7.4818727) q[0];
rz(1.1341473) q[1];
sx q[1];
rz(-2.3448047) q[1];
sx q[1];
rz(-0.30153433) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59936675) q[0];
sx q[0];
rz(-0.24655534) q[0];
sx q[0];
rz(0.85041468) q[0];
x q[1];
rz(2.1211336) q[2];
sx q[2];
rz(-0.34780234) q[2];
sx q[2];
rz(-2.2344799) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0092714) q[1];
sx q[1];
rz(-1.6385211) q[1];
sx q[1];
rz(-0.81791116) q[1];
x q[2];
rz(0.188531) q[3];
sx q[3];
rz(-1.1751231) q[3];
sx q[3];
rz(-2.4238732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.6887168) q[2];
sx q[2];
rz(-2.0417002) q[2];
sx q[2];
rz(-1.1348881) q[2];
rz(3.0196043) q[3];
sx q[3];
rz(-2.205409) q[3];
sx q[3];
rz(-0.056099135) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6720471) q[0];
sx q[0];
rz(-2.8944954) q[0];
sx q[0];
rz(-3.1392414) q[0];
rz(-2.7408842) q[1];
sx q[1];
rz(-1.3688764) q[1];
sx q[1];
rz(-0.8561264) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99994217) q[0];
sx q[0];
rz(-2.460833) q[0];
sx q[0];
rz(-1.1589684) q[0];
rz(1.0547423) q[2];
sx q[2];
rz(-1.3203838) q[2];
sx q[2];
rz(2.5977787) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5634671) q[1];
sx q[1];
rz(-1.620786) q[1];
sx q[1];
rz(1.0432788) q[1];
x q[2];
rz(2.5754588) q[3];
sx q[3];
rz(-1.2581035) q[3];
sx q[3];
rz(1.7305525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1043642) q[2];
sx q[2];
rz(-1.9375485) q[2];
sx q[2];
rz(-0.47002235) q[2];
rz(2.5445599) q[3];
sx q[3];
rz(-3.0266914) q[3];
sx q[3];
rz(-1.6980096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3700579) q[0];
sx q[0];
rz(-0.4929339) q[0];
sx q[0];
rz(-0.22802995) q[0];
rz(-2.6920964) q[1];
sx q[1];
rz(-1.0003961) q[1];
sx q[1];
rz(-2.1953348) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3920604) q[0];
sx q[0];
rz(-0.48198732) q[0];
sx q[0];
rz(-2.1493692) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.94257943) q[2];
sx q[2];
rz(-2.1295071) q[2];
sx q[2];
rz(0.0086580833) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6201278) q[1];
sx q[1];
rz(-1.9072272) q[1];
sx q[1];
rz(-1.7640339) q[1];
rz(1.854293) q[3];
sx q[3];
rz(-1.4186267) q[3];
sx q[3];
rz(-2.7622472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9049282) q[2];
sx q[2];
rz(-1.4651352) q[2];
sx q[2];
rz(1.4795335) q[2];
rz(0.18105257) q[3];
sx q[3];
rz(-0.79020399) q[3];
sx q[3];
rz(2.2553867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7774696) q[0];
sx q[0];
rz(-1.5793707) q[0];
sx q[0];
rz(2.2430578) q[0];
rz(-2.6354375) q[1];
sx q[1];
rz(-1.2944784) q[1];
sx q[1];
rz(1.6815965) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4130044) q[0];
sx q[0];
rz(-2.6046136) q[0];
sx q[0];
rz(-0.12852886) q[0];
rz(-0.74331778) q[2];
sx q[2];
rz(-1.905059) q[2];
sx q[2];
rz(2.0845795) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0713363) q[1];
sx q[1];
rz(-1.500529) q[1];
sx q[1];
rz(1.6855168) q[1];
rz(-1.6494895) q[3];
sx q[3];
rz(-0.65641145) q[3];
sx q[3];
rz(2.479913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8755181) q[2];
sx q[2];
rz(-2.6660599) q[2];
sx q[2];
rz(1.5024705) q[2];
rz(-1.8858887) q[3];
sx q[3];
rz(-1.7399104) q[3];
sx q[3];
rz(-0.046253117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.2332377) q[0];
sx q[0];
rz(-2.4112356) q[0];
sx q[0];
rz(2.9610942) q[0];
rz(1.3795229) q[1];
sx q[1];
rz(-0.5677529) q[1];
sx q[1];
rz(2.0464121) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9126606) q[0];
sx q[0];
rz(-1.4353509) q[0];
sx q[0];
rz(1.5082466) q[0];
rz(-1.4873452) q[2];
sx q[2];
rz(-0.7755643) q[2];
sx q[2];
rz(0.47093876) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4597561) q[1];
sx q[1];
rz(-2.3533258) q[1];
sx q[1];
rz(2.9440109) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7322695) q[3];
sx q[3];
rz(-0.69500178) q[3];
sx q[3];
rz(1.7308066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3695662) q[2];
sx q[2];
rz(-2.2942746) q[2];
sx q[2];
rz(-1.0687211) q[2];
rz(0.42701834) q[3];
sx q[3];
rz(-2.2618099) q[3];
sx q[3];
rz(1.8900185) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0291979) q[0];
sx q[0];
rz(-2.2258832) q[0];
sx q[0];
rz(0.29119626) q[0];
rz(1.683782) q[1];
sx q[1];
rz(-2.1144512) q[1];
sx q[1];
rz(-2.2529032) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55808961) q[0];
sx q[0];
rz(-2.0348685) q[0];
sx q[0];
rz(-2.2232673) q[0];
rz(-pi) q[1];
rz(-3.1142919) q[2];
sx q[2];
rz(-2.3841287) q[2];
sx q[2];
rz(-0.94099076) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7900438) q[1];
sx q[1];
rz(-2.1218532) q[1];
sx q[1];
rz(2.181777) q[1];
rz(-pi) q[2];
rz(-2.3360152) q[3];
sx q[3];
rz(-0.42308261) q[3];
sx q[3];
rz(-2.7132636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7500744) q[2];
sx q[2];
rz(-1.8374279) q[2];
sx q[2];
rz(0.94775003) q[2];
rz(-2.9446972) q[3];
sx q[3];
rz(-0.24053776) q[3];
sx q[3];
rz(-0.31015629) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31203684) q[0];
sx q[0];
rz(-2.0289679) q[0];
sx q[0];
rz(0.46992508) q[0];
rz(1.6795109) q[1];
sx q[1];
rz(-1.8406248) q[1];
sx q[1];
rz(1.7376815) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3199661) q[0];
sx q[0];
rz(-1.2228773) q[0];
sx q[0];
rz(0.97396429) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.092488) q[2];
sx q[2];
rz(-1.7564536) q[2];
sx q[2];
rz(0.80362475) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.39170489) q[1];
sx q[1];
rz(-2.0948695) q[1];
sx q[1];
rz(0.4197555) q[1];
rz(0.49755712) q[3];
sx q[3];
rz(-2.4973739) q[3];
sx q[3];
rz(1.2608755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5002084) q[2];
sx q[2];
rz(-2.310014) q[2];
sx q[2];
rz(0.10759648) q[2];
rz(0.61947668) q[3];
sx q[3];
rz(-1.8457125) q[3];
sx q[3];
rz(-1.3639601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7463995) q[0];
sx q[0];
rz(-2.0893607) q[0];
sx q[0];
rz(2.8885544) q[0];
rz(-3.053275) q[1];
sx q[1];
rz(-2.2679236) q[1];
sx q[1];
rz(0.94892445) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63213962) q[0];
sx q[0];
rz(-1.720819) q[0];
sx q[0];
rz(1.8709183) q[0];
rz(-pi) q[1];
rz(-1.9986834) q[2];
sx q[2];
rz(-0.44369953) q[2];
sx q[2];
rz(0.86715172) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1469377) q[1];
sx q[1];
rz(-2.6587464) q[1];
sx q[1];
rz(-1.6864683) q[1];
x q[2];
rz(1.003951) q[3];
sx q[3];
rz(-2.0210055) q[3];
sx q[3];
rz(2.4382044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.39763149) q[2];
sx q[2];
rz(-2.2546015) q[2];
sx q[2];
rz(2.8308947) q[2];
rz(-2.9001696) q[3];
sx q[3];
rz(-1.2025236) q[3];
sx q[3];
rz(2.0588622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0017515) q[0];
sx q[0];
rz(-0.40626353) q[0];
sx q[0];
rz(-1.235442) q[0];
rz(0.48108092) q[1];
sx q[1];
rz(-2.1886539) q[1];
sx q[1];
rz(-1.4280041) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7803376) q[0];
sx q[0];
rz(-1.5819823) q[0];
sx q[0];
rz(2.3821324) q[0];
rz(-pi) q[1];
rz(-2.9129006) q[2];
sx q[2];
rz(-1.8426101) q[2];
sx q[2];
rz(-2.6120409) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.42983152) q[1];
sx q[1];
rz(-1.3138384) q[1];
sx q[1];
rz(0.25196948) q[1];
rz(2.4190524) q[3];
sx q[3];
rz(-1.819639) q[3];
sx q[3];
rz(-1.794508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1741751) q[2];
sx q[2];
rz(-0.82505161) q[2];
sx q[2];
rz(3.1136801) q[2];
rz(-0.86999718) q[3];
sx q[3];
rz(-0.78012192) q[3];
sx q[3];
rz(1.1869259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1160527) q[0];
sx q[0];
rz(-1.5929796) q[0];
sx q[0];
rz(2.0822339) q[0];
rz(1.4147883) q[1];
sx q[1];
rz(-1.4780412) q[1];
sx q[1];
rz(0.47763225) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85738289) q[0];
sx q[0];
rz(-2.3071831) q[0];
sx q[0];
rz(-1.3697903) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9066493) q[2];
sx q[2];
rz(-1.5726461) q[2];
sx q[2];
rz(0.84650485) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.76948386) q[1];
sx q[1];
rz(-1.6830576) q[1];
sx q[1];
rz(-0.054241971) q[1];
rz(-pi) q[2];
rz(-2.3197993) q[3];
sx q[3];
rz(-2.0374277) q[3];
sx q[3];
rz(2.1958283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.16964218) q[2];
sx q[2];
rz(-2.1596238) q[2];
sx q[2];
rz(-0.36894813) q[2];
rz(1.87489) q[3];
sx q[3];
rz(-2.4167175) q[3];
sx q[3];
rz(0.96316159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0030768) q[0];
sx q[0];
rz(-1.2110447) q[0];
sx q[0];
rz(-1.3883653) q[0];
rz(0.44454642) q[1];
sx q[1];
rz(-2.3311756) q[1];
sx q[1];
rz(-0.12745007) q[1];
rz(1.2426022) q[2];
sx q[2];
rz(-1.9933619) q[2];
sx q[2];
rz(0.91771916) q[2];
rz(2.6740549) q[3];
sx q[3];
rz(-1.5880006) q[3];
sx q[3];
rz(1.4505462) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
