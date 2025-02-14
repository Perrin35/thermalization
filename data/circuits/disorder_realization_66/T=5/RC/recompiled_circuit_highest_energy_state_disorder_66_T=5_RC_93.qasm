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
rz(-0.13770667) q[0];
sx q[0];
rz(-0.51578155) q[0];
sx q[0];
rz(-3.1006324) q[0];
rz(0.75086683) q[1];
sx q[1];
rz(2.6584396) q[1];
sx q[1];
rz(9.2106342) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5222408) q[0];
sx q[0];
rz(-2.9532814) q[0];
sx q[0];
rz(-0.49647402) q[0];
rz(0.08835596) q[2];
sx q[2];
rz(-1.7680829) q[2];
sx q[2];
rz(-0.046869761) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1119207) q[1];
sx q[1];
rz(-0.73603928) q[1];
sx q[1];
rz(-1.8007832) q[1];
rz(-2.0748795) q[3];
sx q[3];
rz(-1.136054) q[3];
sx q[3];
rz(-2.6047516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9302049) q[2];
sx q[2];
rz(-1.7519506) q[2];
sx q[2];
rz(-2.1432121) q[2];
rz(1.8006648) q[3];
sx q[3];
rz(-2.0710335) q[3];
sx q[3];
rz(-2.4836331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.9407161) q[0];
sx q[0];
rz(-0.91280341) q[0];
sx q[0];
rz(0.87769133) q[0];
rz(-0.25582036) q[1];
sx q[1];
rz(-2.1820549) q[1];
sx q[1];
rz(-2.6928601) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2923778) q[0];
sx q[0];
rz(-1.6265944) q[0];
sx q[0];
rz(2.8120338) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0322689) q[2];
sx q[2];
rz(-2.903679) q[2];
sx q[2];
rz(-0.21060056) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.06801735) q[1];
sx q[1];
rz(-2.3367391) q[1];
sx q[1];
rz(-2.6642403) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99440672) q[3];
sx q[3];
rz(-0.72184104) q[3];
sx q[3];
rz(2.5979258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9933219) q[2];
sx q[2];
rz(-1.4485056) q[2];
sx q[2];
rz(1.5987965) q[2];
rz(2.350542) q[3];
sx q[3];
rz(-1.6796716) q[3];
sx q[3];
rz(-2.3492474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1087082) q[0];
sx q[0];
rz(-0.59297639) q[0];
sx q[0];
rz(1.2805043) q[0];
rz(-0.20901021) q[1];
sx q[1];
rz(-2.0031395) q[1];
sx q[1];
rz(1.7144405) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2036594) q[0];
sx q[0];
rz(-1.0190932) q[0];
sx q[0];
rz(-2.0996086) q[0];
x q[1];
rz(1.7383606) q[2];
sx q[2];
rz(-1.6779876) q[2];
sx q[2];
rz(0.78533781) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5722981) q[1];
sx q[1];
rz(-1.3009909) q[1];
sx q[1];
rz(2.5266179) q[1];
x q[2];
rz(-2.9455037) q[3];
sx q[3];
rz(-0.59299378) q[3];
sx q[3];
rz(0.65341572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7872494) q[2];
sx q[2];
rz(-0.83293438) q[2];
sx q[2];
rz(0.99861097) q[2];
rz(-0.35501114) q[3];
sx q[3];
rz(-1.3839) q[3];
sx q[3];
rz(1.2242873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.3806279) q[0];
sx q[0];
rz(-0.17292085) q[0];
sx q[0];
rz(0.91249102) q[0];
rz(-1.6418705) q[1];
sx q[1];
rz(-1.8955756) q[1];
sx q[1];
rz(-1.46371) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3974103) q[0];
sx q[0];
rz(-1.46755) q[0];
sx q[0];
rz(0.29241568) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6557515) q[2];
sx q[2];
rz(-0.96829295) q[2];
sx q[2];
rz(2.9311644) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.50747516) q[1];
sx q[1];
rz(-1.1833131) q[1];
sx q[1];
rz(-1.3204367) q[1];
rz(-pi) q[2];
rz(1.2578418) q[3];
sx q[3];
rz(-2.3944562) q[3];
sx q[3];
rz(1.136387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0616167) q[2];
sx q[2];
rz(-1.4392263) q[2];
sx q[2];
rz(1.3920791) q[2];
rz(-1.1213087) q[3];
sx q[3];
rz(-1.0871355) q[3];
sx q[3];
rz(-0.054718941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-3.1295117) q[0];
sx q[0];
rz(-1.5665781) q[0];
sx q[0];
rz(2.6636301) q[0];
rz(1.6732008) q[1];
sx q[1];
rz(-1.5696399) q[1];
sx q[1];
rz(-2.0782616) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7655575) q[0];
sx q[0];
rz(-2.4004694) q[0];
sx q[0];
rz(-2.3149109) q[0];
rz(-pi) q[1];
rz(-1.0139017) q[2];
sx q[2];
rz(-1.4443358) q[2];
sx q[2];
rz(-1.2173956) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0681615) q[1];
sx q[1];
rz(-0.85660645) q[1];
sx q[1];
rz(2.5769907) q[1];
x q[2];
rz(1.8724779) q[3];
sx q[3];
rz(-1.2815803) q[3];
sx q[3];
rz(2.7986643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8308476) q[2];
sx q[2];
rz(-1.9644535) q[2];
sx q[2];
rz(-1.8929405) q[2];
rz(2.2511075) q[3];
sx q[3];
rz(-1.9200446) q[3];
sx q[3];
rz(-2.7142767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(1.2763727) q[0];
sx q[0];
rz(-2.7557912) q[0];
sx q[0];
rz(-1.2060097) q[0];
rz(-1.457006) q[1];
sx q[1];
rz(-0.99452368) q[1];
sx q[1];
rz(1.0135894) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.834965) q[0];
sx q[0];
rz(-1.6205015) q[0];
sx q[0];
rz(1.6139362) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5382928) q[2];
sx q[2];
rz(-0.68317709) q[2];
sx q[2];
rz(-0.3521093) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4800548) q[1];
sx q[1];
rz(-1.3549003) q[1];
sx q[1];
rz(0.56592249) q[1];
x q[2];
rz(-1.9928855) q[3];
sx q[3];
rz(-1.6339025) q[3];
sx q[3];
rz(-2.1247739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4509585) q[2];
sx q[2];
rz(-0.75814191) q[2];
sx q[2];
rz(-0.071619384) q[2];
rz(-3.0028263) q[3];
sx q[3];
rz(-1.7641188) q[3];
sx q[3];
rz(0.57851401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-3.0305158) q[0];
sx q[0];
rz(-1.0826305) q[0];
sx q[0];
rz(-1.507933) q[0];
rz(-2.0207479) q[1];
sx q[1];
rz(-1.0034674) q[1];
sx q[1];
rz(2.9337163) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6217529) q[0];
sx q[0];
rz(-0.43448453) q[0];
sx q[0];
rz(2.1293917) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6034749) q[2];
sx q[2];
rz(-0.30024116) q[2];
sx q[2];
rz(-3.0150692) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0130335) q[1];
sx q[1];
rz(-2.7454174) q[1];
sx q[1];
rz(-2.106656) q[1];
rz(-pi) q[2];
rz(-0.052938264) q[3];
sx q[3];
rz(-1.8013305) q[3];
sx q[3];
rz(-2.9028877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8230744) q[2];
sx q[2];
rz(-1.069331) q[2];
sx q[2];
rz(11/(7*pi)) q[2];
rz(-0.21627608) q[3];
sx q[3];
rz(-1.4803156) q[3];
sx q[3];
rz(1.1094619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0529387) q[0];
sx q[0];
rz(-1.5733938) q[0];
sx q[0];
rz(1.1608359) q[0];
rz(2.898518) q[1];
sx q[1];
rz(-1.6583775) q[1];
sx q[1];
rz(1.9879139) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50621158) q[0];
sx q[0];
rz(-0.46790174) q[0];
sx q[0];
rz(-2.0718396) q[0];
x q[1];
rz(1.8249545) q[2];
sx q[2];
rz(-1.1593737) q[2];
sx q[2];
rz(-0.14688497) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.28330559) q[1];
sx q[1];
rz(-1.518462) q[1];
sx q[1];
rz(3.1238265) q[1];
rz(2.1716464) q[3];
sx q[3];
rz(-1.5283397) q[3];
sx q[3];
rz(0.24985931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2005177) q[2];
sx q[2];
rz(-2.2633471) q[2];
sx q[2];
rz(2.3556457) q[2];
rz(-0.032622967) q[3];
sx q[3];
rz(-0.91898188) q[3];
sx q[3];
rz(-0.48891625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13474473) q[0];
sx q[0];
rz(-2.5866046) q[0];
sx q[0];
rz(-0.81163374) q[0];
rz(-0.45400485) q[1];
sx q[1];
rz(-1.7894527) q[1];
sx q[1];
rz(0.46617359) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4056835) q[0];
sx q[0];
rz(-2.0529593) q[0];
sx q[0];
rz(-0.10660556) q[0];
x q[1];
rz(0.42731419) q[2];
sx q[2];
rz(-1.1348339) q[2];
sx q[2];
rz(-2.7183926) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5890903) q[1];
sx q[1];
rz(-1.9297761) q[1];
sx q[1];
rz(0.41429452) q[1];
rz(2.3639244) q[3];
sx q[3];
rz(-2.127217) q[3];
sx q[3];
rz(-0.53529763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0635282) q[2];
sx q[2];
rz(-1.177265) q[2];
sx q[2];
rz(1.9130116) q[2];
rz(0.56033963) q[3];
sx q[3];
rz(-1.8599583) q[3];
sx q[3];
rz(1.8165711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90801564) q[0];
sx q[0];
rz(-0.19291872) q[0];
sx q[0];
rz(-0.31780258) q[0];
rz(-2.4494412) q[1];
sx q[1];
rz(-1.0617804) q[1];
sx q[1];
rz(-2.1923776) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91207165) q[0];
sx q[0];
rz(-1.6229652) q[0];
sx q[0];
rz(-2.0020141) q[0];
rz(0.28905343) q[2];
sx q[2];
rz(-0.9346644) q[2];
sx q[2];
rz(-2.4704982) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.91864955) q[1];
sx q[1];
rz(-0.50477305) q[1];
sx q[1];
rz(-0.10279067) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1032221) q[3];
sx q[3];
rz(-0.11603131) q[3];
sx q[3];
rz(-0.86933245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.067351643) q[2];
sx q[2];
rz(-0.94506741) q[2];
sx q[2];
rz(2.2828339) q[2];
rz(-2.7680715) q[3];
sx q[3];
rz(-1.7645323) q[3];
sx q[3];
rz(0.080373272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.2273979) q[0];
sx q[0];
rz(-0.62794958) q[0];
sx q[0];
rz(-3.0218883) q[0];
rz(-1.9627199) q[1];
sx q[1];
rz(-0.96794712) q[1];
sx q[1];
rz(1.6603574) q[1];
rz(-1.0059849) q[2];
sx q[2];
rz(-2.0567675) q[2];
sx q[2];
rz(-1.3153402) q[2];
rz(-0.12283267) q[3];
sx q[3];
rz(-0.67295495) q[3];
sx q[3];
rz(0.72257696) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
