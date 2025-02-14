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
rz(-1.834637) q[0];
sx q[0];
rz(-0.87943465) q[0];
sx q[0];
rz(-0.49361324) q[0];
rz(-1.2797132) q[1];
sx q[1];
rz(3.7682025) q[1];
sx q[1];
rz(11.087853) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0644306) q[0];
sx q[0];
rz(-2.0666984) q[0];
sx q[0];
rz(-2.2999963) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8952719) q[2];
sx q[2];
rz(-2.4518584) q[2];
sx q[2];
rz(1.5472428) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2027758) q[1];
sx q[1];
rz(-0.73828739) q[1];
sx q[1];
rz(-2.8144005) q[1];
x q[2];
rz(-0.93500626) q[3];
sx q[3];
rz(-1.8719487) q[3];
sx q[3];
rz(-0.13553916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0066068) q[2];
sx q[2];
rz(-1.1400433) q[2];
sx q[2];
rz(-2.0331649) q[2];
rz(1.8290352) q[3];
sx q[3];
rz(-2.8969942) q[3];
sx q[3];
rz(-2.826622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7600064) q[0];
sx q[0];
rz(-0.8257603) q[0];
sx q[0];
rz(0.015856892) q[0];
rz(0.99041692) q[1];
sx q[1];
rz(-2.3742193) q[1];
sx q[1];
rz(-2.2784065) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6619157) q[0];
sx q[0];
rz(-0.74115314) q[0];
sx q[0];
rz(2.9984498) q[0];
rz(-pi) q[1];
rz(2.3658889) q[2];
sx q[2];
rz(-0.45368567) q[2];
sx q[2];
rz(-2.016748) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1474804) q[1];
sx q[1];
rz(-1.0801472) q[1];
sx q[1];
rz(1.3730241) q[1];
rz(-pi) q[2];
rz(-0.14948577) q[3];
sx q[3];
rz(-0.85452291) q[3];
sx q[3];
rz(-0.0090741875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.70011675) q[2];
sx q[2];
rz(-0.2773383) q[2];
sx q[2];
rz(-0.46879834) q[2];
rz(0.68764728) q[3];
sx q[3];
rz(-1.5117437) q[3];
sx q[3];
rz(-0.30011737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1402635) q[0];
sx q[0];
rz(-0.92324531) q[0];
sx q[0];
rz(-2.4053251) q[0];
rz(-2.7983792) q[1];
sx q[1];
rz(-0.88756573) q[1];
sx q[1];
rz(-0.18377486) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4030754) q[0];
sx q[0];
rz(-0.74664298) q[0];
sx q[0];
rz(-3.0439506) q[0];
rz(-0.60004514) q[2];
sx q[2];
rz(-0.79216865) q[2];
sx q[2];
rz(-2.5320657) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0818834) q[1];
sx q[1];
rz(-0.78468152) q[1];
sx q[1];
rz(-0.42894026) q[1];
rz(-pi) q[2];
rz(-2.7815357) q[3];
sx q[3];
rz(-2.0814133) q[3];
sx q[3];
rz(-2.1329481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.78202128) q[2];
sx q[2];
rz(-2.6247793) q[2];
sx q[2];
rz(0.51592958) q[2];
rz(-1.1082209) q[3];
sx q[3];
rz(-2.5907232) q[3];
sx q[3];
rz(1.9903323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
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
rz(0.16294031) q[0];
sx q[0];
rz(-1.2115703) q[0];
sx q[0];
rz(-0.51280713) q[0];
rz(0.45979744) q[1];
sx q[1];
rz(-1.5922981) q[1];
sx q[1];
rz(-0.76382452) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5483185) q[0];
sx q[0];
rz(-1.700907) q[0];
sx q[0];
rz(1.4553372) q[0];
x q[1];
rz(-0.41641367) q[2];
sx q[2];
rz(-1.8685307) q[2];
sx q[2];
rz(0.99978757) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0415548) q[1];
sx q[1];
rz(-2.6943992) q[1];
sx q[1];
rz(0.87765043) q[1];
rz(-pi) q[2];
rz(0.34448907) q[3];
sx q[3];
rz(-1.4360582) q[3];
sx q[3];
rz(-1.3400934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1226471) q[2];
sx q[2];
rz(-0.15572369) q[2];
sx q[2];
rz(0.32749185) q[2];
rz(0.28199768) q[3];
sx q[3];
rz(-1.7433386) q[3];
sx q[3];
rz(-1.5871083) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9649488) q[0];
sx q[0];
rz(-2.7556941) q[0];
sx q[0];
rz(-0.97736812) q[0];
rz(-0.36879677) q[1];
sx q[1];
rz(-2.2803523) q[1];
sx q[1];
rz(-1.4126973) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4810886) q[0];
sx q[0];
rz(-1.3810643) q[0];
sx q[0];
rz(2.0697748) q[0];
rz(2.0875889) q[2];
sx q[2];
rz(-1.5587423) q[2];
sx q[2];
rz(-2.4940559) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8846795) q[1];
sx q[1];
rz(-1.2882782) q[1];
sx q[1];
rz(3.1090841) q[1];
x q[2];
rz(1.8983973) q[3];
sx q[3];
rz(-1.0169815) q[3];
sx q[3];
rz(-2.6462951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.29350027) q[2];
sx q[2];
rz(-2.3257181) q[2];
sx q[2];
rz(0.14981848) q[2];
rz(0.39854974) q[3];
sx q[3];
rz(-1.5236676) q[3];
sx q[3];
rz(0.15628763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3715816) q[0];
sx q[0];
rz(-3.0122029) q[0];
sx q[0];
rz(0.55321252) q[0];
rz(0.54168701) q[1];
sx q[1];
rz(-1.0463511) q[1];
sx q[1];
rz(2.6936626) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0174367) q[0];
sx q[0];
rz(-1.7154335) q[0];
sx q[0];
rz(-2.7516624) q[0];
rz(-pi) q[1];
rz(2.7121353) q[2];
sx q[2];
rz(-1.1633368) q[2];
sx q[2];
rz(-2.7384659) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3975911) q[1];
sx q[1];
rz(-0.41644704) q[1];
sx q[1];
rz(-1.9201502) q[1];
rz(2.9782441) q[3];
sx q[3];
rz(-2.8046097) q[3];
sx q[3];
rz(0.28675045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.64122671) q[2];
sx q[2];
rz(-1.9772823) q[2];
sx q[2];
rz(-2.3424303) q[2];
rz(-2.7568119) q[3];
sx q[3];
rz(-0.32618263) q[3];
sx q[3];
rz(-2.1414115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.36050972) q[0];
sx q[0];
rz(-1.0436844) q[0];
sx q[0];
rz(0.59180301) q[0];
rz(-2.7599755) q[1];
sx q[1];
rz(-2.7615669) q[1];
sx q[1];
rz(2.9891678) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7397933) q[0];
sx q[0];
rz(-1.2619979) q[0];
sx q[0];
rz(-1.0303506) q[0];
rz(0.47758684) q[2];
sx q[2];
rz(-1.7268333) q[2];
sx q[2];
rz(-1.0797155) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.85411863) q[1];
sx q[1];
rz(-1.2240613) q[1];
sx q[1];
rz(-1.5513913) q[1];
x q[2];
rz(-2.3197738) q[3];
sx q[3];
rz(-0.83335175) q[3];
sx q[3];
rz(2.5488473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3064208) q[2];
sx q[2];
rz(-0.99205899) q[2];
sx q[2];
rz(1.4527808) q[2];
rz(0.49550223) q[3];
sx q[3];
rz(-0.82363868) q[3];
sx q[3];
rz(2.5775094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38967663) q[0];
sx q[0];
rz(-0.037411995) q[0];
sx q[0];
rz(-2.9801242) q[0];
rz(-3.1096733) q[1];
sx q[1];
rz(-2.4999764) q[1];
sx q[1];
rz(-1.8643103) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.227018) q[0];
sx q[0];
rz(-0.76138568) q[0];
sx q[0];
rz(-0.37933357) q[0];
x q[1];
rz(-2.9109898) q[2];
sx q[2];
rz(-1.1586894) q[2];
sx q[2];
rz(1.3889564) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.016170597) q[1];
sx q[1];
rz(-1.6787663) q[1];
sx q[1];
rz(-0.30534621) q[1];
rz(3.1031514) q[3];
sx q[3];
rz(-1.2769967) q[3];
sx q[3];
rz(0.97803365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.45234597) q[2];
sx q[2];
rz(-2.2301058) q[2];
sx q[2];
rz(-0.39608836) q[2];
rz(2.7416157) q[3];
sx q[3];
rz(-2.5466099) q[3];
sx q[3];
rz(-2.2649435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24092995) q[0];
sx q[0];
rz(-0.018095896) q[0];
sx q[0];
rz(2.990429) q[0];
rz(0.97686544) q[1];
sx q[1];
rz(-2.7604389) q[1];
sx q[1];
rz(0.74473286) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1158651) q[0];
sx q[0];
rz(-0.34894279) q[0];
sx q[0];
rz(0.92744382) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9989002) q[2];
sx q[2];
rz(-1.1617134) q[2];
sx q[2];
rz(-0.2689223) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.1972629) q[1];
sx q[1];
rz(-1.5722317) q[1];
sx q[1];
rz(2.4932842) q[1];
rz(-pi) q[2];
rz(1.8119393) q[3];
sx q[3];
rz(-1.9134812) q[3];
sx q[3];
rz(-1.8543275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.62290827) q[2];
sx q[2];
rz(-1.238995) q[2];
sx q[2];
rz(2.1935479) q[2];
rz(-1.0575804) q[3];
sx q[3];
rz(-1.4761997) q[3];
sx q[3];
rz(1.074056) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0483765) q[0];
sx q[0];
rz(-2.8792448) q[0];
sx q[0];
rz(2.9076305) q[0];
rz(1.9562862) q[1];
sx q[1];
rz(-0.92820853) q[1];
sx q[1];
rz(-1.8918461) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.043086655) q[0];
sx q[0];
rz(-1.0191425) q[0];
sx q[0];
rz(1.9123402) q[0];
x q[1];
rz(0.62057497) q[2];
sx q[2];
rz(-1.7399551) q[2];
sx q[2];
rz(-2.754654) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9712914) q[1];
sx q[1];
rz(-0.93496041) q[1];
sx q[1];
rz(-2.2595992) q[1];
rz(-pi) q[2];
rz(1.3951282) q[3];
sx q[3];
rz(-2.2445956) q[3];
sx q[3];
rz(-1.8927946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.11448161) q[2];
sx q[2];
rz(-1.1610843) q[2];
sx q[2];
rz(-1.7005881) q[2];
rz(3.0103185) q[3];
sx q[3];
rz(-2.6797397) q[3];
sx q[3];
rz(2.89768) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.048653614) q[0];
sx q[0];
rz(-2.1760512) q[0];
sx q[0];
rz(-2.0078134) q[0];
rz(-2.1009905) q[1];
sx q[1];
rz(-0.84747172) q[1];
sx q[1];
rz(-0.52451959) q[1];
rz(1.9807627) q[2];
sx q[2];
rz(-0.66902918) q[2];
sx q[2];
rz(1.4794028) q[2];
rz(-2.5083422) q[3];
sx q[3];
rz(-1.6926822) q[3];
sx q[3];
rz(2.5183986) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
