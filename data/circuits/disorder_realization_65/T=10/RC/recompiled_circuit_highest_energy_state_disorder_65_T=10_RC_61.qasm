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
rz(-2.6429472) q[0];
sx q[0];
rz(-0.54583025) q[0];
sx q[0];
rz(-3.0232865) q[0];
rz(-2.2220597) q[1];
sx q[1];
rz(-1.3032721) q[1];
sx q[1];
rz(1.2143171) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3271546) q[0];
sx q[0];
rz(-1.1390349) q[0];
sx q[0];
rz(-2.9080954) q[0];
rz(-2.1636398) q[2];
sx q[2];
rz(-1.0314157) q[2];
sx q[2];
rz(3.116089) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0561021) q[1];
sx q[1];
rz(-2.506665) q[1];
sx q[1];
rz(1.3391499) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8439775) q[3];
sx q[3];
rz(-1.5915378) q[3];
sx q[3];
rz(-1.0047877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.97293568) q[2];
sx q[2];
rz(-1.1309705) q[2];
sx q[2];
rz(-2.0449779) q[2];
rz(-0.24122572) q[3];
sx q[3];
rz(-1.3696407) q[3];
sx q[3];
rz(0.84478861) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089652561) q[0];
sx q[0];
rz(-1.0900499) q[0];
sx q[0];
rz(0.3845149) q[0];
rz(2.8334726) q[1];
sx q[1];
rz(-0.48079753) q[1];
sx q[1];
rz(0.62517977) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4368432) q[0];
sx q[0];
rz(-1.6024717) q[0];
sx q[0];
rz(0.28625536) q[0];
x q[1];
rz(-0.7169508) q[2];
sx q[2];
rz(-1.4745108) q[2];
sx q[2];
rz(-1.3614298) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2634266) q[1];
sx q[1];
rz(-1.485865) q[1];
sx q[1];
rz(-0.64966605) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1891648) q[3];
sx q[3];
rz(-1.980972) q[3];
sx q[3];
rz(2.4252714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.008931) q[2];
sx q[2];
rz(-2.2098358) q[2];
sx q[2];
rz(-2.8399732) q[2];
rz(0.53145069) q[3];
sx q[3];
rz(-2.2558432) q[3];
sx q[3];
rz(2.0187995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5136435) q[0];
sx q[0];
rz(-2.1301837) q[0];
sx q[0];
rz(-1.9568141) q[0];
rz(-2.9966677) q[1];
sx q[1];
rz(-1.8791608) q[1];
sx q[1];
rz(-1.5473993) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7892706) q[0];
sx q[0];
rz(-1.0836067) q[0];
sx q[0];
rz(-0.67729034) q[0];
rz(-pi) q[1];
rz(2.9582439) q[2];
sx q[2];
rz(-1.1969222) q[2];
sx q[2];
rz(0.093106769) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5963187) q[1];
sx q[1];
rz(-0.68877586) q[1];
sx q[1];
rz(0.94249814) q[1];
rz(-pi) q[2];
rz(1.2482485) q[3];
sx q[3];
rz(-1.9574021) q[3];
sx q[3];
rz(-2.4462388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3323815) q[2];
sx q[2];
rz(-0.96497649) q[2];
sx q[2];
rz(-1.8154209) q[2];
rz(2.2281846) q[3];
sx q[3];
rz(-1.923442) q[3];
sx q[3];
rz(0.53328812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7186385) q[0];
sx q[0];
rz(-2.664743) q[0];
sx q[0];
rz(-1.3612716) q[0];
rz(-0.71296972) q[1];
sx q[1];
rz(-2.0013335) q[1];
sx q[1];
rz(-2.4580809) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26091247) q[0];
sx q[0];
rz(-2.144038) q[0];
sx q[0];
rz(3.1295958) q[0];
rz(2.322898) q[2];
sx q[2];
rz(-0.3017738) q[2];
sx q[2];
rz(-1.2692598) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6660712) q[1];
sx q[1];
rz(-2.254535) q[1];
sx q[1];
rz(0.68902512) q[1];
x q[2];
rz(-2.8972647) q[3];
sx q[3];
rz(-1.667898) q[3];
sx q[3];
rz(1.3299143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6354562) q[2];
sx q[2];
rz(-2.7244302) q[2];
sx q[2];
rz(-0.48640856) q[2];
rz(-0.5168612) q[3];
sx q[3];
rz(-1.9603399) q[3];
sx q[3];
rz(-1.0335056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0559693) q[0];
sx q[0];
rz(-1.6317246) q[0];
sx q[0];
rz(1.3729209) q[0];
rz(1.4068475) q[1];
sx q[1];
rz(-0.56766784) q[1];
sx q[1];
rz(2.6645606) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6876887) q[0];
sx q[0];
rz(-0.46594187) q[0];
sx q[0];
rz(2.0201319) q[0];
rz(-pi) q[1];
rz(0.3933752) q[2];
sx q[2];
rz(-2.057339) q[2];
sx q[2];
rz(-1.7560619) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.20471599) q[1];
sx q[1];
rz(-0.5936247) q[1];
sx q[1];
rz(-2.1267664) q[1];
rz(-pi) q[2];
rz(-0.78272696) q[3];
sx q[3];
rz(-2.1313792) q[3];
sx q[3];
rz(-0.70817876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.10739747) q[2];
sx q[2];
rz(-0.98661462) q[2];
sx q[2];
rz(-2.6743215) q[2];
rz(2.9723736) q[3];
sx q[3];
rz(-1.7305814) q[3];
sx q[3];
rz(-0.30538487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6071103) q[0];
sx q[0];
rz(-0.87962532) q[0];
sx q[0];
rz(0.18950263) q[0];
rz(0.27319187) q[1];
sx q[1];
rz(-2.0676282) q[1];
sx q[1];
rz(-2.3409519) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9439745) q[0];
sx q[0];
rz(-1.5276143) q[0];
sx q[0];
rz(1.8042685) q[0];
rz(-0.24836274) q[2];
sx q[2];
rz(-1.6585322) q[2];
sx q[2];
rz(-1.9875634) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2792017) q[1];
sx q[1];
rz(-1.6995879) q[1];
sx q[1];
rz(-1.3736082) q[1];
x q[2];
rz(-0.54038911) q[3];
sx q[3];
rz(-2.0107993) q[3];
sx q[3];
rz(0.42127452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.93929401) q[2];
sx q[2];
rz(-0.17387986) q[2];
sx q[2];
rz(0.536971) q[2];
rz(1.8592853) q[3];
sx q[3];
rz(-1.3064281) q[3];
sx q[3];
rz(1.6622701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.031161664) q[0];
sx q[0];
rz(-0.81542504) q[0];
sx q[0];
rz(1.331331) q[0];
rz(1.4415461) q[1];
sx q[1];
rz(-1.4794289) q[1];
sx q[1];
rz(1.2225245) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5957907) q[0];
sx q[0];
rz(-1.8078601) q[0];
sx q[0];
rz(-1.9997826) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.020633) q[2];
sx q[2];
rz(-2.5865062) q[2];
sx q[2];
rz(0.70155376) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4536087) q[1];
sx q[1];
rz(-2.2609432) q[1];
sx q[1];
rz(-0.28120561) q[1];
x q[2];
rz(-1.8172713) q[3];
sx q[3];
rz(-0.45650864) q[3];
sx q[3];
rz(-2.498561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.20087251) q[2];
sx q[2];
rz(-0.49786374) q[2];
sx q[2];
rz(-0.3328003) q[2];
rz(0.29916549) q[3];
sx q[3];
rz(-1.9009512) q[3];
sx q[3];
rz(3.1201194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12896319) q[0];
sx q[0];
rz(-2.3733932) q[0];
sx q[0];
rz(2.8666038) q[0];
rz(-3.0601652) q[1];
sx q[1];
rz(-2.3402201) q[1];
sx q[1];
rz(1.8191232) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0819495) q[0];
sx q[0];
rz(-2.9026051) q[0];
sx q[0];
rz(1.8247402) q[0];
rz(-0.88365474) q[2];
sx q[2];
rz(-1.2809086) q[2];
sx q[2];
rz(-1.7366126) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9857707) q[1];
sx q[1];
rz(-0.81822526) q[1];
sx q[1];
rz(-2.6236412) q[1];
rz(-2.4149618) q[3];
sx q[3];
rz(-2.0803578) q[3];
sx q[3];
rz(-1.0973615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.72121173) q[2];
sx q[2];
rz(-3.0215441) q[2];
sx q[2];
rz(-1.7397286) q[2];
rz(-1.8574572) q[3];
sx q[3];
rz(-1.9522342) q[3];
sx q[3];
rz(0.0349667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5597124) q[0];
sx q[0];
rz(-1.5927097) q[0];
sx q[0];
rz(-2.6522719) q[0];
rz(3.1210461) q[1];
sx q[1];
rz(-1.9053883) q[1];
sx q[1];
rz(-0.70125854) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1398292) q[0];
sx q[0];
rz(-2.1368623) q[0];
sx q[0];
rz(2.8553814) q[0];
rz(-pi) q[1];
rz(2.2743938) q[2];
sx q[2];
rz(-1.5343621) q[2];
sx q[2];
rz(1.6619722) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.012246) q[1];
sx q[1];
rz(-1.80629) q[1];
sx q[1];
rz(-0.78462623) q[1];
x q[2];
rz(2.0370767) q[3];
sx q[3];
rz(-0.93887431) q[3];
sx q[3];
rz(1.5620269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4150841) q[2];
sx q[2];
rz(-2.5231611) q[2];
sx q[2];
rz(2.6611967) q[2];
rz(-1.4091617) q[3];
sx q[3];
rz(-1.3536072) q[3];
sx q[3];
rz(-2.8044146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82520634) q[0];
sx q[0];
rz(-1.9012863) q[0];
sx q[0];
rz(-2.7301042) q[0];
rz(1.8972137) q[1];
sx q[1];
rz(-1.139541) q[1];
sx q[1];
rz(0.32434514) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56670226) q[0];
sx q[0];
rz(-1.6234723) q[0];
sx q[0];
rz(-1.943669) q[0];
rz(0.34180832) q[2];
sx q[2];
rz(-1.9035619) q[2];
sx q[2];
rz(1.3403128) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.082338) q[1];
sx q[1];
rz(-1.6538701) q[1];
sx q[1];
rz(-2.1359813) q[1];
x q[2];
rz(0.12121157) q[3];
sx q[3];
rz(-1.5577661) q[3];
sx q[3];
rz(-2.7708294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9559429) q[2];
sx q[2];
rz(-2.0273384) q[2];
sx q[2];
rz(-3.0246217) q[2];
rz(-0.95514917) q[3];
sx q[3];
rz(-0.17894608) q[3];
sx q[3];
rz(0.54168934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87880001) q[0];
sx q[0];
rz(-0.87651064) q[0];
sx q[0];
rz(0.87690092) q[0];
rz(1.1204489) q[1];
sx q[1];
rz(-0.81605492) q[1];
sx q[1];
rz(1.1647404) q[1];
rz(1.9399613) q[2];
sx q[2];
rz(-2.6045447) q[2];
sx q[2];
rz(2.4640026) q[2];
rz(-3.1163327) q[3];
sx q[3];
rz(-1.8413481) q[3];
sx q[3];
rz(-3.1006364) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
