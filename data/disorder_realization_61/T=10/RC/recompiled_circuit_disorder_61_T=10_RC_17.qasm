OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7282038) q[0];
sx q[0];
rz(-2.0079186) q[0];
sx q[0];
rz(-1.5925621) q[0];
rz(1.6917317) q[1];
sx q[1];
rz(5.6258968) q[1];
sx q[1];
rz(13.110553) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.005851) q[0];
sx q[0];
rz(-3.0648181) q[0];
sx q[0];
rz(-2.6430921) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.10138301) q[2];
sx q[2];
rz(-1.5268491) q[2];
sx q[2];
rz(-1.4506538) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.7218329) q[1];
sx q[1];
rz(-0.83336035) q[1];
sx q[1];
rz(-2.1268197) q[1];
rz(0.55239001) q[3];
sx q[3];
rz(-3.1079709) q[3];
sx q[3];
rz(0.53510964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.071775285) q[2];
sx q[2];
rz(-1.2640307) q[2];
sx q[2];
rz(1.3624181) q[2];
rz(-3.1133364) q[3];
sx q[3];
rz(-1.7621721) q[3];
sx q[3];
rz(0.79022592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64489275) q[0];
sx q[0];
rz(-0.70514482) q[0];
sx q[0];
rz(-1.171296) q[0];
rz(0.21121875) q[1];
sx q[1];
rz(-2.6995081) q[1];
sx q[1];
rz(-1.404095) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0140186) q[0];
sx q[0];
rz(-1.9919792) q[0];
sx q[0];
rz(0.76571) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8318769) q[2];
sx q[2];
rz(-0.5030015) q[2];
sx q[2];
rz(2.3878857) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1120226) q[1];
sx q[1];
rz(-0.28407541) q[1];
sx q[1];
rz(3.0594538) q[1];
rz(-pi) q[2];
rz(2.6614463) q[3];
sx q[3];
rz(-1.4501791) q[3];
sx q[3];
rz(0.13531019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.30963787) q[2];
sx q[2];
rz(-1.6654623) q[2];
sx q[2];
rz(-2.1976166) q[2];
rz(-2.5850463) q[3];
sx q[3];
rz(-2.6440933) q[3];
sx q[3];
rz(1.336162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.29207644) q[0];
sx q[0];
rz(-0.76382604) q[0];
sx q[0];
rz(1.4235494) q[0];
rz(-0.81958333) q[1];
sx q[1];
rz(-0.81033605) q[1];
sx q[1];
rz(-0.56366411) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1492796) q[0];
sx q[0];
rz(-0.57846071) q[0];
sx q[0];
rz(3.0119386) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3450422) q[2];
sx q[2];
rz(-0.98276897) q[2];
sx q[2];
rz(-3.0622481) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7970265) q[1];
sx q[1];
rz(-2.1246506) q[1];
sx q[1];
rz(-2.7558541) q[1];
x q[2];
rz(-2.5205069) q[3];
sx q[3];
rz(-1.9120145) q[3];
sx q[3];
rz(0.065545557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6391969) q[2];
sx q[2];
rz(-2.0194619) q[2];
sx q[2];
rz(-1.0602661) q[2];
rz(-0.075573102) q[3];
sx q[3];
rz(-2.2836756) q[3];
sx q[3];
rz(-3.0269572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.32325) q[0];
sx q[0];
rz(-0.30775726) q[0];
sx q[0];
rz(-0.19317214) q[0];
rz(-1.5441783) q[1];
sx q[1];
rz(-1.9924106) q[1];
sx q[1];
rz(-2.4838122) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5415216) q[0];
sx q[0];
rz(-2.9642448) q[0];
sx q[0];
rz(0.17139165) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31099702) q[2];
sx q[2];
rz(-1.1592835) q[2];
sx q[2];
rz(-2.3222773) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.12280497) q[1];
sx q[1];
rz(-2.0969166) q[1];
sx q[1];
rz(0.29463525) q[1];
rz(-pi) q[2];
x q[2];
rz(0.85122078) q[3];
sx q[3];
rz(-2.7606066) q[3];
sx q[3];
rz(1.3850497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0458935) q[2];
sx q[2];
rz(-0.4102439) q[2];
sx q[2];
rz(2.0945385) q[2];
rz(-1.0632769) q[3];
sx q[3];
rz(-1.1626817) q[3];
sx q[3];
rz(-1.261238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3445774) q[0];
sx q[0];
rz(-2.8224967) q[0];
sx q[0];
rz(-2.1176594) q[0];
rz(2.3539885) q[1];
sx q[1];
rz(-1.5658295) q[1];
sx q[1];
rz(-0.62686282) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.244827) q[0];
sx q[0];
rz(-1.5540431) q[0];
sx q[0];
rz(-1.6003952) q[0];
rz(-1.2704029) q[2];
sx q[2];
rz(-1.8744178) q[2];
sx q[2];
rz(2.6189569) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7879466) q[1];
sx q[1];
rz(-1.8783356) q[1];
sx q[1];
rz(-2.8278082) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9634788) q[3];
sx q[3];
rz(-1.6851478) q[3];
sx q[3];
rz(-1.8161731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2287067) q[2];
sx q[2];
rz(-1.9084946) q[2];
sx q[2];
rz(-2.1234925) q[2];
rz(0.61156887) q[3];
sx q[3];
rz(-1.8194149) q[3];
sx q[3];
rz(-0.95156804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(1.6927032) q[0];
sx q[0];
rz(-1.3494116) q[0];
sx q[0];
rz(-1.942379) q[0];
rz(1.9723643) q[1];
sx q[1];
rz(-1.0511845) q[1];
sx q[1];
rz(-2.8170524) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6231461) q[0];
sx q[0];
rz(-1.070797) q[0];
sx q[0];
rz(0.48796939) q[0];
rz(2.2424477) q[2];
sx q[2];
rz(-0.88485826) q[2];
sx q[2];
rz(-1.1555954) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.74450103) q[1];
sx q[1];
rz(-0.53452864) q[1];
sx q[1];
rz(1.6194653) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8683897) q[3];
sx q[3];
rz(-1.9801567) q[3];
sx q[3];
rz(-0.62664947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.67363182) q[2];
sx q[2];
rz(-1.8362074) q[2];
sx q[2];
rz(1.8661873) q[2];
rz(1.0547137) q[3];
sx q[3];
rz(-0.79189363) q[3];
sx q[3];
rz(-1.595343) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4642898) q[0];
sx q[0];
rz(-1.807656) q[0];
sx q[0];
rz(1.4087079) q[0];
rz(2.6858221) q[1];
sx q[1];
rz(-0.20142889) q[1];
sx q[1];
rz(1.2021525) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81284886) q[0];
sx q[0];
rz(-2.3658731) q[0];
sx q[0];
rz(0.98529718) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9583148) q[2];
sx q[2];
rz(-0.85867184) q[2];
sx q[2];
rz(2.9930263) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.89029965) q[1];
sx q[1];
rz(-1.6541462) q[1];
sx q[1];
rz(0.092373089) q[1];
rz(-2.5382302) q[3];
sx q[3];
rz(-2.8238736) q[3];
sx q[3];
rz(2.7226504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1414286) q[2];
sx q[2];
rz(-1.8464073) q[2];
sx q[2];
rz(2.2953575) q[2];
rz(2.6464461) q[3];
sx q[3];
rz(-2.4954002) q[3];
sx q[3];
rz(-1.600986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0705567) q[0];
sx q[0];
rz(-1.8510011) q[0];
sx q[0];
rz(-1.113168) q[0];
rz(2.44599) q[1];
sx q[1];
rz(-0.38989392) q[1];
sx q[1];
rz(-1.5323458) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0198062) q[0];
sx q[0];
rz(-1.1018503) q[0];
sx q[0];
rz(-1.8437587) q[0];
rz(2.3959827) q[2];
sx q[2];
rz(-2.5157305) q[2];
sx q[2];
rz(1.9288043) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6541877) q[1];
sx q[1];
rz(-2.020917) q[1];
sx q[1];
rz(-0.022189157) q[1];
rz(2.5395457) q[3];
sx q[3];
rz(-1.4730519) q[3];
sx q[3];
rz(2.5412113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9939076) q[2];
sx q[2];
rz(-2.9116178) q[2];
sx q[2];
rz(0.85419401) q[2];
rz(-1.2285852) q[3];
sx q[3];
rz(-1.5649786) q[3];
sx q[3];
rz(-1.6555697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5937186) q[0];
sx q[0];
rz(-0.49867189) q[0];
sx q[0];
rz(0.6643995) q[0];
rz(1.1876855) q[1];
sx q[1];
rz(-2.4408051) q[1];
sx q[1];
rz(-1.2845576) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8878471) q[0];
sx q[0];
rz(-0.77639025) q[0];
sx q[0];
rz(-1.9842158) q[0];
rz(0.2662439) q[2];
sx q[2];
rz(-1.2028482) q[2];
sx q[2];
rz(-3.0014696) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5018651) q[1];
sx q[1];
rz(-2.0149391) q[1];
sx q[1];
rz(0.72413866) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7515718) q[3];
sx q[3];
rz(-1.5961002) q[3];
sx q[3];
rz(0.43061531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5294042) q[2];
sx q[2];
rz(-2.2470784) q[2];
sx q[2];
rz(-2.5908296) q[2];
rz(2.4380056) q[3];
sx q[3];
rz(-2.1313322) q[3];
sx q[3];
rz(-1.5065058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.7228912) q[0];
sx q[0];
rz(-0.35690618) q[0];
sx q[0];
rz(-0.044145949) q[0];
rz(-1.528953) q[1];
sx q[1];
rz(-1.2395369) q[1];
sx q[1];
rz(-0.77967656) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89833242) q[0];
sx q[0];
rz(-1.7120449) q[0];
sx q[0];
rz(2.6000573) q[0];
rz(-pi) q[1];
rz(2.0585971) q[2];
sx q[2];
rz(-2.3239845) q[2];
sx q[2];
rz(-1.4873193) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.838678) q[1];
sx q[1];
rz(-2.089114) q[1];
sx q[1];
rz(1.360449) q[1];
rz(-2.036282) q[3];
sx q[3];
rz(-0.84565425) q[3];
sx q[3];
rz(0.5474962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6471275) q[2];
sx q[2];
rz(-0.72138849) q[2];
sx q[2];
rz(-1.54281) q[2];
rz(0.70458448) q[3];
sx q[3];
rz(-1.5141809) q[3];
sx q[3];
rz(-3.1295479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71173944) q[0];
sx q[0];
rz(-1.585351) q[0];
sx q[0];
rz(2.4899695) q[0];
rz(-1.343887) q[1];
sx q[1];
rz(-1.4812891) q[1];
sx q[1];
rz(-0.67566009) q[1];
rz(-2.1716739) q[2];
sx q[2];
rz(-0.51545943) q[2];
sx q[2];
rz(-3.1209844) q[2];
rz(-0.17776168) q[3];
sx q[3];
rz(-2.4958785) q[3];
sx q[3];
rz(2.6444825) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
