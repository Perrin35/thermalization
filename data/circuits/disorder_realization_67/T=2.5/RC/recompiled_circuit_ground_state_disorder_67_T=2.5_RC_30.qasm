OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.32970348) q[0];
sx q[0];
rz(-2.831037) q[0];
sx q[0];
rz(1.9429053) q[0];
rz(-2.0074453) q[1];
sx q[1];
rz(2.3448047) q[1];
sx q[1];
rz(12.264836) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3350007) q[0];
sx q[0];
rz(-1.7552688) q[0];
sx q[0];
rz(-2.9770699) q[0];
rz(2.9542175) q[2];
sx q[2];
rz(-1.865554) q[2];
sx q[2];
rz(2.8127138) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1323213) q[1];
sx q[1];
rz(-1.5030716) q[1];
sx q[1];
rz(0.81791116) q[1];
x q[2];
rz(-1.9925549) q[3];
sx q[3];
rz(-2.7054477) q[3];
sx q[3];
rz(-0.25806043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.6887168) q[2];
sx q[2];
rz(-1.0998925) q[2];
sx q[2];
rz(-2.0067046) q[2];
rz(-3.0196043) q[3];
sx q[3];
rz(-2.205409) q[3];
sx q[3];
rz(0.056099135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6720471) q[0];
sx q[0];
rz(-2.8944954) q[0];
sx q[0];
rz(0.0023512996) q[0];
rz(2.7408842) q[1];
sx q[1];
rz(-1.3688764) q[1];
sx q[1];
rz(0.8561264) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1416505) q[0];
sx q[0];
rz(-0.68075962) q[0];
sx q[0];
rz(1.9826243) q[0];
x q[1];
rz(2.8555774) q[2];
sx q[2];
rz(-2.0692424) q[2];
sx q[2];
rz(1.1666544) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2344831) q[1];
sx q[1];
rz(-0.52965783) q[1];
sx q[1];
rz(1.6698599) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9366185) q[3];
sx q[3];
rz(-1.0351887) q[3];
sx q[3];
rz(0.35282319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.037228435) q[2];
sx q[2];
rz(-1.9375485) q[2];
sx q[2];
rz(2.6715703) q[2];
rz(2.5445599) q[3];
sx q[3];
rz(-3.0266914) q[3];
sx q[3];
rz(1.443583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7715348) q[0];
sx q[0];
rz(-0.4929339) q[0];
sx q[0];
rz(2.9135627) q[0];
rz(-2.6920964) q[1];
sx q[1];
rz(-2.1411965) q[1];
sx q[1];
rz(2.1953348) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3920604) q[0];
sx q[0];
rz(-0.48198732) q[0];
sx q[0];
rz(-0.9922235) q[0];
rz(-pi) q[1];
rz(-2.4837355) q[2];
sx q[2];
rz(-1.0491174) q[2];
sx q[2];
rz(1.9470095) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1552004) q[1];
sx q[1];
rz(-2.7554535) q[1];
sx q[1];
rz(-2.6394352) q[1];
rz(1.069293) q[3];
sx q[3];
rz(-2.8208102) q[3];
sx q[3];
rz(-0.71170413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2366644) q[2];
sx q[2];
rz(-1.6764574) q[2];
sx q[2];
rz(-1.4795335) q[2];
rz(0.18105257) q[3];
sx q[3];
rz(-2.3513887) q[3];
sx q[3];
rz(-2.2553867) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7774696) q[0];
sx q[0];
rz(-1.5793707) q[0];
sx q[0];
rz(0.89853483) q[0];
rz(-0.50615519) q[1];
sx q[1];
rz(-1.8471142) q[1];
sx q[1];
rz(1.6815965) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5622957) q[0];
sx q[0];
rz(-2.102872) q[0];
sx q[0];
rz(1.4946372) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3982749) q[2];
sx q[2];
rz(-1.905059) q[2];
sx q[2];
rz(-1.0570132) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0713363) q[1];
sx q[1];
rz(-1.6410636) q[1];
sx q[1];
rz(1.6855168) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.91588275) q[3];
sx q[3];
rz(-1.5228027) q[3];
sx q[3];
rz(2.1700883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8755181) q[2];
sx q[2];
rz(-2.6660599) q[2];
sx q[2];
rz(1.5024705) q[2];
rz(1.8858887) q[3];
sx q[3];
rz(-1.7399104) q[3];
sx q[3];
rz(0.046253117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
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
rz(-1.9083549) q[0];
sx q[0];
rz(-2.4112356) q[0];
sx q[0];
rz(2.9610942) q[0];
rz(-1.3795229) q[1];
sx q[1];
rz(-0.5677529) q[1];
sx q[1];
rz(-2.0464121) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3334072) q[0];
sx q[0];
rz(-1.6327724) q[0];
sx q[0];
rz(0.13570762) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6542475) q[2];
sx q[2];
rz(-2.3660283) q[2];
sx q[2];
rz(-0.47093876) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4597561) q[1];
sx q[1];
rz(-0.78826681) q[1];
sx q[1];
rz(-0.19758176) q[1];
rz(1.4093231) q[3];
sx q[3];
rz(-0.69500178) q[3];
sx q[3];
rz(-1.410786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3695662) q[2];
sx q[2];
rz(-0.84731805) q[2];
sx q[2];
rz(-2.0728716) q[2];
rz(-2.7145743) q[3];
sx q[3];
rz(-0.87978274) q[3];
sx q[3];
rz(-1.8900185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1123947) q[0];
sx q[0];
rz(-2.2258832) q[0];
sx q[0];
rz(2.8503964) q[0];
rz(-1.4578106) q[1];
sx q[1];
rz(-1.0271415) q[1];
sx q[1];
rz(-0.88868946) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5989953) q[0];
sx q[0];
rz(-0.7805191) q[0];
sx q[0];
rz(2.2602343) q[0];
rz(-0.75727792) q[2];
sx q[2];
rz(-1.5520397) q[2];
sx q[2];
rz(-2.4919486) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2803964) q[1];
sx q[1];
rz(-0.79842192) q[1];
sx q[1];
rz(2.3905928) q[1];
rz(2.8392516) q[3];
sx q[3];
rz(-1.8714222) q[3];
sx q[3];
rz(-1.9019097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7500744) q[2];
sx q[2];
rz(-1.3041648) q[2];
sx q[2];
rz(2.1938426) q[2];
rz(-2.9446972) q[3];
sx q[3];
rz(-2.9010549) q[3];
sx q[3];
rz(-2.8314364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8295558) q[0];
sx q[0];
rz(-1.1126248) q[0];
sx q[0];
rz(-0.46992508) q[0];
rz(1.6795109) q[1];
sx q[1];
rz(-1.3009678) q[1];
sx q[1];
rz(-1.7376815) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3199661) q[0];
sx q[0];
rz(-1.2228773) q[0];
sx q[0];
rz(-2.1676284) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0491046) q[2];
sx q[2];
rz(-1.7564536) q[2];
sx q[2];
rz(0.80362475) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0216071) q[1];
sx q[1];
rz(-0.65896265) q[1];
sx q[1];
rz(0.95665415) q[1];
rz(-pi) q[2];
rz(-2.5581854) q[3];
sx q[3];
rz(-1.2800763) q[3];
sx q[3];
rz(-0.099778508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5002084) q[2];
sx q[2];
rz(-2.310014) q[2];
sx q[2];
rz(-0.10759648) q[2];
rz(-0.61947668) q[3];
sx q[3];
rz(-1.2958801) q[3];
sx q[3];
rz(1.7776325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7463995) q[0];
sx q[0];
rz(-1.0522319) q[0];
sx q[0];
rz(-2.8885544) q[0];
rz(-3.053275) q[1];
sx q[1];
rz(-2.2679236) q[1];
sx q[1];
rz(0.94892445) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2491566) q[0];
sx q[0];
rz(-1.2741486) q[0];
sx q[0];
rz(2.9846624) q[0];
rz(-pi) q[1];
rz(1.1626271) q[2];
sx q[2];
rz(-1.3917096) q[2];
sx q[2];
rz(-0.31291747) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9946549) q[1];
sx q[1];
rz(-0.48284621) q[1];
sx q[1];
rz(1.6864683) q[1];
rz(-2.6213245) q[3];
sx q[3];
rz(-2.0754077) q[3];
sx q[3];
rz(-2.5443973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.39763149) q[2];
sx q[2];
rz(-2.2546015) q[2];
sx q[2];
rz(-2.8308947) q[2];
rz(-2.9001696) q[3];
sx q[3];
rz(-1.2025236) q[3];
sx q[3];
rz(-1.0827304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13984118) q[0];
sx q[0];
rz(-2.7353291) q[0];
sx q[0];
rz(-1.9061506) q[0];
rz(2.6605117) q[1];
sx q[1];
rz(-2.1886539) q[1];
sx q[1];
rz(-1.7135886) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22132231) q[0];
sx q[0];
rz(-2.3820665) q[0];
sx q[0];
rz(-0.016245441) q[0];
x q[1];
rz(1.2920824) q[2];
sx q[2];
rz(-1.3506417) q[2];
sx q[2];
rz(-2.1627592) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7117611) q[1];
sx q[1];
rz(-1.3138384) q[1];
sx q[1];
rz(2.8896232) q[1];
rz(-pi) q[2];
rz(0.72254027) q[3];
sx q[3];
rz(-1.819639) q[3];
sx q[3];
rz(1.794508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1741751) q[2];
sx q[2];
rz(-0.82505161) q[2];
sx q[2];
rz(3.1136801) q[2];
rz(0.86999718) q[3];
sx q[3];
rz(-0.78012192) q[3];
sx q[3];
rz(1.9546668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.025539909) q[0];
sx q[0];
rz(-1.5929796) q[0];
sx q[0];
rz(-1.0593587) q[0];
rz(-1.7268044) q[1];
sx q[1];
rz(-1.4780412) q[1];
sx q[1];
rz(-2.6639604) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57741149) q[0];
sx q[0];
rz(-1.7192657) q[0];
sx q[0];
rz(-0.74651697) q[0];
x q[1];
rz(-1.5651836) q[2];
sx q[2];
rz(-0.33585784) q[2];
sx q[2];
rz(0.72959057) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2207909) q[1];
sx q[1];
rz(-0.12462908) q[1];
sx q[1];
rz(-1.1225379) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.82179339) q[3];
sx q[3];
rz(-1.104165) q[3];
sx q[3];
rz(2.1958283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.16964218) q[2];
sx q[2];
rz(-0.98196882) q[2];
sx q[2];
rz(-2.7726445) q[2];
rz(-1.2667027) q[3];
sx q[3];
rz(-0.72487512) q[3];
sx q[3];
rz(-0.96316159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13851588) q[0];
sx q[0];
rz(-1.930548) q[0];
sx q[0];
rz(1.7532274) q[0];
rz(-2.6970462) q[1];
sx q[1];
rz(-2.3311756) q[1];
sx q[1];
rz(-0.12745007) q[1];
rz(1.2426022) q[2];
sx q[2];
rz(-1.9933619) q[2];
sx q[2];
rz(0.91771916) q[2];
rz(1.5515242) q[3];
sx q[3];
rz(-1.1033333) q[3];
sx q[3];
rz(3.0126572) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
