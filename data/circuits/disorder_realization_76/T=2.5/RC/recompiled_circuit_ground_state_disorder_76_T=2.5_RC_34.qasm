OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0903836) q[0];
sx q[0];
rz(-2.3530355) q[0];
sx q[0];
rz(0.30858827) q[0];
rz(-2.7170972) q[1];
sx q[1];
rz(-1.0841882) q[1];
sx q[1];
rz(0.37556136) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9449995) q[0];
sx q[0];
rz(-0.29175943) q[0];
sx q[0];
rz(-1.8367405) q[0];
x q[1];
rz(-0.1027561) q[2];
sx q[2];
rz(-0.30347541) q[2];
sx q[2];
rz(0.048852531) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6172609) q[1];
sx q[1];
rz(-1.4325805) q[1];
sx q[1];
rz(2.6928965) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4272903) q[3];
sx q[3];
rz(-1.5745899) q[3];
sx q[3];
rz(3.0285545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3382807) q[2];
sx q[2];
rz(-2.8901926) q[2];
sx q[2];
rz(-2.5982762) q[2];
rz(-1.4079037) q[3];
sx q[3];
rz(-1.7270154) q[3];
sx q[3];
rz(2.2442815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59178418) q[0];
sx q[0];
rz(-1.1487288) q[0];
sx q[0];
rz(2.5445004) q[0];
rz(1.0626556) q[1];
sx q[1];
rz(-2.5025044) q[1];
sx q[1];
rz(-2.7426681) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8455578) q[0];
sx q[0];
rz(-0.71656967) q[0];
sx q[0];
rz(-0.96591732) q[0];
x q[1];
rz(0.024256134) q[2];
sx q[2];
rz(-2.127671) q[2];
sx q[2];
rz(-0.66489894) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0831757) q[1];
sx q[1];
rz(-0.99715573) q[1];
sx q[1];
rz(0.87223335) q[1];
x q[2];
rz(0.53360231) q[3];
sx q[3];
rz(-2.1416365) q[3];
sx q[3];
rz(-0.10122964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4661633) q[2];
sx q[2];
rz(-1.9025981) q[2];
sx q[2];
rz(-2.0002401) q[2];
rz(-2.2232248) q[3];
sx q[3];
rz(-0.71056241) q[3];
sx q[3];
rz(0.58951283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.454575) q[0];
sx q[0];
rz(-2.7641251) q[0];
sx q[0];
rz(-0.17534176) q[0];
rz(-1.8270127) q[1];
sx q[1];
rz(-1.8914696) q[1];
sx q[1];
rz(0.66741991) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8274067) q[0];
sx q[0];
rz(-1.6860322) q[0];
sx q[0];
rz(-2.0699571) q[0];
rz(-1.0834789) q[2];
sx q[2];
rz(-2.1082775) q[2];
sx q[2];
rz(-1.4097241) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.35386723) q[1];
sx q[1];
rz(-1.8866061) q[1];
sx q[1];
rz(-1.5701957) q[1];
rz(0.36881558) q[3];
sx q[3];
rz(-2.494156) q[3];
sx q[3];
rz(-1.8819229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1158925) q[2];
sx q[2];
rz(-0.95558715) q[2];
sx q[2];
rz(-1.0828177) q[2];
rz(1.8358021) q[3];
sx q[3];
rz(-2.3256153) q[3];
sx q[3];
rz(-0.79735565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1378655) q[0];
sx q[0];
rz(-0.19345134) q[0];
sx q[0];
rz(2.5306012) q[0];
rz(-2.5773279) q[1];
sx q[1];
rz(-1.0287501) q[1];
sx q[1];
rz(0.70704031) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63227292) q[0];
sx q[0];
rz(-1.5406666) q[0];
sx q[0];
rz(2.3448461) q[0];
x q[1];
rz(-0.41019812) q[2];
sx q[2];
rz(-1.8767523) q[2];
sx q[2];
rz(-3.0434639) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9663716) q[1];
sx q[1];
rz(-2.4680158) q[1];
sx q[1];
rz(-0.74209706) q[1];
rz(-pi) q[2];
rz(-0.50102542) q[3];
sx q[3];
rz(-1.820537) q[3];
sx q[3];
rz(-1.4755032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.45346144) q[2];
sx q[2];
rz(-0.42669272) q[2];
sx q[2];
rz(1.0180265) q[2];
rz(0.38797837) q[3];
sx q[3];
rz(-1.6951025) q[3];
sx q[3];
rz(2.80262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6729048) q[0];
sx q[0];
rz(-1.1274811) q[0];
sx q[0];
rz(-0.52184033) q[0];
rz(-0.28779596) q[1];
sx q[1];
rz(-2.5591873) q[1];
sx q[1];
rz(-0.60307455) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48846196) q[0];
sx q[0];
rz(-1.483027) q[0];
sx q[0];
rz(1.5215389) q[0];
rz(-2.1289044) q[2];
sx q[2];
rz(-0.96622045) q[2];
sx q[2];
rz(-0.1451491) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9085616) q[1];
sx q[1];
rz(-2.7082293) q[1];
sx q[1];
rz(-2.9410081) q[1];
rz(-pi) q[2];
rz(-1.9136004) q[3];
sx q[3];
rz(-1.8209642) q[3];
sx q[3];
rz(-0.34473333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7120984) q[2];
sx q[2];
rz(-2.3488729) q[2];
sx q[2];
rz(0.73954868) q[2];
rz(-1.0150917) q[3];
sx q[3];
rz(-0.52777094) q[3];
sx q[3];
rz(-0.019006193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5571112) q[0];
sx q[0];
rz(-2.8040573) q[0];
sx q[0];
rz(0.8648411) q[0];
rz(-2.5014014) q[1];
sx q[1];
rz(-0.65018153) q[1];
sx q[1];
rz(2.6472299) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.664864) q[0];
sx q[0];
rz(-2.6458178) q[0];
sx q[0];
rz(1.134111) q[0];
x q[1];
rz(0.96832595) q[2];
sx q[2];
rz(-1.5366388) q[2];
sx q[2];
rz(-2.5631225) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3013005) q[1];
sx q[1];
rz(-1.6052142) q[1];
sx q[1];
rz(3.084206) q[1];
rz(0.5036656) q[3];
sx q[3];
rz(-2.0450838) q[3];
sx q[3];
rz(-2.4573457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.87259) q[2];
sx q[2];
rz(-2.5073017) q[2];
sx q[2];
rz(-2.1799901) q[2];
rz(-1.5087992) q[3];
sx q[3];
rz(-1.6274933) q[3];
sx q[3];
rz(-3.0646724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2831869) q[0];
sx q[0];
rz(-0.34342331) q[0];
sx q[0];
rz(-1.0431694) q[0];
rz(-2.4647602) q[1];
sx q[1];
rz(-2.4969641) q[1];
sx q[1];
rz(2.4136995) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2400144) q[0];
sx q[0];
rz(-1.9838617) q[0];
sx q[0];
rz(-0.1211285) q[0];
rz(-pi) q[1];
rz(3.137693) q[2];
sx q[2];
rz(-0.74206381) q[2];
sx q[2];
rz(-2.8067547) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.7866061) q[1];
sx q[1];
rz(-2.1713272) q[1];
sx q[1];
rz(2.7342466) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5520127) q[3];
sx q[3];
rz(-1.2834719) q[3];
sx q[3];
rz(0.048487566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4595043) q[2];
sx q[2];
rz(-2.3900718) q[2];
sx q[2];
rz(-0.53306836) q[2];
rz(-1.8367977) q[3];
sx q[3];
rz(-2.6663836) q[3];
sx q[3];
rz(-2.3615725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.099139) q[0];
sx q[0];
rz(-3.070153) q[0];
sx q[0];
rz(-0.024600994) q[0];
rz(-0.050361659) q[1];
sx q[1];
rz(-0.55614007) q[1];
sx q[1];
rz(0.76739001) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8302916) q[0];
sx q[0];
rz(-2.6008252) q[0];
sx q[0];
rz(-3.0917653) q[0];
rz(-pi) q[1];
rz(0.32979687) q[2];
sx q[2];
rz(-1.7539795) q[2];
sx q[2];
rz(-0.51815301) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.9995867) q[1];
sx q[1];
rz(-2.1697727) q[1];
sx q[1];
rz(1.0890278) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5663451) q[3];
sx q[3];
rz(-1.2367555) q[3];
sx q[3];
rz(0.0041242139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.47160992) q[2];
sx q[2];
rz(-2.6675384) q[2];
sx q[2];
rz(1.6101884) q[2];
rz(-2.1286185) q[3];
sx q[3];
rz(-1.6771202) q[3];
sx q[3];
rz(2.5392635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7811964) q[0];
sx q[0];
rz(-0.66806) q[0];
sx q[0];
rz(-0.49041954) q[0];
rz(2.7410653) q[1];
sx q[1];
rz(-2.782395) q[1];
sx q[1];
rz(-1.5706971) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2928481) q[0];
sx q[0];
rz(-1.2146753) q[0];
sx q[0];
rz(-2.5625021) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2987842) q[2];
sx q[2];
rz(-2.669816) q[2];
sx q[2];
rz(-2.9695784) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3048858) q[1];
sx q[1];
rz(-1.2212824) q[1];
sx q[1];
rz(-3.034449) q[1];
rz(-pi) q[2];
rz(2.3113046) q[3];
sx q[3];
rz(-2.2133996) q[3];
sx q[3];
rz(2.682529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.57764292) q[2];
sx q[2];
rz(-2.9190013) q[2];
sx q[2];
rz(-0.27601784) q[2];
rz(-0.66463071) q[3];
sx q[3];
rz(-2.1641927) q[3];
sx q[3];
rz(2.9450534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4589602) q[0];
sx q[0];
rz(-2.9660872) q[0];
sx q[0];
rz(-1.5631787) q[0];
rz(0.75421929) q[1];
sx q[1];
rz(-1.027532) q[1];
sx q[1];
rz(0.064090699) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0383324) q[0];
sx q[0];
rz(-1.1128367) q[0];
sx q[0];
rz(0.026491212) q[0];
rz(2.5446266) q[2];
sx q[2];
rz(-0.709049) q[2];
sx q[2];
rz(-2.2500791) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8386993) q[1];
sx q[1];
rz(-2.6461227) q[1];
sx q[1];
rz(1.3957146) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.80497165) q[3];
sx q[3];
rz(-2.0197778) q[3];
sx q[3];
rz(1.4624553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.85612193) q[2];
sx q[2];
rz(-0.28588122) q[2];
sx q[2];
rz(0.90606436) q[2];
rz(1.2178577) q[3];
sx q[3];
rz(-1.8702312) q[3];
sx q[3];
rz(-2.0876032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.333552) q[0];
sx q[0];
rz(-1.6186436) q[0];
sx q[0];
rz(2.9343395) q[0];
rz(2.8605657) q[1];
sx q[1];
rz(-1.3916176) q[1];
sx q[1];
rz(2.4236046) q[1];
rz(2.9502921) q[2];
sx q[2];
rz(-0.40861599) q[2];
sx q[2];
rz(1.7392639) q[2];
rz(-1.5850375) q[3];
sx q[3];
rz(-1.0938627) q[3];
sx q[3];
rz(0.59577019) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
