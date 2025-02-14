OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.1104133) q[0];
sx q[0];
rz(-2.2194982) q[0];
sx q[0];
rz(-2.3715012) q[0];
rz(0.71647477) q[1];
sx q[1];
rz(4.0840277) q[1];
sx q[1];
rz(10.066636) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7866369) q[0];
sx q[0];
rz(-1.4572976) q[0];
sx q[0];
rz(-2.9391975) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0014512295) q[2];
sx q[2];
rz(-1.5719169) q[2];
sx q[2];
rz(-1.4933153) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7388307) q[1];
sx q[1];
rz(-1.2114727) q[1];
sx q[1];
rz(1.4043384) q[1];
x q[2];
rz(-2.558488) q[3];
sx q[3];
rz(-2.2197086) q[3];
sx q[3];
rz(2.976604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.98051071) q[2];
sx q[2];
rz(-2.0445721) q[2];
sx q[2];
rz(0.80992997) q[2];
rz(0.03446456) q[3];
sx q[3];
rz(-2.4813215) q[3];
sx q[3];
rz(-3.0380429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63475364) q[0];
sx q[0];
rz(-2.6320808) q[0];
sx q[0];
rz(-2.4822045) q[0];
rz(-1.4393282) q[1];
sx q[1];
rz(-1.5123475) q[1];
sx q[1];
rz(-2.6606681) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.07671) q[0];
sx q[0];
rz(-2.0079029) q[0];
sx q[0];
rz(2.9124898) q[0];
x q[1];
rz(2.4056817) q[2];
sx q[2];
rz(-2.1440268) q[2];
sx q[2];
rz(-0.57511759) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3673076) q[1];
sx q[1];
rz(-1.3626422) q[1];
sx q[1];
rz(-0.03117301) q[1];
rz(-1.102785) q[3];
sx q[3];
rz(-2.0986522) q[3];
sx q[3];
rz(0.46609391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3769569) q[2];
sx q[2];
rz(-2.3048293) q[2];
sx q[2];
rz(-2.5118206) q[2];
rz(1.1791641) q[3];
sx q[3];
rz(-0.7032913) q[3];
sx q[3];
rz(1.111697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4513007) q[0];
sx q[0];
rz(-3.1307104) q[0];
sx q[0];
rz(-2.1731398) q[0];
rz(-3.0015266) q[1];
sx q[1];
rz(-1.7882971) q[1];
sx q[1];
rz(2.5586939) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.934865) q[0];
sx q[0];
rz(-1.7375653) q[0];
sx q[0];
rz(1.616448) q[0];
rz(2.3034322) q[2];
sx q[2];
rz(-1.7485022) q[2];
sx q[2];
rz(0.017824307) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0576833) q[1];
sx q[1];
rz(-1.0094125) q[1];
sx q[1];
rz(-2.3777005) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2887313) q[3];
sx q[3];
rz(-0.4133458) q[3];
sx q[3];
rz(1.0173544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.53050238) q[2];
sx q[2];
rz(-0.077433057) q[2];
sx q[2];
rz(2.9275242) q[2];
rz(0.30238447) q[3];
sx q[3];
rz(-0.74367911) q[3];
sx q[3];
rz(2.7157057) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18836235) q[0];
sx q[0];
rz(-2.132405) q[0];
sx q[0];
rz(-2.0971712) q[0];
rz(2.2638679) q[1];
sx q[1];
rz(-1.5889771) q[1];
sx q[1];
rz(2.9728319) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.343086) q[0];
sx q[0];
rz(-1.0765392) q[0];
sx q[0];
rz(2.8841208) q[0];
x q[1];
rz(1.243882) q[2];
sx q[2];
rz(-1.765328) q[2];
sx q[2];
rz(-0.036613001) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5280594) q[1];
sx q[1];
rz(-2.4984887) q[1];
sx q[1];
rz(-0.36524857) q[1];
rz(-1.5215643) q[3];
sx q[3];
rz(-2.4353046) q[3];
sx q[3];
rz(-0.8362706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7793444) q[2];
sx q[2];
rz(-1.2835953) q[2];
sx q[2];
rz(2.0812422) q[2];
rz(-2.849071) q[3];
sx q[3];
rz(-0.7258324) q[3];
sx q[3];
rz(-0.084107548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58920687) q[0];
sx q[0];
rz(-2.4942106) q[0];
sx q[0];
rz(2.2733083) q[0];
rz(-2.5153416) q[1];
sx q[1];
rz(-1.3887082) q[1];
sx q[1];
rz(-1.3583604) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5007965) q[0];
sx q[0];
rz(-1.8587041) q[0];
sx q[0];
rz(-1.663371) q[0];
x q[1];
rz(-0.38154885) q[2];
sx q[2];
rz(-1.3693491) q[2];
sx q[2];
rz(-0.82009456) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.89857453) q[1];
sx q[1];
rz(-2.822091) q[1];
sx q[1];
rz(-0.77844365) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.503848) q[3];
sx q[3];
rz(-1.9336091) q[3];
sx q[3];
rz(2.2408298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2061578) q[2];
sx q[2];
rz(-0.81659603) q[2];
sx q[2];
rz(-2.6182776) q[2];
rz(-2.7715136) q[3];
sx q[3];
rz(-0.78032929) q[3];
sx q[3];
rz(1.7962598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0321781) q[0];
sx q[0];
rz(-2.9258756) q[0];
sx q[0];
rz(0.35414645) q[0];
rz(0.94193637) q[1];
sx q[1];
rz(-1.4553921) q[1];
sx q[1];
rz(-1.3166434) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9958651) q[0];
sx q[0];
rz(-2.6168129) q[0];
sx q[0];
rz(-2.9324233) q[0];
rz(-0.022158547) q[2];
sx q[2];
rz(-2.5922734) q[2];
sx q[2];
rz(1.4469128) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2159913) q[1];
sx q[1];
rz(-2.9570438) q[1];
sx q[1];
rz(-0.33949236) q[1];
x q[2];
rz(1.8527669) q[3];
sx q[3];
rz(-2.1957955) q[3];
sx q[3];
rz(2.9491912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3256623) q[2];
sx q[2];
rz(-0.38620913) q[2];
sx q[2];
rz(0.33829921) q[2];
rz(0.47652388) q[3];
sx q[3];
rz(-2.3881113) q[3];
sx q[3];
rz(-2.8009955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4395831) q[0];
sx q[0];
rz(-1.5832573) q[0];
sx q[0];
rz(2.6038792) q[0];
rz(-1.733755) q[1];
sx q[1];
rz(-2.6815963) q[1];
sx q[1];
rz(0.62526155) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1921034) q[0];
sx q[0];
rz(-1.92983) q[0];
sx q[0];
rz(-0.60199317) q[0];
rz(0.12282108) q[2];
sx q[2];
rz(-0.54219699) q[2];
sx q[2];
rz(-1.922883) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.96237732) q[1];
sx q[1];
rz(-1.7853702) q[1];
sx q[1];
rz(-2.5291165) q[1];
rz(2.3654371) q[3];
sx q[3];
rz(-1.8381834) q[3];
sx q[3];
rz(-0.54477967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9739146) q[2];
sx q[2];
rz(-2.6748071) q[2];
sx q[2];
rz(-1.3803049) q[2];
rz(0.4365094) q[3];
sx q[3];
rz(-2.1197539) q[3];
sx q[3];
rz(0.71353394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1768782) q[0];
sx q[0];
rz(-0.62129337) q[0];
sx q[0];
rz(2.6253413) q[0];
rz(-0.77975726) q[1];
sx q[1];
rz(-2.1596491) q[1];
sx q[1];
rz(-1.0468743) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16653331) q[0];
sx q[0];
rz(-1.8691081) q[0];
sx q[0];
rz(1.4201643) q[0];
x q[1];
rz(-2.5604232) q[2];
sx q[2];
rz(-2.0072674) q[2];
sx q[2];
rz(1.798686) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6987933) q[1];
sx q[1];
rz(-1.6644786) q[1];
sx q[1];
rz(-3.0832861) q[1];
rz(-pi) q[2];
rz(0.67892142) q[3];
sx q[3];
rz(-2.553294) q[3];
sx q[3];
rz(-0.029840851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9376935) q[2];
sx q[2];
rz(-2.9475309) q[2];
sx q[2];
rz(-2.1843809) q[2];
rz(-0.29414487) q[3];
sx q[3];
rz(-1.1295986) q[3];
sx q[3];
rz(2.8637776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1421563) q[0];
sx q[0];
rz(-0.23040982) q[0];
sx q[0];
rz(-0.18705046) q[0];
rz(2.710178) q[1];
sx q[1];
rz(-0.43771935) q[1];
sx q[1];
rz(1.4923219) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2628415) q[0];
sx q[0];
rz(-3.0410112) q[0];
sx q[0];
rz(-0.79128964) q[0];
x q[1];
rz(-1.0643105) q[2];
sx q[2];
rz(-2.7632755) q[2];
sx q[2];
rz(-1.0418721) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8111251) q[1];
sx q[1];
rz(-2.6973675) q[1];
sx q[1];
rz(-1.1796239) q[1];
rz(-1.727333) q[3];
sx q[3];
rz(-1.6633342) q[3];
sx q[3];
rz(-2.375556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.46818587) q[2];
sx q[2];
rz(-2.2298721) q[2];
sx q[2];
rz(-0.98463303) q[2];
rz(-0.51472384) q[3];
sx q[3];
rz(-0.52557164) q[3];
sx q[3];
rz(1.908186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.89727) q[0];
sx q[0];
rz(-1.4800625) q[0];
sx q[0];
rz(2.3829714) q[0];
rz(1.1933391) q[1];
sx q[1];
rz(-1.1921644) q[1];
sx q[1];
rz(-1.4512482) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0139931) q[0];
sx q[0];
rz(-1.7993449) q[0];
sx q[0];
rz(0.20968135) q[0];
rz(-pi) q[1];
rz(-1.7733712) q[2];
sx q[2];
rz(-1.612101) q[2];
sx q[2];
rz(-1.5243204) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.38498653) q[1];
sx q[1];
rz(-2.1861417) q[1];
sx q[1];
rz(-0.032757515) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8192807) q[3];
sx q[3];
rz(-1.1813643) q[3];
sx q[3];
rz(-3.1403613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5997368) q[2];
sx q[2];
rz(-0.92076045) q[2];
sx q[2];
rz(-0.80121458) q[2];
rz(2.2090705) q[3];
sx q[3];
rz(-1.9263809) q[3];
sx q[3];
rz(2.7264989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.5634609) q[0];
sx q[0];
rz(-1.7628071) q[0];
sx q[0];
rz(2.5264869) q[0];
rz(0.32147944) q[1];
sx q[1];
rz(-2.165806) q[1];
sx q[1];
rz(1.4478366) q[1];
rz(0.64340016) q[2];
sx q[2];
rz(-0.20558029) q[2];
sx q[2];
rz(-2.554648) q[2];
rz(-1.8516171) q[3];
sx q[3];
rz(-1.349189) q[3];
sx q[3];
rz(1.0536581) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
