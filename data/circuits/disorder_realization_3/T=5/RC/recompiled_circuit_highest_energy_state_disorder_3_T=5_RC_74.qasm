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
rz(0.02778223) q[0];
sx q[0];
rz(-1.2454998) q[0];
sx q[0];
rz(-2.0552638) q[0];
rz(-1.0983941) q[1];
sx q[1];
rz(-1.2135222) q[1];
sx q[1];
rz(-1.9953802) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6667696) q[0];
sx q[0];
rz(-0.36628977) q[0];
sx q[0];
rz(1.1046871) q[0];
rz(-pi) q[1];
x q[1];
rz(1.97922) q[2];
sx q[2];
rz(-2.1222489) q[2];
sx q[2];
rz(1.1609032) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.15153989) q[1];
sx q[1];
rz(-1.0672945) q[1];
sx q[1];
rz(-0.997418) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9678461) q[3];
sx q[3];
rz(-2.6045817) q[3];
sx q[3];
rz(2.4853693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9132797) q[2];
sx q[2];
rz(-1.5756807) q[2];
sx q[2];
rz(-0.3621462) q[2];
rz(0.95237887) q[3];
sx q[3];
rz(-2.5824472) q[3];
sx q[3];
rz(-0.70498103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1700965) q[0];
sx q[0];
rz(-2.9269452) q[0];
sx q[0];
rz(-1.5317408) q[0];
rz(1.2486628) q[1];
sx q[1];
rz(-1.6207393) q[1];
sx q[1];
rz(-2.2889287) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.296761) q[0];
sx q[0];
rz(-2.0728025) q[0];
sx q[0];
rz(0.39433483) q[0];
rz(-pi) q[1];
x q[1];
rz(3.018969) q[2];
sx q[2];
rz(-1.8458335) q[2];
sx q[2];
rz(-3.1240535) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.87369253) q[1];
sx q[1];
rz(-2.5136569) q[1];
sx q[1];
rz(-1.2242641) q[1];
rz(-pi) q[2];
rz(2.1308822) q[3];
sx q[3];
rz(-1.9017856) q[3];
sx q[3];
rz(-2.6382642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.93069211) q[2];
sx q[2];
rz(-0.96052581) q[2];
sx q[2];
rz(-1.9453913) q[2];
rz(-0.026737468) q[3];
sx q[3];
rz(-0.49632159) q[3];
sx q[3];
rz(-1.668821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4836327) q[0];
sx q[0];
rz(-2.8909029) q[0];
sx q[0];
rz(2.2657917) q[0];
rz(2.7818413) q[1];
sx q[1];
rz(-1.3641554) q[1];
sx q[1];
rz(2.13805) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7210081) q[0];
sx q[0];
rz(-1.4996753) q[0];
sx q[0];
rz(0.55983243) q[0];
rz(-pi) q[1];
x q[1];
rz(0.92873145) q[2];
sx q[2];
rz(-1.8936529) q[2];
sx q[2];
rz(-1.0869458) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.39611227) q[1];
sx q[1];
rz(-0.89876938) q[1];
sx q[1];
rz(-1.387783) q[1];
rz(-1.0295792) q[3];
sx q[3];
rz(-2.1710058) q[3];
sx q[3];
rz(2.5222491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.935219) q[2];
sx q[2];
rz(-2.424365) q[2];
sx q[2];
rz(0.75500542) q[2];
rz(-0.098091789) q[3];
sx q[3];
rz(-0.57049975) q[3];
sx q[3];
rz(-2.7437955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0574684) q[0];
sx q[0];
rz(-1.5115154) q[0];
sx q[0];
rz(-2.496642) q[0];
rz(1.0954674) q[1];
sx q[1];
rz(-2.7332833) q[1];
sx q[1];
rz(-1.1114978) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7020525) q[0];
sx q[0];
rz(-1.6344392) q[0];
sx q[0];
rz(1.511277) q[0];
x q[1];
rz(2.9745151) q[2];
sx q[2];
rz(-1.4896059) q[2];
sx q[2];
rz(-2.2600391) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9605691) q[1];
sx q[1];
rz(-1.4235068) q[1];
sx q[1];
rz(0.21206124) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74097775) q[3];
sx q[3];
rz(-1.7347226) q[3];
sx q[3];
rz(-1.0865097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.633454) q[2];
sx q[2];
rz(-2.3544669) q[2];
sx q[2];
rz(-2.656929) q[2];
rz(-0.44114068) q[3];
sx q[3];
rz(-1.4128679) q[3];
sx q[3];
rz(-1.2539366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50705528) q[0];
sx q[0];
rz(-2.2978954) q[0];
sx q[0];
rz(-1.7150568) q[0];
rz(-0.85975319) q[1];
sx q[1];
rz(-2.8883002) q[1];
sx q[1];
rz(0.8304798) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62522307) q[0];
sx q[0];
rz(-0.88623669) q[0];
sx q[0];
rz(-2.3806118) q[0];
rz(-pi) q[1];
rz(-2.203622) q[2];
sx q[2];
rz(-1.6170849) q[2];
sx q[2];
rz(2.7673495) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.66248686) q[1];
sx q[1];
rz(-1.7104407) q[1];
sx q[1];
rz(-2.4100811) q[1];
x q[2];
rz(1.7304585) q[3];
sx q[3];
rz(-1.9267907) q[3];
sx q[3];
rz(-1.0245263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6972203) q[2];
sx q[2];
rz(-1.0151007) q[2];
sx q[2];
rz(2.3946136) q[2];
rz(2.8367481) q[3];
sx q[3];
rz(-1.3957142) q[3];
sx q[3];
rz(1.2386809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1759258) q[0];
sx q[0];
rz(-1.8400064) q[0];
sx q[0];
rz(0.23342215) q[0];
rz(0.4625136) q[1];
sx q[1];
rz(-1.9513444) q[1];
sx q[1];
rz(-0.40406427) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2499122) q[0];
sx q[0];
rz(-0.95522987) q[0];
sx q[0];
rz(-2.8717442) q[0];
rz(-pi) q[1];
x q[1];
rz(1.839371) q[2];
sx q[2];
rz(-1.982455) q[2];
sx q[2];
rz(-1.1546763) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0928468) q[1];
sx q[1];
rz(-1.4629852) q[1];
sx q[1];
rz(-2.5786933) q[1];
rz(-pi) q[2];
rz(2.8430953) q[3];
sx q[3];
rz(-0.65234557) q[3];
sx q[3];
rz(2.9708095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.81249753) q[2];
sx q[2];
rz(-1.7915244) q[2];
sx q[2];
rz(-2.4883032) q[2];
rz(-1.6434044) q[3];
sx q[3];
rz(-0.61763063) q[3];
sx q[3];
rz(2.5352855) q[3];
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
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46635982) q[0];
sx q[0];
rz(-1.2132069) q[0];
sx q[0];
rz(2.8208222) q[0];
rz(-2.4877211) q[1];
sx q[1];
rz(-2.8120698) q[1];
sx q[1];
rz(-3.0810862) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9515647) q[0];
sx q[0];
rz(-1.7378283) q[0];
sx q[0];
rz(-1.8442979) q[0];
rz(-pi) q[1];
rz(0.67530151) q[2];
sx q[2];
rz(-2.2847698) q[2];
sx q[2];
rz(1.7973822) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4478895) q[1];
sx q[1];
rz(-2.1809019) q[1];
sx q[1];
rz(-0.59887664) q[1];
x q[2];
rz(-1.9774505) q[3];
sx q[3];
rz(-1.2429626) q[3];
sx q[3];
rz(1.9964807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.74762809) q[2];
sx q[2];
rz(-1.4614033) q[2];
sx q[2];
rz(2.4192269) q[2];
rz(0.28313053) q[3];
sx q[3];
rz(-2.2524998) q[3];
sx q[3];
rz(2.6319671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3515781) q[0];
sx q[0];
rz(-0.55792648) q[0];
sx q[0];
rz(0.18490069) q[0];
rz(2.6376873) q[1];
sx q[1];
rz(-1.7920707) q[1];
sx q[1];
rz(-0.46461836) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5123925) q[0];
sx q[0];
rz(-0.2357395) q[0];
sx q[0];
rz(-2.3605973) q[0];
rz(-pi) q[1];
rz(2.9275168) q[2];
sx q[2];
rz(-0.42595562) q[2];
sx q[2];
rz(1.730012) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.415472) q[1];
sx q[1];
rz(-1.1562043) q[1];
sx q[1];
rz(1.1555671) q[1];
x q[2];
rz(2.4219651) q[3];
sx q[3];
rz(-1.7208281) q[3];
sx q[3];
rz(3.0887136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.23286143) q[2];
sx q[2];
rz(-1.6249012) q[2];
sx q[2];
rz(-2.0880584) q[2];
rz(2.4145224) q[3];
sx q[3];
rz(-1.9126242) q[3];
sx q[3];
rz(-0.54522902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9497249) q[0];
sx q[0];
rz(-1.1911012) q[0];
sx q[0];
rz(-0.18549347) q[0];
rz(-0.56974757) q[1];
sx q[1];
rz(-0.92420095) q[1];
sx q[1];
rz(1.1572908) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13979736) q[0];
sx q[0];
rz(-1.1560062) q[0];
sx q[0];
rz(-1.4152566) q[0];
rz(0.93803381) q[2];
sx q[2];
rz(-1.530458) q[2];
sx q[2];
rz(1.4799839) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.19543513) q[1];
sx q[1];
rz(-1.7845673) q[1];
sx q[1];
rz(2.5740088) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1284105) q[3];
sx q[3];
rz(-2.3250322) q[3];
sx q[3];
rz(0.26536322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.40413228) q[2];
sx q[2];
rz(-0.22455939) q[2];
sx q[2];
rz(-1.7411211) q[2];
rz(-0.36462668) q[3];
sx q[3];
rz(-2.3783763) q[3];
sx q[3];
rz(2.3188685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51065651) q[0];
sx q[0];
rz(-0.83844227) q[0];
sx q[0];
rz(-3.1070218) q[0];
rz(0.336054) q[1];
sx q[1];
rz(-0.91468179) q[1];
sx q[1];
rz(-0.58759442) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60025763) q[0];
sx q[0];
rz(-0.44352874) q[0];
sx q[0];
rz(-2.8044273) q[0];
rz(-pi) q[1];
rz(-1.8950225) q[2];
sx q[2];
rz(-2.2750008) q[2];
sx q[2];
rz(-2.2262822) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1732483) q[1];
sx q[1];
rz(-1.6128596) q[1];
sx q[1];
rz(1.3384046) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.67983277) q[3];
sx q[3];
rz(-1.9380016) q[3];
sx q[3];
rz(-2.1757371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.33362886) q[2];
sx q[2];
rz(-0.89666349) q[2];
sx q[2];
rz(0.72820747) q[2];
rz(2.7306469) q[3];
sx q[3];
rz(-2.9185037) q[3];
sx q[3];
rz(-1.2724916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4651466) q[0];
sx q[0];
rz(-1.5855753) q[0];
sx q[0];
rz(1.3517071) q[0];
rz(0.43988718) q[1];
sx q[1];
rz(-0.37275795) q[1];
sx q[1];
rz(2.8511924) q[1];
rz(2.4287379) q[2];
sx q[2];
rz(-0.44217449) q[2];
sx q[2];
rz(0.90978734) q[2];
rz(-0.41070134) q[3];
sx q[3];
rz(-0.11631415) q[3];
sx q[3];
rz(0.61396413) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
