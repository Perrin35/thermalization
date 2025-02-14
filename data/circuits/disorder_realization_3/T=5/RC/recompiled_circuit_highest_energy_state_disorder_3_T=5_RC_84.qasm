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
rz(-3.1138104) q[0];
sx q[0];
rz(-1.8960928) q[0];
sx q[0];
rz(-1.0863289) q[0];
rz(-1.0983941) q[1];
sx q[1];
rz(-1.2135222) q[1];
sx q[1];
rz(-1.9953802) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47482309) q[0];
sx q[0];
rz(-2.7753029) q[0];
sx q[0];
rz(-2.0369056) q[0];
rz(-2.5511247) q[2];
sx q[2];
rz(-1.9158949) q[2];
sx q[2];
rz(-2.9546628) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9900528) q[1];
sx q[1];
rz(-2.0742982) q[1];
sx q[1];
rz(0.997418) q[1];
rz(-pi) q[2];
x q[2];
rz(0.2262874) q[3];
sx q[3];
rz(-2.0620966) q[3];
sx q[3];
rz(2.0313583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9132797) q[2];
sx q[2];
rz(-1.5756807) q[2];
sx q[2];
rz(-2.7794465) q[2];
rz(-2.1892138) q[3];
sx q[3];
rz(-0.55914545) q[3];
sx q[3];
rz(0.70498103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1700965) q[0];
sx q[0];
rz(-0.2146475) q[0];
sx q[0];
rz(1.6098518) q[0];
rz(-1.8929298) q[1];
sx q[1];
rz(-1.5208533) q[1];
sx q[1];
rz(-0.85266399) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1318052) q[0];
sx q[0];
rz(-0.62776792) q[0];
sx q[0];
rz(0.96012299) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2937828) q[2];
sx q[2];
rz(-1.6887892) q[2];
sx q[2];
rz(-1.5867151) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1601847) q[1];
sx q[1];
rz(-1.369919) q[1];
sx q[1];
rz(-0.9717047) q[1];
rz(-pi) q[2];
rz(1.0107105) q[3];
sx q[3];
rz(-1.9017856) q[3];
sx q[3];
rz(2.6382642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.93069211) q[2];
sx q[2];
rz(-2.1810668) q[2];
sx q[2];
rz(-1.9453913) q[2];
rz(-3.1148552) q[3];
sx q[3];
rz(-2.6452711) q[3];
sx q[3];
rz(-1.668821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4836327) q[0];
sx q[0];
rz(-0.25068972) q[0];
sx q[0];
rz(0.87580097) q[0];
rz(0.35975131) q[1];
sx q[1];
rz(-1.3641554) q[1];
sx q[1];
rz(-2.13805) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7210081) q[0];
sx q[0];
rz(-1.6419174) q[0];
sx q[0];
rz(-0.55983243) q[0];
rz(-0.3957171) q[2];
sx q[2];
rz(-0.96683244) q[2];
sx q[2];
rz(0.25091083) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7454804) q[1];
sx q[1];
rz(-0.89876938) q[1];
sx q[1];
rz(1.7538096) q[1];
x q[2];
rz(1.0295792) q[3];
sx q[3];
rz(-2.1710058) q[3];
sx q[3];
rz(0.61934352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.935219) q[2];
sx q[2];
rz(-0.71722764) q[2];
sx q[2];
rz(-0.75500542) q[2];
rz(-3.0435009) q[3];
sx q[3];
rz(-2.5710929) q[3];
sx q[3];
rz(-2.7437955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0841242) q[0];
sx q[0];
rz(-1.6300772) q[0];
sx q[0];
rz(0.64495069) q[0];
rz(-1.0954674) q[1];
sx q[1];
rz(-0.40830937) q[1];
sx q[1];
rz(-1.1114978) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1350461) q[0];
sx q[0];
rz(-1.5113976) q[0];
sx q[0];
rz(3.0778372) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.686551) q[2];
sx q[2];
rz(-2.9559982) q[2];
sx q[2];
rz(2.0041103) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.42135887) q[1];
sx q[1];
rz(-1.780527) q[1];
sx q[1];
rz(1.7214107) q[1];
rz(-2.4006149) q[3];
sx q[3];
rz(-1.40687) q[3];
sx q[3];
rz(-1.0865097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.633454) q[2];
sx q[2];
rz(-0.78712574) q[2];
sx q[2];
rz(0.48466361) q[2];
rz(-0.44114068) q[3];
sx q[3];
rz(-1.7287247) q[3];
sx q[3];
rz(-1.887656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6345374) q[0];
sx q[0];
rz(-2.2978954) q[0];
sx q[0];
rz(-1.4265358) q[0];
rz(-0.85975319) q[1];
sx q[1];
rz(-0.25329241) q[1];
sx q[1];
rz(2.3111129) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5317212) q[0];
sx q[0];
rz(-0.97519704) q[0];
sx q[0];
rz(-2.2723212) q[0];
rz(-1.6489588) q[2];
sx q[2];
rz(-2.507308) q[2];
sx q[2];
rz(-1.8820349) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3871206) q[1];
sx q[1];
rz(-0.74229147) q[1];
sx q[1];
rz(-2.9341988) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7304585) q[3];
sx q[3];
rz(-1.214802) q[3];
sx q[3];
rz(-2.1170664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.44437235) q[2];
sx q[2];
rz(-2.1264919) q[2];
sx q[2];
rz(2.3946136) q[2];
rz(0.30484453) q[3];
sx q[3];
rz(-1.3957142) q[3];
sx q[3];
rz(-1.2386809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1759258) q[0];
sx q[0];
rz(-1.3015863) q[0];
sx q[0];
rz(-0.23342215) q[0];
rz(-2.6790791) q[1];
sx q[1];
rz(-1.9513444) q[1];
sx q[1];
rz(-0.40406427) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8916805) q[0];
sx q[0];
rz(-0.95522987) q[0];
sx q[0];
rz(2.8717442) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.42521221) q[2];
sx q[2];
rz(-1.325144) q[2];
sx q[2];
rz(-0.52580357) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.048745884) q[1];
sx q[1];
rz(-1.4629852) q[1];
sx q[1];
rz(-0.56289937) q[1];
rz(1.791782) q[3];
sx q[3];
rz(-0.95181247) q[3];
sx q[3];
rz(0.19863752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.81249753) q[2];
sx q[2];
rz(-1.7915244) q[2];
sx q[2];
rz(2.4883032) q[2];
rz(1.4981883) q[3];
sx q[3];
rz(-0.61763063) q[3];
sx q[3];
rz(2.5352855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6752328) q[0];
sx q[0];
rz(-1.9283858) q[0];
sx q[0];
rz(0.32077041) q[0];
rz(-0.65387154) q[1];
sx q[1];
rz(-2.8120698) q[1];
sx q[1];
rz(3.0810862) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084235926) q[0];
sx q[0];
rz(-2.8222146) q[0];
sx q[0];
rz(-2.1288298) q[0];
x q[1];
rz(0.94580067) q[2];
sx q[2];
rz(-0.93987465) q[2];
sx q[2];
rz(0.45880246) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4478895) q[1];
sx q[1];
rz(-2.1809019) q[1];
sx q[1];
rz(-2.542716) q[1];
rz(-pi) q[2];
rz(-0.860608) q[3];
sx q[3];
rz(-0.51653701) q[3];
sx q[3];
rz(0.21632695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3939646) q[2];
sx q[2];
rz(-1.4614033) q[2];
sx q[2];
rz(-2.4192269) q[2];
rz(-0.28313053) q[3];
sx q[3];
rz(-2.2524998) q[3];
sx q[3];
rz(-2.6319671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3515781) q[0];
sx q[0];
rz(-2.5836662) q[0];
sx q[0];
rz(-2.956692) q[0];
rz(0.50390538) q[1];
sx q[1];
rz(-1.3495219) q[1];
sx q[1];
rz(2.6769743) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4330209) q[0];
sx q[0];
rz(-1.4056217) q[0];
sx q[0];
rz(-2.9726221) q[0];
x q[1];
rz(0.41733317) q[2];
sx q[2];
rz(-1.6586896) q[2];
sx q[2];
rz(0.35466874) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.415472) q[1];
sx q[1];
rz(-1.1562043) q[1];
sx q[1];
rz(1.9860255) q[1];
x q[2];
rz(2.9161386) q[3];
sx q[3];
rz(-0.73235336) q[3];
sx q[3];
rz(-1.6868397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.23286143) q[2];
sx q[2];
rz(-1.6249012) q[2];
sx q[2];
rz(2.0880584) q[2];
rz(-2.4145224) q[3];
sx q[3];
rz(-1.2289685) q[3];
sx q[3];
rz(2.5963636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9497249) q[0];
sx q[0];
rz(-1.1911012) q[0];
sx q[0];
rz(2.9560992) q[0];
rz(0.56974757) q[1];
sx q[1];
rz(-0.92420095) q[1];
sx q[1];
rz(-1.1572908) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.630721) q[0];
sx q[0];
rz(-0.44140418) q[0];
sx q[0];
rz(0.33828232) q[0];
rz(-pi) q[1];
rz(3.0915843) q[2];
sx q[2];
rz(-2.2029624) q[2];
sx q[2];
rz(-0.061246733) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2409199) q[1];
sx q[1];
rz(-1.0176588) q[1];
sx q[1];
rz(1.8227804) q[1];
x q[2];
rz(2.7140541) q[3];
sx q[3];
rz(-0.85179177) q[3];
sx q[3];
rz(2.2710272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7374604) q[2];
sx q[2];
rz(-0.22455939) q[2];
sx q[2];
rz(1.7411211) q[2];
rz(-2.776966) q[3];
sx q[3];
rz(-2.3783763) q[3];
sx q[3];
rz(0.82272416) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6309361) q[0];
sx q[0];
rz(-0.83844227) q[0];
sx q[0];
rz(3.1070218) q[0];
rz(-2.8055387) q[1];
sx q[1];
rz(-0.91468179) q[1];
sx q[1];
rz(-0.58759442) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66389342) q[0];
sx q[0];
rz(-1.7132393) q[0];
sx q[0];
rz(-2.7201128) q[0];
x q[1];
rz(2.4108836) q[2];
sx q[2];
rz(-1.3255505) q[2];
sx q[2];
rz(2.2718549) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7149844) q[1];
sx q[1];
rz(-2.9054925) q[1];
sx q[1];
rz(1.7515504) q[1];
rz(-pi) q[2];
rz(0.54909191) q[3];
sx q[3];
rz(-0.7585511) q[3];
sx q[3];
rz(-0.18699924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.33362886) q[2];
sx q[2];
rz(-0.89666349) q[2];
sx q[2];
rz(-0.72820747) q[2];
rz(-0.41094574) q[3];
sx q[3];
rz(-0.22308895) q[3];
sx q[3];
rz(-1.8691011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4651466) q[0];
sx q[0];
rz(-1.5855753) q[0];
sx q[0];
rz(1.3517071) q[0];
rz(-2.7017055) q[1];
sx q[1];
rz(-0.37275795) q[1];
sx q[1];
rz(2.8511924) q[1];
rz(-1.8710623) q[2];
sx q[2];
rz(-1.9004442) q[2];
sx q[2];
rz(1.6729542) q[2];
rz(-3.0348745) q[3];
sx q[3];
rz(-1.617147) q[3];
sx q[3];
rz(-0.54860687) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
