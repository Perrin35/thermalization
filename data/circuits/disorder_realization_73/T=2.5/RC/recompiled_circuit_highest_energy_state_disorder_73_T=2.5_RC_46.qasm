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
rz(1.3923378) q[0];
sx q[0];
rz(-2.7633986) q[0];
sx q[0];
rz(2.7910772) q[0];
rz(-2.2502083) q[1];
sx q[1];
rz(-0.79217029) q[1];
sx q[1];
rz(-1.1729191) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55124768) q[0];
sx q[0];
rz(-2.8714193) q[0];
sx q[0];
rz(-0.7858289) q[0];
rz(-pi) q[1];
rz(0.95845285) q[2];
sx q[2];
rz(-2.8757189) q[2];
sx q[2];
rz(1.0622417) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.71095548) q[1];
sx q[1];
rz(-1.4078377) q[1];
sx q[1];
rz(-2.8204172) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.32781847) q[3];
sx q[3];
rz(-0.49428764) q[3];
sx q[3];
rz(-2.0959299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.17060908) q[2];
sx q[2];
rz(-1.2759408) q[2];
sx q[2];
rz(0.21273908) q[2];
rz(-3.0563266) q[3];
sx q[3];
rz(-1.250896) q[3];
sx q[3];
rz(-1.2712449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3474715) q[0];
sx q[0];
rz(-0.81231064) q[0];
sx q[0];
rz(-1.043327) q[0];
rz(3.0087545) q[1];
sx q[1];
rz(-0.21505198) q[1];
sx q[1];
rz(-1.9827693) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4702328) q[0];
sx q[0];
rz(-1.4263819) q[0];
sx q[0];
rz(-2.7412422) q[0];
x q[1];
rz(1.4857331) q[2];
sx q[2];
rz(-0.16189215) q[2];
sx q[2];
rz(-0.4134824) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0333799) q[1];
sx q[1];
rz(-1.3054818) q[1];
sx q[1];
rz(-2.6379774) q[1];
rz(2.2099233) q[3];
sx q[3];
rz(-0.51995819) q[3];
sx q[3];
rz(0.66321841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1966689) q[2];
sx q[2];
rz(-2.479574) q[2];
sx q[2];
rz(-0.72466737) q[2];
rz(2.479018) q[3];
sx q[3];
rz(-2.1734838) q[3];
sx q[3];
rz(2.841943) q[3];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1835943) q[0];
sx q[0];
rz(-1.8738926) q[0];
sx q[0];
rz(-2.2268353) q[0];
rz(-2.5547408) q[1];
sx q[1];
rz(-1.4205168) q[1];
sx q[1];
rz(0.14952001) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11598524) q[0];
sx q[0];
rz(-1.7113026) q[0];
sx q[0];
rz(2.0951659) q[0];
rz(-pi) q[1];
rz(2.3517866) q[2];
sx q[2];
rz(-2.6457204) q[2];
sx q[2];
rz(-0.42988955) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4853277) q[1];
sx q[1];
rz(-2.1187003) q[1];
sx q[1];
rz(0.62959558) q[1];
rz(-0.71107421) q[3];
sx q[3];
rz(-2.1319763) q[3];
sx q[3];
rz(-0.64519889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.84825039) q[2];
sx q[2];
rz(-1.9363656) q[2];
sx q[2];
rz(-2.3411574) q[2];
rz(2.6750001) q[3];
sx q[3];
rz(-1.1060017) q[3];
sx q[3];
rz(0.65565562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7993497) q[0];
sx q[0];
rz(-1.290134) q[0];
sx q[0];
rz(0.85860646) q[0];
rz(-0.41060064) q[1];
sx q[1];
rz(-2.938439) q[1];
sx q[1];
rz(2.1536749) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84774524) q[0];
sx q[0];
rz(-1.721707) q[0];
sx q[0];
rz(-1.9035089) q[0];
rz(-3.0201077) q[2];
sx q[2];
rz(-2.0628235) q[2];
sx q[2];
rz(-0.58771261) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1204099) q[1];
sx q[1];
rz(-1.2721918) q[1];
sx q[1];
rz(2.8575767) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.32295042) q[3];
sx q[3];
rz(-1.7569555) q[3];
sx q[3];
rz(-1.8017839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9504488) q[2];
sx q[2];
rz(-1.7473651) q[2];
sx q[2];
rz(0.44551715) q[2];
rz(-2.4281003) q[3];
sx q[3];
rz(-0.73234171) q[3];
sx q[3];
rz(0.92002404) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54818654) q[0];
sx q[0];
rz(-1.0819409) q[0];
sx q[0];
rz(1.5997546) q[0];
rz(-1.2708739) q[1];
sx q[1];
rz(-1.4022695) q[1];
sx q[1];
rz(-0.52938968) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.423817) q[0];
sx q[0];
rz(-2.7090906) q[0];
sx q[0];
rz(0.94571094) q[0];
x q[1];
rz(-0.28044706) q[2];
sx q[2];
rz(-2.5051162) q[2];
sx q[2];
rz(-2.150879) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3965949) q[1];
sx q[1];
rz(-2.0740119) q[1];
sx q[1];
rz(1.8580336) q[1];
x q[2];
rz(-2.9934817) q[3];
sx q[3];
rz(-1.9198196) q[3];
sx q[3];
rz(2.7452041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9993837) q[2];
sx q[2];
rz(-1.1915519) q[2];
sx q[2];
rz(-1.8184398) q[2];
rz(-2.7600944) q[3];
sx q[3];
rz(-2.0254717) q[3];
sx q[3];
rz(-0.39314666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8489654) q[0];
sx q[0];
rz(-0.75858527) q[0];
sx q[0];
rz(1.0211771) q[0];
rz(1.609833) q[1];
sx q[1];
rz(-0.94667089) q[1];
sx q[1];
rz(-0.80345947) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5860503) q[0];
sx q[0];
rz(-1.7305505) q[0];
sx q[0];
rz(2.8303888) q[0];
rz(0.41219903) q[2];
sx q[2];
rz(-1.6481109) q[2];
sx q[2];
rz(0.22082034) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.66113055) q[1];
sx q[1];
rz(-1.4509038) q[1];
sx q[1];
rz(-0.57173034) q[1];
rz(-pi) q[2];
rz(1.4705974) q[3];
sx q[3];
rz(-0.9462983) q[3];
sx q[3];
rz(1.1794832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6959186) q[2];
sx q[2];
rz(-0.58738223) q[2];
sx q[2];
rz(-2.7109801) q[2];
rz(2.5066091) q[3];
sx q[3];
rz(-0.018714232) q[3];
sx q[3];
rz(-1.8733321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1465313) q[0];
sx q[0];
rz(-2.596031) q[0];
sx q[0];
rz(1.5437641) q[0];
rz(-2.457288) q[1];
sx q[1];
rz(-1.4477718) q[1];
sx q[1];
rz(-0.79944557) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7921126) q[0];
sx q[0];
rz(-1.6474012) q[0];
sx q[0];
rz(-1.5653866) q[0];
x q[1];
rz(-0.075887738) q[2];
sx q[2];
rz(-0.85691707) q[2];
sx q[2];
rz(-0.80684987) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7974248) q[1];
sx q[1];
rz(-2.0312211) q[1];
sx q[1];
rz(-1.5820222) q[1];
x q[2];
rz(-1.0897899) q[3];
sx q[3];
rz(-2.3481927) q[3];
sx q[3];
rz(2.167986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.53576175) q[2];
sx q[2];
rz(-2.3039218) q[2];
sx q[2];
rz(-2.3568995) q[2];
rz(-3.0253518) q[3];
sx q[3];
rz(-1.616547) q[3];
sx q[3];
rz(1.2522662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5410974) q[0];
sx q[0];
rz(-1.9116115) q[0];
sx q[0];
rz(2.1121209) q[0];
rz(0.12044278) q[1];
sx q[1];
rz(-1.30013) q[1];
sx q[1];
rz(-2.5708503) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.127205) q[0];
sx q[0];
rz(-2.528101) q[0];
sx q[0];
rz(0.86742371) q[0];
rz(-pi) q[1];
rz(-0.12783862) q[2];
sx q[2];
rz(-0.72261506) q[2];
sx q[2];
rz(1.7651209) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.49739079) q[1];
sx q[1];
rz(-0.56193627) q[1];
sx q[1];
rz(-2.2119658) q[1];
rz(-pi) q[2];
rz(-1.9071155) q[3];
sx q[3];
rz(-1.3988931) q[3];
sx q[3];
rz(3.086103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.268078) q[2];
sx q[2];
rz(-2.2669078) q[2];
sx q[2];
rz(2.6178005) q[2];
rz(1.1526456) q[3];
sx q[3];
rz(-2.5609784) q[3];
sx q[3];
rz(-1.4279648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48083392) q[0];
sx q[0];
rz(-1.9402215) q[0];
sx q[0];
rz(2.7177366) q[0];
rz(-2.2747874) q[1];
sx q[1];
rz(-2.1173756) q[1];
sx q[1];
rz(-2.9383235) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.568726) q[0];
sx q[0];
rz(-2.4857268) q[0];
sx q[0];
rz(-3.0650861) q[0];
rz(-pi) q[1];
rz(-1.5166984) q[2];
sx q[2];
rz(-1.8803758) q[2];
sx q[2];
rz(0.54540173) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.79938526) q[1];
sx q[1];
rz(-0.3179271) q[1];
sx q[1];
rz(2.6683776) q[1];
rz(-2.2668082) q[3];
sx q[3];
rz(-2.751707) q[3];
sx q[3];
rz(-0.75857754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.086143494) q[2];
sx q[2];
rz(-1.92675) q[2];
sx q[2];
rz(2.8622368) q[2];
rz(1.2934359) q[3];
sx q[3];
rz(-1.3118298) q[3];
sx q[3];
rz(1.9259341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9781037) q[0];
sx q[0];
rz(-0.80253974) q[0];
sx q[0];
rz(-1.7247024) q[0];
rz(-0.92719999) q[1];
sx q[1];
rz(-1.5617153) q[1];
sx q[1];
rz(-1.4097811) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3438691) q[0];
sx q[0];
rz(-2.2316405) q[0];
sx q[0];
rz(-0.12253527) q[0];
rz(-pi) q[1];
rz(0.3566202) q[2];
sx q[2];
rz(-2.0276514) q[2];
sx q[2];
rz(-2.8518554) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7076787) q[1];
sx q[1];
rz(-2.0053889) q[1];
sx q[1];
rz(0.7264002) q[1];
x q[2];
rz(-1.8133598) q[3];
sx q[3];
rz(-1.6117192) q[3];
sx q[3];
rz(2.248482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7659144) q[2];
sx q[2];
rz(-1.8199074) q[2];
sx q[2];
rz(-2.1709757) q[2];
rz(2.4961903) q[3];
sx q[3];
rz(-0.97964764) q[3];
sx q[3];
rz(1.0055044) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29393016) q[0];
sx q[0];
rz(-1.4959338) q[0];
sx q[0];
rz(-1.3069859) q[0];
rz(2.7597799) q[1];
sx q[1];
rz(-1.3419071) q[1];
sx q[1];
rz(0.76795427) q[1];
rz(0.81224281) q[2];
sx q[2];
rz(-1.1803738) q[2];
sx q[2];
rz(1.3132172) q[2];
rz(0.44245023) q[3];
sx q[3];
rz(-1.9922602) q[3];
sx q[3];
rz(-0.69666399) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
