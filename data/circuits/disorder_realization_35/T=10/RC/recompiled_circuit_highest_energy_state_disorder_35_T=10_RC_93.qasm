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
rz(-2.45911) q[0];
sx q[0];
rz(8.7764813) q[0];
sx q[0];
rz(6.2483151) q[0];
rz(1.3999445) q[1];
sx q[1];
rz(-0.26672426) q[1];
sx q[1];
rz(-1.9884225) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49738231) q[0];
sx q[0];
rz(-0.35150042) q[0];
sx q[0];
rz(-1.6194109) q[0];
rz(0.75208374) q[2];
sx q[2];
rz(-2.8643069) q[2];
sx q[2];
rz(0.54518965) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0583078) q[1];
sx q[1];
rz(-0.95049203) q[1];
sx q[1];
rz(0.22432595) q[1];
x q[2];
rz(-1.6006013) q[3];
sx q[3];
rz(-2.6905224) q[3];
sx q[3];
rz(2.9357437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.030688588) q[2];
sx q[2];
rz(-1.2653753) q[2];
sx q[2];
rz(3.0008924) q[2];
rz(1.5594907) q[3];
sx q[3];
rz(-2.9086106) q[3];
sx q[3];
rz(-1.337602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2577308) q[0];
sx q[0];
rz(-0.95153874) q[0];
sx q[0];
rz(-1.1269493) q[0];
rz(-1.1832712) q[1];
sx q[1];
rz(-2.0014346) q[1];
sx q[1];
rz(2.8347051) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2005655) q[0];
sx q[0];
rz(-0.38492003) q[0];
sx q[0];
rz(-3.0207915) q[0];
x q[1];
rz(0.65019239) q[2];
sx q[2];
rz(-2.1474194) q[2];
sx q[2];
rz(-1.2675831) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3813016) q[1];
sx q[1];
rz(-0.25549251) q[1];
sx q[1];
rz(-0.5762655) q[1];
rz(-pi) q[2];
rz(3.0314802) q[3];
sx q[3];
rz(-1.0854774) q[3];
sx q[3];
rz(1.3934324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.41944501) q[2];
sx q[2];
rz(-2.751613) q[2];
sx q[2];
rz(-1.5004213) q[2];
rz(-1.3106583) q[3];
sx q[3];
rz(-2.1597517) q[3];
sx q[3];
rz(2.380044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0797043) q[0];
sx q[0];
rz(-2.5766928) q[0];
sx q[0];
rz(-1.0796219) q[0];
rz(-2.9234746) q[1];
sx q[1];
rz(-1.7354542) q[1];
sx q[1];
rz(-1.9535779) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69258286) q[0];
sx q[0];
rz(-2.1931358) q[0];
sx q[0];
rz(-1.5908817) q[0];
x q[1];
rz(-2.8883557) q[2];
sx q[2];
rz(-2.5498769) q[2];
sx q[2];
rz(-0.62555056) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6453533) q[1];
sx q[1];
rz(-1.5669973) q[1];
sx q[1];
rz(-3.0423052) q[1];
x q[2];
rz(-3.1165666) q[3];
sx q[3];
rz(-0.2215759) q[3];
sx q[3];
rz(3.0794249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5334566) q[2];
sx q[2];
rz(-2.6965105) q[2];
sx q[2];
rz(2.377887) q[2];
rz(-2.8869827) q[3];
sx q[3];
rz(-2.2201241) q[3];
sx q[3];
rz(-2.7931131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3449645) q[0];
sx q[0];
rz(-2.1944955) q[0];
sx q[0];
rz(3.1297041) q[0];
rz(0.97870007) q[1];
sx q[1];
rz(-1.782676) q[1];
sx q[1];
rz(-1.1241283) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55292144) q[0];
sx q[0];
rz(-1.3359937) q[0];
sx q[0];
rz(2.7166883) q[0];
rz(-pi) q[1];
x q[1];
rz(0.57968906) q[2];
sx q[2];
rz(-1.1986898) q[2];
sx q[2];
rz(-1.6339782) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2024514) q[1];
sx q[1];
rz(-1.141046) q[1];
sx q[1];
rz(0.38265574) q[1];
x q[2];
rz(2.5683764) q[3];
sx q[3];
rz(-2.3251136) q[3];
sx q[3];
rz(3.1117718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0594788) q[2];
sx q[2];
rz(-0.79136005) q[2];
sx q[2];
rz(0.7938844) q[2];
rz(-1.7700899) q[3];
sx q[3];
rz(-0.81727782) q[3];
sx q[3];
rz(0.98384682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7612979) q[0];
sx q[0];
rz(-1.5341606) q[0];
sx q[0];
rz(-2.5979331) q[0];
rz(-1.1477973) q[1];
sx q[1];
rz(-0.59998435) q[1];
sx q[1];
rz(1.071484) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6010671) q[0];
sx q[0];
rz(-2.6357542) q[0];
sx q[0];
rz(0.95212014) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4161827) q[2];
sx q[2];
rz(-1.8616759) q[2];
sx q[2];
rz(2.9335528) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.33794854) q[1];
sx q[1];
rz(-2.0211086) q[1];
sx q[1];
rz(1.922035) q[1];
rz(-pi) q[2];
rz(-2.0918455) q[3];
sx q[3];
rz(-2.30671) q[3];
sx q[3];
rz(-0.3835333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.063283) q[2];
sx q[2];
rz(-2.6368243) q[2];
sx q[2];
rz(-0.90739352) q[2];
rz(3.1388969) q[3];
sx q[3];
rz(-1.4009652) q[3];
sx q[3];
rz(-1.9122972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9232848) q[0];
sx q[0];
rz(-1.2907499) q[0];
sx q[0];
rz(-1.1894591) q[0];
rz(-2.9749191) q[1];
sx q[1];
rz(-2.1107626) q[1];
sx q[1];
rz(1.6698242) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92291245) q[0];
sx q[0];
rz(-2.3793202) q[0];
sx q[0];
rz(2.5566468) q[0];
x q[1];
rz(1.5277063) q[2];
sx q[2];
rz(-0.22858563) q[2];
sx q[2];
rz(0.061554043) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.38684026) q[1];
sx q[1];
rz(-2.0943145) q[1];
sx q[1];
rz(-2.0493815) q[1];
x q[2];
rz(-0.2324403) q[3];
sx q[3];
rz(-2.3101984) q[3];
sx q[3];
rz(1.1375994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2344096) q[2];
sx q[2];
rz(-3.0295591) q[2];
sx q[2];
rz(2.8049133) q[2];
rz(-1.178721) q[3];
sx q[3];
rz(-1.7095292) q[3];
sx q[3];
rz(0.36693507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1692899) q[0];
sx q[0];
rz(-0.65880913) q[0];
sx q[0];
rz(1.7377874) q[0];
rz(-2.9804969) q[1];
sx q[1];
rz(-2.1889071) q[1];
sx q[1];
rz(-2.5162627) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.071043) q[0];
sx q[0];
rz(-0.82938507) q[0];
sx q[0];
rz(0.87096386) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1834945) q[2];
sx q[2];
rz(-2.4110458) q[2];
sx q[2];
rz(-1.1047266) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4703964) q[1];
sx q[1];
rz(-2.0866864) q[1];
sx q[1];
rz(0.080027894) q[1];
x q[2];
rz(-1.9670385) q[3];
sx q[3];
rz(-1.4449287) q[3];
sx q[3];
rz(-2.270442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.27955851) q[2];
sx q[2];
rz(-1.6509849) q[2];
sx q[2];
rz(0.030869182) q[2];
rz(-0.75012642) q[3];
sx q[3];
rz(-2.153219) q[3];
sx q[3];
rz(1.3039024) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52791643) q[0];
sx q[0];
rz(-0.91687098) q[0];
sx q[0];
rz(-0.11012878) q[0];
rz(2.2573709) q[1];
sx q[1];
rz(-1.607211) q[1];
sx q[1];
rz(-1.0482739) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1048909) q[0];
sx q[0];
rz(-2.740777) q[0];
sx q[0];
rz(2.4359279) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3051304) q[2];
sx q[2];
rz(-1.507826) q[2];
sx q[2];
rz(-2.8612325) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6148664) q[1];
sx q[1];
rz(-1.6294565) q[1];
sx q[1];
rz(-0.23414302) q[1];
rz(1.8258207) q[3];
sx q[3];
rz(-2.9035845) q[3];
sx q[3];
rz(1.2824392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5683762) q[2];
sx q[2];
rz(-1.5244851) q[2];
sx q[2];
rz(-3.1061843) q[2];
rz(-3.0951485) q[3];
sx q[3];
rz(-1.2848264) q[3];
sx q[3];
rz(0.057028381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-0.34920084) q[0];
sx q[0];
rz(-2.2430895) q[0];
sx q[0];
rz(-3.0080556) q[0];
rz(-1.7298493) q[1];
sx q[1];
rz(-1.1213419) q[1];
sx q[1];
rz(2.0917361) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8658757) q[0];
sx q[0];
rz(-1.9483074) q[0];
sx q[0];
rz(-1.6927682) q[0];
rz(-3.1341928) q[2];
sx q[2];
rz(-1.0303633) q[2];
sx q[2];
rz(-2.2280047) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0120318) q[1];
sx q[1];
rz(-2.311918) q[1];
sx q[1];
rz(3.0864703) q[1];
rz(-pi) q[2];
rz(1.2796822) q[3];
sx q[3];
rz(-1.9587224) q[3];
sx q[3];
rz(-2.6879626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2213846) q[2];
sx q[2];
rz(-2.460026) q[2];
sx q[2];
rz(-1.3725012) q[2];
rz(-1.3051322) q[3];
sx q[3];
rz(-1.2499481) q[3];
sx q[3];
rz(0.60177747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16804279) q[0];
sx q[0];
rz(-1.16007) q[0];
sx q[0];
rz(-2.2450182) q[0];
rz(-0.98094455) q[1];
sx q[1];
rz(-2.4368821) q[1];
sx q[1];
rz(-0.14095813) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8887863) q[0];
sx q[0];
rz(-1.3355635) q[0];
sx q[0];
rz(1.0381519) q[0];
rz(-1.5522573) q[2];
sx q[2];
rz(-2.435212) q[2];
sx q[2];
rz(-3.0289087) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.861921) q[1];
sx q[1];
rz(-1.035861) q[1];
sx q[1];
rz(-1.5831854) q[1];
rz(-pi) q[2];
rz(2.6550547) q[3];
sx q[3];
rz(-1.3462023) q[3];
sx q[3];
rz(2.3096245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7732546) q[2];
sx q[2];
rz(-0.95803037) q[2];
sx q[2];
rz(2.9216596) q[2];
rz(-2.1851165) q[3];
sx q[3];
rz(-2.317231) q[3];
sx q[3];
rz(-2.1970356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3750951) q[0];
sx q[0];
rz(-1.5909593) q[0];
sx q[0];
rz(2.5084556) q[0];
rz(1.7924869) q[1];
sx q[1];
rz(-0.92020412) q[1];
sx q[1];
rz(-0.79457582) q[1];
rz(-0.53403833) q[2];
sx q[2];
rz(-1.8587458) q[2];
sx q[2];
rz(0.15482422) q[2];
rz(-1.7344162) q[3];
sx q[3];
rz(-0.72973482) q[3];
sx q[3];
rz(-0.90009298) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
