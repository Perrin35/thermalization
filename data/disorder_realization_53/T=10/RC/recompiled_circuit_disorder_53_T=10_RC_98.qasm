OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.18908137) q[0];
sx q[0];
rz(-0.12726769) q[0];
sx q[0];
rz(-2.1531818) q[0];
rz(2.610511) q[1];
sx q[1];
rz(-0.34871066) q[1];
sx q[1];
rz(-1.7928064) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3492337) q[0];
sx q[0];
rz(-1.6867189) q[0];
sx q[0];
rz(2.9754292) q[0];
x q[1];
rz(2.8785107) q[2];
sx q[2];
rz(-1.6533268) q[2];
sx q[2];
rz(-2.6334327) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8197644) q[1];
sx q[1];
rz(-2.1121703) q[1];
sx q[1];
rz(-2.2775047) q[1];
x q[2];
rz(0.19282135) q[3];
sx q[3];
rz(-1.3746327) q[3];
sx q[3];
rz(1.7455268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7621883) q[2];
sx q[2];
rz(-0.72221243) q[2];
sx q[2];
rz(2.7963426) q[2];
rz(2.9521862) q[3];
sx q[3];
rz(-2.6510986) q[3];
sx q[3];
rz(2.1729443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17125601) q[0];
sx q[0];
rz(-2.4617564) q[0];
sx q[0];
rz(-2.7217857) q[0];
rz(-0.35821113) q[1];
sx q[1];
rz(-2.5102291) q[1];
sx q[1];
rz(-2.125724) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43243877) q[0];
sx q[0];
rz(-2.0560987) q[0];
sx q[0];
rz(1.1164467) q[0];
rz(-pi) q[1];
x q[1];
rz(2.303896) q[2];
sx q[2];
rz(-2.838755) q[2];
sx q[2];
rz(2.8267415) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0527555) q[1];
sx q[1];
rz(-1.4375327) q[1];
sx q[1];
rz(2.7620176) q[1];
rz(2.696051) q[3];
sx q[3];
rz(-2.0317151) q[3];
sx q[3];
rz(1.2156435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.35280716) q[2];
sx q[2];
rz(-2.8217227) q[2];
sx q[2];
rz(1.1091728) q[2];
rz(-2.7099113) q[3];
sx q[3];
rz(-2.7717398) q[3];
sx q[3];
rz(3.0811908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1576841) q[0];
sx q[0];
rz(-1.2976054) q[0];
sx q[0];
rz(0.88009673) q[0];
rz(-1.2616715) q[1];
sx q[1];
rz(-0.59363669) q[1];
sx q[1];
rz(0.18149158) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49849579) q[0];
sx q[0];
rz(-0.65319121) q[0];
sx q[0];
rz(0.80723395) q[0];
x q[1];
rz(-0.65012765) q[2];
sx q[2];
rz(-1.0368772) q[2];
sx q[2];
rz(0.18661737) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.1101802) q[1];
sx q[1];
rz(-1.2978147) q[1];
sx q[1];
rz(-3.1168373) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7865903) q[3];
sx q[3];
rz(-2.0424235) q[3];
sx q[3];
rz(2.140688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.024753831) q[2];
sx q[2];
rz(-2.0334058) q[2];
sx q[2];
rz(-1.3428358) q[2];
rz(-1.3252307) q[3];
sx q[3];
rz(-0.67458761) q[3];
sx q[3];
rz(1.0935812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4745859) q[0];
sx q[0];
rz(-2.6285567) q[0];
sx q[0];
rz(2.4530607) q[0];
rz(0.21903285) q[1];
sx q[1];
rz(-2.7929247) q[1];
sx q[1];
rz(-0.15904388) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6200136) q[0];
sx q[0];
rz(-2.0984681) q[0];
sx q[0];
rz(2.2844727) q[0];
rz(2.364047) q[2];
sx q[2];
rz(-0.50707196) q[2];
sx q[2];
rz(0.61447243) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.98475458) q[1];
sx q[1];
rz(-1.595669) q[1];
sx q[1];
rz(-1.3590042) q[1];
rz(-pi) q[2];
rz(0.36985107) q[3];
sx q[3];
rz(-2.5858013) q[3];
sx q[3];
rz(-0.37214798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.29545414) q[2];
sx q[2];
rz(-1.4825772) q[2];
sx q[2];
rz(0.42567483) q[2];
rz(-0.019429026) q[3];
sx q[3];
rz(-2.9255376) q[3];
sx q[3];
rz(-0.69452906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.6762125) q[0];
sx q[0];
rz(-3.1349482) q[0];
sx q[0];
rz(-3.0131969) q[0];
rz(2.45576) q[1];
sx q[1];
rz(-2.1452955) q[1];
sx q[1];
rz(-0.54840666) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054664139) q[0];
sx q[0];
rz(-1.6097277) q[0];
sx q[0];
rz(2.1462847) q[0];
rz(-pi) q[1];
rz(1.7565281) q[2];
sx q[2];
rz(-0.81132946) q[2];
sx q[2];
rz(-0.059543691) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.027443176) q[1];
sx q[1];
rz(-0.47361923) q[1];
sx q[1];
rz(-1.4373535) q[1];
rz(1.8091651) q[3];
sx q[3];
rz(-2.7276464) q[3];
sx q[3];
rz(-0.64869374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2942865) q[2];
sx q[2];
rz(-2.8319478) q[2];
sx q[2];
rz(1.6748641) q[2];
rz(-1.0855801) q[3];
sx q[3];
rz(-1.3054566) q[3];
sx q[3];
rz(-2.2646358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10228957) q[0];
sx q[0];
rz(-1.1941432) q[0];
sx q[0];
rz(-2.0423245) q[0];
rz(0.33637235) q[1];
sx q[1];
rz(-2.4772494) q[1];
sx q[1];
rz(2.1469595) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54906819) q[0];
sx q[0];
rz(-1.7905092) q[0];
sx q[0];
rz(-0.99411221) q[0];
rz(-pi) q[1];
rz(1.5744393) q[2];
sx q[2];
rz(-2.074303) q[2];
sx q[2];
rz(-0.41479455) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1514046) q[1];
sx q[1];
rz(-1.908704) q[1];
sx q[1];
rz(0.0068454725) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46835132) q[3];
sx q[3];
rz(-0.75707179) q[3];
sx q[3];
rz(-0.10145951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2991335) q[2];
sx q[2];
rz(-0.9023388) q[2];
sx q[2];
rz(-2.6966406) q[2];
rz(-0.80727243) q[3];
sx q[3];
rz(-3.0326796) q[3];
sx q[3];
rz(2.3256433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31535661) q[0];
sx q[0];
rz(-0.5744136) q[0];
sx q[0];
rz(0.89609599) q[0];
rz(2.232146) q[1];
sx q[1];
rz(-0.25032955) q[1];
sx q[1];
rz(3.0665841) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79338193) q[0];
sx q[0];
rz(-2.777034) q[0];
sx q[0];
rz(-1.3937464) q[0];
rz(2.1298112) q[2];
sx q[2];
rz(-2.2885239) q[2];
sx q[2];
rz(-0.43874028) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.311104) q[1];
sx q[1];
rz(-1.6357592) q[1];
sx q[1];
rz(-0.91598367) q[1];
rz(-pi) q[2];
rz(-2.1771031) q[3];
sx q[3];
rz(-0.74092591) q[3];
sx q[3];
rz(2.9348506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.16671495) q[2];
sx q[2];
rz(-1.1324984) q[2];
sx q[2];
rz(-1.9141076) q[2];
rz(3.098439) q[3];
sx q[3];
rz(-1.4814601) q[3];
sx q[3];
rz(-0.22668049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1535004) q[0];
sx q[0];
rz(-0.38689125) q[0];
sx q[0];
rz(-0.41326997) q[0];
rz(2.5832672) q[1];
sx q[1];
rz(-2.3007326) q[1];
sx q[1];
rz(-2.4024898) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6571592) q[0];
sx q[0];
rz(-1.6201577) q[0];
sx q[0];
rz(0.56364787) q[0];
x q[1];
rz(1.5584281) q[2];
sx q[2];
rz(-2.6855199) q[2];
sx q[2];
rz(2.9197444) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6369789) q[1];
sx q[1];
rz(-0.74090545) q[1];
sx q[1];
rz(-1.8085338) q[1];
rz(-pi) q[2];
rz(0.056034485) q[3];
sx q[3];
rz(-1.3067712) q[3];
sx q[3];
rz(-0.92428401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.83453137) q[2];
sx q[2];
rz(-2.8246911) q[2];
sx q[2];
rz(2.0397662) q[2];
rz(-0.39673355) q[3];
sx q[3];
rz(-1.705403) q[3];
sx q[3];
rz(2.326899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6535764) q[0];
sx q[0];
rz(-2.5119913) q[0];
sx q[0];
rz(-2.5226412) q[0];
rz(1.0194107) q[1];
sx q[1];
rz(-1.5588201) q[1];
sx q[1];
rz(-0.5272665) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4370678) q[0];
sx q[0];
rz(-1.8509812) q[0];
sx q[0];
rz(-1.8591451) q[0];
rz(2.451532) q[2];
sx q[2];
rz(-1.0933211) q[2];
sx q[2];
rz(0.32327393) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0638709) q[1];
sx q[1];
rz(-1.1472902) q[1];
sx q[1];
rz(1.0788171) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2654952) q[3];
sx q[3];
rz(-0.68317181) q[3];
sx q[3];
rz(-0.20753577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.28425372) q[2];
sx q[2];
rz(-1.3353835) q[2];
sx q[2];
rz(2.2100892) q[2];
rz(2.3305317) q[3];
sx q[3];
rz(-0.61412007) q[3];
sx q[3];
rz(-1.567599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5321524) q[0];
sx q[0];
rz(-2.092822) q[0];
sx q[0];
rz(-0.47927454) q[0];
rz(2.2562064) q[1];
sx q[1];
rz(-0.65941864) q[1];
sx q[1];
rz(2.9634109) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6945267) q[0];
sx q[0];
rz(-0.94213001) q[0];
sx q[0];
rz(-0.64283128) q[0];
rz(0.68797942) q[2];
sx q[2];
rz(-1.2974206) q[2];
sx q[2];
rz(-2.2945987) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0900967) q[1];
sx q[1];
rz(-1.3428543) q[1];
sx q[1];
rz(1.0432613) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8799805) q[3];
sx q[3];
rz(-2.6978921) q[3];
sx q[3];
rz(-1.4902247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.46897727) q[2];
sx q[2];
rz(-2.4185833) q[2];
sx q[2];
rz(-2.318312) q[2];
rz(-3.0205884) q[3];
sx q[3];
rz(-2.3779317) q[3];
sx q[3];
rz(2.1774489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9085893) q[0];
sx q[0];
rz(-2.067726) q[0];
sx q[0];
rz(2.4698972) q[0];
rz(-1.9227149) q[1];
sx q[1];
rz(-1.3296483) q[1];
sx q[1];
rz(1.7760361) q[1];
rz(1.490996) q[2];
sx q[2];
rz(-2.7814293) q[2];
sx q[2];
rz(-0.89520892) q[2];
rz(-0.94869074) q[3];
sx q[3];
rz(-1.8507675) q[3];
sx q[3];
rz(-0.40678195) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
