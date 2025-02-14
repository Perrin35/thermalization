OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.81340462) q[0];
sx q[0];
rz(-0.60941154) q[0];
sx q[0];
rz(3.1031026) q[0];
rz(2.6961532) q[1];
sx q[1];
rz(-0.69094509) q[1];
sx q[1];
rz(-0.12707392) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.049517) q[0];
sx q[0];
rz(-2.319515) q[0];
sx q[0];
rz(2.6263133) q[0];
x q[1];
rz(-1.5655925) q[2];
sx q[2];
rz(-1.5706816) q[2];
sx q[2];
rz(-0.031129908) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2300079) q[1];
sx q[1];
rz(-2.717321) q[1];
sx q[1];
rz(-2.2390635) q[1];
rz(-pi) q[2];
rz(-2.2329406) q[3];
sx q[3];
rz(-2.4342038) q[3];
sx q[3];
rz(0.01568757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.27158633) q[2];
sx q[2];
rz(-0.25124696) q[2];
sx q[2];
rz(1.1146438) q[2];
rz(-2.8031269) q[3];
sx q[3];
rz(-1.4128127) q[3];
sx q[3];
rz(2.4131405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.912643) q[0];
sx q[0];
rz(-2.8680608) q[0];
sx q[0];
rz(1.1852513) q[0];
rz(1.395902) q[1];
sx q[1];
rz(-1.6604559) q[1];
sx q[1];
rz(-3.0286068) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0300354) q[0];
sx q[0];
rz(-1.5148692) q[0];
sx q[0];
rz(-1.5636754) q[0];
rz(2.0143896) q[2];
sx q[2];
rz(-1.3451066) q[2];
sx q[2];
rz(1.1301825) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.13571067) q[1];
sx q[1];
rz(-2.9067079) q[1];
sx q[1];
rz(1.2967111) q[1];
rz(-pi) q[2];
rz(-0.91199888) q[3];
sx q[3];
rz(-2.5781401) q[3];
sx q[3];
rz(-0.4215301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9129755) q[2];
sx q[2];
rz(-1.3591707) q[2];
sx q[2];
rz(-2.1302946) q[2];
rz(-0.060981123) q[3];
sx q[3];
rz(-1.1119548) q[3];
sx q[3];
rz(2.8673867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7762452) q[0];
sx q[0];
rz(-1.341935) q[0];
sx q[0];
rz(1.9976529) q[0];
rz(-1.9260319) q[1];
sx q[1];
rz(-0.85920119) q[1];
sx q[1];
rz(2.5124195) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42581707) q[0];
sx q[0];
rz(-1.5377147) q[0];
sx q[0];
rz(-2.1685409) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3843834) q[2];
sx q[2];
rz(-0.51472419) q[2];
sx q[2];
rz(1.2806614) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6110826) q[1];
sx q[1];
rz(-1.2092918) q[1];
sx q[1];
rz(1.481249) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5190303) q[3];
sx q[3];
rz(-1.3115777) q[3];
sx q[3];
rz(-1.0362877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4048142) q[2];
sx q[2];
rz(-1.3457158) q[2];
sx q[2];
rz(0.2573615) q[2];
rz(1.8836053) q[3];
sx q[3];
rz(-1.4895997) q[3];
sx q[3];
rz(-2.6902698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.314972) q[0];
sx q[0];
rz(-2.4007512) q[0];
sx q[0];
rz(1.524087) q[0];
rz(-1.357088) q[1];
sx q[1];
rz(-2.4391104) q[1];
sx q[1];
rz(1.2290899) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9953908) q[0];
sx q[0];
rz(-2.6553391) q[0];
sx q[0];
rz(2.3381837) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76009373) q[2];
sx q[2];
rz(-2.1017535) q[2];
sx q[2];
rz(-2.107055) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9214529) q[1];
sx q[1];
rz(-1.8760257) q[1];
sx q[1];
rz(-1.6294999) q[1];
x q[2];
rz(-0.48713116) q[3];
sx q[3];
rz(-1.0193977) q[3];
sx q[3];
rz(2.7745807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7288397) q[2];
sx q[2];
rz(-1.7985901) q[2];
sx q[2];
rz(-0.084224852) q[2];
rz(-1.6526875) q[3];
sx q[3];
rz(-3.0907349) q[3];
sx q[3];
rz(-0.97676718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65234891) q[0];
sx q[0];
rz(-0.69559613) q[0];
sx q[0];
rz(3.1074281) q[0];
rz(1.4802406) q[1];
sx q[1];
rz(-2.7378597) q[1];
sx q[1];
rz(1.2108948) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2171777) q[0];
sx q[0];
rz(-2.6519288) q[0];
sx q[0];
rz(-0.90474604) q[0];
x q[1];
rz(-2.6043774) q[2];
sx q[2];
rz(-0.59241931) q[2];
sx q[2];
rz(2.7835961) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.75386) q[1];
sx q[1];
rz(-1.4974623) q[1];
sx q[1];
rz(-2.0078986) q[1];
x q[2];
rz(-1.7922395) q[3];
sx q[3];
rz(-1.4998271) q[3];
sx q[3];
rz(0.25191316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5898798) q[2];
sx q[2];
rz(-2.5248542) q[2];
sx q[2];
rz(1.7390772) q[2];
rz(-1.2276522) q[3];
sx q[3];
rz(-2.472885) q[3];
sx q[3];
rz(-0.68784586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27125204) q[0];
sx q[0];
rz(-1.4372062) q[0];
sx q[0];
rz(1.2649076) q[0];
rz(-1.2875617) q[1];
sx q[1];
rz(-2.2076905) q[1];
sx q[1];
rz(2.5087779) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7408376) q[0];
sx q[0];
rz(-1.3627865) q[0];
sx q[0];
rz(1.237414) q[0];
x q[1];
rz(-0.70437141) q[2];
sx q[2];
rz(-1.272311) q[2];
sx q[2];
rz(2.5260012) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5092689) q[1];
sx q[1];
rz(-0.41328584) q[1];
sx q[1];
rz(-2.762441) q[1];
rz(-pi) q[2];
rz(2.1329844) q[3];
sx q[3];
rz(-0.62963943) q[3];
sx q[3];
rz(-1.8860157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5878933) q[2];
sx q[2];
rz(-2.9554695) q[2];
sx q[2];
rz(0.57559377) q[2];
rz(-2.0268188) q[3];
sx q[3];
rz(-1.2131178) q[3];
sx q[3];
rz(-2.7093844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-1.9734398) q[0];
sx q[0];
rz(-1.7519209) q[0];
sx q[0];
rz(0.51976505) q[0];
rz(-1.686056) q[1];
sx q[1];
rz(-0.74570233) q[1];
sx q[1];
rz(-2.0533452) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80931833) q[0];
sx q[0];
rz(-1.921321) q[0];
sx q[0];
rz(-2.4615898) q[0];
rz(0.23521346) q[2];
sx q[2];
rz(-2.1005582) q[2];
sx q[2];
rz(1.0465682) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.714915) q[1];
sx q[1];
rz(-1.928387) q[1];
sx q[1];
rz(-2.5044051) q[1];
x q[2];
rz(-1.3081) q[3];
sx q[3];
rz(-0.98219959) q[3];
sx q[3];
rz(-1.3610507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9600642) q[2];
sx q[2];
rz(-1.1390511) q[2];
sx q[2];
rz(-0.10614928) q[2];
rz(-1.6541121) q[3];
sx q[3];
rz(-0.57523504) q[3];
sx q[3];
rz(2.9851448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.079025896) q[0];
sx q[0];
rz(-1.8667969) q[0];
sx q[0];
rz(-1.6612735) q[0];
rz(1.5777499) q[1];
sx q[1];
rz(-0.5636951) q[1];
sx q[1];
rz(-0.021934358) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75087529) q[0];
sx q[0];
rz(-0.97033935) q[0];
sx q[0];
rz(-0.12156528) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0742346) q[2];
sx q[2];
rz(-1.1932917) q[2];
sx q[2];
rz(2.2855656) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.96664658) q[1];
sx q[1];
rz(-1.3783598) q[1];
sx q[1];
rz(-0.92144664) q[1];
rz(-pi) q[2];
rz(3.0472067) q[3];
sx q[3];
rz(-1.1991074) q[3];
sx q[3];
rz(-0.92640141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.03453001) q[2];
sx q[2];
rz(-1.7835534) q[2];
sx q[2];
rz(-3.0800842) q[2];
rz(-2.6574078) q[3];
sx q[3];
rz(-1.8248841) q[3];
sx q[3];
rz(-2.9137602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43559647) q[0];
sx q[0];
rz(-1.0974925) q[0];
sx q[0];
rz(-0.47501269) q[0];
rz(1.8571521) q[1];
sx q[1];
rz(-0.78452763) q[1];
sx q[1];
rz(0.49819836) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6817956) q[0];
sx q[0];
rz(-1.4180935) q[0];
sx q[0];
rz(2.5891749) q[0];
rz(-pi) q[1];
rz(0.022749697) q[2];
sx q[2];
rz(-0.91500797) q[2];
sx q[2];
rz(-1.6237669) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0507237) q[1];
sx q[1];
rz(-1.5453639) q[1];
sx q[1];
rz(0.15943105) q[1];
rz(-pi) q[2];
rz(0.60631246) q[3];
sx q[3];
rz(-2.2289235) q[3];
sx q[3];
rz(-2.6325814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4755134) q[2];
sx q[2];
rz(-0.84045118) q[2];
sx q[2];
rz(-1.6125512) q[2];
rz(0.1651925) q[3];
sx q[3];
rz(-1.8742722) q[3];
sx q[3];
rz(0.40273777) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85544473) q[0];
sx q[0];
rz(-2.5832472) q[0];
sx q[0];
rz(-0.98854351) q[0];
rz(-2.5033011) q[1];
sx q[1];
rz(-1.0811564) q[1];
sx q[1];
rz(-0.036103006) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7571018) q[0];
sx q[0];
rz(-0.57390139) q[0];
sx q[0];
rz(-1.4623653) q[0];
x q[1];
rz(-0.86697198) q[2];
sx q[2];
rz(-0.72153202) q[2];
sx q[2];
rz(-1.7448278) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.83602025) q[1];
sx q[1];
rz(-2.0479093) q[1];
sx q[1];
rz(2.1816207) q[1];
rz(-1.345031) q[3];
sx q[3];
rz(-1.3610971) q[3];
sx q[3];
rz(2.5906627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.710076) q[2];
sx q[2];
rz(-0.71147951) q[2];
sx q[2];
rz(-0.21314387) q[2];
rz(1.824481) q[3];
sx q[3];
rz(-0.20753838) q[3];
sx q[3];
rz(-0.801314) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8217736) q[0];
sx q[0];
rz(-1.736301) q[0];
sx q[0];
rz(2.2176493) q[0];
rz(0.79628235) q[1];
sx q[1];
rz(-2.2932107) q[1];
sx q[1];
rz(-0.36416818) q[1];
rz(2.4305565) q[2];
sx q[2];
rz(-1.0461251) q[2];
sx q[2];
rz(0.82991215) q[2];
rz(0.85354652) q[3];
sx q[3];
rz(-0.86409909) q[3];
sx q[3];
rz(1.7649337) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
