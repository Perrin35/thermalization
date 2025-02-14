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
rz(0.12953144) q[0];
sx q[0];
rz(-0.13106267) q[0];
sx q[0];
rz(-0.60453209) q[0];
rz(-0.52855748) q[1];
sx q[1];
rz(2.8837535) q[1];
sx q[1];
rz(16.937994) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.226868) q[0];
sx q[0];
rz(-2.719923) q[0];
sx q[0];
rz(0.64897169) q[0];
x q[1];
rz(-2.5276353) q[2];
sx q[2];
rz(-1.3836622) q[2];
sx q[2];
rz(-0.89636114) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9166419) q[1];
sx q[1];
rz(-1.0826546) q[1];
sx q[1];
rz(0.78190885) q[1];
rz(-1.3377519) q[3];
sx q[3];
rz(-0.50737689) q[3];
sx q[3];
rz(-2.9180067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0872385) q[2];
sx q[2];
rz(-1.1776935) q[2];
sx q[2];
rz(0.13620201) q[2];
rz(0.44998351) q[3];
sx q[3];
rz(-2.4583702) q[3];
sx q[3];
rz(-2.3581678) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1012652) q[0];
sx q[0];
rz(-0.7985006) q[0];
sx q[0];
rz(0.26327565) q[0];
rz(-2.1030262) q[1];
sx q[1];
rz(-2.7949605) q[1];
sx q[1];
rz(2.3668049) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8717125) q[0];
sx q[0];
rz(-2.1795666) q[0];
sx q[0];
rz(1.3193498) q[0];
x q[1];
rz(-0.69157522) q[2];
sx q[2];
rz(-1.0333158) q[2];
sx q[2];
rz(-0.36851685) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9866219) q[1];
sx q[1];
rz(-0.38119527) q[1];
sx q[1];
rz(-1.4710557) q[1];
rz(-pi) q[2];
rz(-0.70330422) q[3];
sx q[3];
rz(-1.1835795) q[3];
sx q[3];
rz(2.4037698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.420555) q[2];
sx q[2];
rz(-1.5169531) q[2];
sx q[2];
rz(-0.59845412) q[2];
rz(1.362644) q[3];
sx q[3];
rz(-0.88038954) q[3];
sx q[3];
rz(-2.8708598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28629985) q[0];
sx q[0];
rz(-0.4628276) q[0];
sx q[0];
rz(-1.6337974) q[0];
rz(-0.15448054) q[1];
sx q[1];
rz(-1.721761) q[1];
sx q[1];
rz(-2.3522164) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3794147) q[0];
sx q[0];
rz(-1.6573849) q[0];
sx q[0];
rz(-2.6103254) q[0];
rz(-pi) q[1];
rz(0.28114281) q[2];
sx q[2];
rz(-1.4340377) q[2];
sx q[2];
rz(3.1355372) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.22948466) q[1];
sx q[1];
rz(-0.69522017) q[1];
sx q[1];
rz(2.0575614) q[1];
x q[2];
rz(-1.0663197) q[3];
sx q[3];
rz(-0.83791997) q[3];
sx q[3];
rz(1.5056226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.35884759) q[2];
sx q[2];
rz(-0.91155702) q[2];
sx q[2];
rz(1.4770799) q[2];
rz(-1.9278256) q[3];
sx q[3];
rz(-1.8002847) q[3];
sx q[3];
rz(1.414813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7724991) q[0];
sx q[0];
rz(-1.0558244) q[0];
sx q[0];
rz(-0.61071998) q[0];
rz(0.67131132) q[1];
sx q[1];
rz(-1.6339615) q[1];
sx q[1];
rz(1.3549365) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5089534) q[0];
sx q[0];
rz(-1.1221948) q[0];
sx q[0];
rz(-0.53855702) q[0];
rz(-pi) q[1];
rz(-1.942039) q[2];
sx q[2];
rz(-1.4535731) q[2];
sx q[2];
rz(-0.10088149) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0860885) q[1];
sx q[1];
rz(-1.5079867) q[1];
sx q[1];
rz(-1.8163866) q[1];
x q[2];
rz(1.4282966) q[3];
sx q[3];
rz(-0.5915407) q[3];
sx q[3];
rz(1.9672036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9352202) q[2];
sx q[2];
rz(-0.78038961) q[2];
sx q[2];
rz(-2.8727403) q[2];
rz(-2.3937461) q[3];
sx q[3];
rz(-2.3220389) q[3];
sx q[3];
rz(1.6929172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95626962) q[0];
sx q[0];
rz(-1.3893501) q[0];
sx q[0];
rz(-0.67778936) q[0];
rz(-2.9970844) q[1];
sx q[1];
rz(-2.3542207) q[1];
sx q[1];
rz(1.4871303) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22600284) q[0];
sx q[0];
rz(-1.5969513) q[0];
sx q[0];
rz(1.7067065) q[0];
rz(-pi) q[1];
rz(-1.4356639) q[2];
sx q[2];
rz(-2.7362974) q[2];
sx q[2];
rz(-2.7050893) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.9095535) q[1];
sx q[1];
rz(-2.0318446) q[1];
sx q[1];
rz(-2.1559245) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6923589) q[3];
sx q[3];
rz(-2.7633177) q[3];
sx q[3];
rz(2.2540783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3378478) q[2];
sx q[2];
rz(-2.8779112) q[2];
sx q[2];
rz(-0.35550508) q[2];
rz(-2.096094) q[3];
sx q[3];
rz(-1.239536) q[3];
sx q[3];
rz(-2.579328) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1657555) q[0];
sx q[0];
rz(-2.5582357) q[0];
sx q[0];
rz(-3.052886) q[0];
rz(-1.6070131) q[1];
sx q[1];
rz(-1.8314223) q[1];
sx q[1];
rz(-3.065899) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82268084) q[0];
sx q[0];
rz(-1.7633798) q[0];
sx q[0];
rz(0.67355021) q[0];
rz(2.3024302) q[2];
sx q[2];
rz(-1.608709) q[2];
sx q[2];
rz(-0.9818075) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0375605) q[1];
sx q[1];
rz(-1.018803) q[1];
sx q[1];
rz(-2.9747573) q[1];
rz(-pi) q[2];
rz(2.663732) q[3];
sx q[3];
rz(-2.5091672) q[3];
sx q[3];
rz(2.6101108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3169516) q[2];
sx q[2];
rz(-1.757916) q[2];
sx q[2];
rz(2.1023882) q[2];
rz(2.9292987) q[3];
sx q[3];
rz(-2.0391235) q[3];
sx q[3];
rz(-1.4958517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-3.0742663) q[0];
sx q[0];
rz(-2.3889611) q[0];
sx q[0];
rz(2.0972032) q[0];
rz(2.3692865) q[1];
sx q[1];
rz(-0.65985313) q[1];
sx q[1];
rz(-0.73371249) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70673101) q[0];
sx q[0];
rz(-2.0594264) q[0];
sx q[0];
rz(-0.74669331) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9915221) q[2];
sx q[2];
rz(-1.2037306) q[2];
sx q[2];
rz(2.31524) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.83625863) q[1];
sx q[1];
rz(-1.4544669) q[1];
sx q[1];
rz(2.8491324) q[1];
x q[2];
rz(-1.6902518) q[3];
sx q[3];
rz(-2.4032695) q[3];
sx q[3];
rz(-2.8286162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.61651984) q[2];
sx q[2];
rz(-2.6313582) q[2];
sx q[2];
rz(-2.5329242) q[2];
rz(-2.0742553) q[3];
sx q[3];
rz(-0.87456861) q[3];
sx q[3];
rz(-2.5909891) q[3];
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
rz(0.50255018) q[0];
sx q[0];
rz(-0.35824963) q[0];
sx q[0];
rz(1.4962366) q[0];
rz(1.7773588) q[1];
sx q[1];
rz(-0.61449209) q[1];
sx q[1];
rz(2.7552557) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8423415) q[0];
sx q[0];
rz(-2.4066397) q[0];
sx q[0];
rz(2.7432261) q[0];
rz(-pi) q[1];
rz(1.9278717) q[2];
sx q[2];
rz(-2.9437967) q[2];
sx q[2];
rz(-1.3551499) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.92703351) q[1];
sx q[1];
rz(-2.7874569) q[1];
sx q[1];
rz(2.781809) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0829857) q[3];
sx q[3];
rz(-2.5502565) q[3];
sx q[3];
rz(0.014733068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6284457) q[2];
sx q[2];
rz(-1.3725504) q[2];
sx q[2];
rz(0.15303843) q[2];
rz(-2.2506574) q[3];
sx q[3];
rz(-2.1864083) q[3];
sx q[3];
rz(0.74131596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9855758) q[0];
sx q[0];
rz(-0.22933904) q[0];
sx q[0];
rz(-0.34737059) q[0];
rz(-3.103745) q[1];
sx q[1];
rz(-1.9719351) q[1];
sx q[1];
rz(-2.6191424) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1570826) q[0];
sx q[0];
rz(-1.9831796) q[0];
sx q[0];
rz(-3.0441112) q[0];
x q[1];
rz(2.7760997) q[2];
sx q[2];
rz(-1.5710982) q[2];
sx q[2];
rz(3.0916758) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.47679893) q[1];
sx q[1];
rz(-2.0719686) q[1];
sx q[1];
rz(1.0887926) q[1];
rz(2.1736828) q[3];
sx q[3];
rz(-0.65015745) q[3];
sx q[3];
rz(-0.12891836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5106421) q[2];
sx q[2];
rz(-2.6783671) q[2];
sx q[2];
rz(1.2003215) q[2];
rz(-0.55780324) q[3];
sx q[3];
rz(-1.5188981) q[3];
sx q[3];
rz(2.7571078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2931622) q[0];
sx q[0];
rz(-0.98386216) q[0];
sx q[0];
rz(1.7278607) q[0];
rz(-3.0005455) q[1];
sx q[1];
rz(-1.3755362) q[1];
sx q[1];
rz(0.65473762) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6976234) q[0];
sx q[0];
rz(-1.3338517) q[0];
sx q[0];
rz(2.9276642) q[0];
rz(-2.0789924) q[2];
sx q[2];
rz(-2.104665) q[2];
sx q[2];
rz(-1.2956308) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4550528) q[1];
sx q[1];
rz(-1.3994263) q[1];
sx q[1];
rz(1.8167348) q[1];
x q[2];
rz(-0.61956866) q[3];
sx q[3];
rz(-1.1244785) q[3];
sx q[3];
rz(2.5321162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.35141382) q[2];
sx q[2];
rz(-1.979579) q[2];
sx q[2];
rz(2.4164825) q[2];
rz(0.890598) q[3];
sx q[3];
rz(-2.2913439) q[3];
sx q[3];
rz(-2.5543673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9643322) q[0];
sx q[0];
rz(-1.5252508) q[0];
sx q[0];
rz(-1.9510212) q[0];
rz(1.1217077) q[1];
sx q[1];
rz(-1.5737166) q[1];
sx q[1];
rz(2.166688) q[1];
rz(0.66365343) q[2];
sx q[2];
rz(-1.8603696) q[2];
sx q[2];
rz(2.6983668) q[2];
rz(-2.9573351) q[3];
sx q[3];
rz(-0.72132106) q[3];
sx q[3];
rz(-3.0527243) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
