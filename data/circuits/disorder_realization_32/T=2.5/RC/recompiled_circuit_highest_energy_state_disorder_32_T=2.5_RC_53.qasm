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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9204694) q[0];
sx q[0];
rz(-1.2386432) q[0];
sx q[0];
rz(-1.3060547) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8240439) q[2];
sx q[2];
rz(-0.63830909) q[2];
sx q[2];
rz(2.2090863) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9166419) q[1];
sx q[1];
rz(-1.0826546) q[1];
sx q[1];
rz(2.3596838) q[1];
rz(1.0749726) q[3];
sx q[3];
rz(-1.6832441) q[3];
sx q[3];
rz(1.5898286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0872385) q[2];
sx q[2];
rz(-1.1776935) q[2];
sx q[2];
rz(-0.13620201) q[2];
rz(-2.6916091) q[3];
sx q[3];
rz(-0.68322244) q[3];
sx q[3];
rz(-0.78342485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0403274) q[0];
sx q[0];
rz(-2.3430921) q[0];
sx q[0];
rz(-0.26327565) q[0];
rz(2.1030262) q[1];
sx q[1];
rz(-2.7949605) q[1];
sx q[1];
rz(-2.3668049) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84768554) q[0];
sx q[0];
rz(-2.4890702) q[0];
sx q[0];
rz(0.34282617) q[0];
rz(-pi) q[1];
rz(0.75157849) q[2];
sx q[2];
rz(-2.2937991) q[2];
sx q[2];
rz(-2.4931049) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.50844932) q[1];
sx q[1];
rz(-1.6078498) q[1];
sx q[1];
rz(-1.9502742) q[1];
rz(2.5790096) q[3];
sx q[3];
rz(-0.78662164) q[3];
sx q[3];
rz(0.41401573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7210377) q[2];
sx q[2];
rz(-1.6246395) q[2];
sx q[2];
rz(-2.5431385) q[2];
rz(1.362644) q[3];
sx q[3];
rz(-0.88038954) q[3];
sx q[3];
rz(0.27073282) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28629985) q[0];
sx q[0];
rz(-2.6787651) q[0];
sx q[0];
rz(-1.6337974) q[0];
rz(0.15448054) q[1];
sx q[1];
rz(-1.4198317) q[1];
sx q[1];
rz(0.78937626) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1868527) q[0];
sx q[0];
rz(-0.53760872) q[0];
sx q[0];
rz(-0.16969271) q[0];
rz(-pi) q[1];
rz(-0.46044431) q[2];
sx q[2];
rz(-2.8297462) q[2];
sx q[2];
rz(2.0058035) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.22948466) q[1];
sx q[1];
rz(-2.4463725) q[1];
sx q[1];
rz(-2.0575614) q[1];
x q[2];
rz(-0.49279883) q[3];
sx q[3];
rz(-0.86233739) q[3];
sx q[3];
rz(0.81564834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7827451) q[2];
sx q[2];
rz(-0.91155702) q[2];
sx q[2];
rz(-1.6645128) q[2];
rz(-1.9278256) q[3];
sx q[3];
rz(-1.8002847) q[3];
sx q[3];
rz(1.414813) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7724991) q[0];
sx q[0];
rz(-2.0857683) q[0];
sx q[0];
rz(2.5308727) q[0];
rz(2.4702813) q[1];
sx q[1];
rz(-1.5076312) q[1];
sx q[1];
rz(1.3549365) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9498861) q[0];
sx q[0];
rz(-1.0903795) q[0];
sx q[0];
rz(2.0818162) q[0];
rz(-1.8846832) q[2];
sx q[2];
rz(-0.38849026) q[2];
sx q[2];
rz(1.1780648) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6411533) q[1];
sx q[1];
rz(-1.8158923) q[1];
sx q[1];
rz(-0.064747253) q[1];
rz(-pi) q[2];
x q[2];
rz(0.98396222) q[3];
sx q[3];
rz(-1.6500743) q[3];
sx q[3];
rz(-0.51494277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2063724) q[2];
sx q[2];
rz(-0.78038961) q[2];
sx q[2];
rz(2.8727403) q[2];
rz(-0.74784652) q[3];
sx q[3];
rz(-0.81955376) q[3];
sx q[3];
rz(1.6929172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95626962) q[0];
sx q[0];
rz(-1.7522426) q[0];
sx q[0];
rz(2.4638033) q[0];
rz(-0.14450821) q[1];
sx q[1];
rz(-0.78737193) q[1];
sx q[1];
rz(1.4871303) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22600284) q[0];
sx q[0];
rz(-1.5969513) q[0];
sx q[0];
rz(1.7067065) q[0];
rz(-3.0838548) q[2];
sx q[2];
rz(-1.9721834) q[2];
sx q[2];
rz(-0.58338469) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.37461149) q[1];
sx q[1];
rz(-2.0882147) q[1];
sx q[1];
rz(-0.53738885) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0934382) q[3];
sx q[3];
rz(-1.9461402) q[3];
sx q[3];
rz(2.1233692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3378478) q[2];
sx q[2];
rz(-2.8779112) q[2];
sx q[2];
rz(-2.7860876) q[2];
rz(-2.096094) q[3];
sx q[3];
rz(-1.239536) q[3];
sx q[3];
rz(-2.579328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1657555) q[0];
sx q[0];
rz(-2.5582357) q[0];
sx q[0];
rz(3.052886) q[0];
rz(1.5345796) q[1];
sx q[1];
rz(-1.3101703) q[1];
sx q[1];
rz(-0.075693695) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3189118) q[0];
sx q[0];
rz(-1.7633798) q[0];
sx q[0];
rz(-2.4680424) q[0];
rz(-1.6275109) q[2];
sx q[2];
rz(-0.73243388) q[2];
sx q[2];
rz(-2.5947941) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5867501) q[1];
sx q[1];
rz(-1.7126709) q[1];
sx q[1];
rz(-2.1290522) q[1];
x q[2];
rz(2.5646943) q[3];
sx q[3];
rz(-1.2954953) q[3];
sx q[3];
rz(-0.64361611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3169516) q[2];
sx q[2];
rz(-1.757916) q[2];
sx q[2];
rz(-1.0392044) q[2];
rz(0.21229395) q[3];
sx q[3];
rz(-2.0391235) q[3];
sx q[3];
rz(-1.645741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0742663) q[0];
sx q[0];
rz(-2.3889611) q[0];
sx q[0];
rz(-2.0972032) q[0];
rz(2.3692865) q[1];
sx q[1];
rz(-0.65985313) q[1];
sx q[1];
rz(2.4078802) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45424592) q[0];
sx q[0];
rz(-0.92763072) q[0];
sx q[0];
rz(-0.9439133) q[0];
x q[1];
rz(1.9416503) q[2];
sx q[2];
rz(-0.39526734) q[2];
sx q[2];
rz(-0.42759174) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0391239) q[1];
sx q[1];
rz(-0.3141292) q[1];
sx q[1];
rz(2.7565095) q[1];
x q[2];
rz(-2.3055607) q[3];
sx q[3];
rz(-1.6510909) q[3];
sx q[3];
rz(1.7952331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.61651984) q[2];
sx q[2];
rz(-0.51023444) q[2];
sx q[2];
rz(-0.60866848) q[2];
rz(-2.0742553) q[3];
sx q[3];
rz(-2.267024) q[3];
sx q[3];
rz(-0.55060351) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6390425) q[0];
sx q[0];
rz(-2.783343) q[0];
sx q[0];
rz(-1.4962366) q[0];
rz(-1.7773588) q[1];
sx q[1];
rz(-0.61449209) q[1];
sx q[1];
rz(0.38633698) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57418121) q[0];
sx q[0];
rz(-1.3076539) q[0];
sx q[0];
rz(-2.4469923) q[0];
x q[1];
rz(-1.7564098) q[2];
sx q[2];
rz(-1.5020554) q[2];
sx q[2];
rz(0.5663213) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.98294241) q[1];
sx q[1];
rz(-1.4483995) q[1];
sx q[1];
rz(-2.8084408) q[1];
rz(-pi) q[2];
rz(-1.531485) q[3];
sx q[3];
rz(-2.1609801) q[3];
sx q[3];
rz(-0.055824669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51314696) q[2];
sx q[2];
rz(-1.3725504) q[2];
sx q[2];
rz(2.9885542) q[2];
rz(0.89093527) q[3];
sx q[3];
rz(-2.1864083) q[3];
sx q[3];
rz(-2.4002767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9855758) q[0];
sx q[0];
rz(-0.22933904) q[0];
sx q[0];
rz(0.34737059) q[0];
rz(3.103745) q[1];
sx q[1];
rz(-1.1696576) q[1];
sx q[1];
rz(0.52245021) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.98451) q[0];
sx q[0];
rz(-1.1584131) q[0];
sx q[0];
rz(-3.0441112) q[0];
x q[1];
rz(-0.36549296) q[2];
sx q[2];
rz(-1.5710982) q[2];
sx q[2];
rz(3.0916758) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3047239) q[1];
sx q[1];
rz(-2.4608399) q[1];
sx q[1];
rz(0.70226837) q[1];
rz(-pi) q[2];
rz(0.96790989) q[3];
sx q[3];
rz(-0.65015745) q[3];
sx q[3];
rz(-3.0126743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5106421) q[2];
sx q[2];
rz(-0.46322552) q[2];
sx q[2];
rz(-1.9412712) q[2];
rz(-2.5837894) q[3];
sx q[3];
rz(-1.5188981) q[3];
sx q[3];
rz(-2.7571078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8484304) q[0];
sx q[0];
rz(-0.98386216) q[0];
sx q[0];
rz(1.7278607) q[0];
rz(-0.14104715) q[1];
sx q[1];
rz(-1.3755362) q[1];
sx q[1];
rz(2.486855) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.065718) q[0];
sx q[0];
rz(-1.362934) q[0];
sx q[0];
rz(-1.8130568) q[0];
rz(-0.59487307) q[2];
sx q[2];
rz(-2.0030177) q[2];
sx q[2];
rz(-3.1405666) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.480994) q[1];
sx q[1];
rz(-2.8428322) q[1];
sx q[1];
rz(-2.1887652) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.522024) q[3];
sx q[3];
rz(-1.1244785) q[3];
sx q[3];
rz(0.60947641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.35141382) q[2];
sx q[2];
rz(-1.979579) q[2];
sx q[2];
rz(2.4164825) q[2];
rz(2.2509947) q[3];
sx q[3];
rz(-2.2913439) q[3];
sx q[3];
rz(2.5543673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17726041) q[0];
sx q[0];
rz(-1.6163419) q[0];
sx q[0];
rz(1.1905715) q[0];
rz(-1.1217077) q[1];
sx q[1];
rz(-1.5678761) q[1];
sx q[1];
rz(-0.97490464) q[1];
rz(1.9323942) q[2];
sx q[2];
rz(-0.93930106) q[2];
sx q[2];
rz(1.3472547) q[2];
rz(1.411047) q[3];
sx q[3];
rz(-0.86426576) q[3];
sx q[3];
rz(-0.15440253) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
