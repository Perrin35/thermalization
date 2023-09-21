OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.71198553) q[0];
sx q[0];
rz(-2.7349732) q[0];
sx q[0];
rz(-0.24917319) q[0];
rz(-0.063440032) q[1];
sx q[1];
rz(-2.1698706) q[1];
sx q[1];
rz(0.5501774) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8367856) q[0];
sx q[0];
rz(-1.3512304) q[0];
sx q[0];
rz(-0.027713393) q[0];
rz(-pi) q[1];
rz(-2.7031541) q[2];
sx q[2];
rz(-2.0619832) q[2];
sx q[2];
rz(0.05471281) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.3783592) q[1];
sx q[1];
rz(-2.2765056) q[1];
sx q[1];
rz(0.50904973) q[1];
x q[2];
rz(-0.95991858) q[3];
sx q[3];
rz(-1.5023408) q[3];
sx q[3];
rz(-1.5271036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7258519) q[2];
sx q[2];
rz(-2.6927413) q[2];
sx q[2];
rz(0.63981167) q[2];
rz(-0.85302991) q[3];
sx q[3];
rz(-2.5363688) q[3];
sx q[3];
rz(0.38133347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2663651) q[0];
sx q[0];
rz(-1.1667644) q[0];
sx q[0];
rz(-2.8711328) q[0];
rz(-0.71331435) q[1];
sx q[1];
rz(-2.1062873) q[1];
sx q[1];
rz(1.6289904) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4166116) q[0];
sx q[0];
rz(-2.1583301) q[0];
sx q[0];
rz(2.7532817) q[0];
rz(-pi) q[1];
rz(-3.1299635) q[2];
sx q[2];
rz(-2.8577869) q[2];
sx q[2];
rz(2.0276045) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1158893) q[1];
sx q[1];
rz(-1.70417) q[1];
sx q[1];
rz(-2.2504836) q[1];
x q[2];
rz(0.87725957) q[3];
sx q[3];
rz(-0.3650529) q[3];
sx q[3];
rz(0.045543268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2018532) q[2];
sx q[2];
rz(-0.24213232) q[2];
sx q[2];
rz(2.3382323) q[2];
rz(-1.057829) q[3];
sx q[3];
rz(-1.4927031) q[3];
sx q[3];
rz(0.025618205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.2597044) q[0];
sx q[0];
rz(-2.0443679) q[0];
sx q[0];
rz(-0.29552466) q[0];
rz(-2.9064536) q[1];
sx q[1];
rz(-1.4328911) q[1];
sx q[1];
rz(2.3957516) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64622067) q[0];
sx q[0];
rz(-1.6618068) q[0];
sx q[0];
rz(1.3489086) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4065811) q[2];
sx q[2];
rz(-2.0914372) q[2];
sx q[2];
rz(-2.8525713) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5513788) q[1];
sx q[1];
rz(-2.5775902) q[1];
sx q[1];
rz(2.6022634) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.48084534) q[3];
sx q[3];
rz(-0.50243176) q[3];
sx q[3];
rz(2.0744827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0535584) q[2];
sx q[2];
rz(-3.1159248) q[2];
sx q[2];
rz(-2.4528465) q[2];
rz(-3.0912494) q[3];
sx q[3];
rz(-0.91209948) q[3];
sx q[3];
rz(-1.6200199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.922309) q[0];
sx q[0];
rz(-2.121121) q[0];
sx q[0];
rz(0.10198378) q[0];
rz(-0.12022262) q[1];
sx q[1];
rz(-0.52821237) q[1];
sx q[1];
rz(-2.8682958) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7057987) q[0];
sx q[0];
rz(-1.4819643) q[0];
sx q[0];
rz(1.8949132) q[0];
x q[1];
rz(-1.0068514) q[2];
sx q[2];
rz(-1.0014921) q[2];
sx q[2];
rz(-3.0727) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.71972388) q[1];
sx q[1];
rz(-2.8641652) q[1];
sx q[1];
rz(-2.1348907) q[1];
x q[2];
rz(-1.7498383) q[3];
sx q[3];
rz(-0.92902196) q[3];
sx q[3];
rz(-0.87953506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3499202) q[2];
sx q[2];
rz(-2.4035954) q[2];
sx q[2];
rz(2.6376574) q[2];
rz(0.079581633) q[3];
sx q[3];
rz(-1.1527529) q[3];
sx q[3];
rz(-0.3751522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2824317) q[0];
sx q[0];
rz(-1.0520042) q[0];
sx q[0];
rz(0.082745634) q[0];
rz(0.67963183) q[1];
sx q[1];
rz(-1.6487164) q[1];
sx q[1];
rz(0.98714978) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25492451) q[0];
sx q[0];
rz(-0.75037557) q[0];
sx q[0];
rz(2.1819941) q[0];
x q[1];
rz(1.5802025) q[2];
sx q[2];
rz(-2.214589) q[2];
sx q[2];
rz(-0.26088342) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.134569) q[1];
sx q[1];
rz(-0.94788523) q[1];
sx q[1];
rz(-1.182215) q[1];
x q[2];
rz(-1.1779285) q[3];
sx q[3];
rz(-1.53189) q[3];
sx q[3];
rz(1.2272569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8905028) q[2];
sx q[2];
rz(-2.3513998) q[2];
sx q[2];
rz(-2.9303072) q[2];
rz(-0.42090297) q[3];
sx q[3];
rz(-2.588721) q[3];
sx q[3];
rz(-0.063407272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6565276) q[0];
sx q[0];
rz(-1.8873029) q[0];
sx q[0];
rz(-0.791839) q[0];
rz(0.99545288) q[1];
sx q[1];
rz(-0.96770006) q[1];
sx q[1];
rz(-1.8621559) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0148221) q[0];
sx q[0];
rz(-1.556067) q[0];
sx q[0];
rz(1.2889839) q[0];
rz(-2.5336669) q[2];
sx q[2];
rz(-1.8106917) q[2];
sx q[2];
rz(-1.8976256) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.8705604) q[1];
sx q[1];
rz(-0.522627) q[1];
sx q[1];
rz(-1.6305627) q[1];
rz(-pi) q[2];
rz(-1.389723) q[3];
sx q[3];
rz(-1.7098134) q[3];
sx q[3];
rz(0.74336038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8906158) q[2];
sx q[2];
rz(-1.0303409) q[2];
sx q[2];
rz(-2.3708169) q[2];
rz(1.6714913) q[3];
sx q[3];
rz(-0.41430587) q[3];
sx q[3];
rz(0.18338403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32271785) q[0];
sx q[0];
rz(-0.6495496) q[0];
sx q[0];
rz(0.055667002) q[0];
rz(0.21559134) q[1];
sx q[1];
rz(-2.3781653) q[1];
sx q[1];
rz(0.0035704426) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63821793) q[0];
sx q[0];
rz(-2.5711381) q[0];
sx q[0];
rz(2.4128782) q[0];
rz(-pi) q[1];
rz(1.1793169) q[2];
sx q[2];
rz(-0.85748312) q[2];
sx q[2];
rz(-1.6187402) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.53686245) q[1];
sx q[1];
rz(-2.4559048) q[1];
sx q[1];
rz(2.0525949) q[1];
x q[2];
rz(1.7690897) q[3];
sx q[3];
rz(-1.6453504) q[3];
sx q[3];
rz(1.9004746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5443762) q[2];
sx q[2];
rz(-0.95149136) q[2];
sx q[2];
rz(-0.92010951) q[2];
rz(-2.9428633) q[3];
sx q[3];
rz(-1.2343497) q[3];
sx q[3];
rz(-2.7238817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51628095) q[0];
sx q[0];
rz(-1.6315062) q[0];
sx q[0];
rz(-1.4177119) q[0];
rz(0.40813804) q[1];
sx q[1];
rz(-1.0393655) q[1];
sx q[1];
rz(-2.4628941) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5355276) q[0];
sx q[0];
rz(-1.8341656) q[0];
sx q[0];
rz(0.33959629) q[0];
x q[1];
rz(2.6762814) q[2];
sx q[2];
rz(-0.88854549) q[2];
sx q[2];
rz(-2.0899783) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7486836) q[1];
sx q[1];
rz(-0.52782413) q[1];
sx q[1];
rz(-1.6960359) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.33903665) q[3];
sx q[3];
rz(-2.6936274) q[3];
sx q[3];
rz(2.8807004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.17710182) q[2];
sx q[2];
rz(-2.047838) q[2];
sx q[2];
rz(2.2237681) q[2];
rz(1.9994036) q[3];
sx q[3];
rz(-0.85421383) q[3];
sx q[3];
rz(-2.2043622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1552102) q[0];
sx q[0];
rz(-0.6739524) q[0];
sx q[0];
rz(-2.881799) q[0];
rz(2.4329176) q[1];
sx q[1];
rz(-2.8627113) q[1];
sx q[1];
rz(-0.52694595) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5848815) q[0];
sx q[0];
rz(-1.7877794) q[0];
sx q[0];
rz(-2.7816539) q[0];
rz(-2.2286345) q[2];
sx q[2];
rz(-1.7857988) q[2];
sx q[2];
rz(-1.022033) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.770551) q[1];
sx q[1];
rz(-2.7275804) q[1];
sx q[1];
rz(2.4310334) q[1];
rz(2.5832289) q[3];
sx q[3];
rz(-0.51279587) q[3];
sx q[3];
rz(2.5496897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.32593411) q[2];
sx q[2];
rz(-2.2951173) q[2];
sx q[2];
rz(0.79088598) q[2];
rz(-0.38665006) q[3];
sx q[3];
rz(-1.0145885) q[3];
sx q[3];
rz(0.33716831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7336422) q[0];
sx q[0];
rz(-0.17245094) q[0];
sx q[0];
rz(-2.1561484) q[0];
rz(-2.573029) q[1];
sx q[1];
rz(-2.0303346) q[1];
sx q[1];
rz(-2.7808166) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0246968) q[0];
sx q[0];
rz(-1.403406) q[0];
sx q[0];
rz(0.025339729) q[0];
rz(-pi) q[1];
rz(-1.0656409) q[2];
sx q[2];
rz(-2.897122) q[2];
sx q[2];
rz(-0.62015647) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1211179) q[1];
sx q[1];
rz(-1.4282244) q[1];
sx q[1];
rz(0.74417443) q[1];
rz(-pi) q[2];
rz(-2.9721857) q[3];
sx q[3];
rz(-1.425404) q[3];
sx q[3];
rz(1.8190847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0976022) q[2];
sx q[2];
rz(-0.68949914) q[2];
sx q[2];
rz(2.1981751) q[2];
rz(-2.5676981) q[3];
sx q[3];
rz(-0.4807764) q[3];
sx q[3];
rz(-1.3963612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5744793) q[0];
sx q[0];
rz(-1.4470826) q[0];
sx q[0];
rz(-0.8599109) q[0];
rz(-1.3600596) q[1];
sx q[1];
rz(-2.139745) q[1];
sx q[1];
rz(2.3812961) q[1];
rz(3.0439539) q[2];
sx q[2];
rz(-1.3664445) q[2];
sx q[2];
rz(0.22899361) q[2];
rz(-0.75136649) q[3];
sx q[3];
rz(-2.1204227) q[3];
sx q[3];
rz(-1.942996) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
