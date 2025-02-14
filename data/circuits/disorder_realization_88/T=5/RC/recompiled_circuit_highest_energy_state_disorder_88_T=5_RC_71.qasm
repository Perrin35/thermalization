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
rz(-2.5992844) q[0];
sx q[0];
rz(-3.0071654) q[0];
sx q[0];
rz(-2.0943213) q[0];
rz(2.8174227) q[1];
sx q[1];
rz(-2.9919762) q[1];
sx q[1];
rz(-0.9447929) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1172663) q[0];
sx q[0];
rz(-2.7136682) q[0];
sx q[0];
rz(-0.86366349) q[0];
rz(-pi) q[1];
rz(0.9247577) q[2];
sx q[2];
rz(-1.7584718) q[2];
sx q[2];
rz(-0.061522324) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.45356645) q[1];
sx q[1];
rz(-2.0658148) q[1];
sx q[1];
rz(2.0903477) q[1];
rz(-pi) q[2];
rz(-3.0202623) q[3];
sx q[3];
rz(-2.0503732) q[3];
sx q[3];
rz(0.12313719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2904539) q[2];
sx q[2];
rz(-1.6893427) q[2];
sx q[2];
rz(1.5770844) q[2];
rz(-2.5751298) q[3];
sx q[3];
rz(-0.62140673) q[3];
sx q[3];
rz(-1.8433146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31747776) q[0];
sx q[0];
rz(-0.7190187) q[0];
sx q[0];
rz(0.7199921) q[0];
rz(-2.8672245) q[1];
sx q[1];
rz(-0.73148483) q[1];
sx q[1];
rz(-0.55363399) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8047129) q[0];
sx q[0];
rz(-1.2186482) q[0];
sx q[0];
rz(0.12356213) q[0];
x q[1];
rz(2.738094) q[2];
sx q[2];
rz(-1.2834233) q[2];
sx q[2];
rz(-1.0822051) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8588176) q[1];
sx q[1];
rz(-1.2056279) q[1];
sx q[1];
rz(-2.3479985) q[1];
rz(-pi) q[2];
x q[2];
rz(0.52218826) q[3];
sx q[3];
rz(-2.0117674) q[3];
sx q[3];
rz(2.9984634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1729892) q[2];
sx q[2];
rz(-1.2022688) q[2];
sx q[2];
rz(2.7776862) q[2];
rz(3.0299752) q[3];
sx q[3];
rz(-2.2026187) q[3];
sx q[3];
rz(0.76044559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97742057) q[0];
sx q[0];
rz(-2.3168679) q[0];
sx q[0];
rz(3.136193) q[0];
rz(-2.3258356) q[1];
sx q[1];
rz(-1.1382256) q[1];
sx q[1];
rz(-2.7556748) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2337795) q[0];
sx q[0];
rz(-2.695437) q[0];
sx q[0];
rz(-2.4491007) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0306401) q[2];
sx q[2];
rz(-0.6331501) q[2];
sx q[2];
rz(1.6748705) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.5792011) q[1];
sx q[1];
rz(-1.1913956) q[1];
sx q[1];
rz(-3.0048278) q[1];
x q[2];
rz(-1.98514) q[3];
sx q[3];
rz(-1.4477535) q[3];
sx q[3];
rz(-2.4407354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.58891121) q[2];
sx q[2];
rz(-2.3083355) q[2];
sx q[2];
rz(1.7595278) q[2];
rz(-2.9803993) q[3];
sx q[3];
rz(-1.7101945) q[3];
sx q[3];
rz(-0.62895697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2564119) q[0];
sx q[0];
rz(-1.2762524) q[0];
sx q[0];
rz(1.3389583) q[0];
rz(-1.3171875) q[1];
sx q[1];
rz(-0.97511292) q[1];
sx q[1];
rz(0.29979527) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9343963) q[0];
sx q[0];
rz(-2.3863433) q[0];
sx q[0];
rz(-0.51934262) q[0];
x q[1];
rz(-2.2251769) q[2];
sx q[2];
rz(-2.4407226) q[2];
sx q[2];
rz(-3.007421) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7924031) q[1];
sx q[1];
rz(-1.2196671) q[1];
sx q[1];
rz(-2.9955088) q[1];
rz(0.93500497) q[3];
sx q[3];
rz(-1.6775369) q[3];
sx q[3];
rz(0.028117953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6224711) q[2];
sx q[2];
rz(-0.93834472) q[2];
sx q[2];
rz(-0.31912121) q[2];
rz(-0.090911344) q[3];
sx q[3];
rz(-0.14486434) q[3];
sx q[3];
rz(-1.7849785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7262481) q[0];
sx q[0];
rz(-2.321796) q[0];
sx q[0];
rz(0.0023284624) q[0];
rz(-1.5709411) q[1];
sx q[1];
rz(-2.3410773) q[1];
sx q[1];
rz(1.9470107) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1462458) q[0];
sx q[0];
rz(-1.6888535) q[0];
sx q[0];
rz(-2.9857568) q[0];
rz(-pi) q[1];
rz(-0.04045893) q[2];
sx q[2];
rz(-2.2881857) q[2];
sx q[2];
rz(-0.77740146) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3227685) q[1];
sx q[1];
rz(-1.4684825) q[1];
sx q[1];
rz(-1.2248071) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8525328) q[3];
sx q[3];
rz(-1.0322555) q[3];
sx q[3];
rz(2.0592214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.883256) q[2];
sx q[2];
rz(-0.77815762) q[2];
sx q[2];
rz(0.077433132) q[2];
rz(2.1954913) q[3];
sx q[3];
rz(-0.48282048) q[3];
sx q[3];
rz(-2.6879123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2133863) q[0];
sx q[0];
rz(-3.0551857) q[0];
sx q[0];
rz(1.4148096) q[0];
rz(0.66697031) q[1];
sx q[1];
rz(-1.357115) q[1];
sx q[1];
rz(-0.40474969) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0694097) q[0];
sx q[0];
rz(-1.2889922) q[0];
sx q[0];
rz(2.826716) q[0];
rz(-pi) q[1];
rz(2.2282529) q[2];
sx q[2];
rz(-1.3619251) q[2];
sx q[2];
rz(0.058519017) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3500028) q[1];
sx q[1];
rz(-1.5096997) q[1];
sx q[1];
rz(1.1076124) q[1];
rz(2.1920106) q[3];
sx q[3];
rz(-2.5646665) q[3];
sx q[3];
rz(0.85437894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.82379597) q[2];
sx q[2];
rz(-2.3522289) q[2];
sx q[2];
rz(-2.6981567) q[2];
rz(-2.7277842) q[3];
sx q[3];
rz(-2.8518854) q[3];
sx q[3];
rz(-2.2558291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39171788) q[0];
sx q[0];
rz(-1.7789919) q[0];
sx q[0];
rz(0.66746563) q[0];
rz(-3.093847) q[1];
sx q[1];
rz(-2.0011438) q[1];
sx q[1];
rz(1.1183636) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55505468) q[0];
sx q[0];
rz(-1.0411052) q[0];
sx q[0];
rz(2.0425969) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7553545) q[2];
sx q[2];
rz(-2.5552351) q[2];
sx q[2];
rz(1.1589915) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0065919) q[1];
sx q[1];
rz(-1.9891095) q[1];
sx q[1];
rz(-2.5043284) q[1];
x q[2];
rz(1.4468071) q[3];
sx q[3];
rz(-0.50278864) q[3];
sx q[3];
rz(-3.0320252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.045791693) q[2];
sx q[2];
rz(-1.8818776) q[2];
sx q[2];
rz(0.0087139159) q[2];
rz(1.2039315) q[3];
sx q[3];
rz(-2.7976076) q[3];
sx q[3];
rz(3.1194527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.1190204) q[0];
sx q[0];
rz(-2.75596) q[0];
sx q[0];
rz(-0.88985306) q[0];
rz(-1.285137) q[1];
sx q[1];
rz(-1.0533918) q[1];
sx q[1];
rz(3.0056675) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1401299) q[0];
sx q[0];
rz(-0.23389947) q[0];
sx q[0];
rz(0.20568307) q[0];
rz(-pi) q[1];
x q[1];
rz(2.459002) q[2];
sx q[2];
rz(-0.75645489) q[2];
sx q[2];
rz(-2.3025049) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4299851) q[1];
sx q[1];
rz(-1.9031591) q[1];
sx q[1];
rz(-1.7794987) q[1];
x q[2];
rz(1.1952885) q[3];
sx q[3];
rz(-2.7673878) q[3];
sx q[3];
rz(-2.1458613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4385684) q[2];
sx q[2];
rz(-2.2615304) q[2];
sx q[2];
rz(2.0795889) q[2];
rz(0.92987531) q[3];
sx q[3];
rz(-2.3558741) q[3];
sx q[3];
rz(0.78833956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.7787665) q[0];
sx q[0];
rz(-0.39149785) q[0];
sx q[0];
rz(-0.40089259) q[0];
rz(-2.3800384) q[1];
sx q[1];
rz(-1.4925894) q[1];
sx q[1];
rz(-0.78899312) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3112199) q[0];
sx q[0];
rz(-1.5867751) q[0];
sx q[0];
rz(-0.1275087) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3775741) q[2];
sx q[2];
rz(-0.68422645) q[2];
sx q[2];
rz(0.88175899) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.096319936) q[1];
sx q[1];
rz(-2.5402656) q[1];
sx q[1];
rz(2.9296404) q[1];
rz(1.6203568) q[3];
sx q[3];
rz(-1.0465099) q[3];
sx q[3];
rz(1.6848642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.57364982) q[2];
sx q[2];
rz(-1.2929792) q[2];
sx q[2];
rz(2.4207777) q[2];
rz(1.4912262) q[3];
sx q[3];
rz(-2.3590915) q[3];
sx q[3];
rz(1.3097552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1184621) q[0];
sx q[0];
rz(-0.11496249) q[0];
sx q[0];
rz(0.7777099) q[0];
rz(-2.481781) q[1];
sx q[1];
rz(-0.86645627) q[1];
sx q[1];
rz(2.7515817) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4068425) q[0];
sx q[0];
rz(-1.6079385) q[0];
sx q[0];
rz(-1.7057306) q[0];
x q[1];
rz(2.397178) q[2];
sx q[2];
rz(-1.0420024) q[2];
sx q[2];
rz(0.86133682) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6811878) q[1];
sx q[1];
rz(-1.8680267) q[1];
sx q[1];
rz(1.1497208) q[1];
x q[2];
rz(2.635119) q[3];
sx q[3];
rz(-0.9829384) q[3];
sx q[3];
rz(-1.7755058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.68221349) q[2];
sx q[2];
rz(-0.43872681) q[2];
sx q[2];
rz(-2.6518346) q[2];
rz(-0.30719906) q[3];
sx q[3];
rz(-0.92370737) q[3];
sx q[3];
rz(-1.3908305) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3601892) q[0];
sx q[0];
rz(-1.6506945) q[0];
sx q[0];
rz(-1.6682464) q[0];
rz(-2.0669943) q[1];
sx q[1];
rz(-2.2886724) q[1];
sx q[1];
rz(1.2235175) q[1];
rz(-2.7424781) q[2];
sx q[2];
rz(-1.21798) q[2];
sx q[2];
rz(-1.1026364) q[2];
rz(2.4018039) q[3];
sx q[3];
rz(-2.3935912) q[3];
sx q[3];
rz(-1.145515) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
