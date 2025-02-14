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
rz(0.4515689) q[0];
sx q[0];
rz(1.4181925) q[0];
sx q[0];
rz(12.137492) q[0];
rz(2.9984527) q[1];
sx q[1];
rz(-1.0909456) q[1];
sx q[1];
rz(-0.080088869) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9432491) q[0];
sx q[0];
rz(-3.0613073) q[0];
sx q[0];
rz(-2.9757418) q[0];
rz(-pi) q[1];
rz(1.3037138) q[2];
sx q[2];
rz(-0.3729698) q[2];
sx q[2];
rz(2.3139062) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.45177191) q[1];
sx q[1];
rz(-2.2597888) q[1];
sx q[1];
rz(-0.34386872) q[1];
rz(-pi) q[2];
rz(1.7350041) q[3];
sx q[3];
rz(-2.1421711) q[3];
sx q[3];
rz(-2.3624275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.074778883) q[2];
sx q[2];
rz(-1.0513693) q[2];
sx q[2];
rz(-1.4220413) q[2];
rz(1.413013) q[3];
sx q[3];
rz(-1.541411) q[3];
sx q[3];
rz(-1.6933256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0195352) q[0];
sx q[0];
rz(-1.3966565) q[0];
sx q[0];
rz(-2.8080217) q[0];
rz(2.9756599) q[1];
sx q[1];
rz(-0.3438147) q[1];
sx q[1];
rz(0.75210345) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.89276) q[0];
sx q[0];
rz(-2.016883) q[0];
sx q[0];
rz(-0.89527881) q[0];
rz(1.5444438) q[2];
sx q[2];
rz(-1.7894723) q[2];
sx q[2];
rz(-1.6353351) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.32622329) q[1];
sx q[1];
rz(-2.0953676) q[1];
sx q[1];
rz(0.12614819) q[1];
x q[2];
rz(-0.26176287) q[3];
sx q[3];
rz(-1.7883375) q[3];
sx q[3];
rz(2.1037773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.077945262) q[2];
sx q[2];
rz(-2.2293978) q[2];
sx q[2];
rz(-2.4845541) q[2];
rz(-2.4296203) q[3];
sx q[3];
rz(-2.4896121) q[3];
sx q[3];
rz(-2.3967192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7065358) q[0];
sx q[0];
rz(-2.5223795) q[0];
sx q[0];
rz(2.1001429) q[0];
rz(-0.374818) q[1];
sx q[1];
rz(-1.638214) q[1];
sx q[1];
rz(0.55694881) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9948363) q[0];
sx q[0];
rz(-0.86849156) q[0];
sx q[0];
rz(2.0392787) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1191423) q[2];
sx q[2];
rz(-1.6385534) q[2];
sx q[2];
rz(-2.6299143) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.16214806) q[1];
sx q[1];
rz(-1.4686161) q[1];
sx q[1];
rz(0.67047944) q[1];
rz(0.30681614) q[3];
sx q[3];
rz(-2.4673438) q[3];
sx q[3];
rz(1.1034654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6046784) q[2];
sx q[2];
rz(-0.92427173) q[2];
sx q[2];
rz(0.43517819) q[2];
rz(0.52470454) q[3];
sx q[3];
rz(-0.33354959) q[3];
sx q[3];
rz(0.13478002) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0404469) q[0];
sx q[0];
rz(-2.0772159) q[0];
sx q[0];
rz(1.6271628) q[0];
rz(-1.0487522) q[1];
sx q[1];
rz(-0.92275134) q[1];
sx q[1];
rz(-1.4357766) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1483106) q[0];
sx q[0];
rz(-1.3956426) q[0];
sx q[0];
rz(-1.6143285) q[0];
x q[1];
rz(2.5975512) q[2];
sx q[2];
rz(-0.61464192) q[2];
sx q[2];
rz(-0.34814542) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8183349) q[1];
sx q[1];
rz(-2.0916953) q[1];
sx q[1];
rz(3.0677615) q[1];
x q[2];
rz(2.1952329) q[3];
sx q[3];
rz(-2.588996) q[3];
sx q[3];
rz(0.70878562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.40020308) q[2];
sx q[2];
rz(-1.783739) q[2];
sx q[2];
rz(-3.1164361) q[2];
rz(1.4870421) q[3];
sx q[3];
rz(-0.39792037) q[3];
sx q[3];
rz(-2.3253843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0461034) q[0];
sx q[0];
rz(-2.4920721) q[0];
sx q[0];
rz(0.15897861) q[0];
rz(-1.1051296) q[1];
sx q[1];
rz(-0.72560328) q[1];
sx q[1];
rz(0.22430688) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.255757) q[0];
sx q[0];
rz(-2.1191264) q[0];
sx q[0];
rz(1.095039) q[0];
x q[1];
rz(-3.1087628) q[2];
sx q[2];
rz(-1.0107892) q[2];
sx q[2];
rz(3.1341886) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4558002) q[1];
sx q[1];
rz(-0.53736254) q[1];
sx q[1];
rz(1.3440029) q[1];
rz(-pi) q[2];
rz(-2.9486233) q[3];
sx q[3];
rz(-2.5558839) q[3];
sx q[3];
rz(-1.0568413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.20532456) q[2];
sx q[2];
rz(-1.3621796) q[2];
sx q[2];
rz(2.3929907) q[2];
rz(-3.0311106) q[3];
sx q[3];
rz(-1.9140849) q[3];
sx q[3];
rz(0.41142472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1614138) q[0];
sx q[0];
rz(-2.7730589) q[0];
sx q[0];
rz(-0.73032105) q[0];
rz(-1.5018564) q[1];
sx q[1];
rz(-1.7514936) q[1];
sx q[1];
rz(-2.9072442) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4880331) q[0];
sx q[0];
rz(-0.83439487) q[0];
sx q[0];
rz(-2.5279765) q[0];
rz(-pi) q[1];
rz(-2.2147398) q[2];
sx q[2];
rz(-1.0711462) q[2];
sx q[2];
rz(0.32678963) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4428574) q[1];
sx q[1];
rz(-2.9720829) q[1];
sx q[1];
rz(2.8560266) q[1];
rz(1.1612915) q[3];
sx q[3];
rz(-2.9886768) q[3];
sx q[3];
rz(-2.5482938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.948287) q[2];
sx q[2];
rz(-0.81393465) q[2];
sx q[2];
rz(1.2813655) q[2];
rz(0.60503259) q[3];
sx q[3];
rz(-2.0421959) q[3];
sx q[3];
rz(-1.2325149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.989711) q[0];
sx q[0];
rz(-1.879377) q[0];
sx q[0];
rz(0.61304098) q[0];
rz(0.68060654) q[1];
sx q[1];
rz(-1.1111518) q[1];
sx q[1];
rz(-1.3253164) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2277057) q[0];
sx q[0];
rz(-1.3038346) q[0];
sx q[0];
rz(0.56446979) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.36779495) q[2];
sx q[2];
rz(-0.22794524) q[2];
sx q[2];
rz(-1.4473372) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9195342) q[1];
sx q[1];
rz(-0.57812837) q[1];
sx q[1];
rz(0.23303746) q[1];
rz(-pi) q[2];
rz(-0.57215055) q[3];
sx q[3];
rz(-2.47743) q[3];
sx q[3];
rz(2.6900239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.94720542) q[2];
sx q[2];
rz(-1.9401865) q[2];
sx q[2];
rz(0.80805937) q[2];
rz(0.64588532) q[3];
sx q[3];
rz(-2.8736727) q[3];
sx q[3];
rz(0.3033692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5757669) q[0];
sx q[0];
rz(-2.4149826) q[0];
sx q[0];
rz(-0.44276825) q[0];
rz(0.59204656) q[1];
sx q[1];
rz(-0.7213842) q[1];
sx q[1];
rz(2.9659042) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73978776) q[0];
sx q[0];
rz(-1.2953838) q[0];
sx q[0];
rz(-2.7729183) q[0];
rz(0.08205361) q[2];
sx q[2];
rz(-2.6887165) q[2];
sx q[2];
rz(2.1575515) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1877039) q[1];
sx q[1];
rz(-2.4989681) q[1];
sx q[1];
rz(-3.0695555) q[1];
x q[2];
rz(2.1814303) q[3];
sx q[3];
rz(-1.9879793) q[3];
sx q[3];
rz(2.7593335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.59777933) q[2];
sx q[2];
rz(-0.74395776) q[2];
sx q[2];
rz(2.234484) q[2];
rz(-0.90726566) q[3];
sx q[3];
rz(-1.3943358) q[3];
sx q[3];
rz(-3.0239014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.9573145) q[0];
sx q[0];
rz(-0.90183455) q[0];
sx q[0];
rz(2.9515475) q[0];
rz(-1.5589145) q[1];
sx q[1];
rz(-1.5946486) q[1];
sx q[1];
rz(1.3072026) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.060576749) q[0];
sx q[0];
rz(-1.2071484) q[0];
sx q[0];
rz(0.15369291) q[0];
rz(-1.1747423) q[2];
sx q[2];
rz(-0.73985928) q[2];
sx q[2];
rz(-2.6137874) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4316201) q[1];
sx q[1];
rz(-3.0289972) q[1];
sx q[1];
rz(-1.8667875) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8966859) q[3];
sx q[3];
rz(-2.6953813) q[3];
sx q[3];
rz(-0.135016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.068453161) q[2];
sx q[2];
rz(-1.0871202) q[2];
sx q[2];
rz(2.9404409) q[2];
rz(-1.0226095) q[3];
sx q[3];
rz(-2.4590731) q[3];
sx q[3];
rz(-1.539591) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1744743) q[0];
sx q[0];
rz(-1.0459463) q[0];
sx q[0];
rz(2.2651267) q[0];
rz(0.62249741) q[1];
sx q[1];
rz(-2.9298156) q[1];
sx q[1];
rz(0.34271398) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4349608) q[0];
sx q[0];
rz(-1.6366658) q[0];
sx q[0];
rz(-0.94191282) q[0];
rz(1.9075526) q[2];
sx q[2];
rz(-1.4554664) q[2];
sx q[2];
rz(0.28240909) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5610159) q[1];
sx q[1];
rz(-2.3388712) q[1];
sx q[1];
rz(2.4965041) q[1];
x q[2];
rz(-1.8233772) q[3];
sx q[3];
rz(-1.1942721) q[3];
sx q[3];
rz(-1.8815709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9524625) q[2];
sx q[2];
rz(-1.4418944) q[2];
sx q[2];
rz(-2.6978037) q[2];
rz(1.1187547) q[3];
sx q[3];
rz(-1.5951364) q[3];
sx q[3];
rz(-0.55383033) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8583869) q[0];
sx q[0];
rz(-1.5571742) q[0];
sx q[0];
rz(-1.5207186) q[0];
rz(0.063477909) q[1];
sx q[1];
rz(-1.3451481) q[1];
sx q[1];
rz(1.4048911) q[1];
rz(-0.061894682) q[2];
sx q[2];
rz(-1.2409004) q[2];
sx q[2];
rz(0.093766669) q[2];
rz(-2.1370112) q[3];
sx q[3];
rz(-1.7773284) q[3];
sx q[3];
rz(-0.048684752) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
