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
rz(0.99474466) q[0];
sx q[0];
rz(3.7288546) q[0];
sx q[0];
rz(10.050339) q[0];
rz(0.62043959) q[1];
sx q[1];
rz(-2.151139) q[1];
sx q[1];
rz(0.78392309) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1358937) q[0];
sx q[0];
rz(-1.5672475) q[0];
sx q[0];
rz(0.43412848) q[0];
rz(-pi) q[1];
rz(-2.6033635) q[2];
sx q[2];
rz(-2.9799649) q[2];
sx q[2];
rz(1.5797576) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.16456024) q[1];
sx q[1];
rz(-2.1592618) q[1];
sx q[1];
rz(-0.87314897) q[1];
rz(-2.3638681) q[3];
sx q[3];
rz(-1.6334264) q[3];
sx q[3];
rz(-0.81111139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1938532) q[2];
sx q[2];
rz(-2.096602) q[2];
sx q[2];
rz(-2.4742773) q[2];
rz(2.8269178) q[3];
sx q[3];
rz(-2.5421725) q[3];
sx q[3];
rz(-2.2144894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84745234) q[0];
sx q[0];
rz(-2.2951234) q[0];
sx q[0];
rz(2.8417929) q[0];
rz(1.0046129) q[1];
sx q[1];
rz(-2.424898) q[1];
sx q[1];
rz(0.85003781) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73180985) q[0];
sx q[0];
rz(-0.96854051) q[0];
sx q[0];
rz(-3.0251363) q[0];
rz(-0.43434914) q[2];
sx q[2];
rz(-2.953989) q[2];
sx q[2];
rz(1.4570731) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5483173) q[1];
sx q[1];
rz(-1.2139582) q[1];
sx q[1];
rz(-0.22290454) q[1];
rz(0.032046958) q[3];
sx q[3];
rz(-0.40478313) q[3];
sx q[3];
rz(-1.7657775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5726418) q[2];
sx q[2];
rz(-0.53223842) q[2];
sx q[2];
rz(-2.4786095) q[2];
rz(2.3891383) q[3];
sx q[3];
rz(-1.821725) q[3];
sx q[3];
rz(-0.19645709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3982584) q[0];
sx q[0];
rz(-0.4011811) q[0];
sx q[0];
rz(-2.4617526) q[0];
rz(-0.72194779) q[1];
sx q[1];
rz(-0.55689055) q[1];
sx q[1];
rz(2.360875) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3328898) q[0];
sx q[0];
rz(-1.5784383) q[0];
sx q[0];
rz(-1.9714668) q[0];
rz(-pi) q[1];
rz(-3.137792) q[2];
sx q[2];
rz(-1.4846621) q[2];
sx q[2];
rz(-1.116556) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7739984) q[1];
sx q[1];
rz(-1.9450608) q[1];
sx q[1];
rz(0.67919977) q[1];
x q[2];
rz(1.8119007) q[3];
sx q[3];
rz(-2.3917197) q[3];
sx q[3];
rz(-2.1151224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.22689936) q[2];
sx q[2];
rz(-1.9136027) q[2];
sx q[2];
rz(2.4719888) q[2];
rz(-2.2412444) q[3];
sx q[3];
rz(-1.0122296) q[3];
sx q[3];
rz(0.77696925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8825619) q[0];
sx q[0];
rz(-2.7029523) q[0];
sx q[0];
rz(2.2341109) q[0];
rz(-2.5471845) q[1];
sx q[1];
rz(-0.92955697) q[1];
sx q[1];
rz(2.0118735) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0967207) q[0];
sx q[0];
rz(-1.7611928) q[0];
sx q[0];
rz(2.2119766) q[0];
x q[1];
rz(-2.9920861) q[2];
sx q[2];
rz(-2.0975916) q[2];
sx q[2];
rz(-0.12509987) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.65986827) q[1];
sx q[1];
rz(-2.4818083) q[1];
sx q[1];
rz(2.6901812) q[1];
rz(-1.2775189) q[3];
sx q[3];
rz(-1.3012816) q[3];
sx q[3];
rz(-0.5130918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.10924673) q[2];
sx q[2];
rz(-1.3999908) q[2];
sx q[2];
rz(-0.67231154) q[2];
rz(-0.0191056) q[3];
sx q[3];
rz(-2.8002383) q[3];
sx q[3];
rz(0.13489558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8371209) q[0];
sx q[0];
rz(-2.3470375) q[0];
sx q[0];
rz(0.16803148) q[0];
rz(-1.9145603) q[1];
sx q[1];
rz(-1.2201759) q[1];
sx q[1];
rz(-0.29774818) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7452324) q[0];
sx q[0];
rz(-1.0567291) q[0];
sx q[0];
rz(2.5395509) q[0];
rz(-0.99397583) q[2];
sx q[2];
rz(-2.0703531) q[2];
sx q[2];
rz(-3.1021038) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7682957) q[1];
sx q[1];
rz(-1.776564) q[1];
sx q[1];
rz(-1.5453669) q[1];
rz(1.566542) q[3];
sx q[3];
rz(-0.15878545) q[3];
sx q[3];
rz(-2.6095853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9510522) q[2];
sx q[2];
rz(-2.0420044) q[2];
sx q[2];
rz(-2.5229048) q[2];
rz(-2.3330073) q[3];
sx q[3];
rz(-2.6638668) q[3];
sx q[3];
rz(1.8019069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.488778) q[0];
sx q[0];
rz(-1.5050911) q[0];
sx q[0];
rz(2.6763647) q[0];
rz(2.6564927) q[1];
sx q[1];
rz(-2.7310889) q[1];
sx q[1];
rz(3.0533275) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8387091) q[0];
sx q[0];
rz(-1.1195991) q[0];
sx q[0];
rz(-2.124435) q[0];
rz(1.5563929) q[2];
sx q[2];
rz(-0.72897899) q[2];
sx q[2];
rz(-0.62504238) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3151206) q[1];
sx q[1];
rz(-1.7355647) q[1];
sx q[1];
rz(-0.75854782) q[1];
rz(-pi) q[2];
rz(1.6591822) q[3];
sx q[3];
rz(-1.2840349) q[3];
sx q[3];
rz(-2.4915316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7469067) q[2];
sx q[2];
rz(-2.2862819) q[2];
sx q[2];
rz(2.6439903) q[2];
rz(2.6615214) q[3];
sx q[3];
rz(-2.2205455) q[3];
sx q[3];
rz(-3.0729496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8467872) q[0];
sx q[0];
rz(-0.46600309) q[0];
sx q[0];
rz(1.4303327) q[0];
rz(-1.826674) q[1];
sx q[1];
rz(-0.39350915) q[1];
sx q[1];
rz(-1.5550782) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1436227) q[0];
sx q[0];
rz(-2.7300832) q[0];
sx q[0];
rz(-1.8726148) q[0];
x q[1];
rz(2.2261593) q[2];
sx q[2];
rz(-0.829773) q[2];
sx q[2];
rz(-0.96617736) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5615214) q[1];
sx q[1];
rz(-0.59021705) q[1];
sx q[1];
rz(-1.3886222) q[1];
rz(-pi) q[2];
rz(-1.4659381) q[3];
sx q[3];
rz(-0.63408454) q[3];
sx q[3];
rz(2.255267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.84307182) q[2];
sx q[2];
rz(-2.4303747) q[2];
sx q[2];
rz(0.92740518) q[2];
rz(-1.983042) q[3];
sx q[3];
rz(-2.4535593) q[3];
sx q[3];
rz(-3.0454175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7189099) q[0];
sx q[0];
rz(-0.89253187) q[0];
sx q[0];
rz(2.6655777) q[0];
rz(-2.1503275) q[1];
sx q[1];
rz(-0.74949336) q[1];
sx q[1];
rz(0.083866619) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.059444025) q[0];
sx q[0];
rz(-2.5184439) q[0];
sx q[0];
rz(2.2150458) q[0];
x q[1];
rz(-2.6798579) q[2];
sx q[2];
rz(-1.2424011) q[2];
sx q[2];
rz(-0.4730313) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.36271154) q[1];
sx q[1];
rz(-1.485752) q[1];
sx q[1];
rz(-1.9910732) q[1];
x q[2];
rz(-1.9355385) q[3];
sx q[3];
rz(-1.6514773) q[3];
sx q[3];
rz(-1.2176568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2241761) q[2];
sx q[2];
rz(-2.8755964) q[2];
sx q[2];
rz(-0.64845294) q[2];
rz(-2.2367541) q[3];
sx q[3];
rz(-1.4584533) q[3];
sx q[3];
rz(-2.8308433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7998578) q[0];
sx q[0];
rz(-0.71001995) q[0];
sx q[0];
rz(0.56890666) q[0];
rz(-2.3078602) q[1];
sx q[1];
rz(-1.8582452) q[1];
sx q[1];
rz(-2.1368829) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2153528) q[0];
sx q[0];
rz(-0.66074255) q[0];
sx q[0];
rz(0.31714971) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9453853) q[2];
sx q[2];
rz(-2.5718712) q[2];
sx q[2];
rz(1.8853055) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1596327) q[1];
sx q[1];
rz(-1.4985604) q[1];
sx q[1];
rz(-0.87149974) q[1];
rz(-3.0187739) q[3];
sx q[3];
rz(-2.5505722) q[3];
sx q[3];
rz(2.1048879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4363165) q[2];
sx q[2];
rz(-2.269561) q[2];
sx q[2];
rz(2.7936068) q[2];
rz(-2.3157388) q[3];
sx q[3];
rz(-2.9233942) q[3];
sx q[3];
rz(-2.5364449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0572877) q[0];
sx q[0];
rz(-1.0658406) q[0];
sx q[0];
rz(-2.9344946) q[0];
rz(-0.53889489) q[1];
sx q[1];
rz(-2.5205044) q[1];
sx q[1];
rz(-0.60212392) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83610632) q[0];
sx q[0];
rz(-2.0417111) q[0];
sx q[0];
rz(-0.54022809) q[0];
rz(-pi) q[1];
rz(1.1634356) q[2];
sx q[2];
rz(-2.1199787) q[2];
sx q[2];
rz(-2.801535) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4796743) q[1];
sx q[1];
rz(-1.9768324) q[1];
sx q[1];
rz(-3.0597294) q[1];
rz(-pi) q[2];
rz(2.3942457) q[3];
sx q[3];
rz(-1.3266194) q[3];
sx q[3];
rz(1.7908733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.78486097) q[2];
sx q[2];
rz(-2.1897903) q[2];
sx q[2];
rz(-1.3447364) q[2];
rz(-2.6209659) q[3];
sx q[3];
rz(-2.8713363) q[3];
sx q[3];
rz(-0.80210137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6945334) q[0];
sx q[0];
rz(-0.72493989) q[0];
sx q[0];
rz(-1.3865393) q[0];
rz(2.3583892) q[1];
sx q[1];
rz(-1.7192817) q[1];
sx q[1];
rz(-1.5789938) q[1];
rz(-2.798525) q[2];
sx q[2];
rz(-1.3160327) q[2];
sx q[2];
rz(2.119488) q[2];
rz(-0.047719638) q[3];
sx q[3];
rz(-0.24038355) q[3];
sx q[3];
rz(3.0265831) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
