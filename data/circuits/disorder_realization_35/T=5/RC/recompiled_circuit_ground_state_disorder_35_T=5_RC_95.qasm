OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72260296) q[0];
sx q[0];
rz(-2.498772) q[0];
sx q[0];
rz(8.2052054) q[0];
rz(0.83528432) q[1];
sx q[1];
rz(-1.3568027) q[1];
sx q[1];
rz(0.9019444) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.565958) q[0];
sx q[0];
rz(-2.1778641) q[0];
sx q[0];
rz(-1.7580887) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7612258) q[2];
sx q[2];
rz(-1.3970301) q[2];
sx q[2];
rz(-1.997239) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.29034258) q[1];
sx q[1];
rz(-2.800673) q[1];
sx q[1];
rz(-2.1147637) q[1];
x q[2];
rz(2.4113184) q[3];
sx q[3];
rz(-1.0515107) q[3];
sx q[3];
rz(1.9744557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.96282643) q[2];
sx q[2];
rz(-2.1361394) q[2];
sx q[2];
rz(-2.315305) q[2];
rz(1.7360342) q[3];
sx q[3];
rz(-2.8345351) q[3];
sx q[3];
rz(0.74339408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5590782) q[0];
sx q[0];
rz(-0.40575108) q[0];
sx q[0];
rz(1.7412809) q[0];
rz(-2.0179613) q[1];
sx q[1];
rz(-0.39389899) q[1];
sx q[1];
rz(2.0434911) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1207635) q[0];
sx q[0];
rz(-1.6801103) q[0];
sx q[0];
rz(0.1569605) q[0];
rz(-1.9015354) q[2];
sx q[2];
rz(-1.2501226) q[2];
sx q[2];
rz(1.641524) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9713116) q[1];
sx q[1];
rz(-2.9990497) q[1];
sx q[1];
rz(2.5851923) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.947883) q[3];
sx q[3];
rz(-1.415731) q[3];
sx q[3];
rz(2.0429275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0850247) q[2];
sx q[2];
rz(-0.3173863) q[2];
sx q[2];
rz(1.8918096) q[2];
rz(1.7791087) q[3];
sx q[3];
rz(-1.0008078) q[3];
sx q[3];
rz(-1.6574297) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63610858) q[0];
sx q[0];
rz(-2.4478069) q[0];
sx q[0];
rz(2.8662477) q[0];
rz(2.6823726) q[1];
sx q[1];
rz(-2.1870435) q[1];
sx q[1];
rz(-3.088248) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5864075) q[0];
sx q[0];
rz(-1.6490899) q[0];
sx q[0];
rz(-0.21348641) q[0];
rz(2.0317475) q[2];
sx q[2];
rz(-1.58733) q[2];
sx q[2];
rz(-1.8422045) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.03618212) q[1];
sx q[1];
rz(-2.1222669) q[1];
sx q[1];
rz(1.5230194) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2417913) q[3];
sx q[3];
rz(-1.7681554) q[3];
sx q[3];
rz(2.1579735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.82103819) q[2];
sx q[2];
rz(-0.36131636) q[2];
sx q[2];
rz(-1.6374755) q[2];
rz(-1.641168) q[3];
sx q[3];
rz(-2.0913561) q[3];
sx q[3];
rz(-0.27568451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42030537) q[0];
sx q[0];
rz(-2.5581701) q[0];
sx q[0];
rz(0.48459184) q[0];
rz(1.8991607) q[1];
sx q[1];
rz(-2.395605) q[1];
sx q[1];
rz(2.7925083) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8139127) q[0];
sx q[0];
rz(-2.591344) q[0];
sx q[0];
rz(1.6525066) q[0];
rz(1.7789087) q[2];
sx q[2];
rz(-2.7914064) q[2];
sx q[2];
rz(1.1454358) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.610747) q[1];
sx q[1];
rz(-1.7602923) q[1];
sx q[1];
rz(-2.9086472) q[1];
x q[2];
rz(-0.99557568) q[3];
sx q[3];
rz(-2.4565426) q[3];
sx q[3];
rz(3.0127061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.35859534) q[2];
sx q[2];
rz(-1.2025669) q[2];
sx q[2];
rz(2.6915468) q[2];
rz(0.83886823) q[3];
sx q[3];
rz(-1.5476371) q[3];
sx q[3];
rz(0.1990327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6435476) q[0];
sx q[0];
rz(-1.6723375) q[0];
sx q[0];
rz(-0.099844649) q[0];
rz(1.8210583) q[1];
sx q[1];
rz(-2.1456199) q[1];
sx q[1];
rz(2.5340396) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12814645) q[0];
sx q[0];
rz(-1.4960086) q[0];
sx q[0];
rz(0.083223968) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5962318) q[2];
sx q[2];
rz(-1.705745) q[2];
sx q[2];
rz(3.1317186) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.93170184) q[1];
sx q[1];
rz(-1.7875515) q[1];
sx q[1];
rz(2.5561129) q[1];
x q[2];
rz(1.2149548) q[3];
sx q[3];
rz(-2.8666593) q[3];
sx q[3];
rz(2.7235018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9472092) q[2];
sx q[2];
rz(-2.3781229) q[2];
sx q[2];
rz(0.53059951) q[2];
rz(-2.7001906) q[3];
sx q[3];
rz(-2.1680021) q[3];
sx q[3];
rz(-0.94943625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43668231) q[0];
sx q[0];
rz(-1.3446151) q[0];
sx q[0];
rz(-0.5740903) q[0];
rz(0.87999815) q[1];
sx q[1];
rz(-1.1209542) q[1];
sx q[1];
rz(1.7990187) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8594151) q[0];
sx q[0];
rz(-2.0382152) q[0];
sx q[0];
rz(2.9168081) q[0];
rz(-1.1583352) q[2];
sx q[2];
rz(-1.9762197) q[2];
sx q[2];
rz(-2.6087922) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.33730971) q[1];
sx q[1];
rz(-1.334467) q[1];
sx q[1];
rz(-1.9419036) q[1];
rz(-pi) q[2];
rz(0.66304147) q[3];
sx q[3];
rz(-1.460307) q[3];
sx q[3];
rz(2.7821531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4673956) q[2];
sx q[2];
rz(-2.8574222) q[2];
sx q[2];
rz(-2.6944842) q[2];
rz(-1.0304662) q[3];
sx q[3];
rz(-1.5117398) q[3];
sx q[3];
rz(0.38465056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36987385) q[0];
sx q[0];
rz(-1.2971224) q[0];
sx q[0];
rz(-0.41279992) q[0];
rz(-1.075047) q[1];
sx q[1];
rz(-1.1378891) q[1];
sx q[1];
rz(-0.80361754) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5301054) q[0];
sx q[0];
rz(-2.3854333) q[0];
sx q[0];
rz(-2.3458781) q[0];
rz(0.92757954) q[2];
sx q[2];
rz(-1.1021492) q[2];
sx q[2];
rz(1.8806632) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.072783006) q[1];
sx q[1];
rz(-0.78662614) q[1];
sx q[1];
rz(-1.5508316) q[1];
x q[2];
rz(0.39866205) q[3];
sx q[3];
rz(-2.9473262) q[3];
sx q[3];
rz(-0.27856058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.86218086) q[2];
sx q[2];
rz(-1.2181686) q[2];
sx q[2];
rz(-1.5009521) q[2];
rz(0.62266478) q[3];
sx q[3];
rz(-0.77087918) q[3];
sx q[3];
rz(-1.4890495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72614661) q[0];
sx q[0];
rz(-1.5203238) q[0];
sx q[0];
rz(2.5836482) q[0];
rz(0.56124148) q[1];
sx q[1];
rz(-1.3713501) q[1];
sx q[1];
rz(-1.4161313) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87420207) q[0];
sx q[0];
rz(-2.6010397) q[0];
sx q[0];
rz(-0.99733277) q[0];
rz(-pi) q[1];
rz(2.2165856) q[2];
sx q[2];
rz(-2.6753798) q[2];
sx q[2];
rz(-3.1104308) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0812644) q[1];
sx q[1];
rz(-2.7980479) q[1];
sx q[1];
rz(1.2765803) q[1];
rz(2.9245695) q[3];
sx q[3];
rz(-0.76730928) q[3];
sx q[3];
rz(1.2485152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.8884362) q[2];
sx q[2];
rz(-0.77740589) q[2];
sx q[2];
rz(1.0003132) q[2];
rz(-3.0194164) q[3];
sx q[3];
rz(-1.6413942) q[3];
sx q[3];
rz(0.15672556) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67824739) q[0];
sx q[0];
rz(-1.954701) q[0];
sx q[0];
rz(-0.004322411) q[0];
rz(2.9777572) q[1];
sx q[1];
rz(-2.7163353) q[1];
sx q[1];
rz(-1.77553) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3177494) q[0];
sx q[0];
rz(-1.5986406) q[0];
sx q[0];
rz(-2.2180024) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.37213426) q[2];
sx q[2];
rz(-1.2805802) q[2];
sx q[2];
rz(-1.8572825) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1613579) q[1];
sx q[1];
rz(-1.5346077) q[1];
sx q[1];
rz(-0.25118942) q[1];
rz(-2.2774599) q[3];
sx q[3];
rz(-2.257405) q[3];
sx q[3];
rz(-3.0326126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0020478) q[2];
sx q[2];
rz(-0.9674955) q[2];
sx q[2];
rz(2.5223993) q[2];
rz(0.55245095) q[3];
sx q[3];
rz(-1.3055472) q[3];
sx q[3];
rz(2.9964871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7904952) q[0];
sx q[0];
rz(-2.9572697) q[0];
sx q[0];
rz(0.47875324) q[0];
rz(1.9888196) q[1];
sx q[1];
rz(-2.8515127) q[1];
sx q[1];
rz(-0.94720381) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1789492) q[0];
sx q[0];
rz(-1.6048204) q[0];
sx q[0];
rz(-1.7555139) q[0];
rz(1.9302773) q[2];
sx q[2];
rz(-0.61865846) q[2];
sx q[2];
rz(1.930069) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3403907) q[1];
sx q[1];
rz(-1.2653192) q[1];
sx q[1];
rz(-2.3190772) q[1];
rz(-pi) q[2];
rz(-0.35263108) q[3];
sx q[3];
rz(-1.8049587) q[3];
sx q[3];
rz(1.682483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7572215) q[2];
sx q[2];
rz(-0.20558509) q[2];
sx q[2];
rz(-0.26665404) q[2];
rz(-1.0673149) q[3];
sx q[3];
rz(-1.2131194) q[3];
sx q[3];
rz(2.1783569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.08854475) q[0];
sx q[0];
rz(-1.6041258) q[0];
sx q[0];
rz(-1.5969101) q[0];
rz(-1.0531986) q[1];
sx q[1];
rz(-1.2146626) q[1];
sx q[1];
rz(0.76520898) q[1];
rz(-2.8996053) q[2];
sx q[2];
rz(-1.568071) q[2];
sx q[2];
rz(-0.65721401) q[2];
rz(-0.33452928) q[3];
sx q[3];
rz(-0.60548895) q[3];
sx q[3];
rz(2.0169712) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
