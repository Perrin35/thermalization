OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3026128) q[0];
sx q[0];
rz(4.8470654) q[0];
sx q[0];
rz(9.6416311) q[0];
rz(2.6537553) q[1];
sx q[1];
rz(-2.2626329) q[1];
sx q[1];
rz(0.15329696) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76591208) q[0];
sx q[0];
rz(-1.3049676) q[0];
sx q[0];
rz(-1.8701118) q[0];
x q[1];
rz(-0.4699444) q[2];
sx q[2];
rz(-2.0389839) q[2];
sx q[2];
rz(0.17344698) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6930506) q[1];
sx q[1];
rz(-1.3227751) q[1];
sx q[1];
rz(-2.6297556) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7342773) q[3];
sx q[3];
rz(-2.370435) q[3];
sx q[3];
rz(0.60265356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.45120254) q[2];
sx q[2];
rz(-2.5610552) q[2];
sx q[2];
rz(-2.2817877) q[2];
rz(-2.9016923) q[3];
sx q[3];
rz(-2.4247215) q[3];
sx q[3];
rz(2.8587604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6567704) q[0];
sx q[0];
rz(-1.7253933) q[0];
sx q[0];
rz(2.2221478) q[0];
rz(0.8473618) q[1];
sx q[1];
rz(-0.58440009) q[1];
sx q[1];
rz(-1.1712317) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6028588) q[0];
sx q[0];
rz(-0.26358381) q[0];
sx q[0];
rz(1.7649505) q[0];
rz(-pi) q[1];
rz(-2.2861346) q[2];
sx q[2];
rz(-1.3163065) q[2];
sx q[2];
rz(1.364691) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1619809) q[1];
sx q[1];
rz(-1.2250326) q[1];
sx q[1];
rz(-2.8678721) q[1];
rz(-pi) q[2];
rz(-0.53816157) q[3];
sx q[3];
rz(-1.5892913) q[3];
sx q[3];
rz(1.5539813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.072711572) q[2];
sx q[2];
rz(-0.88901192) q[2];
sx q[2];
rz(1.9286801) q[2];
rz(-2.9291901) q[3];
sx q[3];
rz(-1.1143149) q[3];
sx q[3];
rz(-2.4685278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4251959) q[0];
sx q[0];
rz(-1.9028417) q[0];
sx q[0];
rz(-0.58309251) q[0];
rz(1.8816226) q[1];
sx q[1];
rz(-2.5446353) q[1];
sx q[1];
rz(-1.3139542) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3667917) q[0];
sx q[0];
rz(-1.9307025) q[0];
sx q[0];
rz(-0.33862305) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5384244) q[2];
sx q[2];
rz(-0.83725032) q[2];
sx q[2];
rz(-2.5493279) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.78216923) q[1];
sx q[1];
rz(-3.0170528) q[1];
sx q[1];
rz(-1.729639) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9745216) q[3];
sx q[3];
rz(-1.9937232) q[3];
sx q[3];
rz(1.058803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6446357) q[2];
sx q[2];
rz(-0.76377112) q[2];
sx q[2];
rz(0.53488564) q[2];
rz(-0.48318091) q[3];
sx q[3];
rz(-1.6202241) q[3];
sx q[3];
rz(3.0218637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0508761) q[0];
sx q[0];
rz(-0.058598761) q[0];
sx q[0];
rz(-1.4709877) q[0];
rz(-1.927467) q[1];
sx q[1];
rz(-1.43707) q[1];
sx q[1];
rz(0.4462744) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8295815) q[0];
sx q[0];
rz(-2.952842) q[0];
sx q[0];
rz(-1.944456) q[0];
rz(-3.1094279) q[2];
sx q[2];
rz(-2.6125561) q[2];
sx q[2];
rz(-0.47495237) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7301) q[1];
sx q[1];
rz(-0.92918438) q[1];
sx q[1];
rz(1.8909251) q[1];
rz(-pi) q[2];
rz(-2.8988188) q[3];
sx q[3];
rz(-1.5847237) q[3];
sx q[3];
rz(-0.33399912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.92397592) q[2];
sx q[2];
rz(-2.6805704) q[2];
sx q[2];
rz(2.8154742) q[2];
rz(0.41607949) q[3];
sx q[3];
rz(-1.6385498) q[3];
sx q[3];
rz(-2.3638341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.071082696) q[0];
sx q[0];
rz(-1.5639045) q[0];
sx q[0];
rz(0.34410205) q[0];
rz(-0.35585078) q[1];
sx q[1];
rz(-0.75332037) q[1];
sx q[1];
rz(2.8867302) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1321642) q[0];
sx q[0];
rz(-0.71064204) q[0];
sx q[0];
rz(-0.406213) q[0];
rz(-pi) q[1];
rz(-3.0414025) q[2];
sx q[2];
rz(-1.1104353) q[2];
sx q[2];
rz(3.0222297) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.566046) q[1];
sx q[1];
rz(-1.1351764) q[1];
sx q[1];
rz(-2.2244885) q[1];
x q[2];
rz(-1.9790824) q[3];
sx q[3];
rz(-0.84813839) q[3];
sx q[3];
rz(1.4681787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7039589) q[2];
sx q[2];
rz(-2.2767229) q[2];
sx q[2];
rz(-0.09058365) q[2];
rz(-2.3352052) q[3];
sx q[3];
rz(-1.839829) q[3];
sx q[3];
rz(-1.3069794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.059747132) q[0];
sx q[0];
rz(-1.2223926) q[0];
sx q[0];
rz(-2.5229689) q[0];
rz(-0.6126569) q[1];
sx q[1];
rz(-2.5109992) q[1];
sx q[1];
rz(0.15288606) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4893462) q[0];
sx q[0];
rz(-1.8033327) q[0];
sx q[0];
rz(1.3287956) q[0];
rz(0.15494187) q[2];
sx q[2];
rz(-0.58523387) q[2];
sx q[2];
rz(2.3335148) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5093435) q[1];
sx q[1];
rz(-1.9652838) q[1];
sx q[1];
rz(0.057827397) q[1];
rz(-pi) q[2];
rz(1.6083999) q[3];
sx q[3];
rz(-2.561224) q[3];
sx q[3];
rz(2.6036724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9360518) q[2];
sx q[2];
rz(-0.80465332) q[2];
sx q[2];
rz(-2.0118227) q[2];
rz(-1.0063082) q[3];
sx q[3];
rz(-1.8215424) q[3];
sx q[3];
rz(-1.6442851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0928918) q[0];
sx q[0];
rz(-2.0755656) q[0];
sx q[0];
rz(0.14376465) q[0];
rz(-0.20214209) q[1];
sx q[1];
rz(-0.527924) q[1];
sx q[1];
rz(-2.0595097) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6444958) q[0];
sx q[0];
rz(-1.93303) q[0];
sx q[0];
rz(-0.80279704) q[0];
rz(-pi) q[1];
rz(0.73418243) q[2];
sx q[2];
rz(-0.95962822) q[2];
sx q[2];
rz(-0.85684674) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.52817594) q[1];
sx q[1];
rz(-2.3966463) q[1];
sx q[1];
rz(1.2553535) q[1];
rz(-pi) q[2];
rz(2.4025926) q[3];
sx q[3];
rz(-1.0987875) q[3];
sx q[3];
rz(-1.7508185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9988592) q[2];
sx q[2];
rz(-1.1799246) q[2];
sx q[2];
rz(-0.72706968) q[2];
rz(1.8266228) q[3];
sx q[3];
rz(-1.246779) q[3];
sx q[3];
rz(-1.4962014) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8975248) q[0];
sx q[0];
rz(-2.0638564) q[0];
sx q[0];
rz(-1.9986073) q[0];
rz(-0.22363981) q[1];
sx q[1];
rz(-0.63839212) q[1];
sx q[1];
rz(1.9167831) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95668225) q[0];
sx q[0];
rz(-1.3089797) q[0];
sx q[0];
rz(-1.2774521) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4203577) q[2];
sx q[2];
rz(-1.2528361) q[2];
sx q[2];
rz(2.6810886) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8401007) q[1];
sx q[1];
rz(-2.9291541) q[1];
sx q[1];
rz(1.849624) q[1];
rz(-pi) q[2];
rz(2.5800621) q[3];
sx q[3];
rz(-1.6609471) q[3];
sx q[3];
rz(-0.60240373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5582922) q[2];
sx q[2];
rz(-1.4536828) q[2];
sx q[2];
rz(-2.3395786) q[2];
rz(1.924104) q[3];
sx q[3];
rz(-0.13974443) q[3];
sx q[3];
rz(1.2453992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8121346) q[0];
sx q[0];
rz(-1.3690925) q[0];
sx q[0];
rz(-1.580397) q[0];
rz(2.2691057) q[1];
sx q[1];
rz(-0.28631887) q[1];
sx q[1];
rz(-2.7365541) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1096538) q[0];
sx q[0];
rz(-2.4068916) q[0];
sx q[0];
rz(-0.74560179) q[0];
rz(-pi) q[1];
rz(-0.59258266) q[2];
sx q[2];
rz(-1.4854476) q[2];
sx q[2];
rz(2.6421955) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8499706) q[1];
sx q[1];
rz(-1.0128145) q[1];
sx q[1];
rz(-0.72429742) q[1];
rz(-pi) q[2];
rz(3.0108589) q[3];
sx q[3];
rz(-2.9933813) q[3];
sx q[3];
rz(3.1053901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8672455) q[2];
sx q[2];
rz(-0.87928191) q[2];
sx q[2];
rz(-0.61152968) q[2];
rz(1.3036171) q[3];
sx q[3];
rz(-2.7795064) q[3];
sx q[3];
rz(1.1360137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63854727) q[0];
sx q[0];
rz(-1.2940116) q[0];
sx q[0];
rz(-0.053330388) q[0];
rz(2.5697925) q[1];
sx q[1];
rz(-1.6245533) q[1];
sx q[1];
rz(1.2791876) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2032971) q[0];
sx q[0];
rz(-0.38559993) q[0];
sx q[0];
rz(2.4176265) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2917742) q[2];
sx q[2];
rz(-0.49683647) q[2];
sx q[2];
rz(1.2087087) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.51579976) q[1];
sx q[1];
rz(-0.98550057) q[1];
sx q[1];
rz(-0.093823508) q[1];
rz(-pi) q[2];
rz(-0.21802302) q[3];
sx q[3];
rz(-1.9338528) q[3];
sx q[3];
rz(-1.3040964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.46644396) q[2];
sx q[2];
rz(-1.6089336) q[2];
sx q[2];
rz(-0.68501985) q[2];
rz(2.3576665) q[3];
sx q[3];
rz(-0.42015606) q[3];
sx q[3];
rz(-0.73721957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6017799) q[0];
sx q[0];
rz(-1.89986) q[0];
sx q[0];
rz(1.4358406) q[0];
rz(0.80541366) q[1];
sx q[1];
rz(-2.6782811) q[1];
sx q[1];
rz(-2.4334999) q[1];
rz(2.5067301) q[2];
sx q[2];
rz(-0.35463968) q[2];
sx q[2];
rz(-1.126033) q[2];
rz(-2.7266034) q[3];
sx q[3];
rz(-1.4616953) q[3];
sx q[3];
rz(-2.2898522) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
