OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8389799) q[0];
sx q[0];
rz(-1.7054727) q[0];
sx q[0];
rz(-0.21685313) q[0];
rz(-3.6294301) q[1];
sx q[1];
rz(4.0205524) q[1];
sx q[1];
rz(9.5780749) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4175999) q[0];
sx q[0];
rz(-1.2823021) q[0];
sx q[0];
rz(2.8640076) q[0];
rz(-pi) q[1];
x q[1];
rz(0.4699444) q[2];
sx q[2];
rz(-2.0389839) q[2];
sx q[2];
rz(-0.17344698) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5342321) q[1];
sx q[1];
rz(-2.5776754) q[1];
sx q[1];
rz(-2.6644071) q[1];
rz(-pi) q[2];
rz(1.4073154) q[3];
sx q[3];
rz(-0.77115763) q[3];
sx q[3];
rz(-0.60265356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6903901) q[2];
sx q[2];
rz(-2.5610552) q[2];
sx q[2];
rz(-0.85980493) q[2];
rz(-0.23990038) q[3];
sx q[3];
rz(-2.4247215) q[3];
sx q[3];
rz(-2.8587604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6567704) q[0];
sx q[0];
rz(-1.4161994) q[0];
sx q[0];
rz(2.2221478) q[0];
rz(-0.8473618) q[1];
sx q[1];
rz(-0.58440009) q[1];
sx q[1];
rz(1.1712317) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5387339) q[0];
sx q[0];
rz(-2.8780088) q[0];
sx q[0];
rz(1.7649505) q[0];
x q[1];
rz(-0.85545808) q[2];
sx q[2];
rz(-1.8252862) q[2];
sx q[2];
rz(1.364691) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9796118) q[1];
sx q[1];
rz(-1.2250326) q[1];
sx q[1];
rz(0.27372056) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53816157) q[3];
sx q[3];
rz(-1.5892913) q[3];
sx q[3];
rz(1.5876113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0688811) q[2];
sx q[2];
rz(-2.2525807) q[2];
sx q[2];
rz(-1.9286801) q[2];
rz(2.9291901) q[3];
sx q[3];
rz(-2.0272777) q[3];
sx q[3];
rz(0.67306486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4251959) q[0];
sx q[0];
rz(-1.238751) q[0];
sx q[0];
rz(-2.5585001) q[0];
rz(-1.8816226) q[1];
sx q[1];
rz(-2.5446353) q[1];
sx q[1];
rz(-1.8276385) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2221916) q[0];
sx q[0];
rz(-1.254651) q[0];
sx q[0];
rz(-1.950398) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.035894265) q[2];
sx q[2];
rz(-2.4074656) q[2];
sx q[2];
rz(0.64059577) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3594234) q[1];
sx q[1];
rz(-0.12453989) q[1];
sx q[1];
rz(-1.729639) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1425972) q[3];
sx q[3];
rz(-1.4185662) q[3];
sx q[3];
rz(-0.58109944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6446357) q[2];
sx q[2];
rz(-0.76377112) q[2];
sx q[2];
rz(0.53488564) q[2];
rz(-0.48318091) q[3];
sx q[3];
rz(-1.5213685) q[3];
sx q[3];
rz(-3.0218637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0508761) q[0];
sx q[0];
rz(-0.058598761) q[0];
sx q[0];
rz(1.6706049) q[0];
rz(1.927467) q[1];
sx q[1];
rz(-1.43707) q[1];
sx q[1];
rz(2.6953183) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31201115) q[0];
sx q[0];
rz(-0.18875067) q[0];
sx q[0];
rz(-1.944456) q[0];
rz(0.52881119) q[2];
sx q[2];
rz(-1.5545648) q[2];
sx q[2];
rz(2.0179786) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2242291) q[1];
sx q[1];
rz(-0.70682061) q[1];
sx q[1];
rz(-0.39868211) q[1];
rz(-2/(11*pi)) q[3];
sx q[3];
rz(-0.24316517) q[3];
sx q[3];
rz(-1.9609708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.92397592) q[2];
sx q[2];
rz(-2.6805704) q[2];
sx q[2];
rz(-2.8154742) q[2];
rz(0.41607949) q[3];
sx q[3];
rz(-1.6385498) q[3];
sx q[3];
rz(0.7777586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.07051) q[0];
sx q[0];
rz(-1.5776881) q[0];
sx q[0];
rz(2.7974906) q[0];
rz(-2.7857419) q[1];
sx q[1];
rz(-0.75332037) q[1];
sx q[1];
rz(-2.8867302) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6159812) q[0];
sx q[0];
rz(-0.92824844) q[0];
sx q[0];
rz(1.8985974) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0331614) q[2];
sx q[2];
rz(-1.6605262) q[2];
sx q[2];
rz(-1.7347908) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5755467) q[1];
sx q[1];
rz(-2.0064163) q[1];
sx q[1];
rz(-0.91710414) q[1];
rz(-1.9790824) q[3];
sx q[3];
rz(-0.84813839) q[3];
sx q[3];
rz(-1.6734139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.43763375) q[2];
sx q[2];
rz(-2.2767229) q[2];
sx q[2];
rz(0.09058365) q[2];
rz(-2.3352052) q[3];
sx q[3];
rz(-1.3017637) q[3];
sx q[3];
rz(1.3069794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(3.0818455) q[0];
sx q[0];
rz(-1.2223926) q[0];
sx q[0];
rz(-0.61862373) q[0];
rz(2.5289358) q[1];
sx q[1];
rz(-2.5109992) q[1];
sx q[1];
rz(-2.9887066) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3090749) q[0];
sx q[0];
rz(-2.8075657) q[0];
sx q[0];
rz(-0.79132737) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4688849) q[2];
sx q[2];
rz(-2.1481124) q[2];
sx q[2];
rz(0.99330639) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.63224918) q[1];
sx q[1];
rz(-1.1763089) q[1];
sx q[1];
rz(3.0837653) q[1];
x q[2];
rz(2.1508407) q[3];
sx q[3];
rz(-1.5914122) q[3];
sx q[3];
rz(1.0643268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9360518) q[2];
sx q[2];
rz(-2.3369393) q[2];
sx q[2];
rz(-1.12977) q[2];
rz(2.1352844) q[3];
sx q[3];
rz(-1.3200503) q[3];
sx q[3];
rz(1.6442851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0928918) q[0];
sx q[0];
rz(-1.0660271) q[0];
sx q[0];
rz(2.997828) q[0];
rz(2.9394506) q[1];
sx q[1];
rz(-0.527924) q[1];
sx q[1];
rz(1.082083) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3976404) q[0];
sx q[0];
rz(-0.86371585) q[0];
sx q[0];
rz(-0.48488626) q[0];
rz(0.80777709) q[2];
sx q[2];
rz(-0.91726724) q[2];
sx q[2];
rz(1.2803004) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.94596568) q[1];
sx q[1];
rz(-2.2711427) q[1];
sx q[1];
rz(0.27865748) q[1];
x q[2];
rz(0.64855002) q[3];
sx q[3];
rz(-2.289384) q[3];
sx q[3];
rz(-0.64283338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.14273345) q[2];
sx q[2];
rz(-1.1799246) q[2];
sx q[2];
rz(0.72706968) q[2];
rz(1.8266228) q[3];
sx q[3];
rz(-1.246779) q[3];
sx q[3];
rz(-1.4962014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24406782) q[0];
sx q[0];
rz(-2.0638564) q[0];
sx q[0];
rz(1.1429853) q[0];
rz(-0.22363981) q[1];
sx q[1];
rz(-2.5032005) q[1];
sx q[1];
rz(-1.9167831) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95668225) q[0];
sx q[0];
rz(-1.3089797) q[0];
sx q[0];
rz(-1.8641406) q[0];
x q[1];
rz(0.72123493) q[2];
sx q[2];
rz(-1.2528361) q[2];
sx q[2];
rz(0.46050408) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.30149192) q[1];
sx q[1];
rz(-0.21243851) q[1];
sx q[1];
rz(1.2919687) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9734328) q[3];
sx q[3];
rz(-2.5736397) q[3];
sx q[3];
rz(1.1105383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.58330047) q[2];
sx q[2];
rz(-1.6879098) q[2];
sx q[2];
rz(2.3395786) q[2];
rz(-1.2174886) q[3];
sx q[3];
rz(-3.0018482) q[3];
sx q[3];
rz(-1.2453992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32945803) q[0];
sx q[0];
rz(-1.3690925) q[0];
sx q[0];
rz(1.580397) q[0];
rz(-2.2691057) q[1];
sx q[1];
rz(-0.28631887) q[1];
sx q[1];
rz(-0.40503851) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1096538) q[0];
sx q[0];
rz(-0.73470107) q[0];
sx q[0];
rz(0.74560179) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4680176) q[2];
sx q[2];
rz(-0.98066247) q[2];
sx q[2];
rz(-2.0128606) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.71724568) q[1];
sx q[1];
rz(-0.97386347) q[1];
sx q[1];
rz(2.2655377) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13073374) q[3];
sx q[3];
rz(-0.14821136) q[3];
sx q[3];
rz(-3.1053901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.27434719) q[2];
sx q[2];
rz(-0.87928191) q[2];
sx q[2];
rz(-0.61152968) q[2];
rz(1.8379755) q[3];
sx q[3];
rz(-0.36208624) q[3];
sx q[3];
rz(-2.005579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5030454) q[0];
sx q[0];
rz(-1.847581) q[0];
sx q[0];
rz(-0.053330388) q[0];
rz(0.57180014) q[1];
sx q[1];
rz(-1.6245533) q[1];
sx q[1];
rz(1.8624051) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.965248) q[0];
sx q[0];
rz(-1.8564447) q[0];
sx q[0];
rz(-1.3081416) q[0];
rz(-2.0512852) q[2];
sx q[2];
rz(-1.7024524) q[2];
sx q[2];
rz(-0.60881329) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4570791) q[1];
sx q[1];
rz(-0.59189905) q[1];
sx q[1];
rz(-1.4303703) q[1];
rz(1.9418777) q[3];
sx q[3];
rz(-1.7744007) q[3];
sx q[3];
rz(-2.9534087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6751487) q[2];
sx q[2];
rz(-1.5326591) q[2];
sx q[2];
rz(-2.4565728) q[2];
rz(2.3576665) q[3];
sx q[3];
rz(-0.42015606) q[3];
sx q[3];
rz(2.4043731) q[3];
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
rz(0.28975363) q[2];
sx q[2];
rz(-1.7782246) q[2];
sx q[2];
rz(-2.092337) q[2];
rz(-0.41498923) q[3];
sx q[3];
rz(-1.6798974) q[3];
sx q[3];
rz(0.85174042) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
