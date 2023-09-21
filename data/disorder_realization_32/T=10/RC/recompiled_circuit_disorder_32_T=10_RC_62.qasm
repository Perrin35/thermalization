OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.10652868) q[0];
sx q[0];
rz(2.0523235) q[0];
sx q[0];
rz(9.2637445) q[0];
rz(1.610202) q[1];
sx q[1];
rz(2.664497) q[1];
sx q[1];
rz(8.9283979) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0954477) q[0];
sx q[0];
rz(-2.5382222) q[0];
sx q[0];
rz(1.3674111) q[0];
rz(-pi) q[1];
rz(2.8018087) q[2];
sx q[2];
rz(-1.3801563) q[2];
sx q[2];
rz(-1.9356188) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1094184) q[1];
sx q[1];
rz(-0.86583455) q[1];
sx q[1];
rz(1.4515619) q[1];
rz(0.22827893) q[3];
sx q[3];
rz(-0.41729673) q[3];
sx q[3];
rz(-1.8018064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6266142) q[2];
sx q[2];
rz(-2.2815621) q[2];
sx q[2];
rz(2.853945) q[2];
rz(1.7488165) q[3];
sx q[3];
rz(-2.1578622) q[3];
sx q[3];
rz(-1.1028517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2709133) q[0];
sx q[0];
rz(-0.7551071) q[0];
sx q[0];
rz(-0.57759181) q[0];
rz(-1.6060991) q[1];
sx q[1];
rz(-2.0237193) q[1];
sx q[1];
rz(1.577852) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0664415) q[0];
sx q[0];
rz(-1.6834714) q[0];
sx q[0];
rz(-0.22002797) q[0];
rz(-pi) q[1];
rz(2.8385542) q[2];
sx q[2];
rz(-1.8327456) q[2];
sx q[2];
rz(2.9532202) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.3040846) q[1];
sx q[1];
rz(-1.4834852) q[1];
sx q[1];
rz(1.1156032) q[1];
rz(-pi) q[2];
rz(1.7277667) q[3];
sx q[3];
rz(-1.7236712) q[3];
sx q[3];
rz(0.65519858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0343798) q[2];
sx q[2];
rz(-2.1652174) q[2];
sx q[2];
rz(-1.0822901) q[2];
rz(-1.1632495) q[3];
sx q[3];
rz(-1.2000368) q[3];
sx q[3];
rz(-2.5015586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0796233) q[0];
sx q[0];
rz(-2.847023) q[0];
sx q[0];
rz(0.18297718) q[0];
rz(3.0984763) q[1];
sx q[1];
rz(-0.95298195) q[1];
sx q[1];
rz(1.144369) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9336752) q[0];
sx q[0];
rz(-1.1753923) q[0];
sx q[0];
rz(2.140722) q[0];
x q[1];
rz(-1.794533) q[2];
sx q[2];
rz(-1.1615331) q[2];
sx q[2];
rz(2.4350016) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.80013393) q[1];
sx q[1];
rz(-2.027249) q[1];
sx q[1];
rz(1.9564499) q[1];
rz(-pi) q[2];
rz(-1.5738437) q[3];
sx q[3];
rz(-1.8584195) q[3];
sx q[3];
rz(2.0730719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.019505067) q[2];
sx q[2];
rz(-1.959789) q[2];
sx q[2];
rz(-2.2369177) q[2];
rz(0.71980643) q[3];
sx q[3];
rz(-0.42680877) q[3];
sx q[3];
rz(-0.75511801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4633789) q[0];
sx q[0];
rz(-0.97020522) q[0];
sx q[0];
rz(-2.4257207) q[0];
rz(-0.32575682) q[1];
sx q[1];
rz(-2.082086) q[1];
sx q[1];
rz(0.99266565) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9216114) q[0];
sx q[0];
rz(-2.5489759) q[0];
sx q[0];
rz(-0.48595925) q[0];
x q[1];
rz(-0.56447345) q[2];
sx q[2];
rz(-2.1285004) q[2];
sx q[2];
rz(1.0587453) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0992972) q[1];
sx q[1];
rz(-2.9534833) q[1];
sx q[1];
rz(-0.52109615) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2654815) q[3];
sx q[3];
rz(-0.70385859) q[3];
sx q[3];
rz(1.0775523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0585534) q[2];
sx q[2];
rz(-0.71648592) q[2];
sx q[2];
rz(-1.7129718) q[2];
rz(0.84609091) q[3];
sx q[3];
rz(-1.5139791) q[3];
sx q[3];
rz(-1.9184453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-2.5457526) q[0];
sx q[0];
rz(-1.4056453) q[0];
sx q[0];
rz(-2.3748421) q[0];
rz(-1.2402361) q[1];
sx q[1];
rz(-2.2940472) q[1];
sx q[1];
rz(0.3516745) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92847604) q[0];
sx q[0];
rz(-0.94841829) q[0];
sx q[0];
rz(-1.2147796) q[0];
x q[1];
rz(-0.93514498) q[2];
sx q[2];
rz(-1.1543373) q[2];
sx q[2];
rz(-2.2098429) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.60492086) q[1];
sx q[1];
rz(-0.39502883) q[1];
sx q[1];
rz(-0.037996304) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0293803) q[3];
sx q[3];
rz(-2.5684528) q[3];
sx q[3];
rz(1.0669607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9436283) q[2];
sx q[2];
rz(-1.6683656) q[2];
sx q[2];
rz(0.53331214) q[2];
rz(0.42896459) q[3];
sx q[3];
rz(-2.6140489) q[3];
sx q[3];
rz(-2.0146446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11739843) q[0];
sx q[0];
rz(-2.0485931) q[0];
sx q[0];
rz(-0.54164106) q[0];
rz(-2.533124) q[1];
sx q[1];
rz(-2.922373) q[1];
sx q[1];
rz(-1.0995964) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.106364) q[0];
sx q[0];
rz(-1.117825) q[0];
sx q[0];
rz(-2.7307672) q[0];
rz(0.18892388) q[2];
sx q[2];
rz(-2.0940229) q[2];
sx q[2];
rz(-2.1649233) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6096566) q[1];
sx q[1];
rz(-1.3927288) q[1];
sx q[1];
rz(0.68105662) q[1];
rz(1.908329) q[3];
sx q[3];
rz(-2.3300155) q[3];
sx q[3];
rz(-0.28901643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5114484) q[2];
sx q[2];
rz(-1.6161329) q[2];
sx q[2];
rz(-0.28298322) q[2];
rz(-0.7061559) q[3];
sx q[3];
rz(-1.5420087) q[3];
sx q[3];
rz(2.506536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0866078) q[0];
sx q[0];
rz(-2.1540756) q[0];
sx q[0];
rz(1.6116066) q[0];
rz(-0.64487547) q[1];
sx q[1];
rz(-2.6480643) q[1];
sx q[1];
rz(-2.9842916) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29042127) q[0];
sx q[0];
rz(-2.3010845) q[0];
sx q[0];
rz(-1.69676) q[0];
x q[1];
rz(0.83432014) q[2];
sx q[2];
rz(-1.7322707) q[2];
sx q[2];
rz(1.8404567) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.09307043) q[1];
sx q[1];
rz(-1.7659566) q[1];
sx q[1];
rz(2.1919769) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.39237202) q[3];
sx q[3];
rz(-0.6654159) q[3];
sx q[3];
rz(2.3371646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1823696) q[2];
sx q[2];
rz(-0.59252858) q[2];
sx q[2];
rz(-2.3106993) q[2];
rz(0.043878555) q[3];
sx q[3];
rz(-1.2336122) q[3];
sx q[3];
rz(0.20251814) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6458994) q[0];
sx q[0];
rz(-2.2024787) q[0];
sx q[0];
rz(-3.1307401) q[0];
rz(-2.8649578) q[1];
sx q[1];
rz(-2.3697772) q[1];
sx q[1];
rz(-2.8483134) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0115259) q[0];
sx q[0];
rz(-1.5648601) q[0];
sx q[0];
rz(-0.54242155) q[0];
rz(-2.2145055) q[2];
sx q[2];
rz(-1.4378528) q[2];
sx q[2];
rz(-1.5333652) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7640904) q[1];
sx q[1];
rz(-1.9891771) q[1];
sx q[1];
rz(2.3368895) q[1];
x q[2];
rz(-2.4268552) q[3];
sx q[3];
rz(-0.5849896) q[3];
sx q[3];
rz(1.294567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.098112) q[2];
sx q[2];
rz(-3.0467693) q[2];
sx q[2];
rz(3.0528255) q[2];
rz(0.25012112) q[3];
sx q[3];
rz(-1.1238031) q[3];
sx q[3];
rz(-0.5789825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1373238) q[0];
sx q[0];
rz(-2.9172638) q[0];
sx q[0];
rz(-1.0639169) q[0];
rz(2.6990199) q[1];
sx q[1];
rz(-1.7174218) q[1];
sx q[1];
rz(-1.3508505) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8352855) q[0];
sx q[0];
rz(-1.6500236) q[0];
sx q[0];
rz(-1.4863192) q[0];
rz(2.4814018) q[2];
sx q[2];
rz(-1.5056416) q[2];
sx q[2];
rz(0.24812631) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0722479) q[1];
sx q[1];
rz(-2.7069271) q[1];
sx q[1];
rz(1.6094631) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.001858) q[3];
sx q[3];
rz(-0.24340478) q[3];
sx q[3];
rz(2.7800625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0315447) q[2];
sx q[2];
rz(-1.8507379) q[2];
sx q[2];
rz(-0.44719493) q[2];
rz(1.3859008) q[3];
sx q[3];
rz(-2.8409676) q[3];
sx q[3];
rz(-0.51469222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8266325) q[0];
sx q[0];
rz(-1.6199912) q[0];
sx q[0];
rz(-2.3610624) q[0];
rz(-2.7138117) q[1];
sx q[1];
rz(-0.3573187) q[1];
sx q[1];
rz(-0.16960493) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0302089) q[0];
sx q[0];
rz(-1.2134524) q[0];
sx q[0];
rz(2.0765199) q[0];
x q[1];
rz(0.42189235) q[2];
sx q[2];
rz(-1.5483529) q[2];
sx q[2];
rz(2.9107712) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5369397) q[1];
sx q[1];
rz(-0.61652684) q[1];
sx q[1];
rz(-1.6517247) q[1];
rz(-1.2492368) q[3];
sx q[3];
rz(-1.7064629) q[3];
sx q[3];
rz(-0.78007573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2910989) q[2];
sx q[2];
rz(-1.1824181) q[2];
sx q[2];
rz(0.69236857) q[2];
rz(-0.46323562) q[3];
sx q[3];
rz(-1.654947) q[3];
sx q[3];
rz(-1.108981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.2904084) q[0];
sx q[0];
rz(-1.52607) q[0];
sx q[0];
rz(-2.4551256) q[0];
rz(0.51207536) q[1];
sx q[1];
rz(-0.33437406) q[1];
sx q[1];
rz(-2.1127111) q[1];
rz(-0.84531534) q[2];
sx q[2];
rz(-0.39388638) q[2];
sx q[2];
rz(-2.3245927) q[2];
rz(-3.0242596) q[3];
sx q[3];
rz(-1.3640079) q[3];
sx q[3];
rz(-0.97466536) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
