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
rz(-1.0892692) q[0];
sx q[0];
rz(-2.9805592) q[0];
rz(1.610202) q[1];
sx q[1];
rz(2.664497) q[1];
sx q[1];
rz(8.9283979) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30713233) q[0];
sx q[0];
rz(-1.4559329) q[0];
sx q[0];
rz(0.97712028) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8018087) q[2];
sx q[2];
rz(-1.7614363) q[2];
sx q[2];
rz(-1.2059739) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6160994) q[1];
sx q[1];
rz(-1.4800737) q[1];
sx q[1];
rz(0.70848042) q[1];
rz(-pi) q[2];
rz(-1.6707889) q[3];
sx q[3];
rz(-1.9766207) q[3];
sx q[3];
rz(2.0506746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6266142) q[2];
sx q[2];
rz(-0.86003059) q[2];
sx q[2];
rz(-2.853945) q[2];
rz(-1.7488165) q[3];
sx q[3];
rz(-0.98373047) q[3];
sx q[3];
rz(2.0387409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2709133) q[0];
sx q[0];
rz(-2.3864855) q[0];
sx q[0];
rz(-2.5640008) q[0];
rz(1.5354935) q[1];
sx q[1];
rz(-2.0237193) q[1];
sx q[1];
rz(1.577852) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6710885) q[0];
sx q[0];
rz(-1.3521863) q[0];
sx q[0];
rz(1.6862306) q[0];
rz(-0.30303843) q[2];
sx q[2];
rz(-1.3088471) q[2];
sx q[2];
rz(-2.9532202) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4430122) q[1];
sx q[1];
rz(-2.6786782) q[1];
sx q[1];
rz(1.3742616) q[1];
rz(-2.9868449) q[3];
sx q[3];
rz(-1.7259211) q[3];
sx q[3];
rz(0.89150067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0343798) q[2];
sx q[2];
rz(-0.97637525) q[2];
sx q[2];
rz(-2.0593026) q[2];
rz(-1.9783431) q[3];
sx q[3];
rz(-1.9415559) q[3];
sx q[3];
rz(-2.5015586) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0619693) q[0];
sx q[0];
rz(-0.29456961) q[0];
sx q[0];
rz(0.18297718) q[0];
rz(-3.0984763) q[1];
sx q[1];
rz(-0.95298195) q[1];
sx q[1];
rz(-1.144369) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5367033) q[0];
sx q[0];
rz(-1.0495782) q[0];
sx q[0];
rz(-2.6813566) q[0];
x q[1];
rz(-0.41855721) q[2];
sx q[2];
rz(-1.7757799) q[2];
sx q[2];
rz(-2.1870854) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.548) q[1];
sx q[1];
rz(-1.2263733) q[1];
sx q[1];
rz(2.6542632) q[1];
rz(-1.5738437) q[3];
sx q[3];
rz(-1.2831732) q[3];
sx q[3];
rz(-2.0730719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.019505067) q[2];
sx q[2];
rz(-1.1818037) q[2];
sx q[2];
rz(-2.2369177) q[2];
rz(2.4217862) q[3];
sx q[3];
rz(-0.42680877) q[3];
sx q[3];
rz(-2.3864746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6782137) q[0];
sx q[0];
rz(-2.1713874) q[0];
sx q[0];
rz(-2.4257207) q[0];
rz(-0.32575682) q[1];
sx q[1];
rz(-2.082086) q[1];
sx q[1];
rz(0.99266565) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0621322) q[0];
sx q[0];
rz(-1.3068763) q[0];
sx q[0];
rz(-2.604565) q[0];
rz(0.86187141) q[2];
sx q[2];
rz(-0.77152354) q[2];
sx q[2];
rz(-2.9574403) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0992972) q[1];
sx q[1];
rz(-2.9534833) q[1];
sx q[1];
rz(-2.6204965) q[1];
rz(-pi) q[2];
rz(-2.2654815) q[3];
sx q[3];
rz(-2.4377341) q[3];
sx q[3];
rz(-1.0775523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0585534) q[2];
sx q[2];
rz(-0.71648592) q[2];
sx q[2];
rz(1.4286208) q[2];
rz(2.2955017) q[3];
sx q[3];
rz(-1.5139791) q[3];
sx q[3];
rz(1.9184453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.59584004) q[0];
sx q[0];
rz(-1.4056453) q[0];
sx q[0];
rz(0.76675057) q[0];
rz(-1.2402361) q[1];
sx q[1];
rz(-0.84754544) q[1];
sx q[1];
rz(2.7899182) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36067097) q[0];
sx q[0];
rz(-2.43649) q[0];
sx q[0];
rz(-0.4522001) q[0];
rz(-pi) q[1];
x q[1];
rz(0.93514498) q[2];
sx q[2];
rz(-1.9872553) q[2];
sx q[2];
rz(0.93174975) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.000948) q[1];
sx q[1];
rz(-1.556177) q[1];
sx q[1];
rz(-2.7468202) q[1];
rz(-2.8204927) q[3];
sx q[3];
rz(-2.0541111) q[3];
sx q[3];
rz(-1.6881642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.19796431) q[2];
sx q[2];
rz(-1.4732271) q[2];
sx q[2];
rz(2.6082805) q[2];
rz(2.7126281) q[3];
sx q[3];
rz(-2.6140489) q[3];
sx q[3];
rz(2.0146446) q[3];
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
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0241942) q[0];
sx q[0];
rz(-2.0485931) q[0];
sx q[0];
rz(0.54164106) q[0];
rz(0.60846865) q[1];
sx q[1];
rz(-2.922373) q[1];
sx q[1];
rz(-1.0995964) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8894419) q[0];
sx q[0];
rz(-2.5398206) q[0];
sx q[0];
rz(-2.2579231) q[0];
rz(-2.9526688) q[2];
sx q[2];
rz(-1.0475698) q[2];
sx q[2];
rz(2.1649233) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.53193608) q[1];
sx q[1];
rz(-1.7488639) q[1];
sx q[1];
rz(0.68105662) q[1];
rz(-pi) q[2];
x q[2];
rz(0.78824708) q[3];
sx q[3];
rz(-1.3282093) q[3];
sx q[3];
rz(1.5188252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5114484) q[2];
sx q[2];
rz(-1.6161329) q[2];
sx q[2];
rz(-0.28298322) q[2];
rz(-0.7061559) q[3];
sx q[3];
rz(-1.599584) q[3];
sx q[3];
rz(0.63505665) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054984897) q[0];
sx q[0];
rz(-2.1540756) q[0];
sx q[0];
rz(-1.5299861) q[0];
rz(-2.4967172) q[1];
sx q[1];
rz(-2.6480643) q[1];
sx q[1];
rz(-0.15730102) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10282117) q[0];
sx q[0];
rz(-0.73909315) q[0];
sx q[0];
rz(0.13939136) q[0];
rz(2.3072725) q[2];
sx q[2];
rz(-1.7322707) q[2];
sx q[2];
rz(-1.8404567) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7423319) q[1];
sx q[1];
rz(-2.4943588) q[1];
sx q[1];
rz(1.2433692) q[1];
rz(-pi) q[2];
rz(2.5141684) q[3];
sx q[3];
rz(-1.3324696) q[3];
sx q[3];
rz(-2.0605007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9592231) q[2];
sx q[2];
rz(-2.5490641) q[2];
sx q[2];
rz(2.3106993) q[2];
rz(0.043878555) q[3];
sx q[3];
rz(-1.2336122) q[3];
sx q[3];
rz(-2.9390745) q[3];
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
rz(-pi/2) q[0];
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
rz(-1.6458994) q[0];
sx q[0];
rz(-2.2024787) q[0];
sx q[0];
rz(3.1307401) q[0];
rz(-2.8649578) q[1];
sx q[1];
rz(-0.77181548) q[1];
sx q[1];
rz(-0.29327926) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54942185) q[0];
sx q[0];
rz(-0.54245078) q[0];
sx q[0];
rz(3.1300934) q[0];
rz(-pi) q[1];
rz(2.9759334) q[2];
sx q[2];
rz(-0.93369166) q[2];
sx q[2];
rz(-0.13656244) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7640904) q[1];
sx q[1];
rz(-1.1524156) q[1];
sx q[1];
rz(0.80470316) q[1];
x q[2];
rz(2.6777612) q[3];
sx q[3];
rz(-1.9411191) q[3];
sx q[3];
rz(2.2390389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.098112) q[2];
sx q[2];
rz(-0.094823368) q[2];
sx q[2];
rz(-0.088767178) q[2];
rz(-2.8914715) q[3];
sx q[3];
rz(-2.0177896) q[3];
sx q[3];
rz(-2.5626101) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1373238) q[0];
sx q[0];
rz(-0.22432888) q[0];
sx q[0];
rz(1.0639169) q[0];
rz(2.6990199) q[1];
sx q[1];
rz(-1.7174218) q[1];
sx q[1];
rz(1.7907422) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3063072) q[0];
sx q[0];
rz(-1.6500236) q[0];
sx q[0];
rz(-1.6552734) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66019085) q[2];
sx q[2];
rz(-1.5056416) q[2];
sx q[2];
rz(2.8934663) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.46637725) q[1];
sx q[1];
rz(-1.5870759) q[1];
sx q[1];
rz(2.0051763) q[1];
x q[2];
rz(-2.001858) q[3];
sx q[3];
rz(-2.8981879) q[3];
sx q[3];
rz(-2.7800625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0315447) q[2];
sx q[2];
rz(-1.2908547) q[2];
sx q[2];
rz(-2.6943977) q[2];
rz(1.7556919) q[3];
sx q[3];
rz(-0.30062506) q[3];
sx q[3];
rz(-0.51469222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8266325) q[0];
sx q[0];
rz(-1.6199912) q[0];
sx q[0];
rz(-0.78053027) q[0];
rz(-0.42778095) q[1];
sx q[1];
rz(-2.784274) q[1];
sx q[1];
rz(-0.16960493) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1113838) q[0];
sx q[0];
rz(-1.2134524) q[0];
sx q[0];
rz(1.0650728) q[0];
x q[1];
rz(1.5461966) q[2];
sx q[2];
rz(-1.9925756) q[2];
sx q[2];
rz(1.8116902) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1093724) q[1];
sx q[1];
rz(-1.5240372) q[1];
sx q[1];
rz(2.1857775) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2492368) q[3];
sx q[3];
rz(-1.4351298) q[3];
sx q[3];
rz(0.78007573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.85049373) q[2];
sx q[2];
rz(-1.9591745) q[2];
sx q[2];
rz(-2.4492241) q[2];
rz(0.46323562) q[3];
sx q[3];
rz(-1.654947) q[3];
sx q[3];
rz(1.108981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.2904084) q[0];
sx q[0];
rz(-1.52607) q[0];
sx q[0];
rz(-2.4551256) q[0];
rz(-2.6295173) q[1];
sx q[1];
rz(-0.33437406) q[1];
sx q[1];
rz(-2.1127111) q[1];
rz(2.2962773) q[2];
sx q[2];
rz(-0.39388638) q[2];
sx q[2];
rz(-2.3245927) q[2];
rz(3.0242596) q[3];
sx q[3];
rz(-1.7775848) q[3];
sx q[3];
rz(2.1669273) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
