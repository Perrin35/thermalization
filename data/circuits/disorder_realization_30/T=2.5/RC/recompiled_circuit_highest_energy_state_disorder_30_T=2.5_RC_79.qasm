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
rz(0.27424115) q[0];
sx q[0];
rz(-2.7354136) q[0];
sx q[0];
rz(0.22845536) q[0];
rz(0.38415456) q[1];
sx q[1];
rz(1.6698807) q[1];
sx q[1];
rz(8.9759965) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3224027) q[0];
sx q[0];
rz(-2.7459641) q[0];
sx q[0];
rz(2.7725793) q[0];
rz(-pi) q[1];
rz(-3.0182462) q[2];
sx q[2];
rz(-1.1575067) q[2];
sx q[2];
rz(-1.1878428) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8463414) q[1];
sx q[1];
rz(-1.59473) q[1];
sx q[1];
rz(1.6344535) q[1];
x q[2];
rz(-3.1136127) q[3];
sx q[3];
rz(-1.5767808) q[3];
sx q[3];
rz(-0.75388349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7479129) q[2];
sx q[2];
rz(-2.9069275) q[2];
sx q[2];
rz(-2.5133384) q[2];
rz(-2.0755532) q[3];
sx q[3];
rz(-0.99848905) q[3];
sx q[3];
rz(-2.0121241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.5265441) q[0];
sx q[0];
rz(-0.55416179) q[0];
sx q[0];
rz(2.2224485) q[0];
rz(2.9650086) q[1];
sx q[1];
rz(-1.44839) q[1];
sx q[1];
rz(-0.44763705) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8887605) q[0];
sx q[0];
rz(-1.5651784) q[0];
sx q[0];
rz(0.0075713099) q[0];
rz(-1.8488919) q[2];
sx q[2];
rz(-1.9534885) q[2];
sx q[2];
rz(1.1130321) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.094504823) q[1];
sx q[1];
rz(-0.89549235) q[1];
sx q[1];
rz(1.7750791) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0466879) q[3];
sx q[3];
rz(-1.1231622) q[3];
sx q[3];
rz(1.6448878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.24719079) q[2];
sx q[2];
rz(-1.9220507) q[2];
sx q[2];
rz(2.2578237) q[2];
rz(-2.4023458) q[3];
sx q[3];
rz(-1.2284307) q[3];
sx q[3];
rz(1.8993186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7582551) q[0];
sx q[0];
rz(-2.1012335) q[0];
sx q[0];
rz(-3.1357159) q[0];
rz(2.8351496) q[1];
sx q[1];
rz(-2.0649316) q[1];
sx q[1];
rz(1.8773016) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41691859) q[0];
sx q[0];
rz(-1.8643799) q[0];
sx q[0];
rz(0.39263607) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7948958) q[2];
sx q[2];
rz(-0.15264855) q[2];
sx q[2];
rz(-0.10748401) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3539999) q[1];
sx q[1];
rz(-1.3718604) q[1];
sx q[1];
rz(0.56942327) q[1];
rz(-pi) q[2];
x q[2];
rz(0.39671582) q[3];
sx q[3];
rz(-1.0457888) q[3];
sx q[3];
rz(-0.58806149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3197202) q[2];
sx q[2];
rz(-1.3786653) q[2];
sx q[2];
rz(2.482448) q[2];
rz(-1.871073) q[3];
sx q[3];
rz(-0.29265413) q[3];
sx q[3];
rz(1.6856153) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4996516) q[0];
sx q[0];
rz(-2.9872276) q[0];
sx q[0];
rz(0.9374215) q[0];
rz(1.591466) q[1];
sx q[1];
rz(-2.4731686) q[1];
sx q[1];
rz(-0.74916565) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5431108) q[0];
sx q[0];
rz(-2.3078663) q[0];
sx q[0];
rz(2.7373256) q[0];
x q[1];
rz(0.90355166) q[2];
sx q[2];
rz(-1.0507116) q[2];
sx q[2];
rz(-1.5992711) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.73203304) q[1];
sx q[1];
rz(-0.91425843) q[1];
sx q[1];
rz(2.892912) q[1];
rz(-0.76633664) q[3];
sx q[3];
rz(-1.727399) q[3];
sx q[3];
rz(-0.58567724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.59430993) q[2];
sx q[2];
rz(-0.82188934) q[2];
sx q[2];
rz(1.8853356) q[2];
rz(0.78993434) q[3];
sx q[3];
rz(-0.94090763) q[3];
sx q[3];
rz(-1.1676211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.054852) q[0];
sx q[0];
rz(-2.5331443) q[0];
sx q[0];
rz(-0.21859455) q[0];
rz(2.356148) q[1];
sx q[1];
rz(-2.0365069) q[1];
sx q[1];
rz(-1.8484263) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79695082) q[0];
sx q[0];
rz(-1.043129) q[0];
sx q[0];
rz(-2.5996741) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.082581994) q[2];
sx q[2];
rz(-1.3090927) q[2];
sx q[2];
rz(-2.341389) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4689815) q[1];
sx q[1];
rz(-2.1124351) q[1];
sx q[1];
rz(0.85503135) q[1];
rz(-1.7863705) q[3];
sx q[3];
rz(-1.6079992) q[3];
sx q[3];
rz(-2.6389909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.30696067) q[2];
sx q[2];
rz(-1.0647651) q[2];
sx q[2];
rz(1.628423) q[2];
rz(-2.4026134) q[3];
sx q[3];
rz(-1.8578953) q[3];
sx q[3];
rz(1.4449323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4479248) q[0];
sx q[0];
rz(-1.6532927) q[0];
sx q[0];
rz(-2.779575) q[0];
rz(-2.0118878) q[1];
sx q[1];
rz(-1.2672707) q[1];
sx q[1];
rz(0.76090181) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7603455) q[0];
sx q[0];
rz(-1.3804304) q[0];
sx q[0];
rz(2.3668853) q[0];
x q[1];
rz(-2.314662) q[2];
sx q[2];
rz(-1.880135) q[2];
sx q[2];
rz(-0.29238809) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.89448358) q[1];
sx q[1];
rz(-1.037408) q[1];
sx q[1];
rz(-2.8937312) q[1];
rz(-pi) q[2];
rz(2.5229955) q[3];
sx q[3];
rz(-1.2391802) q[3];
sx q[3];
rz(2.3823447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1165983) q[2];
sx q[2];
rz(-0.28433871) q[2];
sx q[2];
rz(1.9653758) q[2];
rz(-1.0233915) q[3];
sx q[3];
rz(-2.0659645) q[3];
sx q[3];
rz(2.7682176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0244331) q[0];
sx q[0];
rz(-0.42851448) q[0];
sx q[0];
rz(1.3065216) q[0];
rz(-3.040763) q[1];
sx q[1];
rz(-1.7497302) q[1];
sx q[1];
rz(-0.94591013) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27594179) q[0];
sx q[0];
rz(-1.8748594) q[0];
sx q[0];
rz(1.7070387) q[0];
rz(-pi) q[1];
rz(2.0761515) q[2];
sx q[2];
rz(-1.1123708) q[2];
sx q[2];
rz(-1.6363143) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.99092274) q[1];
sx q[1];
rz(-1.4353818) q[1];
sx q[1];
rz(1.5110498) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.85695388) q[3];
sx q[3];
rz(-1.3009239) q[3];
sx q[3];
rz(1.2601579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1176318) q[2];
sx q[2];
rz(-1.6098216) q[2];
sx q[2];
rz(0.35759932) q[2];
rz(3.0910953) q[3];
sx q[3];
rz(-2.7183618) q[3];
sx q[3];
rz(-0.38664204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9015775) q[0];
sx q[0];
rz(-1.8223263) q[0];
sx q[0];
rz(0.59558076) q[0];
rz(-1.9161842) q[1];
sx q[1];
rz(-1.775368) q[1];
sx q[1];
rz(-2.7092686) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5744988) q[0];
sx q[0];
rz(-1.4050975) q[0];
sx q[0];
rz(0.49543799) q[0];
rz(-pi) q[1];
rz(1.7244699) q[2];
sx q[2];
rz(-1.5072848) q[2];
sx q[2];
rz(0.85086289) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.22307587) q[1];
sx q[1];
rz(-2.1380414) q[1];
sx q[1];
rz(1.5722797) q[1];
rz(-pi) q[2];
x q[2];
rz(0.44345279) q[3];
sx q[3];
rz(-2.2041774) q[3];
sx q[3];
rz(-1.2284281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7633742) q[2];
sx q[2];
rz(-1.32064) q[2];
sx q[2];
rz(2.1136368) q[2];
rz(-1.6592735) q[3];
sx q[3];
rz(-1.6927398) q[3];
sx q[3];
rz(0.75215522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9529652) q[0];
sx q[0];
rz(-1.9099706) q[0];
sx q[0];
rz(0.46284437) q[0];
rz(1.2340744) q[1];
sx q[1];
rz(-1.9201098) q[1];
sx q[1];
rz(3*pi/8) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5355204) q[0];
sx q[0];
rz(-1.7458445) q[0];
sx q[0];
rz(-3.1158064) q[0];
rz(0.47204702) q[2];
sx q[2];
rz(-0.85523048) q[2];
sx q[2];
rz(0.89821076) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2888711) q[1];
sx q[1];
rz(-2.4203972) q[1];
sx q[1];
rz(-1.3836963) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.450348) q[3];
sx q[3];
rz(-1.9102194) q[3];
sx q[3];
rz(-1.2428997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.50596333) q[2];
sx q[2];
rz(-1.7545981) q[2];
sx q[2];
rz(2.3734234) q[2];
rz(-2.4913037) q[3];
sx q[3];
rz(-2.8322329) q[3];
sx q[3];
rz(-0.55997854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77785093) q[0];
sx q[0];
rz(-0.44054458) q[0];
sx q[0];
rz(1.1458696) q[0];
rz(-0.56529415) q[1];
sx q[1];
rz(-0.71563131) q[1];
sx q[1];
rz(0.40564767) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.250428) q[0];
sx q[0];
rz(-1.9180204) q[0];
sx q[0];
rz(0.77607147) q[0];
x q[1];
rz(2.8641607) q[2];
sx q[2];
rz(-1.9346969) q[2];
sx q[2];
rz(0.24574797) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.542681) q[1];
sx q[1];
rz(-0.57334954) q[1];
sx q[1];
rz(-0.47853761) q[1];
rz(-2.9156014) q[3];
sx q[3];
rz(-2.6922494) q[3];
sx q[3];
rz(-2.9087272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.046935) q[2];
sx q[2];
rz(-1.2311225) q[2];
sx q[2];
rz(2.6222353) q[2];
rz(-2.6804067) q[3];
sx q[3];
rz(-1.0591732) q[3];
sx q[3];
rz(0.69380277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2431348) q[0];
sx q[0];
rz(-1.5512137) q[0];
sx q[0];
rz(2.2358538) q[0];
rz(1.4836736) q[1];
sx q[1];
rz(-1.1653733) q[1];
sx q[1];
rz(1.9001874) q[1];
rz(-2.9979669) q[2];
sx q[2];
rz(-0.96114071) q[2];
sx q[2];
rz(-3.0428934) q[2];
rz(0.82571331) q[3];
sx q[3];
rz(-0.63586797) q[3];
sx q[3];
rz(-0.72365367) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
