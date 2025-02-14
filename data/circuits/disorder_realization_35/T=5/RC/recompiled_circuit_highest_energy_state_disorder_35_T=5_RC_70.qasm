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
rz(0.063989446) q[0];
sx q[0];
rz(4.0539157) q[0];
sx q[0];
rz(11.231448) q[0];
rz(-1.1420684) q[1];
sx q[1];
rz(5.7647709) q[1];
sx q[1];
rz(7.9328598) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8700884) q[0];
sx q[0];
rz(-0.056500204) q[0];
sx q[0];
rz(1.2369878) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4465898) q[2];
sx q[2];
rz(-2.0796513) q[2];
sx q[2];
rz(-0.63224224) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8644997) q[1];
sx q[1];
rz(-2.1756209) q[1];
sx q[1];
rz(-2.7983309) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.434695) q[3];
sx q[3];
rz(-2.8968662) q[3];
sx q[3];
rz(-0.32258209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3492744) q[2];
sx q[2];
rz(-0.26733843) q[2];
sx q[2];
rz(3.086669) q[2];
rz(-1.4548291) q[3];
sx q[3];
rz(-0.39595404) q[3];
sx q[3];
rz(1.6247862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2199663) q[0];
sx q[0];
rz(-1.8333789) q[0];
sx q[0];
rz(0.24800214) q[0];
rz(-1.2451046) q[1];
sx q[1];
rz(-0.73687941) q[1];
sx q[1];
rz(-1.2287593) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22511521) q[0];
sx q[0];
rz(-0.18813293) q[0];
sx q[0];
rz(0.1779025) q[0];
rz(-pi) q[1];
rz(0.51313674) q[2];
sx q[2];
rz(-1.0982795) q[2];
sx q[2];
rz(-3.0233011) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9039624) q[1];
sx q[1];
rz(-1.3591237) q[1];
sx q[1];
rz(-2.4304076) q[1];
x q[2];
rz(3.008083) q[3];
sx q[3];
rz(-2.7179681) q[3];
sx q[3];
rz(-2.6299919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9241141) q[2];
sx q[2];
rz(-0.83196297) q[2];
sx q[2];
rz(2.5577616) q[2];
rz(-0.42932388) q[3];
sx q[3];
rz(-1.9177633) q[3];
sx q[3];
rz(2.1302285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15447021) q[0];
sx q[0];
rz(-2.7543289) q[0];
sx q[0];
rz(-1.6744457) q[0];
rz(1.6636498) q[1];
sx q[1];
rz(-1.7324305) q[1];
sx q[1];
rz(-1.0999701) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3397148) q[0];
sx q[0];
rz(-2.9168282) q[0];
sx q[0];
rz(2.6532214) q[0];
rz(-pi) q[1];
rz(-1.9103545) q[2];
sx q[2];
rz(-1.8446088) q[2];
sx q[2];
rz(1.6664315) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.28274714) q[1];
sx q[1];
rz(-1.9044151) q[1];
sx q[1];
rz(-0.40550123) q[1];
rz(0.44902841) q[3];
sx q[3];
rz(-2.0981023) q[3];
sx q[3];
rz(0.15659404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.35525068) q[2];
sx q[2];
rz(-0.92266551) q[2];
sx q[2];
rz(-1.5938909) q[2];
rz(1.8111546) q[3];
sx q[3];
rz(-1.7682313) q[3];
sx q[3];
rz(1.0734585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-2.7610157) q[0];
sx q[0];
rz(-1.8208193) q[0];
sx q[0];
rz(2.9837578) q[0];
rz(-0.24179587) q[1];
sx q[1];
rz(-2.7561185) q[1];
sx q[1];
rz(0.48113021) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7173955) q[0];
sx q[0];
rz(-1.4823874) q[0];
sx q[0];
rz(2.428695) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60637668) q[2];
sx q[2];
rz(-0.42841347) q[2];
sx q[2];
rz(-2.1386752) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.85192902) q[1];
sx q[1];
rz(-0.73411513) q[1];
sx q[1];
rz(-1.0490225) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2838383) q[3];
sx q[3];
rz(-1.7718414) q[3];
sx q[3];
rz(-1.6359117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2699997) q[2];
sx q[2];
rz(-2.3303878) q[2];
sx q[2];
rz(0.2743741) q[2];
rz(-1.6756049) q[3];
sx q[3];
rz(-1.2803187) q[3];
sx q[3];
rz(2.2082641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6896553) q[0];
sx q[0];
rz(-2.2640197) q[0];
sx q[0];
rz(1.0614606) q[0];
rz(2.6929216) q[1];
sx q[1];
rz(-0.45495382) q[1];
sx q[1];
rz(1.4400858) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0827918) q[0];
sx q[0];
rz(-0.77547204) q[0];
sx q[0];
rz(2.916921) q[0];
x q[1];
rz(1.975073) q[2];
sx q[2];
rz(-1.6875522) q[2];
sx q[2];
rz(-1.5991885) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0764112) q[1];
sx q[1];
rz(-1.5440327) q[1];
sx q[1];
rz(-1.8288495) q[1];
x q[2];
rz(0.10847059) q[3];
sx q[3];
rz(-2.1785695) q[3];
sx q[3];
rz(0.25048192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7777286) q[2];
sx q[2];
rz(-1.2960459) q[2];
sx q[2];
rz(-1.9484005) q[2];
rz(0.39919546) q[3];
sx q[3];
rz(-1.4123071) q[3];
sx q[3];
rz(0.027755888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6804009) q[0];
sx q[0];
rz(-1.2874648) q[0];
sx q[0];
rz(-2.4611018) q[0];
rz(-2.3091799) q[1];
sx q[1];
rz(-1.8871555) q[1];
sx q[1];
rz(2.9852273) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.844961) q[0];
sx q[0];
rz(-0.28066844) q[0];
sx q[0];
rz(-2.9513717) q[0];
rz(-2.2920798) q[2];
sx q[2];
rz(-2.2597183) q[2];
sx q[2];
rz(-2.920814) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.029441) q[1];
sx q[1];
rz(-1.2425155) q[1];
sx q[1];
rz(2.8480457) q[1];
x q[2];
rz(1.0632319) q[3];
sx q[3];
rz(-0.2476736) q[3];
sx q[3];
rz(1.4192691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.51123315) q[2];
sx q[2];
rz(-0.66143051) q[2];
sx q[2];
rz(0.88999256) q[2];
rz(-0.47641274) q[3];
sx q[3];
rz(-2.2178631) q[3];
sx q[3];
rz(2.711003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7156859) q[0];
sx q[0];
rz(-1.3193193) q[0];
sx q[0];
rz(0.61221468) q[0];
rz(-1.7828434) q[1];
sx q[1];
rz(-0.76018676) q[1];
sx q[1];
rz(2.8403958) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4801203) q[0];
sx q[0];
rz(-1.299322) q[0];
sx q[0];
rz(-3.0508811) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2280548) q[2];
sx q[2];
rz(-1.2704256) q[2];
sx q[2];
rz(1.1235088) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0750904) q[1];
sx q[1];
rz(-1.9478223) q[1];
sx q[1];
rz(0.45555631) q[1];
x q[2];
rz(2.4047818) q[3];
sx q[3];
rz(-2.7932146) q[3];
sx q[3];
rz(-0.37955561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.907454) q[2];
sx q[2];
rz(-2.335304) q[2];
sx q[2];
rz(-1.0478728) q[2];
rz(-0.35510865) q[3];
sx q[3];
rz(-2.5706048) q[3];
sx q[3];
rz(-0.6855489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.909914) q[0];
sx q[0];
rz(-0.34220085) q[0];
sx q[0];
rz(-0.32383305) q[0];
rz(-2.553885) q[1];
sx q[1];
rz(-0.67981845) q[1];
sx q[1];
rz(0.30276611) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2229109) q[0];
sx q[0];
rz(-0.41707539) q[0];
sx q[0];
rz(-0.66584754) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5220171) q[2];
sx q[2];
rz(-1.80772) q[2];
sx q[2];
rz(1.1295484) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0628793) q[1];
sx q[1];
rz(-0.67725607) q[1];
sx q[1];
rz(-2.1724764) q[1];
rz(-1.9908526) q[3];
sx q[3];
rz(-1.8057293) q[3];
sx q[3];
rz(0.7701503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6451463) q[2];
sx q[2];
rz(-1.030693) q[2];
sx q[2];
rz(-2.5808064) q[2];
rz(-1.0049817) q[3];
sx q[3];
rz(-1.8533665) q[3];
sx q[3];
rz(-2.0463478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0305369) q[0];
sx q[0];
rz(-0.66937864) q[0];
sx q[0];
rz(0.74991599) q[0];
rz(0.38052446) q[1];
sx q[1];
rz(-2.1546202) q[1];
sx q[1];
rz(2.1662625) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0935287) q[0];
sx q[0];
rz(-2.2649327) q[0];
sx q[0];
rz(-0.82406901) q[0];
x q[1];
rz(0.060021518) q[2];
sx q[2];
rz(-1.2916471) q[2];
sx q[2];
rz(-2.699083) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.9561466) q[1];
sx q[1];
rz(-0.6022033) q[1];
sx q[1];
rz(2.5700997) q[1];
rz(2.9191429) q[3];
sx q[3];
rz(-1.8651903) q[3];
sx q[3];
rz(2.6840212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.19899496) q[2];
sx q[2];
rz(-2.2364605) q[2];
sx q[2];
rz(1.0814166) q[2];
rz(0.069247581) q[3];
sx q[3];
rz(-1.705575) q[3];
sx q[3];
rz(-0.92250219) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0104495) q[0];
sx q[0];
rz(-2.8157225) q[0];
sx q[0];
rz(2.8358054) q[0];
rz(0.30793134) q[1];
sx q[1];
rz(-1.4762069) q[1];
sx q[1];
rz(-1.9889132) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9891805) q[0];
sx q[0];
rz(-0.78762857) q[0];
sx q[0];
rz(-0.37352011) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2981554) q[2];
sx q[2];
rz(-0.66323384) q[2];
sx q[2];
rz(-0.014201268) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4977485) q[1];
sx q[1];
rz(-1.9744999) q[1];
sx q[1];
rz(-0.35404842) q[1];
x q[2];
rz(-1.8042555) q[3];
sx q[3];
rz(-0.065107927) q[3];
sx q[3];
rz(-0.59681915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.61047381) q[2];
sx q[2];
rz(-1.0627154) q[2];
sx q[2];
rz(-1.6884241) q[2];
rz(-2.6609663) q[3];
sx q[3];
rz(-0.56699816) q[3];
sx q[3];
rz(1.3948729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4881445) q[0];
sx q[0];
rz(-1.5339889) q[0];
sx q[0];
rz(1.4657159) q[0];
rz(2.0987971) q[1];
sx q[1];
rz(-1.7386309) q[1];
sx q[1];
rz(2.244619) q[1];
rz(1.8389134) q[2];
sx q[2];
rz(-2.0901879) q[2];
sx q[2];
rz(-1.2218634) q[2];
rz(0.10877175) q[3];
sx q[3];
rz(-0.36009195) q[3];
sx q[3];
rz(2.910955) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
