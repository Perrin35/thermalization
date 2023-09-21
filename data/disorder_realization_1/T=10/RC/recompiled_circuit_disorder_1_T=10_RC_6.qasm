OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.82436615) q[0];
sx q[0];
rz(-1.1146201) q[0];
sx q[0];
rz(-0.00014076509) q[0];
rz(-1.8074942) q[1];
sx q[1];
rz(-0.9642095) q[1];
sx q[1];
rz(1.948184) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6407335) q[0];
sx q[0];
rz(-1.8574323) q[0];
sx q[0];
rz(-2.9564234) q[0];
rz(-0.46618669) q[2];
sx q[2];
rz(-0.59980118) q[2];
sx q[2];
rz(0.28238645) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8661583) q[1];
sx q[1];
rz(-0.94238102) q[1];
sx q[1];
rz(2.1584312) q[1];
x q[2];
rz(-3.0317806) q[3];
sx q[3];
rz(-1.3545274) q[3];
sx q[3];
rz(0.11085489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6821735) q[2];
sx q[2];
rz(-0.023962263) q[2];
sx q[2];
rz(-1.9127282) q[2];
rz(1.7284283) q[3];
sx q[3];
rz(-2.0404405) q[3];
sx q[3];
rz(1.4878954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6035778) q[0];
sx q[0];
rz(-1.5025654) q[0];
sx q[0];
rz(2.1287825) q[0];
rz(0.027659841) q[1];
sx q[1];
rz(-0.67359567) q[1];
sx q[1];
rz(2.0181296) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2259953) q[0];
sx q[0];
rz(-1.628327) q[0];
sx q[0];
rz(2.0321839) q[0];
rz(-pi) q[1];
x q[1];
rz(0.063617184) q[2];
sx q[2];
rz(-0.78084313) q[2];
sx q[2];
rz(2.7240679) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0151129) q[1];
sx q[1];
rz(-1.0547332) q[1];
sx q[1];
rz(0.59241398) q[1];
rz(2.4554159) q[3];
sx q[3];
rz(-2.3585329) q[3];
sx q[3];
rz(-0.99883294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3479487) q[2];
sx q[2];
rz(-1.0898033) q[2];
sx q[2];
rz(2.222555) q[2];
rz(2.4675026) q[3];
sx q[3];
rz(-2.489311) q[3];
sx q[3];
rz(-1.526171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(0.27750257) q[0];
sx q[0];
rz(-0.16177495) q[0];
sx q[0];
rz(-1.2751689) q[0];
rz(-0.69349849) q[1];
sx q[1];
rz(-1.8854515) q[1];
sx q[1];
rz(1.1330053) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1904859) q[0];
sx q[0];
rz(-1.5740984) q[0];
sx q[0];
rz(1.4285054) q[0];
rz(-pi) q[1];
rz(-2.0793545) q[2];
sx q[2];
rz(-0.8478176) q[2];
sx q[2];
rz(-1.6194956) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.36724597) q[1];
sx q[1];
rz(-2.5173752) q[1];
sx q[1];
rz(1.063785) q[1];
rz(-pi) q[2];
x q[2];
rz(0.95136178) q[3];
sx q[3];
rz(-2.1054483) q[3];
sx q[3];
rz(1.149328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.2514078) q[2];
sx q[2];
rz(-2.3501985) q[2];
sx q[2];
rz(-1.2934925) q[2];
rz(3.1022762) q[3];
sx q[3];
rz(-1.2189564) q[3];
sx q[3];
rz(-1.2600651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.8816198) q[0];
sx q[0];
rz(-3.0631174) q[0];
sx q[0];
rz(1.9807293) q[0];
rz(-2.2456031) q[1];
sx q[1];
rz(-1.4410102) q[1];
sx q[1];
rz(-3.0060351) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4281222) q[0];
sx q[0];
rz(-2.2566416) q[0];
sx q[0];
rz(1.9786406) q[0];
rz(-pi) q[1];
rz(-1.1866456) q[2];
sx q[2];
rz(-1.640056) q[2];
sx q[2];
rz(2.4406976) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.23952661) q[1];
sx q[1];
rz(-0.24589989) q[1];
sx q[1];
rz(-0.86218254) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0849864) q[3];
sx q[3];
rz(-2.8286472) q[3];
sx q[3];
rz(-1.7900975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.23665145) q[2];
sx q[2];
rz(-0.94649482) q[2];
sx q[2];
rz(0.87990749) q[2];
rz(3.0974292) q[3];
sx q[3];
rz(-1.6396089) q[3];
sx q[3];
rz(-2.8529609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0376461) q[0];
sx q[0];
rz(-0.3750616) q[0];
sx q[0];
rz(1.0132382) q[0];
rz(3.0918616) q[1];
sx q[1];
rz(-2.2278992) q[1];
sx q[1];
rz(-2.0577046) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.097682) q[0];
sx q[0];
rz(-1.327276) q[0];
sx q[0];
rz(2.9527412) q[0];
x q[1];
rz(-2.3124144) q[2];
sx q[2];
rz(-0.36311705) q[2];
sx q[2];
rz(-1.6104289) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6127527) q[1];
sx q[1];
rz(-1.6788947) q[1];
sx q[1];
rz(-1.0428863) q[1];
rz(-pi) q[2];
rz(1.5203939) q[3];
sx q[3];
rz(-1.0574697) q[3];
sx q[3];
rz(0.36171519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.23285398) q[2];
sx q[2];
rz(-0.32662699) q[2];
sx q[2];
rz(2.8971635) q[2];
rz(-2.7092253) q[3];
sx q[3];
rz(-1.7418539) q[3];
sx q[3];
rz(-2.6385245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8571092) q[0];
sx q[0];
rz(-1.720022) q[0];
sx q[0];
rz(3.0474512) q[0];
rz(-0.17177467) q[1];
sx q[1];
rz(-2.005902) q[1];
sx q[1];
rz(0.89541268) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0953513) q[0];
sx q[0];
rz(-1.5329251) q[0];
sx q[0];
rz(2.7992159) q[0];
rz(-pi) q[1];
rz(-1.3307829) q[2];
sx q[2];
rz(-0.95108205) q[2];
sx q[2];
rz(1.0735219) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9858866) q[1];
sx q[1];
rz(-0.17427467) q[1];
sx q[1];
rz(-1.8272912) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1392519) q[3];
sx q[3];
rz(-1.4961092) q[3];
sx q[3];
rz(0.60126388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.133693) q[2];
sx q[2];
rz(-2.7329625) q[2];
sx q[2];
rz(-2.3383979) q[2];
rz(1.9512272) q[3];
sx q[3];
rz(-1.9093711) q[3];
sx q[3];
rz(-0.41263321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.068709277) q[0];
sx q[0];
rz(-0.16462737) q[0];
sx q[0];
rz(-0.51914006) q[0];
rz(2.5601162) q[1];
sx q[1];
rz(-1.1053718) q[1];
sx q[1];
rz(1.2566459) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4706659) q[0];
sx q[0];
rz(-1.5304655) q[0];
sx q[0];
rz(-1.3192024) q[0];
rz(-pi) q[1];
rz(-1.16876) q[2];
sx q[2];
rz(-0.49905825) q[2];
sx q[2];
rz(-1.4025584) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1756806) q[1];
sx q[1];
rz(-1.6351846) q[1];
sx q[1];
rz(0.28190159) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7729633) q[3];
sx q[3];
rz(-0.7437403) q[3];
sx q[3];
rz(-1.3524692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.98823035) q[2];
sx q[2];
rz(-2.1116657) q[2];
sx q[2];
rz(1.3640277) q[2];
rz(0.91056943) q[3];
sx q[3];
rz(-1.986859) q[3];
sx q[3];
rz(1.5301269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0751188) q[0];
sx q[0];
rz(-0.56448889) q[0];
sx q[0];
rz(-2.8334154) q[0];
rz(0.072487436) q[1];
sx q[1];
rz(-2.1283573) q[1];
sx q[1];
rz(0.38696188) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68650866) q[0];
sx q[0];
rz(-0.63502705) q[0];
sx q[0];
rz(1.9576661) q[0];
x q[1];
rz(2.2690291) q[2];
sx q[2];
rz(-1.9930895) q[2];
sx q[2];
rz(-2.1625105) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.19883979) q[1];
sx q[1];
rz(-1.6457335) q[1];
sx q[1];
rz(-0.67237512) q[1];
x q[2];
rz(1.3513226) q[3];
sx q[3];
rz(-1.0843127) q[3];
sx q[3];
rz(-2.6284077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5179634) q[2];
sx q[2];
rz(-1.364418) q[2];
sx q[2];
rz(-2.7015838) q[2];
rz(0.7157588) q[3];
sx q[3];
rz(-1.7093168) q[3];
sx q[3];
rz(-1.0796775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.14116645) q[0];
sx q[0];
rz(-2.3957802) q[0];
sx q[0];
rz(1.0986885) q[0];
rz(0.72775841) q[1];
sx q[1];
rz(-2.7658503) q[1];
sx q[1];
rz(0.049302014) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7002174) q[0];
sx q[0];
rz(-2.2451631) q[0];
sx q[0];
rz(0.81654878) q[0];
rz(0.98408913) q[2];
sx q[2];
rz(-0.77312914) q[2];
sx q[2];
rz(3.0423321) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.94433632) q[1];
sx q[1];
rz(-2.0347056) q[1];
sx q[1];
rz(-1.0524366) q[1];
rz(-pi) q[2];
rz(-0.74069549) q[3];
sx q[3];
rz(-1.1974088) q[3];
sx q[3];
rz(1.545056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0631642) q[2];
sx q[2];
rz(-0.59946632) q[2];
sx q[2];
rz(-2.4196529) q[2];
rz(-0.94349629) q[3];
sx q[3];
rz(-0.75073457) q[3];
sx q[3];
rz(2.8872484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9938875) q[0];
sx q[0];
rz(-1.9839956) q[0];
sx q[0];
rz(-2.06185) q[0];
rz(1.059277) q[1];
sx q[1];
rz(-0.22288999) q[1];
sx q[1];
rz(-1.4019029) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83090529) q[0];
sx q[0];
rz(-1.7831823) q[0];
sx q[0];
rz(-3.0124245) q[0];
x q[1];
rz(-1.5980814) q[2];
sx q[2];
rz(-1.7575043) q[2];
sx q[2];
rz(0.4208496) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.68114963) q[1];
sx q[1];
rz(-0.76421684) q[1];
sx q[1];
rz(-2.0185673) q[1];
rz(-pi) q[2];
rz(-0.63125061) q[3];
sx q[3];
rz(-1.4025941) q[3];
sx q[3];
rz(-1.2390618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6587276) q[2];
sx q[2];
rz(-1.3663224) q[2];
sx q[2];
rz(1.520291) q[2];
rz(-2.5907717) q[3];
sx q[3];
rz(-2.3362624) q[3];
sx q[3];
rz(-2.4441392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14810066) q[0];
sx q[0];
rz(-1.8363331) q[0];
sx q[0];
rz(1.6114417) q[0];
rz(2.2254754) q[1];
sx q[1];
rz(-2.5506908) q[1];
sx q[1];
rz(2.5509902) q[1];
rz(0.79402906) q[2];
sx q[2];
rz(-2.2326438) q[2];
sx q[2];
rz(2.2868962) q[2];
rz(-1.5394474) q[3];
sx q[3];
rz(-2.6242704) q[3];
sx q[3];
rz(-0.27192413) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
