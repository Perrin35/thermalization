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
rz(-0.10997009) q[0];
sx q[0];
rz(-2.2697544) q[0];
sx q[0];
rz(0.81575552) q[0];
rz(-1.8176796) q[1];
sx q[1];
rz(-1.7426999) q[1];
sx q[1];
rz(0.89371347) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92233673) q[0];
sx q[0];
rz(-1.7667662) q[0];
sx q[0];
rz(0.45536572) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.42025153) q[2];
sx q[2];
rz(-0.21152341) q[2];
sx q[2];
rz(-2.6826721) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1964394) q[1];
sx q[1];
rz(-0.47633874) q[1];
sx q[1];
rz(-2.7193386) q[1];
x q[2];
rz(1.516926) q[3];
sx q[3];
rz(-1.3313378) q[3];
sx q[3];
rz(-0.62908632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.90014234) q[2];
sx q[2];
rz(-0.19108471) q[2];
sx q[2];
rz(-0.07766635) q[2];
rz(-0.11134722) q[3];
sx q[3];
rz(-2.1249378) q[3];
sx q[3];
rz(-2.1587423) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10649189) q[0];
sx q[0];
rz(-2.3973871) q[0];
sx q[0];
rz(2.8700854) q[0];
rz(-0.60596451) q[1];
sx q[1];
rz(-0.4393591) q[1];
sx q[1];
rz(2.4966168) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5640372) q[0];
sx q[0];
rz(-2.5112069) q[0];
sx q[0];
rz(-1.0159929) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.740005) q[2];
sx q[2];
rz(-2.2300917) q[2];
sx q[2];
rz(0.35091695) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.067459) q[1];
sx q[1];
rz(-1.7927884) q[1];
sx q[1];
rz(1.4253474) q[1];
rz(-pi) q[2];
rz(-0.88949253) q[3];
sx q[3];
rz(-0.23372641) q[3];
sx q[3];
rz(-2.1672578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6509387) q[2];
sx q[2];
rz(-1.1488612) q[2];
sx q[2];
rz(-2.8059354) q[2];
rz(-3.007174) q[3];
sx q[3];
rz(-0.69177827) q[3];
sx q[3];
rz(-0.604983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8721039) q[0];
sx q[0];
rz(-1.1815) q[0];
sx q[0];
rz(2.9961725) q[0];
rz(-0.91436404) q[1];
sx q[1];
rz(-1.7040355) q[1];
sx q[1];
rz(0.61633715) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5443104) q[0];
sx q[0];
rz(-0.72120404) q[0];
sx q[0];
rz(-2.3116259) q[0];
rz(3.0923244) q[2];
sx q[2];
rz(-2.081909) q[2];
sx q[2];
rz(-0.12898239) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7246252) q[1];
sx q[1];
rz(-2.6706302) q[1];
sx q[1];
rz(-2.4344786) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5068153) q[3];
sx q[3];
rz(-1.758681) q[3];
sx q[3];
rz(-1.5646936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1472037) q[2];
sx q[2];
rz(-2.2115967) q[2];
sx q[2];
rz(2.8623016) q[2];
rz(0.3768557) q[3];
sx q[3];
rz(-2.2241204) q[3];
sx q[3];
rz(-2.3963624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.299861) q[0];
sx q[0];
rz(-1.0693411) q[0];
sx q[0];
rz(-1.8002864) q[0];
rz(0.74814859) q[1];
sx q[1];
rz(-0.52296269) q[1];
sx q[1];
rz(2.0789007) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5077909) q[0];
sx q[0];
rz(-1.4550545) q[0];
sx q[0];
rz(1.0032038) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2466561) q[2];
sx q[2];
rz(-2.3445355) q[2];
sx q[2];
rz(-2.2618798) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6036883) q[1];
sx q[1];
rz(-0.66715566) q[1];
sx q[1];
rz(1.0546507) q[1];
rz(-pi) q[2];
rz(-0.65819169) q[3];
sx q[3];
rz(-2.255283) q[3];
sx q[3];
rz(1.4475065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.62715379) q[2];
sx q[2];
rz(-1.2320765) q[2];
sx q[2];
rz(-0.76986924) q[2];
rz(-1.7047966) q[3];
sx q[3];
rz(-1.0154279) q[3];
sx q[3];
rz(-2.3427486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1771667) q[0];
sx q[0];
rz(-1.8317969) q[0];
sx q[0];
rz(-2.8231296) q[0];
rz(-0.98694363) q[1];
sx q[1];
rz(-1.3401597) q[1];
sx q[1];
rz(-2.9950704) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40020271) q[0];
sx q[0];
rz(-0.70841778) q[0];
sx q[0];
rz(-1.0111459) q[0];
rz(-pi) q[1];
x q[1];
rz(0.71208565) q[2];
sx q[2];
rz(-2.0566517) q[2];
sx q[2];
rz(-2.8949702) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9695786) q[1];
sx q[1];
rz(-0.41436985) q[1];
sx q[1];
rz(-2.5408391) q[1];
rz(2.745689) q[3];
sx q[3];
rz(-0.77953458) q[3];
sx q[3];
rz(-0.59880873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.58586183) q[2];
sx q[2];
rz(-0.18438688) q[2];
sx q[2];
rz(-1.2510074) q[2];
rz(-0.74293724) q[3];
sx q[3];
rz(-1.6299959) q[3];
sx q[3];
rz(1.9353297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10155216) q[0];
sx q[0];
rz(-2.0843625) q[0];
sx q[0];
rz(1.2054766) q[0];
rz(-0.34791738) q[1];
sx q[1];
rz(-1.1707062) q[1];
sx q[1];
rz(-1.9084557) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3065227) q[0];
sx q[0];
rz(-0.26598922) q[0];
sx q[0];
rz(2.8138922) q[0];
rz(-pi) q[1];
rz(2.6377524) q[2];
sx q[2];
rz(-2.2723393) q[2];
sx q[2];
rz(1.1654677) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5844299) q[1];
sx q[1];
rz(-1.129377) q[1];
sx q[1];
rz(-0.95870383) q[1];
rz(-pi) q[2];
rz(1.2258271) q[3];
sx q[3];
rz(-2.0668092) q[3];
sx q[3];
rz(-1.8535875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7463688) q[2];
sx q[2];
rz(-0.93337983) q[2];
sx q[2];
rz(2.7044435) q[2];
rz(-2.7002667) q[3];
sx q[3];
rz(-0.91259846) q[3];
sx q[3];
rz(-0.92379409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31170347) q[0];
sx q[0];
rz(-1.7376816) q[0];
sx q[0];
rz(0.29790685) q[0];
rz(2.9333313) q[1];
sx q[1];
rz(-0.88349897) q[1];
sx q[1];
rz(0.40685245) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7283467) q[0];
sx q[0];
rz(-1.5399031) q[0];
sx q[0];
rz(1.4824162) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2595739) q[2];
sx q[2];
rz(-1.758496) q[2];
sx q[2];
rz(1.3903914) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1081097) q[1];
sx q[1];
rz(-2.2763866) q[1];
sx q[1];
rz(-0.34107855) q[1];
x q[2];
rz(-0.86866697) q[3];
sx q[3];
rz(-0.6955117) q[3];
sx q[3];
rz(1.2630386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.75140658) q[2];
sx q[2];
rz(-2.4441256) q[2];
sx q[2];
rz(2.3263399) q[2];
rz(-0.32663545) q[3];
sx q[3];
rz(-1.1465466) q[3];
sx q[3];
rz(-0.013464125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99454749) q[0];
sx q[0];
rz(-2.4319686) q[0];
sx q[0];
rz(-0.56353322) q[0];
rz(-1.7067478) q[1];
sx q[1];
rz(-2.1610465) q[1];
sx q[1];
rz(2.837406) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72918159) q[0];
sx q[0];
rz(-0.92386073) q[0];
sx q[0];
rz(-1.6135733) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5553914) q[2];
sx q[2];
rz(-2.1888141) q[2];
sx q[2];
rz(1.5780592) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.51830735) q[1];
sx q[1];
rz(-1.3629106) q[1];
sx q[1];
rz(2.1254366) q[1];
x q[2];
rz(1.0902152) q[3];
sx q[3];
rz(-0.46047415) q[3];
sx q[3];
rz(-2.0802896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4284511) q[2];
sx q[2];
rz(-0.17927543) q[2];
sx q[2];
rz(2.156554) q[2];
rz(1.8200412) q[3];
sx q[3];
rz(-1.9165087) q[3];
sx q[3];
rz(-0.2329181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-1.4317076) q[0];
sx q[0];
rz(-2.6490477) q[0];
sx q[0];
rz(2.7370969) q[0];
rz(0.74199039) q[1];
sx q[1];
rz(-1.208035) q[1];
sx q[1];
rz(-0.56125364) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40101926) q[0];
sx q[0];
rz(-0.74441806) q[0];
sx q[0];
rz(2.1932426) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91902952) q[2];
sx q[2];
rz(-2.5723151) q[2];
sx q[2];
rz(-2.809466) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2798279) q[1];
sx q[1];
rz(-2.7142254) q[1];
sx q[1];
rz(-2.159987) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3901479) q[3];
sx q[3];
rz(-1.2925694) q[3];
sx q[3];
rz(3.1415526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.81885091) q[2];
sx q[2];
rz(-0.35917869) q[2];
sx q[2];
rz(2.2364869) q[2];
rz(0.94308606) q[3];
sx q[3];
rz(-1.2312931) q[3];
sx q[3];
rz(-0.73956195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1484225) q[0];
sx q[0];
rz(-0.44023308) q[0];
sx q[0];
rz(-0.37047186) q[0];
rz(0.91520339) q[1];
sx q[1];
rz(-1.4280495) q[1];
sx q[1];
rz(-2.3781093) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3199056) q[0];
sx q[0];
rz(-2.3662462) q[0];
sx q[0];
rz(0.5628235) q[0];
rz(-0.47637003) q[2];
sx q[2];
rz(-0.61846501) q[2];
sx q[2];
rz(-2.9805984) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.50824805) q[1];
sx q[1];
rz(-2.2142383) q[1];
sx q[1];
rz(0.28336759) q[1];
rz(1.6053469) q[3];
sx q[3];
rz(-1.5175382) q[3];
sx q[3];
rz(1.0798423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8107599) q[2];
sx q[2];
rz(-1.5629038) q[2];
sx q[2];
rz(1.2609437) q[2];
rz(1.2817945) q[3];
sx q[3];
rz(-2.1112879) q[3];
sx q[3];
rz(-1.1873881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24487615) q[0];
sx q[0];
rz(-1.5302932) q[0];
sx q[0];
rz(1.2320919) q[0];
rz(-2.2121519) q[1];
sx q[1];
rz(-0.96439958) q[1];
sx q[1];
rz(0.31698116) q[1];
rz(3.0841903) q[2];
sx q[2];
rz(-1.897503) q[2];
sx q[2];
rz(-0.31384604) q[2];
rz(-1.2019602) q[3];
sx q[3];
rz(-1.0452253) q[3];
sx q[3];
rz(-0.54946231) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
