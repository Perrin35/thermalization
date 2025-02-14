OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.1342993) q[0];
sx q[0];
rz(8.2259895) q[0];
sx q[0];
rz(10.209957) q[0];
rz(2.5685318) q[1];
sx q[1];
rz(-0.95102024) q[1];
sx q[1];
rz(-2.5101488) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60053958) q[0];
sx q[0];
rz(-2.1990393) q[0];
sx q[0];
rz(2.598935) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2996743) q[2];
sx q[2];
rz(-2.2158884) q[2];
sx q[2];
rz(-3.0708714) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.59746233) q[1];
sx q[1];
rz(-1.1409581) q[1];
sx q[1];
rz(0.69263117) q[1];
rz(-pi) q[2];
rz(1.4716161) q[3];
sx q[3];
rz(-2.0137312) q[3];
sx q[3];
rz(-2.8888426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.50420061) q[2];
sx q[2];
rz(-1.7039508) q[2];
sx q[2];
rz(0.27153095) q[2];
rz(-0.046393752) q[3];
sx q[3];
rz(-1.4103187) q[3];
sx q[3];
rz(2.7739649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4162327) q[0];
sx q[0];
rz(-0.027149057) q[0];
sx q[0];
rz(-0.44560462) q[0];
rz(0.36000571) q[1];
sx q[1];
rz(-1.5574734) q[1];
sx q[1];
rz(0.97703591) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.925219) q[0];
sx q[0];
rz(-0.77298236) q[0];
sx q[0];
rz(2.3712914) q[0];
rz(0.72190921) q[2];
sx q[2];
rz(-2.0945661) q[2];
sx q[2];
rz(0.80654752) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8277196) q[1];
sx q[1];
rz(-1.4196175) q[1];
sx q[1];
rz(-0.97444356) q[1];
x q[2];
rz(-0.74574836) q[3];
sx q[3];
rz(-2.2848515) q[3];
sx q[3];
rz(2.8929949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.52686024) q[2];
sx q[2];
rz(-1.3440259) q[2];
sx q[2];
rz(2.1993707) q[2];
rz(-0.49643907) q[3];
sx q[3];
rz(-1.6783291) q[3];
sx q[3];
rz(-2.0739323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0838858) q[0];
sx q[0];
rz(-1.8118462) q[0];
sx q[0];
rz(-2.4853793) q[0];
rz(0.48037848) q[1];
sx q[1];
rz(-2.1800133) q[1];
sx q[1];
rz(-2.5286455) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71343173) q[0];
sx q[0];
rz(-1.1640004) q[0];
sx q[0];
rz(-3.0359546) q[0];
x q[1];
rz(0.7360779) q[2];
sx q[2];
rz(-1.5653584) q[2];
sx q[2];
rz(-2.8746282) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9631098) q[1];
sx q[1];
rz(-1.6043538) q[1];
sx q[1];
rz(1.5419649) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6399691) q[3];
sx q[3];
rz(-1.072848) q[3];
sx q[3];
rz(0.70424265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2059242) q[2];
sx q[2];
rz(-2.084145) q[2];
sx q[2];
rz(-2.617188) q[2];
rz(-0.86535254) q[3];
sx q[3];
rz(-1.0582346) q[3];
sx q[3];
rz(2.481148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48786369) q[0];
sx q[0];
rz(-1.3069343) q[0];
sx q[0];
rz(3.0453239) q[0];
rz(0.013821566) q[1];
sx q[1];
rz(-0.50058573) q[1];
sx q[1];
rz(-2.0892508) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5139363) q[0];
sx q[0];
rz(-0.13051438) q[0];
sx q[0];
rz(1.2016199) q[0];
rz(-2.1209893) q[2];
sx q[2];
rz(-1.3250704) q[2];
sx q[2];
rz(-0.042681132) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9279729) q[1];
sx q[1];
rz(-1.3885522) q[1];
sx q[1];
rz(-0.33267633) q[1];
rz(3.105794) q[3];
sx q[3];
rz(-2.0960663) q[3];
sx q[3];
rz(-0.73204277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.72157613) q[2];
sx q[2];
rz(-2.6105011) q[2];
sx q[2];
rz(-0.51900035) q[2];
rz(-1.0985993) q[3];
sx q[3];
rz(-1.4562621) q[3];
sx q[3];
rz(-2.4115653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62243432) q[0];
sx q[0];
rz(-2.9485478) q[0];
sx q[0];
rz(2.4073507) q[0];
rz(1.2892067) q[1];
sx q[1];
rz(-1.0451008) q[1];
sx q[1];
rz(0.9085013) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35373792) q[0];
sx q[0];
rz(-2.2252924) q[0];
sx q[0];
rz(2.2931173) q[0];
x q[1];
rz(-2.5703326) q[2];
sx q[2];
rz(-1.4584342) q[2];
sx q[2];
rz(1.686357) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.31181609) q[1];
sx q[1];
rz(-0.83925345) q[1];
sx q[1];
rz(0.093454425) q[1];
rz(3.1172006) q[3];
sx q[3];
rz(-1.304226) q[3];
sx q[3];
rz(-2.8238058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0970519) q[2];
sx q[2];
rz(-2.3553607) q[2];
sx q[2];
rz(-2.457974) q[2];
rz(2.0475856) q[3];
sx q[3];
rz(-1.8180327) q[3];
sx q[3];
rz(1.5692284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0878658) q[0];
sx q[0];
rz(-1.448486) q[0];
sx q[0];
rz(0.42561612) q[0];
rz(-0.60676891) q[1];
sx q[1];
rz(-2.4538222) q[1];
sx q[1];
rz(1.2026131) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86798651) q[0];
sx q[0];
rz(-2.3487954) q[0];
sx q[0];
rz(-2.2895564) q[0];
rz(-pi) q[1];
rz(0.065826611) q[2];
sx q[2];
rz(-1.5853241) q[2];
sx q[2];
rz(0.17958497) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9252117) q[1];
sx q[1];
rz(-0.93384508) q[1];
sx q[1];
rz(2.2014307) q[1];
x q[2];
rz(2.3601861) q[3];
sx q[3];
rz(-2.8202116) q[3];
sx q[3];
rz(-1.7226294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.46835029) q[2];
sx q[2];
rz(-2.612256) q[2];
sx q[2];
rz(-2.5128515) q[2];
rz(0.33744234) q[3];
sx q[3];
rz(-1.7818261) q[3];
sx q[3];
rz(0.76702816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5479946) q[0];
sx q[0];
rz(-0.070040528) q[0];
sx q[0];
rz(1.8016169) q[0];
rz(-3.0361259) q[1];
sx q[1];
rz(-0.98140812) q[1];
sx q[1];
rz(0.14428219) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6010885) q[0];
sx q[0];
rz(-1.861689) q[0];
sx q[0];
rz(-0.40368794) q[0];
rz(0.86267306) q[2];
sx q[2];
rz(-1.4221592) q[2];
sx q[2];
rz(-1.3512063) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.49715091) q[1];
sx q[1];
rz(-1.9195917) q[1];
sx q[1];
rz(-1.7921053) q[1];
x q[2];
rz(2.4915817) q[3];
sx q[3];
rz(-2.3571005) q[3];
sx q[3];
rz(-0.9264526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0353126) q[2];
sx q[2];
rz(-1.148618) q[2];
sx q[2];
rz(-0.8379575) q[2];
rz(-0.49977866) q[3];
sx q[3];
rz(-2.3482447) q[3];
sx q[3];
rz(0.70118538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4226294) q[0];
sx q[0];
rz(-1.1088277) q[0];
sx q[0];
rz(2.3903589) q[0];
rz(1.2535837) q[1];
sx q[1];
rz(-1.3026078) q[1];
sx q[1];
rz(-1.7734897) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7984895) q[0];
sx q[0];
rz(-2.0105848) q[0];
sx q[0];
rz(-0.90476157) q[0];
x q[1];
rz(0.74636942) q[2];
sx q[2];
rz(-0.41902768) q[2];
sx q[2];
rz(2.7460436) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.83225314) q[1];
sx q[1];
rz(-2.3414438) q[1];
sx q[1];
rz(-3.0235428) q[1];
rz(-pi) q[2];
rz(0.99686025) q[3];
sx q[3];
rz(-1.2632476) q[3];
sx q[3];
rz(2.779863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3086243) q[2];
sx q[2];
rz(-2.7699351) q[2];
sx q[2];
rz(-0.0011681636) q[2];
rz(-0.33231169) q[3];
sx q[3];
rz(-1.6806335) q[3];
sx q[3];
rz(-0.47962475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7697656) q[0];
sx q[0];
rz(-1.9244939) q[0];
sx q[0];
rz(-1.0960854) q[0];
rz(-1.3990043) q[1];
sx q[1];
rz(-1.5053446) q[1];
sx q[1];
rz(-2.2217264) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8815589) q[0];
sx q[0];
rz(-1.1333916) q[0];
sx q[0];
rz(0.046925515) q[0];
rz(-1.1521173) q[2];
sx q[2];
rz(-1.8279462) q[2];
sx q[2];
rz(-0.60794431) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.80650027) q[1];
sx q[1];
rz(-0.92969172) q[1];
sx q[1];
rz(-1.758746) q[1];
rz(0.14443951) q[3];
sx q[3];
rz(-2.0990934) q[3];
sx q[3];
rz(-2.4424255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2190735) q[2];
sx q[2];
rz(-0.24524958) q[2];
sx q[2];
rz(-1.5101439) q[2];
rz(-0.63589969) q[3];
sx q[3];
rz(-2.4363775) q[3];
sx q[3];
rz(0.51457921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9566327) q[0];
sx q[0];
rz(-0.867046) q[0];
sx q[0];
rz(-1.6459203) q[0];
rz(3.0236859) q[1];
sx q[1];
rz(-1.3943744) q[1];
sx q[1];
rz(2.8848021) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3366616) q[0];
sx q[0];
rz(-1.1380702) q[0];
sx q[0];
rz(3.0800729) q[0];
x q[1];
rz(0.43440993) q[2];
sx q[2];
rz(-0.63077482) q[2];
sx q[2];
rz(-1.8112911) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6424055) q[1];
sx q[1];
rz(-2.7491264) q[1];
sx q[1];
rz(2.3981903) q[1];
rz(-pi) q[2];
x q[2];
rz(0.55744265) q[3];
sx q[3];
rz(-2.5187105) q[3];
sx q[3];
rz(2.0940144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.39766496) q[2];
sx q[2];
rz(-2.3981018) q[2];
sx q[2];
rz(0.94997722) q[2];
rz(-2.5395565) q[3];
sx q[3];
rz(-1.9046015) q[3];
sx q[3];
rz(-0.33815798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6947185) q[0];
sx q[0];
rz(-0.61114408) q[0];
sx q[0];
rz(-1.1993988) q[0];
rz(-2.1318204) q[1];
sx q[1];
rz(-0.92242454) q[1];
sx q[1];
rz(-0.19933272) q[1];
rz(-2.845191) q[2];
sx q[2];
rz(-1.2747073) q[2];
sx q[2];
rz(-2.7327197) q[2];
rz(-2.7446741) q[3];
sx q[3];
rz(-0.81868681) q[3];
sx q[3];
rz(0.80358748) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
