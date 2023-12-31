OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5126098) q[0];
sx q[0];
rz(-1.0273758) q[0];
sx q[0];
rz(-2.7789814) q[0];
rz(0.91712046) q[1];
sx q[1];
rz(2.6511104) q[1];
sx q[1];
rz(9.766415) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4877019) q[0];
sx q[0];
rz(-2.747295) q[0];
sx q[0];
rz(-2.8173692) q[0];
rz(-pi) q[1];
rz(2.5506637) q[2];
sx q[2];
rz(-1.2714296) q[2];
sx q[2];
rz(-2.9654944) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8056148) q[1];
sx q[1];
rz(-1.059638) q[1];
sx q[1];
rz(2.7989945) q[1];
rz(-pi) q[2];
rz(-2.9176641) q[3];
sx q[3];
rz(-1.180598) q[3];
sx q[3];
rz(1.8412561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5477156) q[2];
sx q[2];
rz(-1.9888069) q[2];
sx q[2];
rz(0.95648742) q[2];
rz(2.9603413) q[3];
sx q[3];
rz(-2.5528208) q[3];
sx q[3];
rz(-1.5216924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0797121) q[0];
sx q[0];
rz(-0.065608874) q[0];
sx q[0];
rz(1.1454426) q[0];
rz(2.0727797) q[1];
sx q[1];
rz(-0.83199465) q[1];
sx q[1];
rz(-0.72584814) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94978588) q[0];
sx q[0];
rz(-2.1161785) q[0];
sx q[0];
rz(1.1508862) q[0];
rz(-pi) q[1];
x q[1];
rz(0.59402324) q[2];
sx q[2];
rz(-1.1813287) q[2];
sx q[2];
rz(0.71690744) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5729546) q[1];
sx q[1];
rz(-1.8796088) q[1];
sx q[1];
rz(-0.026334865) q[1];
rz(1.782136) q[3];
sx q[3];
rz(-2.1826715) q[3];
sx q[3];
rz(2.8733727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6590865) q[2];
sx q[2];
rz(-1.3799474) q[2];
sx q[2];
rz(0.99728161) q[2];
rz(1.4128134) q[3];
sx q[3];
rz(-2.0480859) q[3];
sx q[3];
rz(-0.93878448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41855758) q[0];
sx q[0];
rz(-0.96590531) q[0];
sx q[0];
rz(1.1285271) q[0];
rz(1.8380802) q[1];
sx q[1];
rz(-2.4120657) q[1];
sx q[1];
rz(-0.9054786) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0008464) q[0];
sx q[0];
rz(-1.8119537) q[0];
sx q[0];
rz(1.7667889) q[0];
x q[1];
rz(-0.75240527) q[2];
sx q[2];
rz(-1.2453326) q[2];
sx q[2];
rz(-0.67093713) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3382197) q[1];
sx q[1];
rz(-1.6446911) q[1];
sx q[1];
rz(-0.85892962) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.043025322) q[3];
sx q[3];
rz(-2.283841) q[3];
sx q[3];
rz(0.42792861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.7604312) q[2];
sx q[2];
rz(-0.16813777) q[2];
sx q[2];
rz(-0.9872438) q[2];
rz(-0.045771249) q[3];
sx q[3];
rz(-1.2922623) q[3];
sx q[3];
rz(0.18019095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0320597) q[0];
sx q[0];
rz(-2.4592472) q[0];
sx q[0];
rz(0.12606829) q[0];
rz(2.9571422) q[1];
sx q[1];
rz(-1.0686864) q[1];
sx q[1];
rz(-1.9741845) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6895034) q[0];
sx q[0];
rz(-2.3514682) q[0];
sx q[0];
rz(2.8566314) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.5037751) q[2];
sx q[2];
rz(-1.7283895) q[2];
sx q[2];
rz(2.4368844) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3525225) q[1];
sx q[1];
rz(-1.3067424) q[1];
sx q[1];
rz(-0.67184429) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.64658029) q[3];
sx q[3];
rz(-1.0160035) q[3];
sx q[3];
rz(-3.0559029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2771153) q[2];
sx q[2];
rz(-0.41431674) q[2];
sx q[2];
rz(0.63151044) q[2];
rz(0.19550066) q[3];
sx q[3];
rz(-1.2771527) q[3];
sx q[3];
rz(-2.7235532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
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
rz(1.1642078) q[0];
sx q[0];
rz(-3.003484) q[0];
sx q[0];
rz(-2.6389129) q[0];
rz(2.0282822) q[1];
sx q[1];
rz(-2.138442) q[1];
sx q[1];
rz(-2.6745093) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53411667) q[0];
sx q[0];
rz(-0.42629888) q[0];
sx q[0];
rz(2.0712453) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5555243) q[2];
sx q[2];
rz(-1.35891) q[2];
sx q[2];
rz(-0.52473247) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.621532) q[1];
sx q[1];
rz(-0.48081765) q[1];
sx q[1];
rz(-2.0562999) q[1];
x q[2];
rz(2.2561982) q[3];
sx q[3];
rz(-1.6939031) q[3];
sx q[3];
rz(-0.76243329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6749394) q[2];
sx q[2];
rz(-1.8920687) q[2];
sx q[2];
rz(0.49218991) q[2];
rz(2.8594033) q[3];
sx q[3];
rz(-1.7084833) q[3];
sx q[3];
rz(1.5757489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(2.2714587) q[0];
sx q[0];
rz(-2.6014355) q[0];
sx q[0];
rz(-3.0555994) q[0];
rz(1.4280691) q[1];
sx q[1];
rz(-2.1336887) q[1];
sx q[1];
rz(1.6451947) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2841227) q[0];
sx q[0];
rz(-2.5792482) q[0];
sx q[0];
rz(-1.9447295) q[0];
rz(-1.1804579) q[2];
sx q[2];
rz(-0.57363735) q[2];
sx q[2];
rz(0.2804642) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8138258) q[1];
sx q[1];
rz(-1.5944949) q[1];
sx q[1];
rz(-1.9253255) q[1];
x q[2];
rz(-0.014724894) q[3];
sx q[3];
rz(-2.428599) q[3];
sx q[3];
rz(-1.7409489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4795586) q[2];
sx q[2];
rz(-2.5223314) q[2];
sx q[2];
rz(-0.39624828) q[2];
rz(-2.5475492) q[3];
sx q[3];
rz(-1.0915353) q[3];
sx q[3];
rz(-1.5007639) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2632161) q[0];
sx q[0];
rz(-2.5686869) q[0];
sx q[0];
rz(-2.0492045) q[0];
rz(1.3573525) q[1];
sx q[1];
rz(-1.4211632) q[1];
sx q[1];
rz(3.1299652) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7250925) q[0];
sx q[0];
rz(-1.9511127) q[0];
sx q[0];
rz(-1.2854544) q[0];
x q[1];
rz(1.5426703) q[2];
sx q[2];
rz(-1.5962432) q[2];
sx q[2];
rz(-0.59288073) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7203622) q[1];
sx q[1];
rz(-2.1167397) q[1];
sx q[1];
rz(3.1139042) q[1];
rz(0.67354789) q[3];
sx q[3];
rz(-1.3140956) q[3];
sx q[3];
rz(1.7208769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6860883) q[2];
sx q[2];
rz(-0.65767613) q[2];
sx q[2];
rz(0.28798506) q[2];
rz(-1.9912432) q[3];
sx q[3];
rz(-1.3214313) q[3];
sx q[3];
rz(-2.5434125) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8312254) q[0];
sx q[0];
rz(-2.6103525) q[0];
sx q[0];
rz(1.9294552) q[0];
rz(2.7638226) q[1];
sx q[1];
rz(-1.2021844) q[1];
sx q[1];
rz(1.4935965) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9858915) q[0];
sx q[0];
rz(-2.771286) q[0];
sx q[0];
rz(-2.5168688) q[0];
rz(-1.9940894) q[2];
sx q[2];
rz(-0.14483843) q[2];
sx q[2];
rz(-0.3651948) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8131866) q[1];
sx q[1];
rz(-1.9771736) q[1];
sx q[1];
rz(2.3996668) q[1];
rz(-pi) q[2];
rz(1.3573523) q[3];
sx q[3];
rz(-0.5118013) q[3];
sx q[3];
rz(1.7970999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7887855) q[2];
sx q[2];
rz(-1.1151423) q[2];
sx q[2];
rz(-1.2517694) q[2];
rz(2.2552323) q[3];
sx q[3];
rz(-2.3754407) q[3];
sx q[3];
rz(-0.10929508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66824526) q[0];
sx q[0];
rz(-2.1003523) q[0];
sx q[0];
rz(2.9633203) q[0];
rz(3.0687304) q[1];
sx q[1];
rz(-0.59385121) q[1];
sx q[1];
rz(-1.9624306) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3423791) q[0];
sx q[0];
rz(-1.2782492) q[0];
sx q[0];
rz(-0.53036687) q[0];
rz(-pi) q[1];
rz(-0.91782848) q[2];
sx q[2];
rz(-2.4036916) q[2];
sx q[2];
rz(0.81673056) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0852016) q[1];
sx q[1];
rz(-1.3136275) q[1];
sx q[1];
rz(-2.609054) q[1];
rz(-2.7618802) q[3];
sx q[3];
rz(-1.100193) q[3];
sx q[3];
rz(-3.0203366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.69958413) q[2];
sx q[2];
rz(-0.99094355) q[2];
sx q[2];
rz(-2.1098095) q[2];
rz(1.7992841) q[3];
sx q[3];
rz(-1.4167901) q[3];
sx q[3];
rz(2.9163196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.853302) q[0];
sx q[0];
rz(-2.0993711) q[0];
sx q[0];
rz(0.061696079) q[0];
rz(1.8348947) q[1];
sx q[1];
rz(-1.6625762) q[1];
sx q[1];
rz(1.0888938) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13578829) q[0];
sx q[0];
rz(-0.09011589) q[0];
sx q[0];
rz(-2.3401005) q[0];
x q[1];
rz(1.8082952) q[2];
sx q[2];
rz(-2.4727159) q[2];
sx q[2];
rz(-2.3046658) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2838193) q[1];
sx q[1];
rz(-0.76446629) q[1];
sx q[1];
rz(-0.70160265) q[1];
rz(-pi) q[2];
rz(2.0652566) q[3];
sx q[3];
rz(-0.38443243) q[3];
sx q[3];
rz(-1.2151375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5237727) q[2];
sx q[2];
rz(-2.4536295) q[2];
sx q[2];
rz(3.0715959) q[2];
rz(2.2648515) q[3];
sx q[3];
rz(-1.8165959) q[3];
sx q[3];
rz(-2.783412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-0.6488279) q[0];
sx q[0];
rz(-1.6236826) q[0];
sx q[0];
rz(0.59326011) q[0];
rz(-1.6246673) q[1];
sx q[1];
rz(-2.0943874) q[1];
sx q[1];
rz(2.692683) q[1];
rz(1.0629874) q[2];
sx q[2];
rz(-2.4175736) q[2];
sx q[2];
rz(-2.8833817) q[2];
rz(-1.2354479) q[3];
sx q[3];
rz(-1.609758) q[3];
sx q[3];
rz(2.5555425) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
