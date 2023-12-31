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
rz(-2.2244722) q[1];
sx q[1];
rz(-2.6511104) q[1];
sx q[1];
rz(-2.7999556) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61603123) q[0];
sx q[0];
rz(-1.6934868) q[0];
sx q[0];
rz(-2.7659155) q[0];
rz(-pi) q[1];
rz(0.59092893) q[2];
sx q[2];
rz(-1.2714296) q[2];
sx q[2];
rz(2.9654944) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8056148) q[1];
sx q[1];
rz(-1.059638) q[1];
sx q[1];
rz(0.34259818) q[1];
x q[2];
rz(1.9699691) q[3];
sx q[3];
rz(-1.7776383) q[3];
sx q[3];
rz(-2.9575461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.59387702) q[2];
sx q[2];
rz(-1.9888069) q[2];
sx q[2];
rz(-2.1851052) q[2];
rz(-0.18125136) q[3];
sx q[3];
rz(-2.5528208) q[3];
sx q[3];
rz(-1.5216924) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0618806) q[0];
sx q[0];
rz(-0.065608874) q[0];
sx q[0];
rz(-1.1454426) q[0];
rz(2.0727797) q[1];
sx q[1];
rz(-2.309598) q[1];
sx q[1];
rz(0.72584814) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9024591) q[0];
sx q[0];
rz(-0.67507889) q[0];
sx q[0];
rz(2.549987) q[0];
x q[1];
rz(0.59402324) q[2];
sx q[2];
rz(-1.1813287) q[2];
sx q[2];
rz(0.71690744) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5686381) q[1];
sx q[1];
rz(-1.8796088) q[1];
sx q[1];
rz(-3.1152578) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5190982) q[3];
sx q[3];
rz(-1.7433634) q[3];
sx q[3];
rz(-1.7164001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.48250616) q[2];
sx q[2];
rz(-1.7616452) q[2];
sx q[2];
rz(-2.144311) q[2];
rz(-1.7287792) q[3];
sx q[3];
rz(-2.0480859) q[3];
sx q[3];
rz(2.2028082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.41855758) q[0];
sx q[0];
rz(-0.96590531) q[0];
sx q[0];
rz(-1.1285271) q[0];
rz(-1.3035125) q[1];
sx q[1];
rz(-0.729527) q[1];
sx q[1];
rz(0.9054786) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3826686) q[0];
sx q[0];
rz(-1.3805458) q[0];
sx q[0];
rz(0.24567901) q[0];
rz(-pi) q[1];
rz(2.0037903) q[2];
sx q[2];
rz(-2.275122) q[2];
sx q[2];
rz(-1.1906884) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3382197) q[1];
sx q[1];
rz(-1.4969016) q[1];
sx q[1];
rz(-0.85892962) q[1];
x q[2];
rz(-0.043025322) q[3];
sx q[3];
rz(-2.283841) q[3];
sx q[3];
rz(-2.713664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.7604312) q[2];
sx q[2];
rz(-0.16813777) q[2];
sx q[2];
rz(-2.1543489) q[2];
rz(3.0958214) q[3];
sx q[3];
rz(-1.8493303) q[3];
sx q[3];
rz(-0.18019095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10953294) q[0];
sx q[0];
rz(-0.68234545) q[0];
sx q[0];
rz(3.0155244) q[0];
rz(2.9571422) q[1];
sx q[1];
rz(-1.0686864) q[1];
sx q[1];
rz(1.1674081) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6895034) q[0];
sx q[0];
rz(-2.3514682) q[0];
sx q[0];
rz(-2.8566314) q[0];
x q[1];
rz(0.5037751) q[2];
sx q[2];
rz(-1.7283895) q[2];
sx q[2];
rz(0.7047082) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5645204) q[1];
sx q[1];
rz(-2.2153691) q[1];
sx q[1];
rz(1.9034027) q[1];
x q[2];
rz(0.79951841) q[3];
sx q[3];
rz(-0.82516731) q[3];
sx q[3];
rz(-1.0471917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.8644774) q[2];
sx q[2];
rz(-2.7272759) q[2];
sx q[2];
rz(2.5100822) q[2];
rz(-0.19550066) q[3];
sx q[3];
rz(-1.2771527) q[3];
sx q[3];
rz(2.7235532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9773848) q[0];
sx q[0];
rz(-0.13810869) q[0];
sx q[0];
rz(2.6389129) q[0];
rz(-2.0282822) q[1];
sx q[1];
rz(-1.0031507) q[1];
sx q[1];
rz(-2.6745093) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5669117) q[0];
sx q[0];
rz(-1.7705288) q[0];
sx q[0];
rz(-1.1916222) q[0];
rz(2.7706625) q[2];
sx q[2];
rz(-0.61894722) q[2];
sx q[2];
rz(-1.3528454) q[2];
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
rz(-pi) q[2];
rz(-1.3777556) q[3];
sx q[3];
rz(-0.69460624) q[3];
sx q[3];
rz(2.4823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6749394) q[2];
sx q[2];
rz(-1.249524) q[2];
sx q[2];
rz(0.49218991) q[2];
rz(-2.8594033) q[3];
sx q[3];
rz(-1.4331093) q[3];
sx q[3];
rz(1.5757489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.870134) q[0];
sx q[0];
rz(-0.54015714) q[0];
sx q[0];
rz(0.085993275) q[0];
rz(1.7135235) q[1];
sx q[1];
rz(-1.0079039) q[1];
sx q[1];
rz(1.6451947) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8498358) q[0];
sx q[0];
rz(-1.0514326) q[0];
sx q[0];
rz(-0.22626466) q[0];
rz(-2.9005269) q[2];
sx q[2];
rz(-1.0450372) q[2];
sx q[2];
rz(-2.966553) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.23425809) q[1];
sx q[1];
rz(-1.9252216) q[1];
sx q[1];
rz(-0.025269421) q[1];
rz(-pi) q[2];
rz(2.4286527) q[3];
sx q[3];
rz(-1.561165) q[3];
sx q[3];
rz(-2.9603017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4795586) q[2];
sx q[2];
rz(-2.5223314) q[2];
sx q[2];
rz(-2.7453444) q[2];
rz(-0.59404343) q[3];
sx q[3];
rz(-2.0500573) q[3];
sx q[3];
rz(1.6408287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2632161) q[0];
sx q[0];
rz(-0.57290572) q[0];
sx q[0];
rz(2.0492045) q[0];
rz(-1.3573525) q[1];
sx q[1];
rz(-1.4211632) q[1];
sx q[1];
rz(-3.1299652) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8893338) q[0];
sx q[0];
rz(-0.47124915) q[0];
sx q[0];
rz(-2.528119) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1161357) q[2];
sx q[2];
rz(-1.5426794) q[2];
sx q[2];
rz(-2.1643929) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9776476) q[1];
sx q[1];
rz(-1.5471336) q[1];
sx q[1];
rz(2.1169099) q[1];
rz(-1.2467975) q[3];
sx q[3];
rz(-2.218459) q[3];
sx q[3];
rz(3.0917633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.45550436) q[2];
sx q[2];
rz(-0.65767613) q[2];
sx q[2];
rz(-2.8536076) q[2];
rz(1.1503495) q[3];
sx q[3];
rz(-1.3214313) q[3];
sx q[3];
rz(0.59818017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31036723) q[0];
sx q[0];
rz(-0.53124017) q[0];
sx q[0];
rz(1.2121375) q[0];
rz(2.7638226) q[1];
sx q[1];
rz(-1.2021844) q[1];
sx q[1];
rz(-1.6479962) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0069107) q[0];
sx q[0];
rz(-1.3575166) q[0];
sx q[0];
rz(-2.8365305) q[0];
rz(0.059842589) q[2];
sx q[2];
rz(-1.4388196) q[2];
sx q[2];
rz(0.062054141) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.89477506) q[1];
sx q[1];
rz(-2.2404039) q[1];
sx q[1];
rz(-1.0424022) q[1];
rz(-pi) q[2];
rz(-3.0231608) q[3];
sx q[3];
rz(-1.0716972) q[3];
sx q[3];
rz(-1.1008319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7887855) q[2];
sx q[2];
rz(-2.0264503) q[2];
sx q[2];
rz(1.2517694) q[2];
rz(-0.88636032) q[3];
sx q[3];
rz(-2.3754407) q[3];
sx q[3];
rz(3.0322976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4733474) q[0];
sx q[0];
rz(-1.0412403) q[0];
sx q[0];
rz(-2.9633203) q[0];
rz(3.0687304) q[1];
sx q[1];
rz(-2.5477414) q[1];
sx q[1];
rz(-1.1791621) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68547738) q[0];
sx q[0];
rz(-0.59887409) q[0];
sx q[0];
rz(-0.53703888) q[0];
rz(0.945325) q[2];
sx q[2];
rz(-1.1497467) q[2];
sx q[2];
rz(2.9025214) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0563911) q[1];
sx q[1];
rz(-1.3136275) q[1];
sx q[1];
rz(-0.53253865) q[1];
rz(1.0696899) q[3];
sx q[3];
rz(-1.9074829) q[3];
sx q[3];
rz(1.6285553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.69958413) q[2];
sx q[2];
rz(-0.99094355) q[2];
sx q[2];
rz(2.1098095) q[2];
rz(1.3423086) q[3];
sx q[3];
rz(-1.4167901) q[3];
sx q[3];
rz(0.22527307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28829065) q[0];
sx q[0];
rz(-1.0422215) q[0];
sx q[0];
rz(3.0798966) q[0];
rz(-1.8348947) q[1];
sx q[1];
rz(-1.4790165) q[1];
sx q[1];
rz(1.0888938) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63554791) q[0];
sx q[0];
rz(-1.6354927) q[0];
sx q[0];
rz(-3.0788172) q[0];
x q[1];
rz(0.18386545) q[2];
sx q[2];
rz(-2.2176748) q[2];
sx q[2];
rz(-1.136214) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1354462) q[1];
sx q[1];
rz(-1.0137614) q[1];
sx q[1];
rz(1.0165434) q[1];
rz(0.18967929) q[3];
sx q[3];
rz(-1.2343725) q[3];
sx q[3];
rz(1.7419025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5237727) q[2];
sx q[2];
rz(-2.4536295) q[2];
sx q[2];
rz(3.0715959) q[2];
rz(-2.2648515) q[3];
sx q[3];
rz(-1.3249967) q[3];
sx q[3];
rz(-2.783412) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4927647) q[0];
sx q[0];
rz(-1.6236826) q[0];
sx q[0];
rz(0.59326011) q[0];
rz(-1.5169253) q[1];
sx q[1];
rz(-1.0472052) q[1];
sx q[1];
rz(-0.44890961) q[1];
rz(2.0786053) q[2];
sx q[2];
rz(-0.72401902) q[2];
sx q[2];
rz(0.25821092) q[2];
rz(-3.1003351) q[3];
sx q[3];
rz(-1.2357124) q[3];
sx q[3];
rz(0.99832051) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
