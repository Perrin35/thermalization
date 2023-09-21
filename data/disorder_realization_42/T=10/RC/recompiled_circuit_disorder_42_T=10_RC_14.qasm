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
rz(-0.4904823) q[1];
sx q[1];
rz(-0.34163707) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0029966) q[0];
sx q[0];
rz(-1.1980822) q[0];
sx q[0];
rz(-1.7025823) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5506637) q[2];
sx q[2];
rz(-1.2714296) q[2];
sx q[2];
rz(2.9654944) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8475854) q[1];
sx q[1];
rz(-0.60677401) q[1];
sx q[1];
rz(-2.1104382) q[1];
rz(-pi) q[2];
rz(2.9176641) q[3];
sx q[3];
rz(-1.180598) q[3];
sx q[3];
rz(1.3003365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5477156) q[2];
sx q[2];
rz(-1.1527858) q[2];
sx q[2];
rz(2.1851052) q[2];
rz(-2.9603413) q[3];
sx q[3];
rz(-2.5528208) q[3];
sx q[3];
rz(1.5216924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.0618806) q[0];
sx q[0];
rz(-0.065608874) q[0];
sx q[0];
rz(-1.99615) q[0];
rz(-1.068813) q[1];
sx q[1];
rz(-2.309598) q[1];
sx q[1];
rz(0.72584814) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9024591) q[0];
sx q[0];
rz(-2.4665138) q[0];
sx q[0];
rz(-0.59160561) q[0];
rz(-0.59402324) q[2];
sx q[2];
rz(-1.1813287) q[2];
sx q[2];
rz(2.4246852) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1314288) q[1];
sx q[1];
rz(-1.5457075) q[1];
sx q[1];
rz(-1.2618834) q[1];
x q[2];
rz(-1.3594567) q[3];
sx q[3];
rz(-0.95892116) q[3];
sx q[3];
rz(0.26821995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.48250616) q[2];
sx q[2];
rz(-1.3799474) q[2];
sx q[2];
rz(-0.99728161) q[2];
rz(-1.4128134) q[3];
sx q[3];
rz(-2.0480859) q[3];
sx q[3];
rz(0.93878448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7230351) q[0];
sx q[0];
rz(-2.1756873) q[0];
sx q[0];
rz(2.0130656) q[0];
rz(1.3035125) q[1];
sx q[1];
rz(-2.4120657) q[1];
sx q[1];
rz(-2.2361141) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.758924) q[0];
sx q[0];
rz(-1.7610468) q[0];
sx q[0];
rz(2.8959136) q[0];
rz(-0.75240527) q[2];
sx q[2];
rz(-1.2453326) q[2];
sx q[2];
rz(-0.67093713) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.29618759) q[1];
sx q[1];
rz(-0.86128174) q[1];
sx q[1];
rz(-3.0441277) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.85729349) q[3];
sx q[3];
rz(-1.5382574) q[3];
sx q[3];
rz(-1.9705704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3811615) q[2];
sx q[2];
rz(-2.9734549) q[2];
sx q[2];
rz(-0.9872438) q[2];
rz(-3.0958214) q[3];
sx q[3];
rz(-1.2922623) q[3];
sx q[3];
rz(-0.18019095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10953294) q[0];
sx q[0];
rz(-0.68234545) q[0];
sx q[0];
rz(3.0155244) q[0];
rz(0.18445045) q[1];
sx q[1];
rz(-1.0686864) q[1];
sx q[1];
rz(-1.1674081) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2950738) q[0];
sx q[0];
rz(-2.3210038) q[0];
sx q[0];
rz(1.8473162) q[0];
rz(-pi) q[1];
x q[1];
rz(0.5037751) q[2];
sx q[2];
rz(-1.7283895) q[2];
sx q[2];
rz(-2.4368844) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5645204) q[1];
sx q[1];
rz(-0.92622354) q[1];
sx q[1];
rz(-1.2381899) q[1];
rz(-pi) q[2];
rz(2.231009) q[3];
sx q[3];
rz(-1.0331717) q[3];
sx q[3];
rz(-1.1066574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2771153) q[2];
sx q[2];
rz(-0.41431674) q[2];
sx q[2];
rz(-2.5100822) q[2];
rz(-2.946092) q[3];
sx q[3];
rz(-1.86444) q[3];
sx q[3];
rz(-0.41803944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1642078) q[0];
sx q[0];
rz(-0.13810869) q[0];
sx q[0];
rz(-0.50267977) q[0];
rz(-2.0282822) q[1];
sx q[1];
rz(-2.138442) q[1];
sx q[1];
rz(-0.46708333) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53411667) q[0];
sx q[0];
rz(-2.7152938) q[0];
sx q[0];
rz(1.0703474) q[0];
x q[1];
rz(-1.8234812) q[2];
sx q[2];
rz(-0.99950302) q[2];
sx q[2];
rz(0.90734446) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.621532) q[1];
sx q[1];
rz(-2.660775) q[1];
sx q[1];
rz(-2.0562999) q[1];
x q[2];
rz(0.15848666) q[3];
sx q[3];
rz(-0.89156686) q[3];
sx q[3];
rz(0.90843006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.46665329) q[2];
sx q[2];
rz(-1.8920687) q[2];
sx q[2];
rz(0.49218991) q[2];
rz(2.8594033) q[3];
sx q[3];
rz(-1.7084833) q[3];
sx q[3];
rz(-1.5658437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(2.2714587) q[0];
sx q[0];
rz(-0.54015714) q[0];
sx q[0];
rz(3.0555994) q[0];
rz(1.7135235) q[1];
sx q[1];
rz(-2.1336887) q[1];
sx q[1];
rz(-1.6451947) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2917568) q[0];
sx q[0];
rz(-2.09016) q[0];
sx q[0];
rz(-2.915328) q[0];
rz(-pi) q[1];
rz(0.24106579) q[2];
sx q[2];
rz(-2.0965555) q[2];
sx q[2];
rz(-0.17503967) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8346356) q[1];
sx q[1];
rz(-0.35528696) q[1];
sx q[1];
rz(1.5026232) q[1];
x q[2];
rz(0.014724894) q[3];
sx q[3];
rz(-0.71299362) q[3];
sx q[3];
rz(-1.7409489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.66203403) q[2];
sx q[2];
rz(-0.61926121) q[2];
sx q[2];
rz(-2.7453444) q[2];
rz(-2.5475492) q[3];
sx q[3];
rz(-2.0500573) q[3];
sx q[3];
rz(-1.6408287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2632161) q[0];
sx q[0];
rz(-2.5686869) q[0];
sx q[0];
rz(-1.0923882) q[0];
rz(1.7842402) q[1];
sx q[1];
rz(-1.7204294) q[1];
sx q[1];
rz(3.1299652) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8788293) q[0];
sx q[0];
rz(-1.3063523) q[0];
sx q[0];
rz(-0.39475616) q[0];
rz(-pi) q[1];
rz(-2.3064012) q[2];
sx q[2];
rz(-0.037926849) q[2];
sx q[2];
rz(1.4284301) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9776476) q[1];
sx q[1];
rz(-1.5944591) q[1];
sx q[1];
rz(2.1169099) q[1];
rz(-pi) q[2];
x q[2];
rz(2.743268) q[3];
sx q[3];
rz(-2.4279804) q[3];
sx q[3];
rz(0.45809612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.45550436) q[2];
sx q[2];
rz(-2.4839165) q[2];
sx q[2];
rz(-0.28798506) q[2];
rz(1.1503495) q[3];
sx q[3];
rz(-1.8201613) q[3];
sx q[3];
rz(-0.59818017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31036723) q[0];
sx q[0];
rz(-0.53124017) q[0];
sx q[0];
rz(-1.2121375) q[0];
rz(2.7638226) q[1];
sx q[1];
rz(-1.9394082) q[1];
sx q[1];
rz(-1.4935965) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0069107) q[0];
sx q[0];
rz(-1.784076) q[0];
sx q[0];
rz(2.8365305) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9940894) q[2];
sx q[2];
rz(-0.14483843) q[2];
sx q[2];
rz(-0.3651948) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8131866) q[1];
sx q[1];
rz(-1.1644191) q[1];
sx q[1];
rz(-2.3996668) q[1];
rz(-pi) q[2];
rz(-2.0728552) q[3];
sx q[3];
rz(-1.6747253) q[3];
sx q[3];
rz(2.7285189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7887855) q[2];
sx q[2];
rz(-1.1151423) q[2];
sx q[2];
rz(1.2517694) q[2];
rz(0.88636032) q[3];
sx q[3];
rz(-0.76615196) q[3];
sx q[3];
rz(3.0322976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66824526) q[0];
sx q[0];
rz(-1.0412403) q[0];
sx q[0];
rz(2.9633203) q[0];
rz(-3.0687304) q[1];
sx q[1];
rz(-0.59385121) q[1];
sx q[1];
rz(1.9624306) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0807063) q[0];
sx q[0];
rz(-1.0651677) q[0];
sx q[0];
rz(1.2348742) q[0];
rz(-pi) q[1];
rz(-2.1962677) q[2];
sx q[2];
rz(-1.9918459) q[2];
sx q[2];
rz(-2.9025214) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.33680962) q[1];
sx q[1];
rz(-2.0840624) q[1];
sx q[1];
rz(1.8670765) q[1];
rz(-pi) q[2];
rz(-2.7618802) q[3];
sx q[3];
rz(-2.0413997) q[3];
sx q[3];
rz(-0.1212561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4420085) q[2];
sx q[2];
rz(-2.1506491) q[2];
sx q[2];
rz(-1.0317831) q[2];
rz(-1.7992841) q[3];
sx q[3];
rz(-1.4167901) q[3];
sx q[3];
rz(0.22527307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.853302) q[0];
sx q[0];
rz(-1.0422215) q[0];
sx q[0];
rz(-0.061696079) q[0];
rz(-1.306698) q[1];
sx q[1];
rz(-1.4790165) q[1];
sx q[1];
rz(-1.0888938) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93931224) q[0];
sx q[0];
rz(-1.6334403) q[0];
sx q[0];
rz(-1.5059727) q[0];
rz(-pi) q[1];
rz(0.18386545) q[2];
sx q[2];
rz(-2.2176748) q[2];
sx q[2];
rz(-1.136214) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1354462) q[1];
sx q[1];
rz(-2.1278312) q[1];
sx q[1];
rz(2.1250493) q[1];
rz(1.0763361) q[3];
sx q[3];
rz(-2.7571602) q[3];
sx q[3];
rz(-1.2151375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5237727) q[2];
sx q[2];
rz(-2.4536295) q[2];
sx q[2];
rz(-3.0715959) q[2];
rz(2.2648515) q[3];
sx q[3];
rz(-1.8165959) q[3];
sx q[3];
rz(0.35818067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4927647) q[0];
sx q[0];
rz(-1.6236826) q[0];
sx q[0];
rz(0.59326011) q[0];
rz(1.5169253) q[1];
sx q[1];
rz(-2.0943874) q[1];
sx q[1];
rz(2.692683) q[1];
rz(-0.91296997) q[2];
sx q[2];
rz(-1.2428478) q[2];
sx q[2];
rz(-1.7075677) q[2];
rz(3.1003351) q[3];
sx q[3];
rz(-1.9058803) q[3];
sx q[3];
rz(-2.1432721) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];