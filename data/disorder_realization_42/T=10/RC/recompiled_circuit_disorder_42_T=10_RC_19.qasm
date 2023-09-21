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
rz(0.36261121) q[0];
rz(0.91712046) q[1];
sx q[1];
rz(-0.4904823) q[1];
sx q[1];
rz(2.7999556) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5255614) q[0];
sx q[0];
rz(-1.4481059) q[0];
sx q[0];
rz(-2.7659155) q[0];
x q[1];
rz(2.5506637) q[2];
sx q[2];
rz(-1.2714296) q[2];
sx q[2];
rz(-2.9654944) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8056148) q[1];
sx q[1];
rz(-2.0819547) q[1];
sx q[1];
rz(-2.7989945) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.2239286) q[3];
sx q[3];
rz(-1.180598) q[3];
sx q[3];
rz(-1.8412561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5477156) q[2];
sx q[2];
rz(-1.1527858) q[2];
sx q[2];
rz(-0.95648742) q[2];
rz(-2.9603413) q[3];
sx q[3];
rz(-0.58877188) q[3];
sx q[3];
rz(-1.5216924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0797121) q[0];
sx q[0];
rz(-3.0759838) q[0];
sx q[0];
rz(1.99615) q[0];
rz(-2.0727797) q[1];
sx q[1];
rz(-2.309598) q[1];
sx q[1];
rz(-0.72584814) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9024591) q[0];
sx q[0];
rz(-2.4665138) q[0];
sx q[0];
rz(2.549987) q[0];
x q[1];
rz(-1.1109353) q[2];
sx q[2];
rz(-2.1150555) q[2];
sx q[2];
rz(2.0366675) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5686381) q[1];
sx q[1];
rz(-1.8796088) q[1];
sx q[1];
rz(-0.026334865) q[1];
rz(-pi) q[2];
rz(-1.3594567) q[3];
sx q[3];
rz(-0.95892116) q[3];
sx q[3];
rz(0.26821995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.48250616) q[2];
sx q[2];
rz(-1.3799474) q[2];
sx q[2];
rz(-0.99728161) q[2];
rz(-1.4128134) q[3];
sx q[3];
rz(-2.0480859) q[3];
sx q[3];
rz(-2.2028082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41855758) q[0];
sx q[0];
rz(-2.1756873) q[0];
sx q[0];
rz(2.0130656) q[0];
rz(-1.3035125) q[1];
sx q[1];
rz(-0.729527) q[1];
sx q[1];
rz(0.9054786) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0008464) q[0];
sx q[0];
rz(-1.329639) q[0];
sx q[0];
rz(1.3748037) q[0];
rz(0.75240527) q[2];
sx q[2];
rz(-1.89626) q[2];
sx q[2];
rz(-0.67093713) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.803373) q[1];
sx q[1];
rz(-1.6446911) q[1];
sx q[1];
rz(-0.85892962) q[1];
rz(-pi) q[2];
x q[2];
rz(1.521103) q[3];
sx q[3];
rz(-2.4274785) q[3];
sx q[3];
rz(-0.36220887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.7604312) q[2];
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
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0320597) q[0];
sx q[0];
rz(-2.4592472) q[0];
sx q[0];
rz(-0.12606829) q[0];
rz(-0.18445045) q[1];
sx q[1];
rz(-2.0729063) q[1];
sx q[1];
rz(-1.1674081) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8465189) q[0];
sx q[0];
rz(-0.82058883) q[0];
sx q[0];
rz(1.2942765) q[0];
rz(1.7502968) q[2];
sx q[2];
rz(-2.0677535) q[2];
sx q[2];
rz(-2.3617982) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3525225) q[1];
sx q[1];
rz(-1.8348502) q[1];
sx q[1];
rz(-2.4697484) q[1];
rz(-pi) q[2];
rz(-2.4950124) q[3];
sx q[3];
rz(-1.0160035) q[3];
sx q[3];
rz(3.0559029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.8644774) q[2];
sx q[2];
rz(-2.7272759) q[2];
sx q[2];
rz(-0.63151044) q[2];
rz(-0.19550066) q[3];
sx q[3];
rz(-1.2771527) q[3];
sx q[3];
rz(2.7235532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.607476) q[0];
sx q[0];
rz(-0.42629888) q[0];
sx q[0];
rz(2.0712453) q[0];
x q[1];
rz(1.3181114) q[2];
sx q[2];
rz(-0.99950302) q[2];
sx q[2];
rz(-2.2342482) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0568697) q[1];
sx q[1];
rz(-1.9922171) q[1];
sx q[1];
rz(0.23878581) q[1];
rz(-0.15848666) q[3];
sx q[3];
rz(-2.2500258) q[3];
sx q[3];
rz(-2.2331626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6749394) q[2];
sx q[2];
rz(-1.249524) q[2];
sx q[2];
rz(-2.6494027) q[2];
rz(2.8594033) q[3];
sx q[3];
rz(-1.4331093) q[3];
sx q[3];
rz(1.5658437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.870134) q[0];
sx q[0];
rz(-2.6014355) q[0];
sx q[0];
rz(-3.0555994) q[0];
rz(1.7135235) q[1];
sx q[1];
rz(-2.1336887) q[1];
sx q[1];
rz(1.496398) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7487885) q[0];
sx q[0];
rz(-1.7668056) q[0];
sx q[0];
rz(-1.0402354) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1804579) q[2];
sx q[2];
rz(-2.5679553) q[2];
sx q[2];
rz(-2.8611285) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.23425809) q[1];
sx q[1];
rz(-1.2163711) q[1];
sx q[1];
rz(-0.025269421) q[1];
rz(-1.5580642) q[3];
sx q[3];
rz(-0.85789645) q[3];
sx q[3];
rz(-1.3811779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4795586) q[2];
sx q[2];
rz(-0.61926121) q[2];
sx q[2];
rz(-2.7453444) q[2];
rz(-0.59404343) q[3];
sx q[3];
rz(-2.0500573) q[3];
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
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8783766) q[0];
sx q[0];
rz(-2.5686869) q[0];
sx q[0];
rz(-1.0923882) q[0];
rz(1.7842402) q[1];
sx q[1];
rz(-1.7204294) q[1];
sx q[1];
rz(-0.011627442) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41650018) q[0];
sx q[0];
rz(-1.9511127) q[0];
sx q[0];
rz(-1.8561383) q[0];
x q[1];
rz(-3.1161357) q[2];
sx q[2];
rz(-1.5989132) q[2];
sx q[2];
rz(0.97719976) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1639451) q[1];
sx q[1];
rz(-1.5471336) q[1];
sx q[1];
rz(1.0246828) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2467975) q[3];
sx q[3];
rz(-0.9231336) q[3];
sx q[3];
rz(-3.0917633) q[3];
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
rz(1.9912432) q[3];
sx q[3];
rz(-1.8201613) q[3];
sx q[3];
rz(-2.5434125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31036723) q[0];
sx q[0];
rz(-2.6103525) q[0];
sx q[0];
rz(-1.9294552) q[0];
rz(-0.37777004) q[1];
sx q[1];
rz(-1.2021844) q[1];
sx q[1];
rz(1.4935965) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50266788) q[0];
sx q[0];
rz(-1.8687316) q[0];
sx q[0];
rz(-1.3475247) q[0];
rz(-pi) q[1];
rz(1.9940894) q[2];
sx q[2];
rz(-2.9967542) q[2];
sx q[2];
rz(-0.3651948) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6495325) q[1];
sx q[1];
rz(-0.82693716) q[1];
sx q[1];
rz(-2.5745113) q[1];
rz(-pi) q[2];
rz(-2.0728552) q[3];
sx q[3];
rz(-1.4668674) q[3];
sx q[3];
rz(-2.7285189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7887855) q[2];
sx q[2];
rz(-2.0264503) q[2];
sx q[2];
rz(-1.8898233) q[2];
rz(2.2552323) q[3];
sx q[3];
rz(-0.76615196) q[3];
sx q[3];
rz(-3.0322976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66824526) q[0];
sx q[0];
rz(-2.1003523) q[0];
sx q[0];
rz(0.1782724) q[0];
rz(3.0687304) q[1];
sx q[1];
rz(-0.59385121) q[1];
sx q[1];
rz(1.1791621) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0807063) q[0];
sx q[0];
rz(-2.076425) q[0];
sx q[0];
rz(-1.2348742) q[0];
x q[1];
rz(0.5046919) q[2];
sx q[2];
rz(-2.1345277) q[2];
sx q[2];
rz(1.5228412) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.33680962) q[1];
sx q[1];
rz(-2.0840624) q[1];
sx q[1];
rz(1.8670765) q[1];
rz(-pi) q[2];
rz(2.7618802) q[3];
sx q[3];
rz(-1.100193) q[3];
sx q[3];
rz(3.0203366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.69958413) q[2];
sx q[2];
rz(-0.99094355) q[2];
sx q[2];
rz(-2.1098095) q[2];
rz(-1.7992841) q[3];
sx q[3];
rz(-1.4167901) q[3];
sx q[3];
rz(-2.9163196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28829065) q[0];
sx q[0];
rz(-2.0993711) q[0];
sx q[0];
rz(0.061696079) q[0];
rz(1.306698) q[1];
sx q[1];
rz(-1.4790165) q[1];
sx q[1];
rz(1.0888938) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13578829) q[0];
sx q[0];
rz(-0.09011589) q[0];
sx q[0];
rz(0.80149217) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9577272) q[2];
sx q[2];
rz(-0.92391787) q[2];
sx q[2];
rz(-2.0053787) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.006146487) q[1];
sx q[1];
rz(-1.0137614) q[1];
sx q[1];
rz(2.1250493) q[1];
x q[2];
rz(-2.9519134) q[3];
sx q[3];
rz(-1.9072201) q[3];
sx q[3];
rz(1.3996901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5237727) q[2];
sx q[2];
rz(-2.4536295) q[2];
sx q[2];
rz(0.069996746) q[2];
rz(2.2648515) q[3];
sx q[3];
rz(-1.8165959) q[3];
sx q[3];
rz(-2.783412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4927647) q[0];
sx q[0];
rz(-1.6236826) q[0];
sx q[0];
rz(0.59326011) q[0];
rz(-1.6246673) q[1];
sx q[1];
rz(-2.0943874) q[1];
sx q[1];
rz(2.692683) q[1];
rz(2.0786053) q[2];
sx q[2];
rz(-0.72401902) q[2];
sx q[2];
rz(0.25821092) q[2];
rz(1.6886961) q[3];
sx q[3];
rz(-2.8040734) q[3];
sx q[3];
rz(0.87344195) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];