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
rz(1.7293575) q[0];
sx q[0];
rz(-2.5492302) q[0];
sx q[0];
rz(1.2047729) q[0];
rz(1.0869979) q[1];
sx q[1];
rz(-2.6462061) q[1];
sx q[1];
rz(-1.9538716) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3760494) q[0];
sx q[0];
rz(-0.57224579) q[0];
sx q[0];
rz(-2.8577198) q[0];
x q[1];
rz(-1.1101025) q[2];
sx q[2];
rz(-1.0506127) q[2];
sx q[2];
rz(2.4536163) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9662142) q[1];
sx q[1];
rz(-1.3982275) q[1];
sx q[1];
rz(0.047974384) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9378565) q[3];
sx q[3];
rz(-1.7008743) q[3];
sx q[3];
rz(0.50103044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5997233) q[2];
sx q[2];
rz(-0.60167998) q[2];
sx q[2];
rz(-2.0987233) q[2];
rz(-0.045182191) q[3];
sx q[3];
rz(-2.9729645) q[3];
sx q[3];
rz(1.5359623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.0586108) q[0];
sx q[0];
rz(-0.11317145) q[0];
sx q[0];
rz(0.97852069) q[0];
rz(1.9784031) q[1];
sx q[1];
rz(-2.1470943) q[1];
sx q[1];
rz(-1.3226343) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90962553) q[0];
sx q[0];
rz(-0.84950209) q[0];
sx q[0];
rz(2.2366174) q[0];
x q[1];
rz(1.2873994) q[2];
sx q[2];
rz(-1.3958418) q[2];
sx q[2];
rz(2.8307479) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6860477) q[1];
sx q[1];
rz(-1.3551222) q[1];
sx q[1];
rz(0.91051813) q[1];
rz(-pi) q[2];
rz(1.421881) q[3];
sx q[3];
rz(-2.2595539) q[3];
sx q[3];
rz(1.0721579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3375552) q[2];
sx q[2];
rz(-0.20122169) q[2];
sx q[2];
rz(-3.1031754) q[2];
rz(-1.6328968) q[3];
sx q[3];
rz(-1.326606) q[3];
sx q[3];
rz(-1.4940777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9134193) q[0];
sx q[0];
rz(-2.4689624) q[0];
sx q[0];
rz(-0.2359373) q[0];
rz(1.963223) q[1];
sx q[1];
rz(-1.447568) q[1];
sx q[1];
rz(-2.6441914) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0097522) q[0];
sx q[0];
rz(-1.3924283) q[0];
sx q[0];
rz(1.7531036) q[0];
x q[1];
rz(-1.986662) q[2];
sx q[2];
rz(-1.207104) q[2];
sx q[2];
rz(1.9726582) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6420501) q[1];
sx q[1];
rz(-1.1772369) q[1];
sx q[1];
rz(1.5469013) q[1];
x q[2];
rz(1.3128993) q[3];
sx q[3];
rz(-1.329943) q[3];
sx q[3];
rz(-1.9633213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1824789) q[2];
sx q[2];
rz(-1.6751869) q[2];
sx q[2];
rz(1.2314931) q[2];
rz(1.69151) q[3];
sx q[3];
rz(-1.8385889) q[3];
sx q[3];
rz(-0.85795295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0602144) q[0];
sx q[0];
rz(-0.84351081) q[0];
sx q[0];
rz(2.9486935) q[0];
rz(1.7861722) q[1];
sx q[1];
rz(-1.5602292) q[1];
sx q[1];
rz(-0.17098175) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20342228) q[0];
sx q[0];
rz(-1.5121542) q[0];
sx q[0];
rz(-2.0873656) q[0];
rz(2.5582983) q[2];
sx q[2];
rz(-0.82628822) q[2];
sx q[2];
rz(2.9982931) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1384004) q[1];
sx q[1];
rz(-0.44282162) q[1];
sx q[1];
rz(-1.3243616) q[1];
x q[2];
rz(1.2661352) q[3];
sx q[3];
rz(-2.4551292) q[3];
sx q[3];
rz(1.0047508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1811447) q[2];
sx q[2];
rz(-1.9713216) q[2];
sx q[2];
rz(-2.81874) q[2];
rz(1.3747831) q[3];
sx q[3];
rz(-1.0389453) q[3];
sx q[3];
rz(0.84097451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5508995) q[0];
sx q[0];
rz(-1.0223848) q[0];
sx q[0];
rz(2.2160231) q[0];
rz(-3.0135221) q[1];
sx q[1];
rz(-1.4984727) q[1];
sx q[1];
rz(-3.0421323) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8680222) q[0];
sx q[0];
rz(-1.5707468) q[0];
sx q[0];
rz(-2.2873061) q[0];
rz(-2.7814513) q[2];
sx q[2];
rz(-1.8966873) q[2];
sx q[2];
rz(-1.858409) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4872553) q[1];
sx q[1];
rz(-2.0981826) q[1];
sx q[1];
rz(0.22589639) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8811536) q[3];
sx q[3];
rz(-2.884293) q[3];
sx q[3];
rz(-2.8043487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9596935) q[2];
sx q[2];
rz(-1.5450031) q[2];
sx q[2];
rz(2.7679475) q[2];
rz(-0.89961189) q[3];
sx q[3];
rz(-2.3989232) q[3];
sx q[3];
rz(-1.1348178) q[3];
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
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9099092) q[0];
sx q[0];
rz(-0.21655701) q[0];
sx q[0];
rz(-2.5496971) q[0];
rz(2.9369211) q[1];
sx q[1];
rz(-2.1494631) q[1];
sx q[1];
rz(0.051102292) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.013553) q[0];
sx q[0];
rz(-2.8951515) q[0];
sx q[0];
rz(0.67009573) q[0];
rz(0.25801664) q[2];
sx q[2];
rz(-0.92890255) q[2];
sx q[2];
rz(-2.8947865) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8569736) q[1];
sx q[1];
rz(-1.2220807) q[1];
sx q[1];
rz(-0.32964175) q[1];
x q[2];
rz(1.072364) q[3];
sx q[3];
rz(-1.0921156) q[3];
sx q[3];
rz(-2.5414027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0245612) q[2];
sx q[2];
rz(-1.1144964) q[2];
sx q[2];
rz(1.0943817) q[2];
rz(0.38703212) q[3];
sx q[3];
rz(-0.47457591) q[3];
sx q[3];
rz(2.835137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36560202) q[0];
sx q[0];
rz(-2.2633573) q[0];
sx q[0];
rz(0.24755092) q[0];
rz(-1.6869102) q[1];
sx q[1];
rz(-1.2431815) q[1];
sx q[1];
rz(-0.036458485) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55900967) q[0];
sx q[0];
rz(-1.2866308) q[0];
sx q[0];
rz(-2.3733632) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.22561947) q[2];
sx q[2];
rz(-0.7157514) q[2];
sx q[2];
rz(0.72474397) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9769638) q[1];
sx q[1];
rz(-1.1067302) q[1];
sx q[1];
rz(1.0010757) q[1];
x q[2];
rz(2.2685675) q[3];
sx q[3];
rz(-0.11618488) q[3];
sx q[3];
rz(-0.98770638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.054691943) q[2];
sx q[2];
rz(-2.1754103) q[2];
sx q[2];
rz(-1.3521693) q[2];
rz(-1.2482721) q[3];
sx q[3];
rz(-1.5779481) q[3];
sx q[3];
rz(1.9498391) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9098814) q[0];
sx q[0];
rz(-2.8928962) q[0];
sx q[0];
rz(-1.6491718) q[0];
rz(2.9630648) q[1];
sx q[1];
rz(-1.2622204) q[1];
sx q[1];
rz(-0.55007225) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72075413) q[0];
sx q[0];
rz(-1.5394569) q[0];
sx q[0];
rz(-0.019958812) q[0];
rz(3.078207) q[2];
sx q[2];
rz(-1.3235385) q[2];
sx q[2];
rz(-0.72615577) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7272721) q[1];
sx q[1];
rz(-1.7309233) q[1];
sx q[1];
rz(-2.5989207) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86217238) q[3];
sx q[3];
rz(-0.58364999) q[3];
sx q[3];
rz(-0.72008163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7940346) q[2];
sx q[2];
rz(-1.7358235) q[2];
sx q[2];
rz(0.75322914) q[2];
rz(2.8473162) q[3];
sx q[3];
rz(-0.080065057) q[3];
sx q[3];
rz(3.0729955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2733961) q[0];
sx q[0];
rz(-2.8559339) q[0];
sx q[0];
rz(-1.4519325) q[0];
rz(-2.946335) q[1];
sx q[1];
rz(-2.0611019) q[1];
sx q[1];
rz(-1.4917699) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0367835) q[0];
sx q[0];
rz(-1.097479) q[0];
sx q[0];
rz(-1.9154873) q[0];
rz(-pi) q[1];
rz(0.765018) q[2];
sx q[2];
rz(-1.6510626) q[2];
sx q[2];
rz(1.2470055) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0780219) q[1];
sx q[1];
rz(-1.6352904) q[1];
sx q[1];
rz(1.7266573) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5641065) q[3];
sx q[3];
rz(-1.5433528) q[3];
sx q[3];
rz(-0.60718482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0407654) q[2];
sx q[2];
rz(-0.46697524) q[2];
sx q[2];
rz(-2.2838498) q[2];
rz(-1.6929251) q[3];
sx q[3];
rz(-1.6489776) q[3];
sx q[3];
rz(-2.9878476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0070852) q[0];
sx q[0];
rz(-0.15468287) q[0];
sx q[0];
rz(2.3585228) q[0];
rz(-2.018441) q[1];
sx q[1];
rz(-1.4479366) q[1];
sx q[1];
rz(-3.0812841) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6864844) q[0];
sx q[0];
rz(-1.9251072) q[0];
sx q[0];
rz(3.0418116) q[0];
rz(-pi) q[1];
x q[1];
rz(0.98452703) q[2];
sx q[2];
rz(-1.482555) q[2];
sx q[2];
rz(2.4172104) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7544514) q[1];
sx q[1];
rz(-1.5031414) q[1];
sx q[1];
rz(2.4715502) q[1];
x q[2];
rz(-0.29197146) q[3];
sx q[3];
rz(-0.37759104) q[3];
sx q[3];
rz(0.99536125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.15410885) q[2];
sx q[2];
rz(-2.2769112) q[2];
sx q[2];
rz(-1.3101428) q[2];
rz(1.8080669) q[3];
sx q[3];
rz(-1.3430877) q[3];
sx q[3];
rz(-1.8346627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6780728) q[0];
sx q[0];
rz(-1.1165883) q[0];
sx q[0];
rz(-0.22620871) q[0];
rz(0.75640596) q[1];
sx q[1];
rz(-2.6466128) q[1];
sx q[1];
rz(-0.80642798) q[1];
rz(3.1210774) q[2];
sx q[2];
rz(-0.43045577) q[2];
sx q[2];
rz(0.53878709) q[2];
rz(1.9202833) q[3];
sx q[3];
rz(-1.3987473) q[3];
sx q[3];
rz(-1.7940298) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
