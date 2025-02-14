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
rz(-2.7491375) q[0];
sx q[0];
rz(-2.8578704) q[0];
sx q[0];
rz(-1.1878045) q[0];
rz(2.6999733) q[1];
sx q[1];
rz(-1.266357) q[1];
sx q[1];
rz(-1.1511572) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1063163) q[0];
sx q[0];
rz(-1.6589249) q[0];
sx q[0];
rz(3.1195537) q[0];
x q[1];
rz(0.63090905) q[2];
sx q[2];
rz(-0.5782457) q[2];
sx q[2];
rz(-1.9672245) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3656039) q[1];
sx q[1];
rz(-1.1465441) q[1];
sx q[1];
rz(1.7107023) q[1];
rz(-1.0488767) q[3];
sx q[3];
rz(-1.8303695) q[3];
sx q[3];
rz(2.5775016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2426131) q[2];
sx q[2];
rz(-0.8967163) q[2];
sx q[2];
rz(0.29956451) q[2];
rz(2.7859531) q[3];
sx q[3];
rz(-0.99323291) q[3];
sx q[3];
rz(-2.7049844) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75293175) q[0];
sx q[0];
rz(-2.5932073) q[0];
sx q[0];
rz(-2.8652628) q[0];
rz(-3.062124) q[1];
sx q[1];
rz(-0.53229585) q[1];
sx q[1];
rz(-0.4963378) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4547098) q[0];
sx q[0];
rz(-2.080603) q[0];
sx q[0];
rz(-2.8532993) q[0];
rz(-2.8032643) q[2];
sx q[2];
rz(-1.4979432) q[2];
sx q[2];
rz(-3.0908527) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4071634) q[1];
sx q[1];
rz(-1.8968655) q[1];
sx q[1];
rz(0.77496518) q[1];
rz(1.4088936) q[3];
sx q[3];
rz(-1.5739958) q[3];
sx q[3];
rz(2.7174983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.77208272) q[2];
sx q[2];
rz(-1.3630867) q[2];
sx q[2];
rz(-2.962964) q[2];
rz(1.5313088) q[3];
sx q[3];
rz(-2.2154) q[3];
sx q[3];
rz(2.3791651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.4302706) q[0];
sx q[0];
rz(-1.0866168) q[0];
sx q[0];
rz(-1.7311199) q[0];
rz(1.3871644) q[1];
sx q[1];
rz(-1.6491363) q[1];
sx q[1];
rz(-14/(3*pi)) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0091772) q[0];
sx q[0];
rz(-0.88736358) q[0];
sx q[0];
rz(1.528549) q[0];
rz(-pi) q[1];
rz(1.1457655) q[2];
sx q[2];
rz(-2.3299937) q[2];
sx q[2];
rz(-1.5852864) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8456257) q[1];
sx q[1];
rz(-2.6085771) q[1];
sx q[1];
rz(0.94818268) q[1];
rz(-2.0250506) q[3];
sx q[3];
rz(-0.71703934) q[3];
sx q[3];
rz(-1.7771135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2014655) q[2];
sx q[2];
rz(-2.1374233) q[2];
sx q[2];
rz(2.6489769) q[2];
rz(0.41871873) q[3];
sx q[3];
rz(-2.1188354) q[3];
sx q[3];
rz(-2.54971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.017460499) q[0];
sx q[0];
rz(-3.0191665) q[0];
sx q[0];
rz(-0.19805743) q[0];
rz(0.65912229) q[1];
sx q[1];
rz(-0.56050572) q[1];
sx q[1];
rz(0.19703279) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8640932) q[0];
sx q[0];
rz(-0.53170337) q[0];
sx q[0];
rz(-3.0752027) q[0];
rz(0.87984271) q[2];
sx q[2];
rz(-2.7759399) q[2];
sx q[2];
rz(-1.8118389) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9777364) q[1];
sx q[1];
rz(-1.9242373) q[1];
sx q[1];
rz(-0.56956915) q[1];
rz(3.1048685) q[3];
sx q[3];
rz(-1.0261327) q[3];
sx q[3];
rz(0.36946083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.78843242) q[2];
sx q[2];
rz(-1.8400695) q[2];
sx q[2];
rz(-1.2737459) q[2];
rz(-1.5431917) q[3];
sx q[3];
rz(-1.9485995) q[3];
sx q[3];
rz(0.25632349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0565599) q[0];
sx q[0];
rz(-1.7340478) q[0];
sx q[0];
rz(-0.10277596) q[0];
rz(-3.00434) q[1];
sx q[1];
rz(-2.6922701) q[1];
sx q[1];
rz(-1.7030254) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53046561) q[0];
sx q[0];
rz(-2.1411863) q[0];
sx q[0];
rz(0.80369759) q[0];
x q[1];
rz(-1.5248143) q[2];
sx q[2];
rz(-1.1080989) q[2];
sx q[2];
rz(1.785194) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.66617353) q[1];
sx q[1];
rz(-1.8295145) q[1];
sx q[1];
rz(-2.9889876) q[1];
rz(0.12326316) q[3];
sx q[3];
rz(-1.4597037) q[3];
sx q[3];
rz(2.9934903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.30456257) q[2];
sx q[2];
rz(-1.6670767) q[2];
sx q[2];
rz(0.70011082) q[2];
rz(-0.89811283) q[3];
sx q[3];
rz(-2.5105295) q[3];
sx q[3];
rz(-1.496544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9333618) q[0];
sx q[0];
rz(-0.58657402) q[0];
sx q[0];
rz(-1.5775648) q[0];
rz(-2.481752) q[1];
sx q[1];
rz(-0.78958646) q[1];
sx q[1];
rz(-0.97445625) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.385372) q[0];
sx q[0];
rz(-1.5819183) q[0];
sx q[0];
rz(-1.5626989) q[0];
rz(-pi) q[1];
rz(1.7436557) q[2];
sx q[2];
rz(-2.0522016) q[2];
sx q[2];
rz(-2.0916033) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8169957) q[1];
sx q[1];
rz(-1.5600509) q[1];
sx q[1];
rz(0.50159494) q[1];
x q[2];
rz(0.10147126) q[3];
sx q[3];
rz(-0.54059282) q[3];
sx q[3];
rz(1.056788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4727891) q[2];
sx q[2];
rz(-0.69837022) q[2];
sx q[2];
rz(-0.93639708) q[2];
rz(0.89013731) q[3];
sx q[3];
rz(-1.3606768) q[3];
sx q[3];
rz(2.9357125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5278006) q[0];
sx q[0];
rz(-2.9449936) q[0];
sx q[0];
rz(0.11372456) q[0];
rz(1.744489) q[1];
sx q[1];
rz(-1.0662268) q[1];
sx q[1];
rz(-0.024959175) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7222116) q[0];
sx q[0];
rz(-1.7891277) q[0];
sx q[0];
rz(-3.0857357) q[0];
rz(-pi) q[1];
x q[1];
rz(0.34100248) q[2];
sx q[2];
rz(-0.67217584) q[2];
sx q[2];
rz(-2.9580046) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.125122) q[1];
sx q[1];
rz(-0.85713398) q[1];
sx q[1];
rz(3.0550856) q[1];
rz(0.27351348) q[3];
sx q[3];
rz(-2.9363147) q[3];
sx q[3];
rz(-2.5375348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.52594841) q[2];
sx q[2];
rz(-1.7082461) q[2];
sx q[2];
rz(2.136266) q[2];
rz(-2.2510236) q[3];
sx q[3];
rz(-1.4054479) q[3];
sx q[3];
rz(0.65517202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6152076) q[0];
sx q[0];
rz(-2.7387775) q[0];
sx q[0];
rz(-1.1108387) q[0];
rz(3.0530744) q[1];
sx q[1];
rz(-2.7169777) q[1];
sx q[1];
rz(1.3540589) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41218573) q[0];
sx q[0];
rz(-2.3581187) q[0];
sx q[0];
rz(-1.4108578) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0646787) q[2];
sx q[2];
rz(-0.30914834) q[2];
sx q[2];
rz(0.34999166) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.718487) q[1];
sx q[1];
rz(-0.25782789) q[1];
sx q[1];
rz(2.8400374) q[1];
rz(-0.60175394) q[3];
sx q[3];
rz(-2.2838998) q[3];
sx q[3];
rz(2.8181638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2985349) q[2];
sx q[2];
rz(-0.83710805) q[2];
sx q[2];
rz(2.7561772) q[2];
rz(-0.75000969) q[3];
sx q[3];
rz(-0.95732006) q[3];
sx q[3];
rz(0.065940417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1239922) q[0];
sx q[0];
rz(-2.0662859) q[0];
sx q[0];
rz(-0.83220926) q[0];
rz(2.0596793) q[1];
sx q[1];
rz(-0.81580201) q[1];
sx q[1];
rz(1.6197416) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6523217) q[0];
sx q[0];
rz(-2.0390563) q[0];
sx q[0];
rz(0.44083198) q[0];
rz(2.887313) q[2];
sx q[2];
rz(-1.2966208) q[2];
sx q[2];
rz(0.57480301) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.79241133) q[1];
sx q[1];
rz(-2.1508625) q[1];
sx q[1];
rz(0.68934727) q[1];
rz(-0.91401871) q[3];
sx q[3];
rz(-1.7631712) q[3];
sx q[3];
rz(0.094948204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1443783) q[2];
sx q[2];
rz(-1.245456) q[2];
sx q[2];
rz(-1.2388371) q[2];
rz(1.291412) q[3];
sx q[3];
rz(-1.2715128) q[3];
sx q[3];
rz(-0.13135697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0111897) q[0];
sx q[0];
rz(-2.1881396) q[0];
sx q[0];
rz(-2.3977872) q[0];
rz(-3.0939843) q[1];
sx q[1];
rz(-2.036939) q[1];
sx q[1];
rz(2.0000752) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6689169) q[0];
sx q[0];
rz(-2.0341691) q[0];
sx q[0];
rz(2.1457304) q[0];
rz(-pi) q[1];
rz(2.6180778) q[2];
sx q[2];
rz(-2.4552629) q[2];
sx q[2];
rz(1.4586142) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9824578) q[1];
sx q[1];
rz(-1.9736814) q[1];
sx q[1];
rz(0.77352685) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.64472736) q[3];
sx q[3];
rz(-1.4253741) q[3];
sx q[3];
rz(2.4032088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6514674) q[2];
sx q[2];
rz(-2.9059124) q[2];
sx q[2];
rz(-0.30533314) q[2];
rz(1.7134282) q[3];
sx q[3];
rz(-1.7686663) q[3];
sx q[3];
rz(1.292424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3931428) q[0];
sx q[0];
rz(-0.23101692) q[0];
sx q[0];
rz(1.1363181) q[0];
rz(-2.582386) q[1];
sx q[1];
rz(-1.8103841) q[1];
sx q[1];
rz(0.24787535) q[1];
rz(0.66721532) q[2];
sx q[2];
rz(-1.2922774) q[2];
sx q[2];
rz(-0.86624672) q[2];
rz(-1.4356267) q[3];
sx q[3];
rz(-1.6886938) q[3];
sx q[3];
rz(2.0791048) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
