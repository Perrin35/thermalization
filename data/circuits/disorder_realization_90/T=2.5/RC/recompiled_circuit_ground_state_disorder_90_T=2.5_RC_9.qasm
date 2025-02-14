OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5034135) q[0];
sx q[0];
rz(-1.2794275) q[0];
sx q[0];
rz(0.77603618) q[0];
rz(-2.4322721) q[1];
sx q[1];
rz(-1.5172989) q[1];
sx q[1];
rz(-2.5425743) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.020928362) q[0];
sx q[0];
rz(-2.0428223) q[0];
sx q[0];
rz(0.13735227) q[0];
x q[1];
rz(-1.9089977) q[2];
sx q[2];
rz(-0.71824348) q[2];
sx q[2];
rz(-2.3423549) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.30556074) q[1];
sx q[1];
rz(-1.9232043) q[1];
sx q[1];
rz(-2.6120017) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1407858) q[3];
sx q[3];
rz(-2.3035192) q[3];
sx q[3];
rz(0.091146745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0029995) q[2];
sx q[2];
rz(-1.1370167) q[2];
sx q[2];
rz(-0.19621672) q[2];
rz(-1.044322) q[3];
sx q[3];
rz(-2.315867) q[3];
sx q[3];
rz(-0.31061068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-1.0394548) q[0];
sx q[0];
rz(-0.65334833) q[0];
sx q[0];
rz(-0.85773221) q[0];
rz(2.6013382) q[1];
sx q[1];
rz(-1.0539571) q[1];
sx q[1];
rz(0.78261715) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1038541) q[0];
sx q[0];
rz(-1.6574934) q[0];
sx q[0];
rz(-2.187665) q[0];
rz(-pi) q[1];
rz(-2.906053) q[2];
sx q[2];
rz(-2.1793834) q[2];
sx q[2];
rz(-2.4289301) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5021584) q[1];
sx q[1];
rz(-2.5089988) q[1];
sx q[1];
rz(-0.32986395) q[1];
x q[2];
rz(-0.55763839) q[3];
sx q[3];
rz(-2.1774315) q[3];
sx q[3];
rz(-2.5733657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.1648078) q[2];
sx q[2];
rz(-0.79443496) q[2];
sx q[2];
rz(2.5505193) q[2];
rz(1.0278541) q[3];
sx q[3];
rz(-2.5058392) q[3];
sx q[3];
rz(-3.0052321) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22981055) q[0];
sx q[0];
rz(-2.6214143) q[0];
sx q[0];
rz(-2.0853364) q[0];
rz(-2.1504869) q[1];
sx q[1];
rz(-1.3737563) q[1];
sx q[1];
rz(-0.80449218) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1209542) q[0];
sx q[0];
rz(-0.94392969) q[0];
sx q[0];
rz(-1.2406209) q[0];
rz(-pi) q[1];
rz(1.7011004) q[2];
sx q[2];
rz(-2.5087803) q[2];
sx q[2];
rz(1.9595343) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.75315216) q[1];
sx q[1];
rz(-0.33666753) q[1];
sx q[1];
rz(0.80516385) q[1];
x q[2];
rz(0.64739703) q[3];
sx q[3];
rz(-1.7248099) q[3];
sx q[3];
rz(1.2274418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.450401) q[2];
sx q[2];
rz(-0.85019008) q[2];
sx q[2];
rz(-2.5737393) q[2];
rz(-1.19207) q[3];
sx q[3];
rz(-0.53693938) q[3];
sx q[3];
rz(0.15538628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
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
rz(1.6464226) q[0];
sx q[0];
rz(-1.3986873) q[0];
sx q[0];
rz(1.9544741) q[0];
rz(0.25539708) q[1];
sx q[1];
rz(-1.4030158) q[1];
sx q[1];
rz(2.752221) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0559278) q[0];
sx q[0];
rz(-1.7048258) q[0];
sx q[0];
rz(2.9197951) q[0];
rz(1.537048) q[2];
sx q[2];
rz(-0.72740388) q[2];
sx q[2];
rz(0.45798618) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.43453211) q[1];
sx q[1];
rz(-1.7599306) q[1];
sx q[1];
rz(-2.2347336) q[1];
x q[2];
rz(0.98020966) q[3];
sx q[3];
rz(-2.906153) q[3];
sx q[3];
rz(0.93577945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7369467) q[2];
sx q[2];
rz(-2.8179822) q[2];
sx q[2];
rz(1.3673937) q[2];
rz(0.5395475) q[3];
sx q[3];
rz(-1.5893785) q[3];
sx q[3];
rz(1.3246983) q[3];
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
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0747727) q[0];
sx q[0];
rz(-2.9434151) q[0];
sx q[0];
rz(2.5812126) q[0];
rz(-0.21891521) q[1];
sx q[1];
rz(-1.0812662) q[1];
sx q[1];
rz(2.970649) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8290299) q[0];
sx q[0];
rz(-1.1033022) q[0];
sx q[0];
rz(0.74689052) q[0];
x q[1];
rz(-2.937491) q[2];
sx q[2];
rz(-1.4320254) q[2];
sx q[2];
rz(-2.6134174) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.584632) q[1];
sx q[1];
rz(-1.0611532) q[1];
sx q[1];
rz(2.4181448) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4724949) q[3];
sx q[3];
rz(-2.519033) q[3];
sx q[3];
rz(-3.0873201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.35974744) q[2];
sx q[2];
rz(-1.1283987) q[2];
sx q[2];
rz(-0.95853364) q[2];
rz(0.16252276) q[3];
sx q[3];
rz(-0.84238094) q[3];
sx q[3];
rz(-1.5490279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.573134) q[0];
sx q[0];
rz(-3.077226) q[0];
sx q[0];
rz(0.74813133) q[0];
rz(2.023078) q[1];
sx q[1];
rz(-0.41901127) q[1];
sx q[1];
rz(-0.31707877) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9592764) q[0];
sx q[0];
rz(-1.8166564) q[0];
sx q[0];
rz(-1.5005174) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.19642475) q[2];
sx q[2];
rz(-1.6187917) q[2];
sx q[2];
rz(2.6065207) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.85641454) q[1];
sx q[1];
rz(-0.76907255) q[1];
sx q[1];
rz(-1.4682707) q[1];
rz(-pi) q[2];
rz(2.3646687) q[3];
sx q[3];
rz(-1.5990997) q[3];
sx q[3];
rz(-1.5212718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.76724425) q[2];
sx q[2];
rz(-2.2129462) q[2];
sx q[2];
rz(1.7283776) q[2];
rz(-3.0610541) q[3];
sx q[3];
rz(-1.3486515) q[3];
sx q[3];
rz(0.22182375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11248511) q[0];
sx q[0];
rz(-1.8928098) q[0];
sx q[0];
rz(-1.2741733) q[0];
rz(1.796465) q[1];
sx q[1];
rz(-0.74256623) q[1];
sx q[1];
rz(1.3412195) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31650922) q[0];
sx q[0];
rz(-2.1232455) q[0];
sx q[0];
rz(2.1567206) q[0];
rz(-pi) q[1];
rz(-0.48168452) q[2];
sx q[2];
rz(-1.8284214) q[2];
sx q[2];
rz(3.1236908) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.37137274) q[1];
sx q[1];
rz(-2.8951982) q[1];
sx q[1];
rz(-2.2486992) q[1];
rz(-pi) q[2];
rz(2.2209211) q[3];
sx q[3];
rz(-1.4772282) q[3];
sx q[3];
rz(-0.9182932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.207927) q[2];
sx q[2];
rz(-1.1873446) q[2];
sx q[2];
rz(2.6742317) q[2];
rz(0.76006877) q[3];
sx q[3];
rz(-1.6642539) q[3];
sx q[3];
rz(1.8925586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46045983) q[0];
sx q[0];
rz(-2.12119) q[0];
sx q[0];
rz(1.6946174) q[0];
rz(2.815822) q[1];
sx q[1];
rz(-2.9278946) q[1];
sx q[1];
rz(1.5209341) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9252773) q[0];
sx q[0];
rz(-2.0774573) q[0];
sx q[0];
rz(-1.265488) q[0];
rz(-pi) q[1];
rz(3.1374819) q[2];
sx q[2];
rz(-0.09164022) q[2];
sx q[2];
rz(0.81697538) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3586761) q[1];
sx q[1];
rz(-2.7442928) q[1];
sx q[1];
rz(2.0372169) q[1];
rz(-pi) q[2];
rz(2.9321947) q[3];
sx q[3];
rz(-2.2639674) q[3];
sx q[3];
rz(1.5537167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.96523607) q[2];
sx q[2];
rz(-0.74702817) q[2];
sx q[2];
rz(0.61526862) q[2];
rz(1.2728914) q[3];
sx q[3];
rz(-2.8764909) q[3];
sx q[3];
rz(-2.7191775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0203005) q[0];
sx q[0];
rz(-1.2739807) q[0];
sx q[0];
rz(0.85025382) q[0];
rz(-1.1212768) q[1];
sx q[1];
rz(-1.3669776) q[1];
sx q[1];
rz(0.77176315) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4798652) q[0];
sx q[0];
rz(-0.87718946) q[0];
sx q[0];
rz(-3.1105009) q[0];
rz(-pi) q[1];
rz(-2.9608126) q[2];
sx q[2];
rz(-1.639699) q[2];
sx q[2];
rz(3.1071752) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1771639) q[1];
sx q[1];
rz(-1.6258952) q[1];
sx q[1];
rz(-1.8146309) q[1];
rz(-pi) q[2];
rz(1.2466627) q[3];
sx q[3];
rz(-1.3132902) q[3];
sx q[3];
rz(2.3070564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0261592) q[2];
sx q[2];
rz(-1.5449646) q[2];
sx q[2];
rz(0.8405295) q[2];
rz(0.080951512) q[3];
sx q[3];
rz(-2.5637124) q[3];
sx q[3];
rz(-2.8988083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7131272) q[0];
sx q[0];
rz(-0.19793887) q[0];
sx q[0];
rz(-1.9848829) q[0];
rz(-1.0276065) q[1];
sx q[1];
rz(-1.7637858) q[1];
sx q[1];
rz(-0.81046945) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58662232) q[0];
sx q[0];
rz(-1.3548618) q[0];
sx q[0];
rz(-0.76540375) q[0];
rz(2.5963411) q[2];
sx q[2];
rz(-1.8601189) q[2];
sx q[2];
rz(2.3985753) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4833853) q[1];
sx q[1];
rz(-1.2772296) q[1];
sx q[1];
rz(-0.29694966) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23328554) q[3];
sx q[3];
rz(-0.57396997) q[3];
sx q[3];
rz(-3.1149394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.345574) q[2];
sx q[2];
rz(-1.8859665) q[2];
sx q[2];
rz(1.6176809) q[2];
rz(-2.0412622) q[3];
sx q[3];
rz(-1.1805781) q[3];
sx q[3];
rz(2.9617917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0236459) q[0];
sx q[0];
rz(-1.7272341) q[0];
sx q[0];
rz(-0.9247307) q[0];
rz(0.042451518) q[1];
sx q[1];
rz(-1.1141384) q[1];
sx q[1];
rz(1.3420807) q[1];
rz(-0.056006147) q[2];
sx q[2];
rz(-1.9765035) q[2];
sx q[2];
rz(-0.34513173) q[2];
rz(0.45741882) q[3];
sx q[3];
rz(-1.6030967) q[3];
sx q[3];
rz(0.23489192) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
