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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65389079) q[0];
sx q[0];
rz(-2.747295) q[0];
sx q[0];
rz(-0.32422347) q[0];
rz(-pi) q[1];
rz(-0.59092893) q[2];
sx q[2];
rz(-1.2714296) q[2];
sx q[2];
rz(0.17609827) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.33597782) q[1];
sx q[1];
rz(-1.059638) q[1];
sx q[1];
rz(0.34259818) q[1];
rz(-pi) q[2];
rz(-0.2239286) q[3];
sx q[3];
rz(-1.180598) q[3];
sx q[3];
rz(-1.8412561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.59387702) q[2];
sx q[2];
rz(-1.1527858) q[2];
sx q[2];
rz(-2.1851052) q[2];
rz(0.18125136) q[3];
sx q[3];
rz(-2.5528208) q[3];
sx q[3];
rz(-1.6199002) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0618806) q[0];
sx q[0];
rz(-0.065608874) q[0];
sx q[0];
rz(1.1454426) q[0];
rz(-2.0727797) q[1];
sx q[1];
rz(-2.309598) q[1];
sx q[1];
rz(-0.72584814) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1918068) q[0];
sx q[0];
rz(-1.0254142) q[0];
sx q[0];
rz(-1.1508862) q[0];
x q[1];
rz(-2.0306573) q[2];
sx q[2];
rz(-1.0265372) q[2];
sx q[2];
rz(2.0366675) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5729546) q[1];
sx q[1];
rz(-1.8796088) q[1];
sx q[1];
rz(-0.026334865) q[1];
rz(-pi) q[2];
rz(-0.29048357) q[3];
sx q[3];
rz(-0.64290128) q[3];
sx q[3];
rz(3.052352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.48250616) q[2];
sx q[2];
rz(-1.3799474) q[2];
sx q[2];
rz(2.144311) q[2];
rz(-1.4128134) q[3];
sx q[3];
rz(-2.0480859) q[3];
sx q[3];
rz(0.93878448) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
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
rz(2.2361141) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3826686) q[0];
sx q[0];
rz(-1.7610468) q[0];
sx q[0];
rz(-0.24567901) q[0];
rz(-pi) q[1];
rz(-2.682914) q[2];
sx q[2];
rz(-0.80692601) q[2];
sx q[2];
rz(2.5708831) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.29618759) q[1];
sx q[1];
rz(-2.2803109) q[1];
sx q[1];
rz(3.0441277) q[1];
rz(-1.521103) q[3];
sx q[3];
rz(-0.71411413) q[3];
sx q[3];
rz(-0.36220887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3811615) q[2];
sx q[2];
rz(-0.16813777) q[2];
sx q[2];
rz(0.9872438) q[2];
rz(-3.0958214) q[3];
sx q[3];
rz(-1.8493303) q[3];
sx q[3];
rz(-2.9614017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.10953294) q[0];
sx q[0];
rz(-0.68234545) q[0];
sx q[0];
rz(0.12606829) q[0];
rz(-2.9571422) q[1];
sx q[1];
rz(-2.0729063) q[1];
sx q[1];
rz(1.1674081) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084598736) q[0];
sx q[0];
rz(-1.3697249) q[0];
sx q[0];
rz(-2.3720471) q[0];
rz(2.6378176) q[2];
sx q[2];
rz(-1.7283895) q[2];
sx q[2];
rz(2.4368844) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0988401) q[1];
sx q[1];
rz(-0.71430695) q[1];
sx q[1];
rz(-2.7318098) q[1];
rz(-2.4950124) q[3];
sx q[3];
rz(-1.0160035) q[3];
sx q[3];
rz(3.0559029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.8644774) q[2];
sx q[2];
rz(-0.41431674) q[2];
sx q[2];
rz(-2.5100822) q[2];
rz(-0.19550066) q[3];
sx q[3];
rz(-1.86444) q[3];
sx q[3];
rz(0.41803944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
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
rz(-1.0031507) q[1];
sx q[1];
rz(0.46708333) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57468092) q[0];
sx q[0];
rz(-1.7705288) q[0];
sx q[0];
rz(1.9499705) q[0];
x q[1];
rz(0.37093016) q[2];
sx q[2];
rz(-0.61894722) q[2];
sx q[2];
rz(1.3528454) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.621532) q[1];
sx q[1];
rz(-0.48081765) q[1];
sx q[1];
rz(-2.0562999) q[1];
rz(-0.88539447) q[3];
sx q[3];
rz(-1.6939031) q[3];
sx q[3];
rz(2.3791594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6749394) q[2];
sx q[2];
rz(-1.249524) q[2];
sx q[2];
rz(-2.6494027) q[2];
rz(2.8594033) q[3];
sx q[3];
rz(-1.4331093) q[3];
sx q[3];
rz(-1.5757489) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2714587) q[0];
sx q[0];
rz(-2.6014355) q[0];
sx q[0];
rz(-3.0555994) q[0];
rz(-1.4280691) q[1];
sx q[1];
rz(-1.0079039) q[1];
sx q[1];
rz(-1.496398) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8498358) q[0];
sx q[0];
rz(-1.0514326) q[0];
sx q[0];
rz(-0.22626466) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9005269) q[2];
sx q[2];
rz(-1.0450372) q[2];
sx q[2];
rz(-0.17503967) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8346356) q[1];
sx q[1];
rz(-2.7863057) q[1];
sx q[1];
rz(1.6389695) q[1];
rz(1.5580642) q[3];
sx q[3];
rz(-0.85789645) q[3];
sx q[3];
rz(1.3811779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.66203403) q[2];
sx q[2];
rz(-2.5223314) q[2];
sx q[2];
rz(0.39624828) q[2];
rz(0.59404343) q[3];
sx q[3];
rz(-1.0915353) q[3];
sx q[3];
rz(-1.5007639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8783766) q[0];
sx q[0];
rz(-0.57290572) q[0];
sx q[0];
rz(2.0492045) q[0];
rz(1.3573525) q[1];
sx q[1];
rz(-1.7204294) q[1];
sx q[1];
rz(0.011627442) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7250925) q[0];
sx q[0];
rz(-1.9511127) q[0];
sx q[0];
rz(-1.2854544) q[0];
rz(-1.5989223) q[2];
sx q[2];
rz(-1.5962432) q[2];
sx q[2];
rz(2.5487119) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.36794084) q[1];
sx q[1];
rz(-0.546574) q[1];
sx q[1];
rz(-1.6163338) q[1];
x q[2];
rz(0.39832468) q[3];
sx q[3];
rz(-0.71361226) q[3];
sx q[3];
rz(0.45809612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.45550436) q[2];
sx q[2];
rz(-0.65767613) q[2];
sx q[2];
rz(2.8536076) q[2];
rz(1.1503495) q[3];
sx q[3];
rz(-1.3214313) q[3];
sx q[3];
rz(-2.5434125) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8312254) q[0];
sx q[0];
rz(-2.6103525) q[0];
sx q[0];
rz(1.9294552) q[0];
rz(-2.7638226) q[1];
sx q[1];
rz(-1.9394082) q[1];
sx q[1];
rz(-1.6479962) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0069107) q[0];
sx q[0];
rz(-1.3575166) q[0];
sx q[0];
rz(2.8365305) q[0];
rz(-pi) q[1];
x q[1];
rz(0.059842589) q[2];
sx q[2];
rz(-1.4388196) q[2];
sx q[2];
rz(-3.0795385) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8131866) q[1];
sx q[1];
rz(-1.9771736) q[1];
sx q[1];
rz(0.74192585) q[1];
x q[2];
rz(-1.7842403) q[3];
sx q[3];
rz(-2.6297914) q[3];
sx q[3];
rz(-1.7970999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3528072) q[2];
sx q[2];
rz(-2.0264503) q[2];
sx q[2];
rz(-1.8898233) q[2];
rz(-2.2552323) q[3];
sx q[3];
rz(-0.76615196) q[3];
sx q[3];
rz(3.0322976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4733474) q[0];
sx q[0];
rz(-2.1003523) q[0];
sx q[0];
rz(-0.1782724) q[0];
rz(3.0687304) q[1];
sx q[1];
rz(-2.5477414) q[1];
sx q[1];
rz(1.9624306) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3423791) q[0];
sx q[0];
rz(-1.2782492) q[0];
sx q[0];
rz(0.53036687) q[0];
rz(-2.6369008) q[2];
sx q[2];
rz(-1.0070649) q[2];
sx q[2];
rz(1.6187514) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.804783) q[1];
sx q[1];
rz(-2.0840624) q[1];
sx q[1];
rz(-1.8670765) q[1];
rz(-pi) q[2];
rz(2.0719028) q[3];
sx q[3];
rz(-1.2341098) q[3];
sx q[3];
rz(-1.5130373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4420085) q[2];
sx q[2];
rz(-2.1506491) q[2];
sx q[2];
rz(1.0317831) q[2];
rz(1.3423086) q[3];
sx q[3];
rz(-1.4167901) q[3];
sx q[3];
rz(-2.9163196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28829065) q[0];
sx q[0];
rz(-1.0422215) q[0];
sx q[0];
rz(-3.0798966) q[0];
rz(-1.8348947) q[1];
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
rz(2.2022804) q[0];
sx q[0];
rz(-1.5081524) q[0];
sx q[0];
rz(1.5059727) q[0];
x q[1];
rz(-1.3332974) q[2];
sx q[2];
rz(-2.4727159) q[2];
sx q[2];
rz(-2.3046658) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.8577733) q[1];
sx q[1];
rz(-0.76446629) q[1];
sx q[1];
rz(-0.70160265) q[1];
x q[2];
rz(2.0652566) q[3];
sx q[3];
rz(-0.38443243) q[3];
sx q[3];
rz(1.9264551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6178199) q[2];
sx q[2];
rz(-2.4536295) q[2];
sx q[2];
rz(-3.0715959) q[2];
rz(2.2648515) q[3];
sx q[3];
rz(-1.3249967) q[3];
sx q[3];
rz(2.783412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6488279) q[0];
sx q[0];
rz(-1.5179101) q[0];
sx q[0];
rz(-2.5483325) q[0];
rz(1.5169253) q[1];
sx q[1];
rz(-2.0943874) q[1];
sx q[1];
rz(2.692683) q[1];
rz(2.2286227) q[2];
sx q[2];
rz(-1.2428478) q[2];
sx q[2];
rz(-1.7075677) q[2];
rz(-1.6886961) q[3];
sx q[3];
rz(-0.33751925) q[3];
sx q[3];
rz(-2.2681507) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
