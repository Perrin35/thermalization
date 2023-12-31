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
rz(-0.34163707) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61603123) q[0];
sx q[0];
rz(-1.6934868) q[0];
sx q[0];
rz(-0.37567715) q[0];
rz(2.5506637) q[2];
sx q[2];
rz(-1.8701631) q[2];
sx q[2];
rz(-0.17609827) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.29400723) q[1];
sx q[1];
rz(-0.60677401) q[1];
sx q[1];
rz(-1.0311544) q[1];
rz(-pi) q[2];
rz(-2.0658675) q[3];
sx q[3];
rz(-0.44701156) q[3];
sx q[3];
rz(-1.8398374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.59387702) q[2];
sx q[2];
rz(-1.1527858) q[2];
sx q[2];
rz(0.95648742) q[2];
rz(2.9603413) q[3];
sx q[3];
rz(-0.58877188) q[3];
sx q[3];
rz(-1.6199002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0797121) q[0];
sx q[0];
rz(-0.065608874) q[0];
sx q[0];
rz(-1.99615) q[0];
rz(-2.0727797) q[1];
sx q[1];
rz(-0.83199465) q[1];
sx q[1];
rz(0.72584814) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2929935) q[0];
sx q[0];
rz(-1.2147875) q[0];
sx q[0];
rz(-2.5550935) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1109353) q[2];
sx q[2];
rz(-1.0265372) q[2];
sx q[2];
rz(1.1049251) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4821856) q[1];
sx q[1];
rz(-2.8316951) q[1];
sx q[1];
rz(-1.4884429) q[1];
rz(2.8511091) q[3];
sx q[3];
rz(-0.64290128) q[3];
sx q[3];
rz(-0.089240616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6590865) q[2];
sx q[2];
rz(-1.7616452) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
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
rz(-1.1285271) q[0];
rz(-1.3035125) q[1];
sx q[1];
rz(-2.4120657) q[1];
sx q[1];
rz(-0.9054786) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3072976) q[0];
sx q[0];
rz(-2.8320438) q[0];
sx q[0];
rz(-2.4718667) q[0];
rz(-pi) q[1];
rz(2.682914) q[2];
sx q[2];
rz(-2.3346666) q[2];
sx q[2];
rz(2.5708831) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9943774) q[1];
sx q[1];
rz(-0.71502393) q[1];
sx q[1];
rz(1.4579525) q[1];
rz(-2.2842992) q[3];
sx q[3];
rz(-1.6033353) q[3];
sx q[3];
rz(-1.9705704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0320597) q[0];
sx q[0];
rz(-0.68234545) q[0];
sx q[0];
rz(-3.0155244) q[0];
rz(0.18445045) q[1];
sx q[1];
rz(-1.0686864) q[1];
sx q[1];
rz(-1.1674081) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6895034) q[0];
sx q[0];
rz(-0.79012442) q[0];
sx q[0];
rz(2.8566314) q[0];
x q[1];
rz(-2.8235769) q[2];
sx q[2];
rz(-2.615775) q[2];
sx q[2];
rz(-1.1434681) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3525225) q[1];
sx q[1];
rz(-1.8348502) q[1];
sx q[1];
rz(-2.4697484) q[1];
x q[2];
rz(-0.79951841) q[3];
sx q[3];
rz(-2.3164253) q[3];
sx q[3];
rz(2.094401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.8644774) q[2];
sx q[2];
rz(-2.7272759) q[2];
sx q[2];
rz(2.5100822) q[2];
rz(2.946092) q[3];
sx q[3];
rz(-1.86444) q[3];
sx q[3];
rz(-2.7235532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9773848) q[0];
sx q[0];
rz(-0.13810869) q[0];
sx q[0];
rz(-0.50267977) q[0];
rz(2.0282822) q[1];
sx q[1];
rz(-2.138442) q[1];
sx q[1];
rz(-2.6745093) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0750077) q[0];
sx q[0];
rz(-1.9420615) q[0];
sx q[0];
rz(2.9270372) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5555243) q[2];
sx q[2];
rz(-1.7826826) q[2];
sx q[2];
rz(0.52473247) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0847229) q[1];
sx q[1];
rz(-1.9922171) q[1];
sx q[1];
rz(0.23878581) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.763837) q[3];
sx q[3];
rz(-2.4469864) q[3];
sx q[3];
rz(2.4823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.46665329) q[2];
sx q[2];
rz(-1.249524) q[2];
sx q[2];
rz(-0.49218991) q[2];
rz(-2.8594033) q[3];
sx q[3];
rz(-1.7084833) q[3];
sx q[3];
rz(1.5658437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.870134) q[0];
sx q[0];
rz(-0.54015714) q[0];
sx q[0];
rz(-3.0555994) q[0];
rz(-1.7135235) q[1];
sx q[1];
rz(-1.0079039) q[1];
sx q[1];
rz(1.496398) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7487885) q[0];
sx q[0];
rz(-1.3747871) q[0];
sx q[0];
rz(-2.1013573) q[0];
x q[1];
rz(0.24106579) q[2];
sx q[2];
rz(-2.0965555) q[2];
sx q[2];
rz(2.966553) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9073346) q[1];
sx q[1];
rz(-1.2163711) q[1];
sx q[1];
rz(3.1163232) q[1];
rz(0.014724894) q[3];
sx q[3];
rz(-2.428599) q[3];
sx q[3];
rz(1.7409489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.66203403) q[2];
sx q[2];
rz(-2.5223314) q[2];
sx q[2];
rz(2.7453444) q[2];
rz(-0.59404343) q[3];
sx q[3];
rz(-1.0915353) q[3];
sx q[3];
rz(-1.6408287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2632161) q[0];
sx q[0];
rz(-0.57290572) q[0];
sx q[0];
rz(2.0492045) q[0];
rz(-1.7842402) q[1];
sx q[1];
rz(-1.4211632) q[1];
sx q[1];
rz(-0.011627442) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8893338) q[0];
sx q[0];
rz(-2.6703435) q[0];
sx q[0];
rz(0.61347368) q[0];
x q[1];
rz(2.3064012) q[2];
sx q[2];
rz(-0.037926849) q[2];
sx q[2];
rz(1.7131625) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1639451) q[1];
sx q[1];
rz(-1.5944591) q[1];
sx q[1];
rz(-1.0246828) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2467975) q[3];
sx q[3];
rz(-2.218459) q[3];
sx q[3];
rz(-3.0917633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6860883) q[2];
sx q[2];
rz(-2.4839165) q[2];
sx q[2];
rz(2.8536076) q[2];
rz(-1.1503495) q[3];
sx q[3];
rz(-1.3214313) q[3];
sx q[3];
rz(2.5434125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-1.9394082) q[1];
sx q[1];
rz(1.6479962) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50266788) q[0];
sx q[0];
rz(-1.272861) q[0];
sx q[0];
rz(-1.3475247) q[0];
rz(-pi) q[1];
rz(1.9940894) q[2];
sx q[2];
rz(-0.14483843) q[2];
sx q[2];
rz(0.3651948) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6495325) q[1];
sx q[1];
rz(-2.3146555) q[1];
sx q[1];
rz(2.5745113) q[1];
rz(-1.7842403) q[3];
sx q[3];
rz(-0.5118013) q[3];
sx q[3];
rz(-1.3444927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3528072) q[2];
sx q[2];
rz(-2.0264503) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66824526) q[0];
sx q[0];
rz(-1.0412403) q[0];
sx q[0];
rz(0.1782724) q[0];
rz(-0.072862236) q[1];
sx q[1];
rz(-2.5477414) q[1];
sx q[1];
rz(1.9624306) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68547738) q[0];
sx q[0];
rz(-2.5427186) q[0];
sx q[0];
rz(-2.6045538) q[0];
rz(-pi) q[1];
rz(-2.6369008) q[2];
sx q[2];
rz(-1.0070649) q[2];
sx q[2];
rz(-1.5228412) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2485679) q[1];
sx q[1];
rz(-2.5556459) q[1];
sx q[1];
rz(-2.663661) q[1];
rz(-pi) q[2];
rz(-2.7618802) q[3];
sx q[3];
rz(-1.100193) q[3];
sx q[3];
rz(-3.0203366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4420085) q[2];
sx q[2];
rz(-0.99094355) q[2];
sx q[2];
rz(-2.1098095) q[2];
rz(1.7992841) q[3];
sx q[3];
rz(-1.7248025) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.853302) q[0];
sx q[0];
rz(-2.0993711) q[0];
sx q[0];
rz(3.0798966) q[0];
rz(-1.8348947) q[1];
sx q[1];
rz(-1.6625762) q[1];
sx q[1];
rz(2.0526989) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63554791) q[0];
sx q[0];
rz(-1.6354927) q[0];
sx q[0];
rz(-0.062775469) q[0];
rz(0.18386545) q[2];
sx q[2];
rz(-2.2176748) q[2];
sx q[2];
rz(2.0053787) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8809001) q[1];
sx q[1];
rz(-2.0339194) q[1];
sx q[1];
rz(-0.63219627) q[1];
rz(1.9128996) q[3];
sx q[3];
rz(-1.7497239) q[3];
sx q[3];
rz(0.10781328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6178199) q[2];
sx q[2];
rz(-2.4536295) q[2];
sx q[2];
rz(-0.069996746) q[2];
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
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
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
rz(-2.2286227) q[2];
sx q[2];
rz(-1.8987449) q[2];
sx q[2];
rz(1.434025) q[2];
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
