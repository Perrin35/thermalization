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
rz(-0.67796081) q[0];
sx q[0];
rz(-0.27422658) q[0];
sx q[0];
rz(1.8507313) q[0];
rz(-3.3554606) q[1];
sx q[1];
rz(5.9670347) q[1];
sx q[1];
rz(11.066758) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16786538) q[0];
sx q[0];
rz(-1.8059547) q[0];
sx q[0];
rz(0.27882378) q[0];
rz(-pi) q[1];
x q[1];
rz(0.29408367) q[2];
sx q[2];
rz(-2.3399647) q[2];
sx q[2];
rz(1.5844567) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5747524) q[1];
sx q[1];
rz(-2.2623203) q[1];
sx q[1];
rz(1.1471466) q[1];
rz(-pi) q[2];
x q[2];
rz(0.95176272) q[3];
sx q[3];
rz(-2.4118399) q[3];
sx q[3];
rz(-2.8918134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9700254) q[2];
sx q[2];
rz(-2.1000523) q[2];
sx q[2];
rz(0.901326) q[2];
rz(2.6775635) q[3];
sx q[3];
rz(-1.8601067) q[3];
sx q[3];
rz(-0.06981167) q[3];
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
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0483911) q[0];
sx q[0];
rz(-0.036660107) q[0];
sx q[0];
rz(1.2777591) q[0];
rz(0.094206421) q[1];
sx q[1];
rz(-2.6779046) q[1];
sx q[1];
rz(-1.5346079) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33714715) q[0];
sx q[0];
rz(-0.60071105) q[0];
sx q[0];
rz(-1.5128193) q[0];
x q[1];
rz(3.0846527) q[2];
sx q[2];
rz(-1.7527662) q[2];
sx q[2];
rz(-1.6873311) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0329735) q[1];
sx q[1];
rz(-1.0796282) q[1];
sx q[1];
rz(-1.2742548) q[1];
rz(-2.3393163) q[3];
sx q[3];
rz(-0.55219383) q[3];
sx q[3];
rz(1.4331762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4216807) q[2];
sx q[2];
rz(-1.2016502) q[2];
sx q[2];
rz(-0.16033944) q[2];
rz(-2.4028589) q[3];
sx q[3];
rz(-2.2952047) q[3];
sx q[3];
rz(1.4275449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5147603) q[0];
sx q[0];
rz(-1.7034096) q[0];
sx q[0];
rz(-0.50977388) q[0];
rz(-1.7104507) q[1];
sx q[1];
rz(-1.1809843) q[1];
sx q[1];
rz(0.74657718) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0261542) q[0];
sx q[0];
rz(-0.55284772) q[0];
sx q[0];
rz(-1.5325969) q[0];
x q[1];
rz(-2.4116729) q[2];
sx q[2];
rz(-2.0509208) q[2];
sx q[2];
rz(-1.4381222) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4164734) q[1];
sx q[1];
rz(-1.9320357) q[1];
sx q[1];
rz(0.84049417) q[1];
rz(-pi) q[2];
rz(-0.38807773) q[3];
sx q[3];
rz(-0.92736926) q[3];
sx q[3];
rz(-2.3449947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9868682) q[2];
sx q[2];
rz(-0.26686033) q[2];
sx q[2];
rz(-1.223684) q[2];
rz(2.0679421) q[3];
sx q[3];
rz(-1.7045538) q[3];
sx q[3];
rz(-0.8980208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8680442) q[0];
sx q[0];
rz(-0.88669625) q[0];
sx q[0];
rz(-2.7622188) q[0];
rz(-0.1768449) q[1];
sx q[1];
rz(-1.6601446) q[1];
sx q[1];
rz(0.79536974) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0220227) q[0];
sx q[0];
rz(-1.2031462) q[0];
sx q[0];
rz(3.0421542) q[0];
rz(-pi) q[1];
rz(-2.3248939) q[2];
sx q[2];
rz(-1.7672605) q[2];
sx q[2];
rz(-2.8981371) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8898192) q[1];
sx q[1];
rz(-1.3135513) q[1];
sx q[1];
rz(-2.1104269) q[1];
rz(-pi) q[2];
rz(2.7164812) q[3];
sx q[3];
rz(-0.13775857) q[3];
sx q[3];
rz(2.3099096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7103601) q[2];
sx q[2];
rz(-1.7962339) q[2];
sx q[2];
rz(1.7809407) q[2];
rz(1.0294754) q[3];
sx q[3];
rz(-1.4808713) q[3];
sx q[3];
rz(-1.8195389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65115702) q[0];
sx q[0];
rz(-1.5999726) q[0];
sx q[0];
rz(1.5486451) q[0];
rz(-2.7410638) q[1];
sx q[1];
rz(-1.488204) q[1];
sx q[1];
rz(-1.4097479) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9072394) q[0];
sx q[0];
rz(-0.84915224) q[0];
sx q[0];
rz(-1.3035151) q[0];
rz(-2.7821252) q[2];
sx q[2];
rz(-1.4233982) q[2];
sx q[2];
rz(1.1254252) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5286136) q[1];
sx q[1];
rz(-2.3720212) q[1];
sx q[1];
rz(0.90958325) q[1];
rz(-pi) q[2];
rz(-1.665776) q[3];
sx q[3];
rz(-0.62415571) q[3];
sx q[3];
rz(2.399866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3132402) q[2];
sx q[2];
rz(-1.7007549) q[2];
sx q[2];
rz(-2.7102615) q[2];
rz(-3.0865772) q[3];
sx q[3];
rz(-0.31156817) q[3];
sx q[3];
rz(-2.5855248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7677652) q[0];
sx q[0];
rz(-2.8887833) q[0];
sx q[0];
rz(-2.1164236) q[0];
rz(-0.45285666) q[1];
sx q[1];
rz(-0.57944524) q[1];
sx q[1];
rz(-2.2377009) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9867803) q[0];
sx q[0];
rz(-0.93789369) q[0];
sx q[0];
rz(0.38582071) q[0];
rz(-pi) q[1];
rz(-0.91458648) q[2];
sx q[2];
rz(-1.0433955) q[2];
sx q[2];
rz(1.5114776) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3749668) q[1];
sx q[1];
rz(-1.300525) q[1];
sx q[1];
rz(2.388776) q[1];
x q[2];
rz(-0.17136161) q[3];
sx q[3];
rz(-1.7336066) q[3];
sx q[3];
rz(-2.5622615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8338833) q[2];
sx q[2];
rz(-1.4557975) q[2];
sx q[2];
rz(-2.9252388) q[2];
rz(0.89573914) q[3];
sx q[3];
rz(-0.68550617) q[3];
sx q[3];
rz(1.640813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1322587) q[0];
sx q[0];
rz(-2.0348771) q[0];
sx q[0];
rz(2.4427781) q[0];
rz(1.1837333) q[1];
sx q[1];
rz(-1.7203169) q[1];
sx q[1];
rz(2.2898477) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18139938) q[0];
sx q[0];
rz(-1.6087247) q[0];
sx q[0];
rz(1.0344124) q[0];
x q[1];
rz(0.57735195) q[2];
sx q[2];
rz(-1.5923579) q[2];
sx q[2];
rz(2.2039764) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7523562) q[1];
sx q[1];
rz(-0.79006486) q[1];
sx q[1];
rz(0.52100928) q[1];
x q[2];
rz(0.61148879) q[3];
sx q[3];
rz(-0.73966714) q[3];
sx q[3];
rz(-2.8304493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.8375887) q[2];
sx q[2];
rz(-1.8303266) q[2];
sx q[2];
rz(-0.50789976) q[2];
rz(-1.0287644) q[3];
sx q[3];
rz(-2.2977836) q[3];
sx q[3];
rz(2.540551) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4429338) q[0];
sx q[0];
rz(-1.826094) q[0];
sx q[0];
rz(-0.0096631924) q[0];
rz(1.5254321) q[1];
sx q[1];
rz(-1.7639672) q[1];
sx q[1];
rz(-2.756871) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98823901) q[0];
sx q[0];
rz(-1.8073449) q[0];
sx q[0];
rz(-2.9896215) q[0];
rz(-pi) q[1];
rz(-1.184395) q[2];
sx q[2];
rz(-1.3829872) q[2];
sx q[2];
rz(-1.1593429) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2710451) q[1];
sx q[1];
rz(-1.3636075) q[1];
sx q[1];
rz(-1.8533587) q[1];
rz(-pi) q[2];
x q[2];
rz(2.557005) q[3];
sx q[3];
rz(-0.86170022) q[3];
sx q[3];
rz(-0.60230909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.559451) q[2];
sx q[2];
rz(-2.1820575) q[2];
sx q[2];
rz(0.01586308) q[2];
rz(1.2265685) q[3];
sx q[3];
rz(-0.80569402) q[3];
sx q[3];
rz(-0.62172833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3847619) q[0];
sx q[0];
rz(-2.4767196) q[0];
sx q[0];
rz(2.050198) q[0];
rz(0.81820828) q[1];
sx q[1];
rz(-1.810377) q[1];
sx q[1];
rz(0.6689201) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92907897) q[0];
sx q[0];
rz(-2.2336279) q[0];
sx q[0];
rz(0.79848358) q[0];
rz(-pi) q[1];
rz(-2.0537655) q[2];
sx q[2];
rz(-2.3685799) q[2];
sx q[2];
rz(2.9376466) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0099677) q[1];
sx q[1];
rz(-2.5796081) q[1];
sx q[1];
rz(1.704292) q[1];
rz(-2.1534641) q[3];
sx q[3];
rz(-1.8496397) q[3];
sx q[3];
rz(-3.0089965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6341256) q[2];
sx q[2];
rz(-0.52906817) q[2];
sx q[2];
rz(-2.1779306) q[2];
rz(-1.4108747) q[3];
sx q[3];
rz(-1.4286634) q[3];
sx q[3];
rz(1.7760407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1249579) q[0];
sx q[0];
rz(-0.5683012) q[0];
sx q[0];
rz(1.4547263) q[0];
rz(-1.1085054) q[1];
sx q[1];
rz(-1.6051555) q[1];
sx q[1];
rz(0.32807168) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5755558) q[0];
sx q[0];
rz(-1.6854316) q[0];
sx q[0];
rz(0.057997142) q[0];
x q[1];
rz(-2.9607707) q[2];
sx q[2];
rz(-1.4563113) q[2];
sx q[2];
rz(2.3608728) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0340424) q[1];
sx q[1];
rz(-0.75354105) q[1];
sx q[1];
rz(0.57489245) q[1];
x q[2];
rz(-1.7441196) q[3];
sx q[3];
rz(-1.4651235) q[3];
sx q[3];
rz(2.8901576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6803117) q[2];
sx q[2];
rz(-1.842247) q[2];
sx q[2];
rz(2.7247562) q[2];
rz(-2.4428115) q[3];
sx q[3];
rz(-2.4650033) q[3];
sx q[3];
rz(0.39066395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1998491) q[0];
sx q[0];
rz(-1.2171634) q[0];
sx q[0];
rz(-1.4456277) q[0];
rz(0.70882123) q[1];
sx q[1];
rz(-0.91239057) q[1];
sx q[1];
rz(-0.053587996) q[1];
rz(1.3016635) q[2];
sx q[2];
rz(-2.4872645) q[2];
sx q[2];
rz(-0.80254868) q[2];
rz(2.3117456) q[3];
sx q[3];
rz(-1.1329069) q[3];
sx q[3];
rz(-2.2923242) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
