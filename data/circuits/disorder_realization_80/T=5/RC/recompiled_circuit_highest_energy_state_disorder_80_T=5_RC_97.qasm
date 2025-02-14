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
rz(5.0168283) q[1];
sx q[1];
rz(14.556806) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53357994) q[0];
sx q[0];
rz(-1.5927497) q[0];
sx q[0];
rz(1.6589462) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5106836) q[2];
sx q[2];
rz(-2.563347) q[2];
sx q[2];
rz(1.9672245) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.77598876) q[1];
sx q[1];
rz(-1.1465441) q[1];
sx q[1];
rz(-1.7107023) q[1];
x q[2];
rz(2.092716) q[3];
sx q[3];
rz(-1.3112231) q[3];
sx q[3];
rz(-2.5775016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2426131) q[2];
sx q[2];
rz(-2.2448764) q[2];
sx q[2];
rz(-2.8420281) q[2];
rz(-2.7859531) q[3];
sx q[3];
rz(-2.1483597) q[3];
sx q[3];
rz(0.43660823) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75293175) q[0];
sx q[0];
rz(-0.54838538) q[0];
sx q[0];
rz(-2.8652628) q[0];
rz(0.07946864) q[1];
sx q[1];
rz(-0.53229585) q[1];
sx q[1];
rz(2.6452549) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4547098) q[0];
sx q[0];
rz(-1.0609896) q[0];
sx q[0];
rz(2.8532993) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4935819) q[2];
sx q[2];
rz(-1.9081915) q[2];
sx q[2];
rz(1.5959306) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.53238324) q[1];
sx q[1];
rz(-2.2953798) q[1];
sx q[1];
rz(2.0128472) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4088936) q[3];
sx q[3];
rz(-1.5739958) q[3];
sx q[3];
rz(2.7174983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.77208272) q[2];
sx q[2];
rz(-1.3630867) q[2];
sx q[2];
rz(-0.17862865) q[2];
rz(-1.6102839) q[3];
sx q[3];
rz(-2.2154) q[3];
sx q[3];
rz(2.3791651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71132201) q[0];
sx q[0];
rz(-2.0549759) q[0];
sx q[0];
rz(-1.4104728) q[0];
rz(1.7544282) q[1];
sx q[1];
rz(-1.6491363) q[1];
sx q[1];
rz(-1.656146) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94233275) q[0];
sx q[0];
rz(-0.6845277) q[0];
sx q[0];
rz(0.051816414) q[0];
rz(2.7316667) q[2];
sx q[2];
rz(-2.2927444) q[2];
sx q[2];
rz(-0.97463911) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.3988603) q[1];
sx q[1];
rz(-1.9963063) q[1];
sx q[1];
rz(0.33136435) q[1];
x q[2];
rz(-0.90624202) q[3];
sx q[3];
rz(-1.2782879) q[3];
sx q[3];
rz(-0.14636235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9401271) q[2];
sx q[2];
rz(-1.0041693) q[2];
sx q[2];
rz(-0.49261576) q[2];
rz(-0.41871873) q[3];
sx q[3];
rz(-1.0227572) q[3];
sx q[3];
rz(0.59188265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1241322) q[0];
sx q[0];
rz(-5/(13*pi)) q[0];
sx q[0];
rz(-2.9435352) q[0];
rz(0.65912229) q[1];
sx q[1];
rz(-2.5810869) q[1];
sx q[1];
rz(2.9445599) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7871088) q[0];
sx q[0];
rz(-1.0403883) q[0];
sx q[0];
rz(-1.6097989) q[0];
x q[1];
rz(2.2617499) q[2];
sx q[2];
rz(-2.7759399) q[2];
sx q[2];
rz(-1.3297538) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2391165) q[1];
sx q[1];
rz(-0.65985876) q[1];
sx q[1];
rz(2.541594) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5102618) q[3];
sx q[3];
rz(-0.54577561) q[3];
sx q[3];
rz(-0.4402524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3531602) q[2];
sx q[2];
rz(-1.3015231) q[2];
sx q[2];
rz(1.8678467) q[2];
rz(1.5431917) q[3];
sx q[3];
rz(-1.9485995) q[3];
sx q[3];
rz(2.8852692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0850328) q[0];
sx q[0];
rz(-1.4075449) q[0];
sx q[0];
rz(0.10277596) q[0];
rz(3.00434) q[1];
sx q[1];
rz(-2.6922701) q[1];
sx q[1];
rz(1.7030254) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5202839) q[0];
sx q[0];
rz(-2.1946593) q[0];
sx q[0];
rz(0.72786967) q[0];
rz(0.091890865) q[2];
sx q[2];
rz(-0.46481195) q[2];
sx q[2];
rz(-1.6824695) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.66617353) q[1];
sx q[1];
rz(-1.8295145) q[1];
sx q[1];
rz(2.9889876) q[1];
rz(-pi) q[2];
rz(-1.4588615) q[3];
sx q[3];
rz(-1.4482968) q[3];
sx q[3];
rz(-1.4089597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.30456257) q[2];
sx q[2];
rz(-1.6670767) q[2];
sx q[2];
rz(-2.4414818) q[2];
rz(-2.2434798) q[3];
sx q[3];
rz(-0.63106314) q[3];
sx q[3];
rz(-1.496544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
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
rz(2.9333618) q[0];
sx q[0];
rz(-2.5550186) q[0];
sx q[0];
rz(-1.5775648) q[0];
rz(-0.6598407) q[1];
sx q[1];
rz(-2.3520062) q[1];
sx q[1];
rz(-0.97445625) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81466573) q[0];
sx q[0];
rz(-1.5788933) q[0];
sx q[0];
rz(3.1304703) q[0];
x q[1];
rz(0.48759385) q[2];
sx q[2];
rz(-1.7238443) q[2];
sx q[2];
rz(0.4401373) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.32459697) q[1];
sx q[1];
rz(-1.5600509) q[1];
sx q[1];
rz(-2.6399977) q[1];
rz(-pi) q[2];
rz(3.0401214) q[3];
sx q[3];
rz(-0.54059282) q[3];
sx q[3];
rz(-1.056788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4727891) q[2];
sx q[2];
rz(-0.69837022) q[2];
sx q[2];
rz(-2.2051956) q[2];
rz(2.2514553) q[3];
sx q[3];
rz(-1.3606768) q[3];
sx q[3];
rz(-2.9357125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1393046) q[0];
sx q[0];
rz(-1.6253259) q[0];
sx q[0];
rz(1.3521347) q[0];
x q[1];
rz(-2.4981047) q[2];
sx q[2];
rz(-1.7805791) q[2];
sx q[2];
rz(1.6580559) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2568367) q[1];
sx q[1];
rz(-2.4236226) q[1];
sx q[1];
rz(-1.4713478) q[1];
rz(-2.9437441) q[3];
sx q[3];
rz(-1.6258844) q[3];
sx q[3];
rz(-0.69869374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6156442) q[2];
sx q[2];
rz(-1.7082461) q[2];
sx q[2];
rz(-2.136266) q[2];
rz(-2.2510236) q[3];
sx q[3];
rz(-1.4054479) q[3];
sx q[3];
rz(0.65517202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6152076) q[0];
sx q[0];
rz(-2.7387775) q[0];
sx q[0];
rz(-1.1108387) q[0];
rz(0.088518294) q[1];
sx q[1];
rz(-0.42461494) q[1];
sx q[1];
rz(-1.7875338) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.272404) q[0];
sx q[0];
rz(-1.6834295) q[0];
sx q[0];
rz(2.3478481) q[0];
rz(-1.2983923) q[2];
sx q[2];
rz(-1.7188311) q[2];
sx q[2];
rz(-1.4349951) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.718487) q[1];
sx q[1];
rz(-2.8837648) q[1];
sx q[1];
rz(0.30155525) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.76137162) q[3];
sx q[3];
rz(-1.1283481) q[3];
sx q[3];
rz(-0.82514742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2985349) q[2];
sx q[2];
rz(-0.83710805) q[2];
sx q[2];
rz(-2.7561772) q[2];
rz(-2.391583) q[3];
sx q[3];
rz(-2.1842726) q[3];
sx q[3];
rz(-3.0756522) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0176004) q[0];
sx q[0];
rz(-1.0753068) q[0];
sx q[0];
rz(-0.83220926) q[0];
rz(-1.0819134) q[1];
sx q[1];
rz(-2.3257906) q[1];
sx q[1];
rz(-1.6197416) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2913301) q[0];
sx q[0];
rz(-1.1801774) q[0];
sx q[0];
rz(1.0608835) q[0];
rz(0.25427962) q[2];
sx q[2];
rz(-1.2966208) q[2];
sx q[2];
rz(2.5667896) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.202692) q[1];
sx q[1];
rz(-1.0098411) q[1];
sx q[1];
rz(-0.86679082) q[1];
x q[2];
rz(1.8796104) q[3];
sx q[3];
rz(-2.4612456) q[3];
sx q[3];
rz(1.4226567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1443783) q[2];
sx q[2];
rz(-1.8961366) q[2];
sx q[2];
rz(1.2388371) q[2];
rz(1.291412) q[3];
sx q[3];
rz(-1.2715128) q[3];
sx q[3];
rz(3.0102357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.130403) q[0];
sx q[0];
rz(-0.953453) q[0];
sx q[0];
rz(2.3977872) q[0];
rz(-0.047608308) q[1];
sx q[1];
rz(-2.036939) q[1];
sx q[1];
rz(-2.0000752) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6358546) q[0];
sx q[0];
rz(-0.72162745) q[0];
sx q[0];
rz(0.82765915) q[0];
rz(-0.6170527) q[2];
sx q[2];
rz(-1.8931555) q[2];
sx q[2];
rz(2.6094391) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9824578) q[1];
sx q[1];
rz(-1.1679113) q[1];
sx q[1];
rz(-2.3680658) q[1];
rz(-pi) q[2];
rz(-1.7520245) q[3];
sx q[3];
rz(-0.93398213) q[3];
sx q[3];
rz(2.2006478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6514674) q[2];
sx q[2];
rz(-0.23568025) q[2];
sx q[2];
rz(0.30533314) q[2];
rz(-1.7134282) q[3];
sx q[3];
rz(-1.3729264) q[3];
sx q[3];
rz(-1.8491687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74844985) q[0];
sx q[0];
rz(-0.23101692) q[0];
sx q[0];
rz(1.1363181) q[0];
rz(-2.582386) q[1];
sx q[1];
rz(-1.8103841) q[1];
sx q[1];
rz(0.24787535) q[1];
rz(-0.43287533) q[2];
sx q[2];
rz(-2.4268711) q[2];
sx q[2];
rz(-2.101244) q[2];
rz(2.2918626) q[3];
sx q[3];
rz(-2.9624675) q[3];
sx q[3];
rz(-0.20476838) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
