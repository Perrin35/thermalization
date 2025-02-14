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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6080127) q[0];
sx q[0];
rz(-1.5488429) q[0];
sx q[0];
rz(1.4826464) q[0];
x q[1];
rz(-0.63090905) q[2];
sx q[2];
rz(-0.5782457) q[2];
sx q[2];
rz(1.9672245) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.44635998) q[1];
sx q[1];
rz(-0.44539136) q[1];
sx q[1];
rz(-2.8421418) q[1];
rz(-2.0602442) q[3];
sx q[3];
rz(-2.5640805) q[3];
sx q[3];
rz(-0.58693991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.89897951) q[2];
sx q[2];
rz(-2.2448764) q[2];
sx q[2];
rz(-0.29956451) q[2];
rz(0.35563955) q[3];
sx q[3];
rz(-2.1483597) q[3];
sx q[3];
rz(-2.7049844) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3886609) q[0];
sx q[0];
rz(-2.5932073) q[0];
sx q[0];
rz(-2.8652628) q[0];
rz(3.062124) q[1];
sx q[1];
rz(-0.53229585) q[1];
sx q[1];
rz(-2.6452549) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68688289) q[0];
sx q[0];
rz(-2.080603) q[0];
sx q[0];
rz(2.8532993) q[0];
rz(-0.33832835) q[2];
sx q[2];
rz(-1.4979432) q[2];
sx q[2];
rz(-0.05073994) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6092094) q[1];
sx q[1];
rz(-0.84621284) q[1];
sx q[1];
rz(-1.1287454) q[1];
rz(-pi) q[2];
rz(-1.4088936) q[3];
sx q[3];
rz(-1.5675968) q[3];
sx q[3];
rz(2.7174983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3695099) q[2];
sx q[2];
rz(-1.778506) q[2];
sx q[2];
rz(-2.962964) q[2];
rz(-1.6102839) q[3];
sx q[3];
rz(-0.9261927) q[3];
sx q[3];
rz(0.76242751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71132201) q[0];
sx q[0];
rz(-1.0866168) q[0];
sx q[0];
rz(1.7311199) q[0];
rz(-1.3871644) q[1];
sx q[1];
rz(-1.6491363) q[1];
sx q[1];
rz(-1.656146) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5532865) q[0];
sx q[0];
rz(-1.5380412) q[0];
sx q[0];
rz(0.68386987) q[0];
rz(2.335821) q[2];
sx q[2];
rz(-1.2670332) q[2];
sx q[2];
rz(2.8250776) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8456257) q[1];
sx q[1];
rz(-0.53301552) q[1];
sx q[1];
rz(0.94818268) q[1];
x q[2];
rz(2.0250506) q[3];
sx q[3];
rz(-0.71703934) q[3];
sx q[3];
rz(-1.3644791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2014655) q[2];
sx q[2];
rz(-1.0041693) q[2];
sx q[2];
rz(0.49261576) q[2];
rz(-0.41871873) q[3];
sx q[3];
rz(-2.1188354) q[3];
sx q[3];
rz(-0.59188265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
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
rz(0.19805743) q[0];
rz(-2.4824704) q[1];
sx q[1];
rz(-0.56050572) q[1];
sx q[1];
rz(-2.9445599) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27749946) q[0];
sx q[0];
rz(-0.53170337) q[0];
sx q[0];
rz(-0.066389991) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23931673) q[2];
sx q[2];
rz(-1.2916358) q[2];
sx q[2];
rz(-2.5365732) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9777364) q[1];
sx q[1];
rz(-1.2173554) q[1];
sx q[1];
rz(2.5720235) q[1];
rz(-pi) q[2];
rz(2.1157589) q[3];
sx q[3];
rz(-1.539388) q[3];
sx q[3];
rz(-1.1823013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3531602) q[2];
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
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0565599) q[0];
sx q[0];
rz(-1.4075449) q[0];
sx q[0];
rz(3.0388167) q[0];
rz(3.00434) q[1];
sx q[1];
rz(-2.6922701) q[1];
sx q[1];
rz(-1.4385673) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.611821) q[0];
sx q[0];
rz(-0.91982929) q[0];
sx q[0];
rz(0.82470993) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.46311997) q[2];
sx q[2];
rz(-1.5296522) q[2];
sx q[2];
rz(0.1938614) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2074429) q[1];
sx q[1];
rz(-0.29950073) q[1];
sx q[1];
rz(-2.092176) q[1];
rz(-pi) q[2];
rz(-3.0183295) q[3];
sx q[3];
rz(-1.6818889) q[3];
sx q[3];
rz(-2.9934903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8370301) q[2];
sx q[2];
rz(-1.6670767) q[2];
sx q[2];
rz(-0.70011082) q[2];
rz(0.89811283) q[3];
sx q[3];
rz(-0.63106314) q[3];
sx q[3];
rz(-1.496544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9333618) q[0];
sx q[0];
rz(-2.5550186) q[0];
sx q[0];
rz(-1.5640278) q[0];
rz(-2.481752) q[1];
sx q[1];
rz(-0.78958646) q[1];
sx q[1];
rz(2.1671364) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.385372) q[0];
sx q[0];
rz(-1.5819183) q[0];
sx q[0];
rz(1.5626989) q[0];
rz(-pi) q[1];
rz(-1.3979369) q[2];
sx q[2];
rz(-2.0522016) q[2];
sx q[2];
rz(1.0499894) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.240307) q[1];
sx q[1];
rz(-1.069233) q[1];
sx q[1];
rz(1.5585414) q[1];
x q[2];
rz(3.0401214) q[3];
sx q[3];
rz(-2.6009998) q[3];
sx q[3];
rz(1.056788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.66880354) q[2];
sx q[2];
rz(-0.69837022) q[2];
sx q[2];
rz(2.2051956) q[2];
rz(0.89013731) q[3];
sx q[3];
rz(-1.3606768) q[3];
sx q[3];
rz(-0.20588017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5278006) q[0];
sx q[0];
rz(-2.9449936) q[0];
sx q[0];
rz(-0.11372456) q[0];
rz(-1.3971036) q[1];
sx q[1];
rz(-1.0662268) q[1];
sx q[1];
rz(3.1166335) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1393046) q[0];
sx q[0];
rz(-1.6253259) q[0];
sx q[0];
rz(1.789458) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4981047) q[2];
sx q[2];
rz(-1.7805791) q[2];
sx q[2];
rz(1.6580559) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2568367) q[1];
sx q[1];
rz(-0.71797007) q[1];
sx q[1];
rz(-1.6702449) q[1];
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
rz(-2.6156442) q[2];
sx q[2];
rz(-1.4333466) q[2];
sx q[2];
rz(-1.0053267) q[2];
rz(0.89056906) q[3];
sx q[3];
rz(-1.4054479) q[3];
sx q[3];
rz(-2.4864206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5263851) q[0];
sx q[0];
rz(-0.40281519) q[0];
sx q[0];
rz(-2.030754) q[0];
rz(0.088518294) q[1];
sx q[1];
rz(-0.42461494) q[1];
sx q[1];
rz(-1.7875338) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8691886) q[0];
sx q[0];
rz(-1.4581632) q[0];
sx q[0];
rz(-0.79374452) q[0];
x q[1];
rz(-2.9879775) q[2];
sx q[2];
rz(-1.8401463) q[2];
sx q[2];
rz(2.9646123) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4231056) q[1];
sx q[1];
rz(-2.8837648) q[1];
sx q[1];
rz(-2.8400374) q[1];
rz(0.99127524) q[3];
sx q[3];
rz(-2.2439828) q[3];
sx q[3];
rz(-2.0085101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.84305772) q[2];
sx q[2];
rz(-2.3044846) q[2];
sx q[2];
rz(2.7561772) q[2];
rz(-2.391583) q[3];
sx q[3];
rz(-2.1842726) q[3];
sx q[3];
rz(-3.0756522) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1239922) q[0];
sx q[0];
rz(-2.0662859) q[0];
sx q[0];
rz(2.3093834) q[0];
rz(1.0819134) q[1];
sx q[1];
rz(-2.3257906) q[1];
sx q[1];
rz(1.6197416) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48927096) q[0];
sx q[0];
rz(-2.0390563) q[0];
sx q[0];
rz(-2.7007607) q[0];
x q[1];
rz(-1.85361) q[2];
sx q[2];
rz(-1.3262124) q[2];
sx q[2];
rz(-2.2158538) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.950022) q[1];
sx q[1];
rz(-2.2723928) q[1];
sx q[1];
rz(0.8002886) q[1];
rz(-0.91401871) q[3];
sx q[3];
rz(-1.3784215) q[3];
sx q[3];
rz(3.0466444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.99721432) q[2];
sx q[2];
rz(-1.245456) q[2];
sx q[2];
rz(-1.2388371) q[2];
rz(-1.291412) q[3];
sx q[3];
rz(-1.8700799) q[3];
sx q[3];
rz(3.0102357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.130403) q[0];
sx q[0];
rz(-2.1881396) q[0];
sx q[0];
rz(0.74380547) q[0];
rz(0.047608308) q[1];
sx q[1];
rz(-2.036939) q[1];
sx q[1];
rz(-1.1415175) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6358546) q[0];
sx q[0];
rz(-2.4199652) q[0];
sx q[0];
rz(-0.82765915) q[0];
x q[1];
rz(0.52351482) q[2];
sx q[2];
rz(-0.68632978) q[2];
sx q[2];
rz(1.4586142) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3642753) q[1];
sx q[1];
rz(-2.2688443) q[1];
sx q[1];
rz(-1.0335405) q[1];
rz(-1.7520245) q[3];
sx q[3];
rz(-2.2076105) q[3];
sx q[3];
rz(0.94094482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6514674) q[2];
sx q[2];
rz(-2.9059124) q[2];
sx q[2];
rz(-2.8362595) q[2];
rz(-1.4281645) q[3];
sx q[3];
rz(-1.7686663) q[3];
sx q[3];
rz(1.292424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74844985) q[0];
sx q[0];
rz(-0.23101692) q[0];
sx q[0];
rz(1.1363181) q[0];
rz(0.5592067) q[1];
sx q[1];
rz(-1.8103841) q[1];
sx q[1];
rz(0.24787535) q[1];
rz(-1.2216907) q[2];
sx q[2];
rz(-2.2080055) q[2];
sx q[2];
rz(0.49127084) q[2];
rz(-2.2918626) q[3];
sx q[3];
rz(-0.17912514) q[3];
sx q[3];
rz(2.9368243) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
