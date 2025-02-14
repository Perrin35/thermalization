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
rz(3.5340478) q[0];
sx q[0];
rz(3.4253149) q[0];
sx q[0];
rz(14.520159) q[0];
rz(2.6999733) q[1];
sx q[1];
rz(-1.266357) q[1];
sx q[1];
rz(-1.1511572) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6080127) q[0];
sx q[0];
rz(-1.5927497) q[0];
sx q[0];
rz(1.6589462) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5106836) q[2];
sx q[2];
rz(-2.563347) q[2];
sx q[2];
rz(-1.9672245) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.77598876) q[1];
sx q[1];
rz(-1.9950486) q[1];
sx q[1];
rz(-1.7107023) q[1];
x q[2];
rz(-2.0602442) q[3];
sx q[3];
rz(-2.5640805) q[3];
sx q[3];
rz(2.5546527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2426131) q[2];
sx q[2];
rz(-2.2448764) q[2];
sx q[2];
rz(-0.29956451) q[2];
rz(0.35563955) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75293175) q[0];
sx q[0];
rz(-2.5932073) q[0];
sx q[0];
rz(0.27632982) q[0];
rz(0.07946864) q[1];
sx q[1];
rz(-0.53229585) q[1];
sx q[1];
rz(-0.4963378) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9086592) q[0];
sx q[0];
rz(-2.56224) q[0];
sx q[0];
rz(-2.0412372) q[0];
rz(-pi) q[1];
rz(-1.4935819) q[2];
sx q[2];
rz(-1.9081915) q[2];
sx q[2];
rz(1.5959306) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.73442921) q[1];
sx q[1];
rz(-1.2447272) q[1];
sx q[1];
rz(2.3666275) q[1];
rz(-pi) q[2];
rz(-1.5509505) q[3];
sx q[3];
rz(-0.16193411) q[3];
sx q[3];
rz(-1.9753044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.77208272) q[2];
sx q[2];
rz(-1.3630867) q[2];
sx q[2];
rz(-2.962964) q[2];
rz(1.5313088) q[3];
sx q[3];
rz(-0.9261927) q[3];
sx q[3];
rz(0.76242751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4302706) q[0];
sx q[0];
rz(-1.0866168) q[0];
sx q[0];
rz(1.4104728) q[0];
rz(-1.7544282) q[1];
sx q[1];
rz(-1.6491363) q[1];
sx q[1];
rz(-14/(3*pi)) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94233275) q[0];
sx q[0];
rz(-0.6845277) q[0];
sx q[0];
rz(0.051816414) q[0];
rz(-2.335821) q[2];
sx q[2];
rz(-1.2670332) q[2];
sx q[2];
rz(-2.8250776) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.29596696) q[1];
sx q[1];
rz(-2.6085771) q[1];
sx q[1];
rz(-0.94818268) q[1];
x q[2];
rz(2.0250506) q[3];
sx q[3];
rz(-2.4245533) q[3];
sx q[3];
rz(1.3644791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9401271) q[2];
sx q[2];
rz(-1.0041693) q[2];
sx q[2];
rz(-2.6489769) q[2];
rz(-2.7228739) q[3];
sx q[3];
rz(-1.0227572) q[3];
sx q[3];
rz(2.54971) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.017460499) q[0];
sx q[0];
rz(-5/(13*pi)) q[0];
sx q[0];
rz(2.9435352) q[0];
rz(0.65912229) q[1];
sx q[1];
rz(-2.5810869) q[1];
sx q[1];
rz(-0.19703279) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2360508) q[0];
sx q[0];
rz(-1.6044378) q[0];
sx q[0];
rz(-2.6108526) q[0];
rz(1.8577099) q[2];
sx q[2];
rz(-1.3409145) q[2];
sx q[2];
rz(-0.89864602) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9024762) q[1];
sx q[1];
rz(-2.4817339) q[1];
sx q[1];
rz(-2.541594) q[1];
rz(3.1048685) q[3];
sx q[3];
rz(-2.1154599) q[3];
sx q[3];
rz(2.7721318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3531602) q[2];
sx q[2];
rz(-1.8400695) q[2];
sx q[2];
rz(1.2737459) q[2];
rz(-1.5984009) q[3];
sx q[3];
rz(-1.1929932) q[3];
sx q[3];
rz(0.25632349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.0850328) q[0];
sx q[0];
rz(-1.4075449) q[0];
sx q[0];
rz(-3.0388167) q[0];
rz(-0.13725266) q[1];
sx q[1];
rz(-0.44932258) q[1];
sx q[1];
rz(-1.7030254) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52977167) q[0];
sx q[0];
rz(-0.91982929) q[0];
sx q[0];
rz(0.82470993) q[0];
x q[1];
rz(-2.6784727) q[2];
sx q[2];
rz(-1.5296522) q[2];
sx q[2];
rz(-0.1938614) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2074429) q[1];
sx q[1];
rz(-2.8420919) q[1];
sx q[1];
rz(2.092176) q[1];
rz(-pi) q[2];
rz(-3.0183295) q[3];
sx q[3];
rz(-1.4597037) q[3];
sx q[3];
rz(2.9934903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.30456257) q[2];
sx q[2];
rz(-1.4745159) q[2];
sx q[2];
rz(0.70011082) q[2];
rz(2.2434798) q[3];
sx q[3];
rz(-0.63106314) q[3];
sx q[3];
rz(1.496544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9333618) q[0];
sx q[0];
rz(-2.5550186) q[0];
sx q[0];
rz(-1.5775648) q[0];
rz(0.6598407) q[1];
sx q[1];
rz(-2.3520062) q[1];
sx q[1];
rz(-2.1671364) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3269269) q[0];
sx q[0];
rz(-1.5788933) q[0];
sx q[0];
rz(0.011122313) q[0];
rz(-pi) q[1];
rz(0.48759385) q[2];
sx q[2];
rz(-1.7238443) q[2];
sx q[2];
rz(0.4401373) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.240307) q[1];
sx q[1];
rz(-1.069233) q[1];
sx q[1];
rz(-1.5830513) q[1];
rz(3.0401214) q[3];
sx q[3];
rz(-0.54059282) q[3];
sx q[3];
rz(-1.056788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.66880354) q[2];
sx q[2];
rz(-2.4432224) q[2];
sx q[2];
rz(-0.93639708) q[2];
rz(2.2514553) q[3];
sx q[3];
rz(-1.3606768) q[3];
sx q[3];
rz(-2.9357125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5278006) q[0];
sx q[0];
rz(-2.9449936) q[0];
sx q[0];
rz(-3.0278681) q[0];
rz(-1.744489) q[1];
sx q[1];
rz(-1.0662268) q[1];
sx q[1];
rz(-3.1166335) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41938101) q[0];
sx q[0];
rz(-1.352465) q[0];
sx q[0];
rz(-3.0857357) q[0];
x q[1];
rz(2.8005902) q[2];
sx q[2];
rz(-2.4694168) q[2];
sx q[2];
rz(0.18358809) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.125122) q[1];
sx q[1];
rz(-2.2844587) q[1];
sx q[1];
rz(3.0550856) q[1];
rz(-pi) q[2];
rz(-2.8680792) q[3];
sx q[3];
rz(-2.9363147) q[3];
sx q[3];
rz(0.60405788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.52594841) q[2];
sx q[2];
rz(-1.4333466) q[2];
sx q[2];
rz(-2.136266) q[2];
rz(0.89056906) q[3];
sx q[3];
rz(-1.7361448) q[3];
sx q[3];
rz(2.4864206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5263851) q[0];
sx q[0];
rz(-0.40281519) q[0];
sx q[0];
rz(-1.1108387) q[0];
rz(-3.0530744) q[1];
sx q[1];
rz(-0.42461494) q[1];
sx q[1];
rz(-1.7875338) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18830472) q[0];
sx q[0];
rz(-2.3416356) q[0];
sx q[0];
rz(2.9842581) q[0];
x q[1];
rz(1.0646787) q[2];
sx q[2];
rz(-0.30914834) q[2];
sx q[2];
rz(-0.34999166) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.14446324) q[1];
sx q[1];
rz(-1.6465997) q[1];
sx q[1];
rz(2.8949225) q[1];
rz(-pi) q[2];
rz(-0.99127524) q[3];
sx q[3];
rz(-2.2439828) q[3];
sx q[3];
rz(2.0085101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.84305772) q[2];
sx q[2];
rz(-0.83710805) q[2];
sx q[2];
rz(0.38541547) q[2];
rz(2.391583) q[3];
sx q[3];
rz(-0.95732006) q[3];
sx q[3];
rz(-3.0756522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.1239922) q[0];
sx q[0];
rz(-2.0662859) q[0];
sx q[0];
rz(0.83220926) q[0];
rz(-1.0819134) q[1];
sx q[1];
rz(-2.3257906) q[1];
sx q[1];
rz(-1.6197416) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8233021) q[0];
sx q[0];
rz(-2.5099237) q[0];
sx q[0];
rz(0.87000991) q[0];
rz(-0.25427962) q[2];
sx q[2];
rz(-1.2966208) q[2];
sx q[2];
rz(0.57480301) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.19157067) q[1];
sx q[1];
rz(-0.86919987) q[1];
sx q[1];
rz(-2.3413041) q[1];
rz(1.8796104) q[3];
sx q[3];
rz(-2.4612456) q[3];
sx q[3];
rz(-1.718936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1443783) q[2];
sx q[2];
rz(-1.8961366) q[2];
sx q[2];
rz(1.2388371) q[2];
rz(-1.8501806) q[3];
sx q[3];
rz(-1.2715128) q[3];
sx q[3];
rz(3.0102357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0111897) q[0];
sx q[0];
rz(-2.1881396) q[0];
sx q[0];
rz(-0.74380547) q[0];
rz(0.047608308) q[1];
sx q[1];
rz(-1.1046537) q[1];
sx q[1];
rz(-2.0000752) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6689169) q[0];
sx q[0];
rz(-2.0341691) q[0];
sx q[0];
rz(2.1457304) q[0];
rz(0.6170527) q[2];
sx q[2];
rz(-1.8931555) q[2];
sx q[2];
rz(0.53215357) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7773173) q[1];
sx q[1];
rz(-2.2688443) q[1];
sx q[1];
rz(-2.1080522) q[1];
rz(0.64472736) q[3];
sx q[3];
rz(-1.4253741) q[3];
sx q[3];
rz(0.73838387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.49012524) q[2];
sx q[2];
rz(-0.23568025) q[2];
sx q[2];
rz(2.8362595) q[2];
rz(1.4281645) q[3];
sx q[3];
rz(-1.3729264) q[3];
sx q[3];
rz(-1.8491687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3931428) q[0];
sx q[0];
rz(-0.23101692) q[0];
sx q[0];
rz(1.1363181) q[0];
rz(2.582386) q[1];
sx q[1];
rz(-1.3312085) q[1];
sx q[1];
rz(-2.8937173) q[1];
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
