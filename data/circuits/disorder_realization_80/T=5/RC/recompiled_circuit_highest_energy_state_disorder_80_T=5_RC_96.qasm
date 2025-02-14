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
rz(0.39245519) q[0];
sx q[0];
rz(-0.28372228) q[0];
sx q[0];
rz(-1.9537881) q[0];
rz(-0.44161931) q[1];
sx q[1];
rz(-1.8752357) q[1];
sx q[1];
rz(1.1511572) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6080127) q[0];
sx q[0];
rz(-1.5927497) q[0];
sx q[0];
rz(-1.4826464) q[0];
rz(0.48502977) q[2];
sx q[2];
rz(-1.2425307) q[2];
sx q[2];
rz(-0.94543823) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3656039) q[1];
sx q[1];
rz(-1.9950486) q[1];
sx q[1];
rz(1.4308903) q[1];
rz(-pi) q[2];
rz(-2.092716) q[3];
sx q[3];
rz(-1.3112231) q[3];
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
rz(-2.2448764) q[2];
sx q[2];
rz(2.8420281) q[2];
rz(-0.35563955) q[3];
sx q[3];
rz(-2.1483597) q[3];
sx q[3];
rz(2.7049844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75293175) q[0];
sx q[0];
rz(-0.54838538) q[0];
sx q[0];
rz(0.27632982) q[0];
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
rz(2.4014028) q[0];
sx q[0];
rz(-1.3200134) q[0];
sx q[0];
rz(2.0987505) q[0];
rz(0.33832835) q[2];
sx q[2];
rz(-1.6436495) q[2];
sx q[2];
rz(-0.05073994) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4071634) q[1];
sx q[1];
rz(-1.8968655) q[1];
sx q[1];
rz(2.3666275) q[1];
rz(-pi) q[2];
rz(1.4088936) q[3];
sx q[3];
rz(-1.5675968) q[3];
sx q[3];
rz(-2.7174983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.77208272) q[2];
sx q[2];
rz(-1.778506) q[2];
sx q[2];
rz(2.962964) q[2];
rz(-1.5313088) q[3];
sx q[3];
rz(-2.2154) q[3];
sx q[3];
rz(-2.3791651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-1.4924563) q[1];
sx q[1];
rz(-1.656146) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0091772) q[0];
sx q[0];
rz(-0.88736358) q[0];
sx q[0];
rz(-1.528549) q[0];
rz(-2.7316667) q[2];
sx q[2];
rz(-2.2927444) q[2];
sx q[2];
rz(0.97463911) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.29596696) q[1];
sx q[1];
rz(-0.53301552) q[1];
sx q[1];
rz(-2.19341) q[1];
rz(-pi) q[2];
rz(2.2353506) q[3];
sx q[3];
rz(-1.8633047) q[3];
sx q[3];
rz(-2.9952303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9401271) q[2];
sx q[2];
rz(-1.0041693) q[2];
sx q[2];
rz(-2.6489769) q[2];
rz(-2.7228739) q[3];
sx q[3];
rz(-1.0227572) q[3];
sx q[3];
rz(-0.59188265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.017460499) q[0];
sx q[0];
rz(-3.0191665) q[0];
sx q[0];
rz(-0.19805743) q[0];
rz(-0.65912229) q[1];
sx q[1];
rz(-2.5810869) q[1];
sx q[1];
rz(-2.9445599) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8640932) q[0];
sx q[0];
rz(-0.53170337) q[0];
sx q[0];
rz(-3.0752027) q[0];
x q[1];
rz(0.87984271) q[2];
sx q[2];
rz(-0.36565271) q[2];
sx q[2];
rz(1.8118389) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9777364) q[1];
sx q[1];
rz(-1.9242373) q[1];
sx q[1];
rz(-0.56956915) q[1];
rz(-pi) q[2];
rz(1.0258338) q[3];
sx q[3];
rz(-1.6022046) q[3];
sx q[3];
rz(1.9592913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.78843242) q[2];
sx q[2];
rz(-1.8400695) q[2];
sx q[2];
rz(-1.8678467) q[2];
rz(1.5984009) q[3];
sx q[3];
rz(-1.1929932) q[3];
sx q[3];
rz(-0.25632349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0850328) q[0];
sx q[0];
rz(-1.7340478) q[0];
sx q[0];
rz(-0.10277596) q[0];
rz(-0.13725266) q[1];
sx q[1];
rz(-2.6922701) q[1];
sx q[1];
rz(-1.4385673) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5202839) q[0];
sx q[0];
rz(-0.94693333) q[0];
sx q[0];
rz(0.72786967) q[0];
x q[1];
rz(-0.091890865) q[2];
sx q[2];
rz(-2.6767807) q[2];
sx q[2];
rz(-1.6824695) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4754191) q[1];
sx q[1];
rz(-1.3120781) q[1];
sx q[1];
rz(-0.15260507) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4047671) q[3];
sx q[3];
rz(-2.9758436) q[3];
sx q[3];
rz(-0.98893569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8370301) q[2];
sx q[2];
rz(-1.4745159) q[2];
sx q[2];
rz(2.4414818) q[2];
rz(2.2434798) q[3];
sx q[3];
rz(-2.5105295) q[3];
sx q[3];
rz(-1.496544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(-0.20823088) q[0];
sx q[0];
rz(-2.5550186) q[0];
sx q[0];
rz(1.5640278) q[0];
rz(-2.481752) q[1];
sx q[1];
rz(-2.3520062) q[1];
sx q[1];
rz(0.97445625) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81466573) q[0];
sx q[0];
rz(-1.5788933) q[0];
sx q[0];
rz(0.011122313) q[0];
x q[1];
rz(1.3979369) q[2];
sx q[2];
rz(-1.0893911) q[2];
sx q[2];
rz(1.0499894) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8169957) q[1];
sx q[1];
rz(-1.5815418) q[1];
sx q[1];
rz(-0.50159494) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53831999) q[3];
sx q[3];
rz(-1.622952) q[3];
sx q[3];
rz(-2.5405034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4727891) q[2];
sx q[2];
rz(-2.4432224) q[2];
sx q[2];
rz(2.2051956) q[2];
rz(-0.89013731) q[3];
sx q[3];
rz(-1.3606768) q[3];
sx q[3];
rz(0.20588017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5278006) q[0];
sx q[0];
rz(-2.9449936) q[0];
sx q[0];
rz(3.0278681) q[0];
rz(1.3971036) q[1];
sx q[1];
rz(-2.0753658) q[1];
sx q[1];
rz(3.1166335) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7222116) q[0];
sx q[0];
rz(-1.352465) q[0];
sx q[0];
rz(0.055856987) q[0];
rz(-0.34100248) q[2];
sx q[2];
rz(-2.4694168) q[2];
sx q[2];
rz(0.18358809) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.61103454) q[1];
sx q[1];
rz(-1.636158) q[1];
sx q[1];
rz(2.2863131) q[1];
rz(-1.5146144) q[3];
sx q[3];
rz(-1.7683408) q[3];
sx q[3];
rz(-0.88314013) q[3];
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
rz(2.2510236) q[3];
sx q[3];
rz(-1.7361448) q[3];
sx q[3];
rz(-2.4864206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6152076) q[0];
sx q[0];
rz(-0.40281519) q[0];
sx q[0];
rz(2.030754) q[0];
rz(0.088518294) q[1];
sx q[1];
rz(-2.7169777) q[1];
sx q[1];
rz(1.7875338) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18830472) q[0];
sx q[0];
rz(-2.3416356) q[0];
sx q[0];
rz(-2.9842581) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9879775) q[2];
sx q[2];
rz(-1.3014463) q[2];
sx q[2];
rz(0.17698032) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4072666) q[1];
sx q[1];
rz(-1.8167434) q[1];
sx q[1];
rz(1.4926364) q[1];
x q[2];
rz(2.5398387) q[3];
sx q[3];
rz(-2.2838998) q[3];
sx q[3];
rz(-0.32342887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2985349) q[2];
sx q[2];
rz(-0.83710805) q[2];
sx q[2];
rz(-0.38541547) q[2];
rz(-0.75000969) q[3];
sx q[3];
rz(-2.1842726) q[3];
sx q[3];
rz(3.0756522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0176004) q[0];
sx q[0];
rz(-2.0662859) q[0];
sx q[0];
rz(2.3093834) q[0];
rz(1.0819134) q[1];
sx q[1];
rz(-2.3257906) q[1];
sx q[1];
rz(-1.5218511) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6523217) q[0];
sx q[0];
rz(-1.1025363) q[0];
sx q[0];
rz(2.7007607) q[0];
rz(-2.887313) q[2];
sx q[2];
rz(-1.2966208) q[2];
sx q[2];
rz(-0.57480301) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.79241133) q[1];
sx q[1];
rz(-2.1508625) q[1];
sx q[1];
rz(2.4522454) q[1];
rz(-pi) q[2];
rz(-0.24116349) q[3];
sx q[3];
rz(-2.2134288) q[3];
sx q[3];
rz(1.8120899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.99721432) q[2];
sx q[2];
rz(-1.8961366) q[2];
sx q[2];
rz(1.9027556) q[2];
rz(1.291412) q[3];
sx q[3];
rz(-1.8700799) q[3];
sx q[3];
rz(-3.0102357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.130403) q[0];
sx q[0];
rz(-2.1881396) q[0];
sx q[0];
rz(-0.74380547) q[0];
rz(3.0939843) q[1];
sx q[1];
rz(-2.036939) q[1];
sx q[1];
rz(-2.0000752) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6358546) q[0];
sx q[0];
rz(-0.72162745) q[0];
sx q[0];
rz(0.82765915) q[0];
rz(-pi) q[1];
x q[1];
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
rz(-pi) q[0];
x q[0];
rz(2.9824578) q[1];
sx q[1];
rz(-1.1679113) q[1];
sx q[1];
rz(-0.77352685) q[1];
rz(1.3895681) q[3];
sx q[3];
rz(-2.2076105) q[3];
sx q[3];
rz(0.94094482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.49012524) q[2];
sx q[2];
rz(-0.23568025) q[2];
sx q[2];
rz(-0.30533314) q[2];
rz(-1.4281645) q[3];
sx q[3];
rz(-1.7686663) q[3];
sx q[3];
rz(-1.8491687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.4743773) q[2];
sx q[2];
rz(-1.2922774) q[2];
sx q[2];
rz(-0.86624672) q[2];
rz(-0.84973002) q[3];
sx q[3];
rz(-2.9624675) q[3];
sx q[3];
rz(-0.20476838) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
