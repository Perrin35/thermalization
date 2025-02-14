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
rz(-2.8615992) q[0];
sx q[0];
rz(-2.1092829) q[0];
sx q[0];
rz(-1.9814459) q[0];
rz(-0.87814826) q[1];
sx q[1];
rz(3.5801764) q[1];
sx q[1];
rz(11.68625) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6683249) q[0];
sx q[0];
rz(-0.62444326) q[0];
sx q[0];
rz(-1.2732328) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9737259) q[2];
sx q[2];
rz(-1.4968728) q[2];
sx q[2];
rz(1.5396724) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1714726) q[1];
sx q[1];
rz(-2.5181541) q[1];
sx q[1];
rz(-1.2894676) q[1];
rz(-2.6875917) q[3];
sx q[3];
rz(-1.4635283) q[3];
sx q[3];
rz(-2.0219181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1151513) q[2];
sx q[2];
rz(-0.5618962) q[2];
sx q[2];
rz(1.119841) q[2];
rz(-1.0037054) q[3];
sx q[3];
rz(-0.34990889) q[3];
sx q[3];
rz(-0.85163918) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.411946) q[0];
sx q[0];
rz(-0.8929407) q[0];
sx q[0];
rz(0.44806421) q[0];
rz(1.0753151) q[1];
sx q[1];
rz(-1.0766462) q[1];
sx q[1];
rz(-3.101128) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1115904) q[0];
sx q[0];
rz(-2.4424822) q[0];
sx q[0];
rz(-1.2707236) q[0];
rz(-1.1574817) q[2];
sx q[2];
rz(-1.4885474) q[2];
sx q[2];
rz(2.8678107) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.76772122) q[1];
sx q[1];
rz(-1.5295856) q[1];
sx q[1];
rz(-0.26880625) q[1];
rz(-0.026626066) q[3];
sx q[3];
rz(-0.31272075) q[3];
sx q[3];
rz(-1.7975765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.86541092) q[2];
sx q[2];
rz(-2.2409596) q[2];
sx q[2];
rz(0.5033699) q[2];
rz(-1.447575) q[3];
sx q[3];
rz(-1.2332799) q[3];
sx q[3];
rz(-0.90731049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4564948) q[0];
sx q[0];
rz(-0.95360294) q[0];
sx q[0];
rz(0.29846919) q[0];
rz(1.5553156) q[1];
sx q[1];
rz(-1.7319873) q[1];
sx q[1];
rz(-1.635199) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28373805) q[0];
sx q[0];
rz(-1.8687792) q[0];
sx q[0];
rz(1.4922754) q[0];
rz(-pi) q[1];
rz(-0.2770642) q[2];
sx q[2];
rz(-2.2383406) q[2];
sx q[2];
rz(2.2570222) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7506444) q[1];
sx q[1];
rz(-1.4381289) q[1];
sx q[1];
rz(-2.4544117) q[1];
rz(-pi) q[2];
rz(-3.0873624) q[3];
sx q[3];
rz(-1.7600865) q[3];
sx q[3];
rz(-2.7638276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.35206587) q[2];
sx q[2];
rz(-0.76849285) q[2];
sx q[2];
rz(-2.3699769) q[2];
rz(1.4113034) q[3];
sx q[3];
rz(-1.227042) q[3];
sx q[3];
rz(2.3119149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4378433) q[0];
sx q[0];
rz(-2.4898536) q[0];
sx q[0];
rz(-1.334345) q[0];
rz(1.5127381) q[1];
sx q[1];
rz(-1.520243) q[1];
sx q[1];
rz(-0.13660647) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2481987) q[0];
sx q[0];
rz(-1.8780439) q[0];
sx q[0];
rz(-2.224438) q[0];
x q[1];
rz(3.0380546) q[2];
sx q[2];
rz(-1.2627201) q[2];
sx q[2];
rz(3.0335226) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.881286) q[1];
sx q[1];
rz(-0.76943123) q[1];
sx q[1];
rz(2.6854482) q[1];
x q[2];
rz(1.2219239) q[3];
sx q[3];
rz(-1.8851981) q[3];
sx q[3];
rz(-1.507198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2887349) q[2];
sx q[2];
rz(-1.6625657) q[2];
sx q[2];
rz(0.23957254) q[2];
rz(-0.95909709) q[3];
sx q[3];
rz(-1.9775016) q[3];
sx q[3];
rz(-1.0378999) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94944823) q[0];
sx q[0];
rz(-2.0862155) q[0];
sx q[0];
rz(1.5437833) q[0];
rz(-1.1100618) q[1];
sx q[1];
rz(-1.9381356) q[1];
sx q[1];
rz(-1.69453) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.869782) q[0];
sx q[0];
rz(-2.8431034) q[0];
sx q[0];
rz(-3.117352) q[0];
rz(1.4571152) q[2];
sx q[2];
rz(-0.34514134) q[2];
sx q[2];
rz(0.98649401) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.93433468) q[1];
sx q[1];
rz(-0.90597766) q[1];
sx q[1];
rz(1.3477172) q[1];
rz(0.54357678) q[3];
sx q[3];
rz(-1.2050306) q[3];
sx q[3];
rz(-0.32777946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1453104) q[2];
sx q[2];
rz(-1.1386394) q[2];
sx q[2];
rz(2.3223257) q[2];
rz(-2.1938358) q[3];
sx q[3];
rz(-0.78815931) q[3];
sx q[3];
rz(-1.2558827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6325697) q[0];
sx q[0];
rz(-0.13795723) q[0];
sx q[0];
rz(0.55009681) q[0];
rz(-2.8944648) q[1];
sx q[1];
rz(-1.1526266) q[1];
sx q[1];
rz(-2.1955042) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45029634) q[0];
sx q[0];
rz(-1.2142203) q[0];
sx q[0];
rz(-2.8647115) q[0];
rz(-1.8515023) q[2];
sx q[2];
rz(-2.9200313) q[2];
sx q[2];
rz(-0.46125476) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4934686) q[1];
sx q[1];
rz(-1.5350684) q[1];
sx q[1];
rz(-0.46983998) q[1];
x q[2];
rz(2.7339659) q[3];
sx q[3];
rz(-0.25926155) q[3];
sx q[3];
rz(1.2490937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.75521022) q[2];
sx q[2];
rz(-0.55042616) q[2];
sx q[2];
rz(2.4605301) q[2];
rz(0.81124535) q[3];
sx q[3];
rz(-2.097605) q[3];
sx q[3];
rz(-2.5913141) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1193467) q[0];
sx q[0];
rz(-2.5211054) q[0];
sx q[0];
rz(2.3959809) q[0];
rz(-1.9781808) q[1];
sx q[1];
rz(-0.62143505) q[1];
sx q[1];
rz(-0.17722873) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39800571) q[0];
sx q[0];
rz(-1.9015802) q[0];
sx q[0];
rz(-2.8547006) q[0];
x q[1];
rz(1.4800328) q[2];
sx q[2];
rz(-1.6338193) q[2];
sx q[2];
rz(1.742996) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.191135) q[1];
sx q[1];
rz(-1.0656351) q[1];
sx q[1];
rz(2.7209366) q[1];
rz(-2.4468462) q[3];
sx q[3];
rz(-1.5383895) q[3];
sx q[3];
rz(-2.29914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5916198) q[2];
sx q[2];
rz(-0.98229304) q[2];
sx q[2];
rz(-1.6922692) q[2];
rz(-2.5189404) q[3];
sx q[3];
rz(-0.52716523) q[3];
sx q[3];
rz(1.3194293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1249238) q[0];
sx q[0];
rz(-1.387384) q[0];
sx q[0];
rz(-1.5119875) q[0];
rz(0.92388693) q[1];
sx q[1];
rz(-1.7216564) q[1];
sx q[1];
rz(-2.5249544) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23417191) q[0];
sx q[0];
rz(-1.7230095) q[0];
sx q[0];
rz(2.175075) q[0];
x q[1];
rz(3.055279) q[2];
sx q[2];
rz(-2.6722135) q[2];
sx q[2];
rz(-0.98938194) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0860446) q[1];
sx q[1];
rz(-1.4851855) q[1];
sx q[1];
rz(0.61998359) q[1];
rz(-0.40469669) q[3];
sx q[3];
rz(-2.5202978) q[3];
sx q[3];
rz(2.9018847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3672678) q[2];
sx q[2];
rz(-2.7820945) q[2];
sx q[2];
rz(0.16505879) q[2];
rz(-0.74602357) q[3];
sx q[3];
rz(-1.6227928) q[3];
sx q[3];
rz(3.0987958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54166722) q[0];
sx q[0];
rz(-0.2921108) q[0];
sx q[0];
rz(0.41411972) q[0];
rz(2.7060624) q[1];
sx q[1];
rz(-2.6256517) q[1];
sx q[1];
rz(-1.2783277) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3992452) q[0];
sx q[0];
rz(-1.3297526) q[0];
sx q[0];
rz(-0.15055116) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6428618) q[2];
sx q[2];
rz(-0.35158508) q[2];
sx q[2];
rz(-1.716223) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8929249) q[1];
sx q[1];
rz(-1.5826514) q[1];
sx q[1];
rz(0.067889386) q[1];
x q[2];
rz(1.4632332) q[3];
sx q[3];
rz(-0.53036571) q[3];
sx q[3];
rz(-1.4211185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.72329676) q[2];
sx q[2];
rz(-1.7593242) q[2];
sx q[2];
rz(-2.3988775) q[2];
rz(-2.3501979) q[3];
sx q[3];
rz(-0.68358889) q[3];
sx q[3];
rz(-1.7323823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0940229) q[0];
sx q[0];
rz(-2.8122734) q[0];
sx q[0];
rz(2.2299715) q[0];
rz(-2.1564663) q[1];
sx q[1];
rz(-2.7148425) q[1];
sx q[1];
rz(-1.3381348) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8050753) q[0];
sx q[0];
rz(-1.0960082) q[0];
sx q[0];
rz(-2.7874448) q[0];
x q[1];
rz(0.31735898) q[2];
sx q[2];
rz(-1.9322763) q[2];
sx q[2];
rz(-1.1731479) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0185408) q[1];
sx q[1];
rz(-1.1710852) q[1];
sx q[1];
rz(-0.91688658) q[1];
rz(-pi) q[2];
x q[2];
rz(0.84566859) q[3];
sx q[3];
rz(-2.659043) q[3];
sx q[3];
rz(3.1180814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7868339) q[2];
sx q[2];
rz(-1.4749196) q[2];
sx q[2];
rz(-1.9762529) q[2];
rz(1.3960086) q[3];
sx q[3];
rz(-1.1895807) q[3];
sx q[3];
rz(-1.5918119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6373445) q[0];
sx q[0];
rz(-1.5550384) q[0];
sx q[0];
rz(0.73873781) q[0];
rz(2.5042116) q[1];
sx q[1];
rz(-1.5370054) q[1];
sx q[1];
rz(-0.76269033) q[1];
rz(1.1470096) q[2];
sx q[2];
rz(-1.4374566) q[2];
sx q[2];
rz(-0.79604436) q[2];
rz(1.1203587) q[3];
sx q[3];
rz(-2.0831345) q[3];
sx q[3];
rz(1.9048445) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
