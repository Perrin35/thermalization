OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6361976) q[0];
sx q[0];
rz(6.0072748) q[0];
sx q[0];
rz(10.732565) q[0];
rz(1.1360599) q[1];
sx q[1];
rz(-0.93568957) q[1];
sx q[1];
rz(-1.5712665) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9309064) q[0];
sx q[0];
rz(-1.9084198) q[0];
sx q[0];
rz(-0.36436413) q[0];
rz(2.2575126) q[2];
sx q[2];
rz(-2.2097217) q[2];
sx q[2];
rz(-2.0761348) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.624144) q[1];
sx q[1];
rz(-1.9074719) q[1];
sx q[1];
rz(-1.8271853) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0516112) q[3];
sx q[3];
rz(-0.40502031) q[3];
sx q[3];
rz(-2.4570176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2661665) q[2];
sx q[2];
rz(-0.29310075) q[2];
sx q[2];
rz(2.0092633) q[2];
rz(1.6752361) q[3];
sx q[3];
rz(-1.3365859) q[3];
sx q[3];
rz(1.0124538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9448626) q[0];
sx q[0];
rz(-2.9319627) q[0];
sx q[0];
rz(-0.18584132) q[0];
rz(-0.56022412) q[1];
sx q[1];
rz(-1.8461684) q[1];
sx q[1];
rz(2.9247608) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29138716) q[0];
sx q[0];
rz(-0.7190401) q[0];
sx q[0];
rz(-1.1262116) q[0];
rz(2.1255323) q[2];
sx q[2];
rz(-1.9811355) q[2];
sx q[2];
rz(1.0085269) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.20570457) q[1];
sx q[1];
rz(-1.3591896) q[1];
sx q[1];
rz(0.88797027) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.64036815) q[3];
sx q[3];
rz(-0.98207563) q[3];
sx q[3];
rz(-1.1191739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8314787) q[2];
sx q[2];
rz(-2.3159413) q[2];
sx q[2];
rz(1.8537834) q[2];
rz(2.3790322) q[3];
sx q[3];
rz(-1.972714) q[3];
sx q[3];
rz(0.30502239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6771616) q[0];
sx q[0];
rz(-2.7966249) q[0];
sx q[0];
rz(-2.537354) q[0];
rz(1.8151981) q[1];
sx q[1];
rz(-1.3605958) q[1];
sx q[1];
rz(-0.93260971) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37768294) q[0];
sx q[0];
rz(-1.8881067) q[0];
sx q[0];
rz(3.0850287) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2432125) q[2];
sx q[2];
rz(-2.8542238) q[2];
sx q[2];
rz(-2.3786366) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.052913594) q[1];
sx q[1];
rz(-1.0148078) q[1];
sx q[1];
rz(1.4312137) q[1];
rz(-pi) q[2];
rz(1.014939) q[3];
sx q[3];
rz(-1.1851289) q[3];
sx q[3];
rz(-2.8876497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9937667) q[2];
sx q[2];
rz(-1.0819165) q[2];
sx q[2];
rz(-2.0489342) q[2];
rz(-2.5993733) q[3];
sx q[3];
rz(-1.0850302) q[3];
sx q[3];
rz(-2.1742163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3595235) q[0];
sx q[0];
rz(-3.0451267) q[0];
sx q[0];
rz(0.50022593) q[0];
rz(2.3362828) q[1];
sx q[1];
rz(-1.1601245) q[1];
sx q[1];
rz(-1.6436228) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7786176) q[0];
sx q[0];
rz(-2.551429) q[0];
sx q[0];
rz(1.0233364) q[0];
x q[1];
rz(-1.285032) q[2];
sx q[2];
rz(-2.9432202) q[2];
sx q[2];
rz(-0.49516585) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.30724635) q[1];
sx q[1];
rz(-1.8028959) q[1];
sx q[1];
rz(2.8744065) q[1];
rz(-pi) q[2];
rz(-2.0235396) q[3];
sx q[3];
rz(-1.1467883) q[3];
sx q[3];
rz(-2.8297569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.74636373) q[2];
sx q[2];
rz(-0.56240288) q[2];
sx q[2];
rz(0.70181075) q[2];
rz(-0.83135215) q[3];
sx q[3];
rz(-2.1777007) q[3];
sx q[3];
rz(-2.519616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9005301) q[0];
sx q[0];
rz(-0.59589544) q[0];
sx q[0];
rz(2.3262614) q[0];
rz(-1.6197846) q[1];
sx q[1];
rz(-0.83414572) q[1];
sx q[1];
rz(1.048208) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4002776) q[0];
sx q[0];
rz(-1.321723) q[0];
sx q[0];
rz(1.8427909) q[0];
rz(-pi) q[1];
rz(-1.5830718) q[2];
sx q[2];
rz(-2.2222812) q[2];
sx q[2];
rz(-0.94142454) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8909) q[1];
sx q[1];
rz(-0.39784583) q[1];
sx q[1];
rz(2.489151) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7311677) q[3];
sx q[3];
rz(-2.1459208) q[3];
sx q[3];
rz(0.23469532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.52577019) q[2];
sx q[2];
rz(-0.56695357) q[2];
sx q[2];
rz(1.099951) q[2];
rz(0.82529092) q[3];
sx q[3];
rz(-1.0422948) q[3];
sx q[3];
rz(0.88551372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0681756) q[0];
sx q[0];
rz(-0.59403479) q[0];
sx q[0];
rz(2.2391879) q[0];
rz(2.1249318) q[1];
sx q[1];
rz(-1.0598176) q[1];
sx q[1];
rz(0.12983233) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3467305) q[0];
sx q[0];
rz(-2.4299893) q[0];
sx q[0];
rz(2.5559588) q[0];
x q[1];
rz(2.7733299) q[2];
sx q[2];
rz(-1.0266745) q[2];
sx q[2];
rz(0.95168176) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2643913) q[1];
sx q[1];
rz(-1.0069205) q[1];
sx q[1];
rz(-2.0045723) q[1];
rz(-0.31452175) q[3];
sx q[3];
rz(-0.57146996) q[3];
sx q[3];
rz(-2.275327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8292024) q[2];
sx q[2];
rz(-0.94909334) q[2];
sx q[2];
rz(2.9373346) q[2];
rz(-1.9355109) q[3];
sx q[3];
rz(-1.6198502) q[3];
sx q[3];
rz(0.23541418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4181353) q[0];
sx q[0];
rz(-1.8122939) q[0];
sx q[0];
rz(1.4468505) q[0];
rz(-1.2591259) q[1];
sx q[1];
rz(-2.1513758) q[1];
sx q[1];
rz(-2.4553305) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8780554) q[0];
sx q[0];
rz(-1.1362846) q[0];
sx q[0];
rz(0.52699071) q[0];
rz(-pi) q[1];
x q[1];
rz(0.94481988) q[2];
sx q[2];
rz(-2.2959024) q[2];
sx q[2];
rz(-3.0996029) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4916617) q[1];
sx q[1];
rz(-1.4964536) q[1];
sx q[1];
rz(1.7564303) q[1];
x q[2];
rz(0.63728441) q[3];
sx q[3];
rz(-0.50989671) q[3];
sx q[3];
rz(-0.84264681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.69616047) q[2];
sx q[2];
rz(-1.7636718) q[2];
sx q[2];
rz(-3.1398204) q[2];
rz(2.5799675) q[3];
sx q[3];
rz(-2.2300945) q[3];
sx q[3];
rz(1.6368438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6034265) q[0];
sx q[0];
rz(-2.4551233) q[0];
sx q[0];
rz(1.6954533) q[0];
rz(-2.360545) q[1];
sx q[1];
rz(-1.3054409) q[1];
sx q[1];
rz(1.5015645) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1743463) q[0];
sx q[0];
rz(-1.5902728) q[0];
sx q[0];
rz(-2.0635701) q[0];
rz(-1.1649706) q[2];
sx q[2];
rz(-0.067194447) q[2];
sx q[2];
rz(-0.42025987) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.20128076) q[1];
sx q[1];
rz(-0.33946013) q[1];
sx q[1];
rz(0.27075726) q[1];
rz(-pi) q[2];
x q[2];
rz(0.6286962) q[3];
sx q[3];
rz(-2.0572212) q[3];
sx q[3];
rz(1.4592255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.35187307) q[2];
sx q[2];
rz(-1.7799653) q[2];
sx q[2];
rz(1.3191351) q[2];
rz(-1.2119279) q[3];
sx q[3];
rz(-1.2865678) q[3];
sx q[3];
rz(-0.31931988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8050352) q[0];
sx q[0];
rz(-2.5890077) q[0];
sx q[0];
rz(-1.9375027) q[0];
rz(0.38326344) q[1];
sx q[1];
rz(-0.52572322) q[1];
sx q[1];
rz(-2.7899172) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8512745) q[0];
sx q[0];
rz(-1.9418678) q[0];
sx q[0];
rz(-2.0593658) q[0];
x q[1];
rz(-2.7792764) q[2];
sx q[2];
rz(-10*pi/13) q[2];
sx q[2];
rz(-2.97646) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6955399) q[1];
sx q[1];
rz(-0.76247588) q[1];
sx q[1];
rz(1.4941077) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6978108) q[3];
sx q[3];
rz(-0.1212596) q[3];
sx q[3];
rz(1.4951984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.3433156) q[2];
sx q[2];
rz(-2.0337992) q[2];
sx q[2];
rz(-1.8593672) q[2];
rz(1.6451689) q[3];
sx q[3];
rz(-1.5346425) q[3];
sx q[3];
rz(1.055868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6431817) q[0];
sx q[0];
rz(-1.8739941) q[0];
sx q[0];
rz(0.19432755) q[0];
rz(-1.0378029) q[1];
sx q[1];
rz(-0.56832814) q[1];
sx q[1];
rz(1.0338354) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1949085) q[0];
sx q[0];
rz(-0.91741981) q[0];
sx q[0];
rz(-2.9772467) q[0];
x q[1];
rz(0.98722234) q[2];
sx q[2];
rz(-0.7910896) q[2];
sx q[2];
rz(2.1949878) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.918805) q[1];
sx q[1];
rz(-2.5606887) q[1];
sx q[1];
rz(1.7564299) q[1];
rz(-0.25063534) q[3];
sx q[3];
rz(-0.72892979) q[3];
sx q[3];
rz(-2.7540516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0795435) q[2];
sx q[2];
rz(-2.1958308) q[2];
sx q[2];
rz(0.6357843) q[2];
rz(0.27030269) q[3];
sx q[3];
rz(-2.342194) q[3];
sx q[3];
rz(1.5283782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4476267) q[0];
sx q[0];
rz(-1.3128558) q[0];
sx q[0];
rz(-2.0679612) q[0];
rz(1.7059965) q[1];
sx q[1];
rz(-1.5626848) q[1];
sx q[1];
rz(-2.3609153) q[1];
rz(-3.1096935) q[2];
sx q[2];
rz(-0.96822856) q[2];
sx q[2];
rz(-0.4005489) q[2];
rz(-1.2407606) q[3];
sx q[3];
rz(-1.5041372) q[3];
sx q[3];
rz(2.1048673) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
