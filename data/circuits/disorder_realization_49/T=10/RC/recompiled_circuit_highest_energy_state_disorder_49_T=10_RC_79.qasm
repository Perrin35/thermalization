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
rz(0.3857412) q[0];
sx q[0];
rz(-0.98303151) q[0];
sx q[0];
rz(-0.55735832) q[0];
rz(1.214667) q[1];
sx q[1];
rz(-0.85433706) q[1];
sx q[1];
rz(2.1995423) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1250114) q[0];
sx q[0];
rz(-1.6894369) q[0];
sx q[0];
rz(-1.7811379) q[0];
rz(0.98359707) q[2];
sx q[2];
rz(-2.1237806) q[2];
sx q[2];
rz(1.3830954) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5696733) q[1];
sx q[1];
rz(-1.5077356) q[1];
sx q[1];
rz(0.38625269) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5317134) q[3];
sx q[3];
rz(-2.0470662) q[3];
sx q[3];
rz(2.7000138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.72630924) q[2];
sx q[2];
rz(-1.7674663) q[2];
sx q[2];
rz(-1.3709925) q[2];
rz(-2.3224984) q[3];
sx q[3];
rz(-0.84688014) q[3];
sx q[3];
rz(0.11025652) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46721989) q[0];
sx q[0];
rz(-1.233036) q[0];
sx q[0];
rz(2.9357173) q[0];
rz(-1.3174093) q[1];
sx q[1];
rz(-1.2818047) q[1];
sx q[1];
rz(-2.6726216) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7531573) q[0];
sx q[0];
rz(-2.4860365) q[0];
sx q[0];
rz(-0.038319328) q[0];
rz(1.9838721) q[2];
sx q[2];
rz(-0.5013572) q[2];
sx q[2];
rz(0.92228466) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9392021) q[1];
sx q[1];
rz(-1.1031045) q[1];
sx q[1];
rz(2.2012996) q[1];
x q[2];
rz(-0.75403611) q[3];
sx q[3];
rz(-2.0181351) q[3];
sx q[3];
rz(0.099180982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.59490243) q[2];
sx q[2];
rz(-2.3089843) q[2];
sx q[2];
rz(-2.5140095) q[2];
rz(-1.0707431) q[3];
sx q[3];
rz(-0.50361931) q[3];
sx q[3];
rz(0.32881769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4048432) q[0];
sx q[0];
rz(-0.81031814) q[0];
sx q[0];
rz(-2.3531083) q[0];
rz(0.57111797) q[1];
sx q[1];
rz(-1.8126789) q[1];
sx q[1];
rz(-1.6956537) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75748945) q[0];
sx q[0];
rz(-1.5575842) q[0];
sx q[0];
rz(1.3020975) q[0];
x q[1];
rz(3.0769602) q[2];
sx q[2];
rz(-2.6026313) q[2];
sx q[2];
rz(0.85899437) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.87679497) q[1];
sx q[1];
rz(-1.670627) q[1];
sx q[1];
rz(-2.8431176) q[1];
rz(-pi) q[2];
rz(2.568073) q[3];
sx q[3];
rz(-1.8692501) q[3];
sx q[3];
rz(-0.87707601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0420405) q[2];
sx q[2];
rz(-2.6991762) q[2];
sx q[2];
rz(-2.7624847) q[2];
rz(1.3742617) q[3];
sx q[3];
rz(-0.97135025) q[3];
sx q[3];
rz(-2.1578535) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9336201) q[0];
sx q[0];
rz(-2.3868028) q[0];
sx q[0];
rz(-2.2291613) q[0];
rz(-1.3937996) q[1];
sx q[1];
rz(-0.50519609) q[1];
sx q[1];
rz(-2.8392653) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49719822) q[0];
sx q[0];
rz(-0.14046758) q[0];
sx q[0];
rz(-0.44151784) q[0];
rz(-pi) q[1];
rz(-2.0197844) q[2];
sx q[2];
rz(-0.76631472) q[2];
sx q[2];
rz(1.7036071) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.15790882) q[1];
sx q[1];
rz(-1.527546) q[1];
sx q[1];
rz(-1.3614348) q[1];
rz(-pi) q[2];
rz(-2.264278) q[3];
sx q[3];
rz(-2.1319509) q[3];
sx q[3];
rz(-2.0356242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.48431524) q[2];
sx q[2];
rz(-2.4399098) q[2];
sx q[2];
rz(-2.5122128) q[2];
rz(1.4194007) q[3];
sx q[3];
rz(-1.8604167) q[3];
sx q[3];
rz(-2.4331376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.092954271) q[0];
sx q[0];
rz(-2.104367) q[0];
sx q[0];
rz(-1.3193489) q[0];
rz(2.0535779) q[1];
sx q[1];
rz(-1.8698147) q[1];
sx q[1];
rz(1.6768657) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3977023) q[0];
sx q[0];
rz(-2.4410998) q[0];
sx q[0];
rz(1.6206436) q[0];
rz(-2.508923) q[2];
sx q[2];
rz(-1.5969689) q[2];
sx q[2];
rz(0.90961344) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8301009) q[1];
sx q[1];
rz(-1.7602008) q[1];
sx q[1];
rz(2.7968391) q[1];
x q[2];
rz(2.5358236) q[3];
sx q[3];
rz(-1.7585227) q[3];
sx q[3];
rz(-0.62357215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4488843) q[2];
sx q[2];
rz(-3.0227737) q[2];
sx q[2];
rz(-2.9046655) q[2];
rz(0.98053011) q[3];
sx q[3];
rz(-1.7122995) q[3];
sx q[3];
rz(0.94394365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-2.9834845) q[0];
sx q[0];
rz(-0.10949245) q[0];
sx q[0];
rz(0.14516251) q[0];
rz(1.6204087) q[1];
sx q[1];
rz(-2.3416134) q[1];
sx q[1];
rz(-0.0072341166) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63165346) q[0];
sx q[0];
rz(-1.8064033) q[0];
sx q[0];
rz(-2.2805236) q[0];
rz(-pi) q[1];
rz(2.1223567) q[2];
sx q[2];
rz(-2.2182815) q[2];
sx q[2];
rz(-1.0873356) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.09649) q[1];
sx q[1];
rz(-1.4670925) q[1];
sx q[1];
rz(1.3357049) q[1];
rz(-0.16611734) q[3];
sx q[3];
rz(-1.5990077) q[3];
sx q[3];
rz(1.3835554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.89473692) q[2];
sx q[2];
rz(-1.9523018) q[2];
sx q[2];
rz(-2.5640633) q[2];
rz(1.1524811) q[3];
sx q[3];
rz(-0.58497506) q[3];
sx q[3];
rz(-1.729689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-2.5685527) q[0];
sx q[0];
rz(-0.19565208) q[0];
sx q[0];
rz(-1.0618807) q[0];
rz(1.3878239) q[1];
sx q[1];
rz(-1.9111218) q[1];
sx q[1];
rz(1.6434297) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47344917) q[0];
sx q[0];
rz(-1.4381583) q[0];
sx q[0];
rz(-2.7154865) q[0];
rz(2.7154244) q[2];
sx q[2];
rz(-0.69103253) q[2];
sx q[2];
rz(2.6994801) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0688144) q[1];
sx q[1];
rz(-1.0176976) q[1];
sx q[1];
rz(-0.26550737) q[1];
rz(-0.56363799) q[3];
sx q[3];
rz(-2.5803714) q[3];
sx q[3];
rz(1.3432251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7428703) q[2];
sx q[2];
rz(-0.12464945) q[2];
sx q[2];
rz(2.2275662) q[2];
rz(2.3877609) q[3];
sx q[3];
rz(-1.2316615) q[3];
sx q[3];
rz(2.6851173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1330426) q[0];
sx q[0];
rz(-0.32519105) q[0];
sx q[0];
rz(-3.131026) q[0];
rz(1.0039302) q[1];
sx q[1];
rz(-1.7251451) q[1];
sx q[1];
rz(0.038987003) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7332135) q[0];
sx q[0];
rz(-0.36535242) q[0];
sx q[0];
rz(1.062458) q[0];
rz(1.580367) q[2];
sx q[2];
rz(-1.2540069) q[2];
sx q[2];
rz(0.94632705) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7379563) q[1];
sx q[1];
rz(-1.1925329) q[1];
sx q[1];
rz(-0.72241108) q[1];
x q[2];
rz(0.66308588) q[3];
sx q[3];
rz(-1.5001825) q[3];
sx q[3];
rz(-0.021573349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5616592) q[2];
sx q[2];
rz(-1.9114405) q[2];
sx q[2];
rz(0.099543355) q[2];
rz(1.8482515) q[3];
sx q[3];
rz(-1.8563396) q[3];
sx q[3];
rz(-2.7660363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.67000166) q[0];
sx q[0];
rz(-1.3674068) q[0];
sx q[0];
rz(-1.312183) q[0];
rz(-2.6147764) q[1];
sx q[1];
rz(-2.3878038) q[1];
sx q[1];
rz(2.3048293) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0055702607) q[0];
sx q[0];
rz(-1.3107398) q[0];
sx q[0];
rz(1.8751161) q[0];
rz(0.5469162) q[2];
sx q[2];
rz(-1.0267057) q[2];
sx q[2];
rz(2.5227199) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4556344) q[1];
sx q[1];
rz(-2.1966329) q[1];
sx q[1];
rz(-1.246093) q[1];
rz(-pi) q[2];
x q[2];
rz(2.760254) q[3];
sx q[3];
rz(-2.0102889) q[3];
sx q[3];
rz(0.56909305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7709363) q[2];
sx q[2];
rz(-0.61345658) q[2];
sx q[2];
rz(1.2072309) q[2];
rz(-1.3197445) q[3];
sx q[3];
rz(-1.5650257) q[3];
sx q[3];
rz(2.3453662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0470444) q[0];
sx q[0];
rz(-0.47484174) q[0];
sx q[0];
rz(0.71665254) q[0];
rz(1.5776177) q[1];
sx q[1];
rz(-1.3037222) q[1];
sx q[1];
rz(-0.61734739) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8113239) q[0];
sx q[0];
rz(-1.399938) q[0];
sx q[0];
rz(0.85612423) q[0];
x q[1];
rz(1.1025589) q[2];
sx q[2];
rz(-1.7044633) q[2];
sx q[2];
rz(0.70883162) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.58891728) q[1];
sx q[1];
rz(-2.008634) q[1];
sx q[1];
rz(2.4312996) q[1];
rz(-pi) q[2];
rz(-2.2111597) q[3];
sx q[3];
rz(-1.798822) q[3];
sx q[3];
rz(1.7785398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3295595) q[2];
sx q[2];
rz(-2.7541408) q[2];
sx q[2];
rz(-0.47920245) q[2];
rz(-1.6241578) q[3];
sx q[3];
rz(-1.7370217) q[3];
sx q[3];
rz(1.4994538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1351521) q[0];
sx q[0];
rz(-1.7870446) q[0];
sx q[0];
rz(2.5743299) q[0];
rz(-0.80815036) q[1];
sx q[1];
rz(-1.6193401) q[1];
sx q[1];
rz(-1.5338939) q[1];
rz(2.0907503) q[2];
sx q[2];
rz(-1.4563917) q[2];
sx q[2];
rz(0.2105486) q[2];
rz(2.3324689) q[3];
sx q[3];
rz(-1.7715142) q[3];
sx q[3];
rz(0.51668744) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
