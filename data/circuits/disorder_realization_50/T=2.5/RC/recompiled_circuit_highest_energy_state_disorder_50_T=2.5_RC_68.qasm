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
rz(1.7293575) q[0];
sx q[0];
rz(3.7339551) q[0];
sx q[0];
rz(10.629551) q[0];
rz(-2.0545948) q[1];
sx q[1];
rz(-0.49538651) q[1];
sx q[1];
rz(1.9538716) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.710085) q[0];
sx q[0];
rz(-2.1174701) q[0];
sx q[0];
rz(-1.7492848) q[0];
rz(-pi) q[1];
rz(-2.0314902) q[2];
sx q[2];
rz(-1.0506127) q[2];
sx q[2];
rz(-2.4536163) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0443327) q[1];
sx q[1];
rz(-9/(16*pi)) q[1];
sx q[1];
rz(-1.3023085) q[1];
rz(-pi) q[2];
rz(3.0023469) q[3];
sx q[3];
rz(-1.2069824) q[3];
sx q[3];
rz(-1.1195967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5997233) q[2];
sx q[2];
rz(-2.5399127) q[2];
sx q[2];
rz(2.0987233) q[2];
rz(3.0964105) q[3];
sx q[3];
rz(-0.16862814) q[3];
sx q[3];
rz(1.6056304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0586108) q[0];
sx q[0];
rz(-3.0284212) q[0];
sx q[0];
rz(-2.163072) q[0];
rz(-1.1631896) q[1];
sx q[1];
rz(-2.1470943) q[1];
sx q[1];
rz(-1.3226343) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.103687) q[0];
sx q[0];
rz(-0.9390489) q[0];
sx q[0];
rz(-0.61236713) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8541932) q[2];
sx q[2];
rz(-1.3958418) q[2];
sx q[2];
rz(-2.8307479) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.95056995) q[1];
sx q[1];
rz(-0.92838597) q[1];
sx q[1];
rz(-0.2705785) q[1];
rz(-pi) q[2];
rz(-2.4473684) q[3];
sx q[3];
rz(-1.4560008) q[3];
sx q[3];
rz(0.59370422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.80403745) q[2];
sx q[2];
rz(-2.940371) q[2];
sx q[2];
rz(3.1031754) q[2];
rz(1.6328968) q[3];
sx q[3];
rz(-1.8149866) q[3];
sx q[3];
rz(1.647515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22817336) q[0];
sx q[0];
rz(-2.4689624) q[0];
sx q[0];
rz(2.9056554) q[0];
rz(1.1783696) q[1];
sx q[1];
rz(-1.6940247) q[1];
sx q[1];
rz(-2.6441914) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5478514) q[0];
sx q[0];
rz(-1.3914131) q[0];
sx q[0];
rz(-0.18130882) q[0];
x q[1];
rz(-1.986662) q[2];
sx q[2];
rz(-1.9344887) q[2];
sx q[2];
rz(1.1689344) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0795035) q[1];
sx q[1];
rz(-1.5487284) q[1];
sx q[1];
rz(2.7479321) q[1];
rz(-pi) q[2];
rz(-1.3128993) q[3];
sx q[3];
rz(-1.329943) q[3];
sx q[3];
rz(-1.1782714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1824789) q[2];
sx q[2];
rz(-1.4664058) q[2];
sx q[2];
rz(1.2314931) q[2];
rz(-1.69151) q[3];
sx q[3];
rz(-1.3030038) q[3];
sx q[3];
rz(2.2836397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0813783) q[0];
sx q[0];
rz(-0.84351081) q[0];
sx q[0];
rz(-2.9486935) q[0];
rz(1.7861722) q[1];
sx q[1];
rz(-1.5813634) q[1];
sx q[1];
rz(-2.9706109) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9381704) q[0];
sx q[0];
rz(-1.5121542) q[0];
sx q[0];
rz(-2.0873656) q[0];
rz(-pi) q[1];
rz(2.5582983) q[2];
sx q[2];
rz(-2.3153044) q[2];
sx q[2];
rz(0.14329958) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.8668756) q[1];
sx q[1];
rz(-1.1422542) q[1];
sx q[1];
rz(-0.11517597) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.234262) q[3];
sx q[3];
rz(-1.7620834) q[3];
sx q[3];
rz(-0.80463791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1811447) q[2];
sx q[2];
rz(-1.9713216) q[2];
sx q[2];
rz(-0.32285264) q[2];
rz(1.3747831) q[3];
sx q[3];
rz(-2.1026473) q[3];
sx q[3];
rz(-0.84097451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5906931) q[0];
sx q[0];
rz(-1.0223848) q[0];
sx q[0];
rz(-2.2160231) q[0];
rz(0.12807056) q[1];
sx q[1];
rz(-1.4984727) q[1];
sx q[1];
rz(0.099460348) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29716897) q[0];
sx q[0];
rz(-0.71650973) q[0];
sx q[0];
rz(1.5708718) q[0];
x q[1];
rz(1.9173309) q[2];
sx q[2];
rz(-1.230403) q[2];
sx q[2];
rz(2.9739591) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.65433733) q[1];
sx q[1];
rz(-2.0981826) q[1];
sx q[1];
rz(2.9156963) q[1];
rz(3.0614047) q[3];
sx q[3];
rz(-1.8155451) q[3];
sx q[3];
rz(0.65746869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1818992) q[2];
sx q[2];
rz(-1.5965896) q[2];
sx q[2];
rz(-0.37364513) q[2];
rz(-0.89961189) q[3];
sx q[3];
rz(-0.74266946) q[3];
sx q[3];
rz(-2.0067748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9099092) q[0];
sx q[0];
rz(-0.21655701) q[0];
sx q[0];
rz(-2.5496971) q[0];
rz(-0.20467155) q[1];
sx q[1];
rz(-2.1494631) q[1];
sx q[1];
rz(-3.0904904) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1280397) q[0];
sx q[0];
rz(-0.2464412) q[0];
sx q[0];
rz(2.4714969) q[0];
rz(-pi) q[1];
rz(-1.899753) q[2];
sx q[2];
rz(-2.4566413) q[2];
sx q[2];
rz(-2.9732413) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.49858677) q[1];
sx q[1];
rz(-2.6664147) q[1];
sx q[1];
rz(2.2982486) q[1];
rz(-pi) q[2];
rz(-1.072364) q[3];
sx q[3];
rz(-1.0921156) q[3];
sx q[3];
rz(-0.60018998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0245612) q[2];
sx q[2];
rz(-2.0270963) q[2];
sx q[2];
rz(-2.047211) q[2];
rz(-0.38703212) q[3];
sx q[3];
rz(-2.6670167) q[3];
sx q[3];
rz(-0.3064557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7759906) q[0];
sx q[0];
rz(-2.2633573) q[0];
sx q[0];
rz(-2.8940417) q[0];
rz(-1.4546825) q[1];
sx q[1];
rz(-1.8984112) q[1];
sx q[1];
rz(-0.036458485) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.847408) q[0];
sx q[0];
rz(-0.80889055) q[0];
sx q[0];
rz(-2.7436867) q[0];
x q[1];
rz(0.70306422) q[2];
sx q[2];
rz(-1.7181239) q[2];
sx q[2];
rz(-0.67455233) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6853906) q[1];
sx q[1];
rz(-1.0674607) q[1];
sx q[1];
rz(-2.6052942) q[1];
rz(-pi) q[2];
rz(0.8730252) q[3];
sx q[3];
rz(-0.11618488) q[3];
sx q[3];
rz(0.98770638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0869007) q[2];
sx q[2];
rz(-0.96618235) q[2];
sx q[2];
rz(1.3521693) q[2];
rz(-1.8933206) q[3];
sx q[3];
rz(-1.5636445) q[3];
sx q[3];
rz(-1.1917535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2317113) q[0];
sx q[0];
rz(-2.8928962) q[0];
sx q[0];
rz(-1.6491718) q[0];
rz(-2.9630648) q[1];
sx q[1];
rz(-1.8793722) q[1];
sx q[1];
rz(2.5915204) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85066768) q[0];
sx q[0];
rz(-1.5508473) q[0];
sx q[0];
rz(-1.5394506) q[0];
rz(-pi) q[1];
rz(-1.3230611) q[2];
sx q[2];
rz(-1.5093409) q[2];
sx q[2];
rz(-2.2814192) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8892563) q[1];
sx q[1];
rz(-1.0358216) q[1];
sx q[1];
rz(1.757213) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7356803) q[3];
sx q[3];
rz(-2.0024869) q[3];
sx q[3];
rz(-1.622705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.34755808) q[2];
sx q[2];
rz(-1.7358235) q[2];
sx q[2];
rz(-0.75322914) q[2];
rz(0.29427648) q[3];
sx q[3];
rz(-0.080065057) q[3];
sx q[3];
rz(0.068597138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.8681965) q[0];
sx q[0];
rz(-2.8559339) q[0];
sx q[0];
rz(1.6896601) q[0];
rz(2.946335) q[1];
sx q[1];
rz(-1.0804907) q[1];
sx q[1];
rz(1.6498227) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8378252) q[0];
sx q[0];
rz(-1.2653102) q[0];
sx q[0];
rz(0.49834337) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.765018) q[2];
sx q[2];
rz(-1.6510626) q[2];
sx q[2];
rz(1.8945872) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0635707) q[1];
sx q[1];
rz(-1.6352904) q[1];
sx q[1];
rz(1.7266573) q[1];
rz(-pi) q[2];
rz(-1.5641065) q[3];
sx q[3];
rz(-1.5982398) q[3];
sx q[3];
rz(2.5344078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1008272) q[2];
sx q[2];
rz(-0.46697524) q[2];
sx q[2];
rz(-2.2838498) q[2];
rz(1.4486676) q[3];
sx q[3];
rz(-1.6489776) q[3];
sx q[3];
rz(0.15374507) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0070852) q[0];
sx q[0];
rz(-2.9869098) q[0];
sx q[0];
rz(-0.78306985) q[0];
rz(1.1231517) q[1];
sx q[1];
rz(-1.4479366) q[1];
sx q[1];
rz(0.060308594) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9911847) q[0];
sx q[0];
rz(-1.6643608) q[0];
sx q[0];
rz(1.9267328) q[0];
x q[1];
rz(3.0357828) q[2];
sx q[2];
rz(-2.1544837) q[2];
sx q[2];
rz(0.90487827) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.13014252) q[1];
sx q[1];
rz(-0.9025652) q[1];
sx q[1];
rz(1.6570309) q[1];
rz(-pi) q[2];
rz(-1.4571244) q[3];
sx q[3];
rz(-1.2099363) q[3];
sx q[3];
rz(-1.833503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9874838) q[2];
sx q[2];
rz(-0.86468148) q[2];
sx q[2];
rz(1.8314499) q[2];
rz(1.8080669) q[3];
sx q[3];
rz(-1.798505) q[3];
sx q[3];
rz(1.8346627) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46351984) q[0];
sx q[0];
rz(-1.1165883) q[0];
sx q[0];
rz(-0.22620871) q[0];
rz(0.75640596) q[1];
sx q[1];
rz(-2.6466128) q[1];
sx q[1];
rz(-0.80642798) q[1];
rz(-0.020515223) q[2];
sx q[2];
rz(-0.43045577) q[2];
sx q[2];
rz(0.53878709) q[2];
rz(-1.2213094) q[3];
sx q[3];
rz(-1.3987473) q[3];
sx q[3];
rz(-1.7940298) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
