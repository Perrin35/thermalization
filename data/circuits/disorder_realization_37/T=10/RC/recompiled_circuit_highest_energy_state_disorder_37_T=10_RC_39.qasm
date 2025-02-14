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
rz(-1.3353424) q[0];
sx q[0];
rz(3.5080533) q[0];
sx q[0];
rz(8.9656497) q[0];
rz(2.4899809) q[1];
sx q[1];
rz(-1.4067283) q[1];
sx q[1];
rz(-0.23174098) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15329862) q[0];
sx q[0];
rz(-2.7911515) q[0];
sx q[0];
rz(0.20163433) q[0];
rz(-1.5135147) q[2];
sx q[2];
rz(-2.1165032) q[2];
sx q[2];
rz(1.6451665) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.3331795) q[1];
sx q[1];
rz(-1.1929379) q[1];
sx q[1];
rz(-2.4757132) q[1];
rz(-1.6123338) q[3];
sx q[3];
rz(-1.2327063) q[3];
sx q[3];
rz(1.6281782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.090652466) q[2];
sx q[2];
rz(-0.53477627) q[2];
sx q[2];
rz(-1.2789307) q[2];
rz(2.2517962) q[3];
sx q[3];
rz(-1.5225478) q[3];
sx q[3];
rz(2.8029627) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5559674) q[0];
sx q[0];
rz(-0.90921679) q[0];
sx q[0];
rz(-2.2157748) q[0];
rz(2.0766808) q[1];
sx q[1];
rz(-1.5644667) q[1];
sx q[1];
rz(3.056114) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8877713) q[0];
sx q[0];
rz(-1.5740875) q[0];
sx q[0];
rz(-0.38981593) q[0];
rz(-pi) q[1];
rz(2.1824942) q[2];
sx q[2];
rz(-1.2818205) q[2];
sx q[2];
rz(-1.1346357) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.036219941) q[1];
sx q[1];
rz(-0.58507996) q[1];
sx q[1];
rz(-1.1452504) q[1];
x q[2];
rz(0.55172975) q[3];
sx q[3];
rz(-0.76653102) q[3];
sx q[3];
rz(2.6006002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.815879) q[2];
sx q[2];
rz(-1.4789944) q[2];
sx q[2];
rz(-0.386664) q[2];
rz(-1.2169085) q[3];
sx q[3];
rz(-2.7972126) q[3];
sx q[3];
rz(0.2200505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3306408) q[0];
sx q[0];
rz(-0.37136677) q[0];
sx q[0];
rz(2.6445342) q[0];
rz(2.0992384) q[1];
sx q[1];
rz(-2.1940239) q[1];
sx q[1];
rz(-0.29464468) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2089068) q[0];
sx q[0];
rz(-2.8021376) q[0];
sx q[0];
rz(0.20821388) q[0];
rz(-pi) q[1];
rz(-2.5025236) q[2];
sx q[2];
rz(-1.7099755) q[2];
sx q[2];
rz(2.7509909) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.68874796) q[1];
sx q[1];
rz(-2.0032681) q[1];
sx q[1];
rz(0.8853064) q[1];
rz(-pi) q[2];
rz(-2.9975843) q[3];
sx q[3];
rz(-2.1668808) q[3];
sx q[3];
rz(-2.2726187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3543388) q[2];
sx q[2];
rz(-2.2806809) q[2];
sx q[2];
rz(-1.5617237) q[2];
rz(1.076237) q[3];
sx q[3];
rz(-1.8393686) q[3];
sx q[3];
rz(0.92897433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-1.314986) q[0];
sx q[0];
rz(-1.6153233) q[0];
sx q[0];
rz(-0.22931799) q[0];
rz(1.6353105) q[1];
sx q[1];
rz(-1.6835667) q[1];
sx q[1];
rz(0.062072676) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45907606) q[0];
sx q[0];
rz(-2.354265) q[0];
sx q[0];
rz(0.83700257) q[0];
x q[1];
rz(0.96458413) q[2];
sx q[2];
rz(-2.6754489) q[2];
sx q[2];
rz(3.0678444) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5531299) q[1];
sx q[1];
rz(-2.3847918) q[1];
sx q[1];
rz(-0.41017814) q[1];
x q[2];
rz(1.4775949) q[3];
sx q[3];
rz(-1.7708994) q[3];
sx q[3];
rz(2.4148108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0486003) q[2];
sx q[2];
rz(-2.2255662) q[2];
sx q[2];
rz(1.830706) q[2];
rz(3.1033031) q[3];
sx q[3];
rz(-1.3524651) q[3];
sx q[3];
rz(-0.37461764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49486092) q[0];
sx q[0];
rz(-2.4202388) q[0];
sx q[0];
rz(2.2485961) q[0];
rz(-1.0294186) q[1];
sx q[1];
rz(-2.4118377) q[1];
sx q[1];
rz(-0.095887862) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9150118) q[0];
sx q[0];
rz(-1.9862729) q[0];
sx q[0];
rz(-1.0092495) q[0];
rz(1.65972) q[2];
sx q[2];
rz(-1.6409573) q[2];
sx q[2];
rz(1.905575) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.694987) q[1];
sx q[1];
rz(-1.0195273) q[1];
sx q[1];
rz(1.1405844) q[1];
rz(-2.5378102) q[3];
sx q[3];
rz(-2.5386435) q[3];
sx q[3];
rz(0.83151885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.15478495) q[2];
sx q[2];
rz(-1.751535) q[2];
sx q[2];
rz(0.42547697) q[2];
rz(2.7191539) q[3];
sx q[3];
rz(-2.0368302) q[3];
sx q[3];
rz(-0.5031684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-1.7285889) q[0];
sx q[0];
rz(-0.63218963) q[0];
sx q[0];
rz(0.11716209) q[0];
rz(1.6475742) q[1];
sx q[1];
rz(-0.90227503) q[1];
sx q[1];
rz(-2.07043) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2334918) q[0];
sx q[0];
rz(-1.0885493) q[0];
sx q[0];
rz(-1.2511176) q[0];
x q[1];
rz(2.6980459) q[2];
sx q[2];
rz(-0.58813349) q[2];
sx q[2];
rz(3.0556222) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.92381645) q[1];
sx q[1];
rz(-2.4144396) q[1];
sx q[1];
rz(-1.9132861) q[1];
x q[2];
rz(0.31757327) q[3];
sx q[3];
rz(-1.7264719) q[3];
sx q[3];
rz(-0.50567108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5693207) q[2];
sx q[2];
rz(-0.57491493) q[2];
sx q[2];
rz(-0.039483698) q[2];
rz(-0.076400541) q[3];
sx q[3];
rz(-1.1698086) q[3];
sx q[3];
rz(0.45342818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70306784) q[0];
sx q[0];
rz(-2.330133) q[0];
sx q[0];
rz(-0.95034289) q[0];
rz(-1.1514661) q[1];
sx q[1];
rz(-1.5299503) q[1];
sx q[1];
rz(1.3075525) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83043419) q[0];
sx q[0];
rz(-1.3518999) q[0];
sx q[0];
rz(-1.9040742) q[0];
rz(-1.4979532) q[2];
sx q[2];
rz(-1.9460287) q[2];
sx q[2];
rz(2.2323687) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2491409) q[1];
sx q[1];
rz(-1.080852) q[1];
sx q[1];
rz(2.9016414) q[1];
x q[2];
rz(-3.0072104) q[3];
sx q[3];
rz(-0.9690949) q[3];
sx q[3];
rz(-0.4225522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3211956) q[2];
sx q[2];
rz(-0.68009192) q[2];
sx q[2];
rz(2.5879228) q[2];
rz(2.4456444) q[3];
sx q[3];
rz(-1.4724933) q[3];
sx q[3];
rz(1.455201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11947908) q[0];
sx q[0];
rz(-1.4520293) q[0];
sx q[0];
rz(2.8346862) q[0];
rz(-1.8652929) q[1];
sx q[1];
rz(-2.6752094) q[1];
sx q[1];
rz(1.2409522) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3281455) q[0];
sx q[0];
rz(-1.4147583) q[0];
sx q[0];
rz(1.1213844) q[0];
rz(-pi) q[1];
x q[1];
rz(0.5979171) q[2];
sx q[2];
rz(-1.2094524) q[2];
sx q[2];
rz(-1.8449699) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.13682374) q[1];
sx q[1];
rz(-1.7589749) q[1];
sx q[1];
rz(-2.8269563) q[1];
rz(-2.5961865) q[3];
sx q[3];
rz(-0.68531236) q[3];
sx q[3];
rz(-2.3092676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3500195) q[2];
sx q[2];
rz(-0.80524033) q[2];
sx q[2];
rz(1.2903068) q[2];
rz(0.6811412) q[3];
sx q[3];
rz(-0.9001503) q[3];
sx q[3];
rz(-1.6397938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27555585) q[0];
sx q[0];
rz(-2.10502) q[0];
sx q[0];
rz(1.1736897) q[0];
rz(-0.58700079) q[1];
sx q[1];
rz(-2.360011) q[1];
sx q[1];
rz(0.079708286) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83277205) q[0];
sx q[0];
rz(-1.3448633) q[0];
sx q[0];
rz(1.9507879) q[0];
x q[1];
rz(-1.4585895) q[2];
sx q[2];
rz(-1.3507028) q[2];
sx q[2];
rz(-0.28973636) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4533055) q[1];
sx q[1];
rz(-0.88677471) q[1];
sx q[1];
rz(-2.3832393) q[1];
rz(1.8675141) q[3];
sx q[3];
rz(-0.78700262) q[3];
sx q[3];
rz(-2.3001461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.54032636) q[2];
sx q[2];
rz(-1.9129632) q[2];
sx q[2];
rz(-1.3339174) q[2];
rz(-1.5506844) q[3];
sx q[3];
rz(-1.0100789) q[3];
sx q[3];
rz(2.9529115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.659336) q[0];
sx q[0];
rz(-1.5784669) q[0];
sx q[0];
rz(-1.851086) q[0];
rz(-3.105063) q[1];
sx q[1];
rz(-1.9772269) q[1];
sx q[1];
rz(2.0972924) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.999015) q[0];
sx q[0];
rz(-1.8685088) q[0];
sx q[0];
rz(1.1243058) q[0];
x q[1];
rz(1.2408592) q[2];
sx q[2];
rz(-1.82429) q[2];
sx q[2];
rz(1.1546749) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.83936939) q[1];
sx q[1];
rz(-1.4166792) q[1];
sx q[1];
rz(-2.3030512) q[1];
rz(-2.812522) q[3];
sx q[3];
rz(-1.4621203) q[3];
sx q[3];
rz(1.4377126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2269939) q[2];
sx q[2];
rz(-1.031216) q[2];
sx q[2];
rz(-1.3247789) q[2];
rz(-0.75657183) q[3];
sx q[3];
rz(-0.888266) q[3];
sx q[3];
rz(-0.85048401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4751547) q[0];
sx q[0];
rz(-1.403724) q[0];
sx q[0];
rz(-2.4028461) q[0];
rz(-1.7186164) q[1];
sx q[1];
rz(-1.1971133) q[1];
sx q[1];
rz(-2.8308629) q[1];
rz(1.995261) q[2];
sx q[2];
rz(-0.88256114) q[2];
sx q[2];
rz(0.12086856) q[2];
rz(0.66154578) q[3];
sx q[3];
rz(-1.1894124) q[3];
sx q[3];
rz(-0.00581707) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
