OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.39188448) q[0];
sx q[0];
rz(-0.19667974) q[0];
sx q[0];
rz(-1.952202) q[0];
rz(0.2285129) q[1];
sx q[1];
rz(-0.84140468) q[1];
sx q[1];
rz(0.37766159) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.022764) q[0];
sx q[0];
rz(-1.6435992) q[0];
sx q[0];
rz(-1.4490436) q[0];
rz(3.112552) q[2];
sx q[2];
rz(-1.8543058) q[2];
sx q[2];
rz(-3.0344506) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9610112) q[1];
sx q[1];
rz(-2.2716224) q[1];
sx q[1];
rz(2.5518774) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7226726) q[3];
sx q[3];
rz(-1.6994611) q[3];
sx q[3];
rz(0.045623771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1315786) q[2];
sx q[2];
rz(-2.6476314) q[2];
sx q[2];
rz(0.67260355) q[2];
rz(-0.16942313) q[3];
sx q[3];
rz(-0.38893458) q[3];
sx q[3];
rz(-1.3385564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.910903) q[0];
sx q[0];
rz(-0.90086532) q[0];
sx q[0];
rz(-0.22856523) q[0];
rz(2.9810492) q[1];
sx q[1];
rz(-1.4385782) q[1];
sx q[1];
rz(0.28796089) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44293091) q[0];
sx q[0];
rz(-1.3804111) q[0];
sx q[0];
rz(2.9605244) q[0];
x q[1];
rz(0.38950133) q[2];
sx q[2];
rz(-2.1283538) q[2];
sx q[2];
rz(-0.26706375) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5353363) q[1];
sx q[1];
rz(-1.9793946) q[1];
sx q[1];
rz(2.4819083) q[1];
x q[2];
rz(0.35034758) q[3];
sx q[3];
rz(-1.5347893) q[3];
sx q[3];
rz(0.85067526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5471389) q[2];
sx q[2];
rz(-1.889166) q[2];
sx q[2];
rz(-1.4734369) q[2];
rz(-0.93747059) q[3];
sx q[3];
rz(-0.44527403) q[3];
sx q[3];
rz(-2.766585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(0.32524747) q[0];
sx q[0];
rz(-2.6419817) q[0];
sx q[0];
rz(-0.26741272) q[0];
rz(1.422241) q[1];
sx q[1];
rz(-2.0098675) q[1];
sx q[1];
rz(2.1898988) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50599498) q[0];
sx q[0];
rz(-1.6615168) q[0];
sx q[0];
rz(-1.3784301) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0551664) q[2];
sx q[2];
rz(-1.7679169) q[2];
sx q[2];
rz(-0.65161639) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3141331) q[1];
sx q[1];
rz(-0.810312) q[1];
sx q[1];
rz(2.00287) q[1];
rz(2.1943201) q[3];
sx q[3];
rz(-1.1215278) q[3];
sx q[3];
rz(-2.505213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0114228) q[2];
sx q[2];
rz(-1.6458076) q[2];
sx q[2];
rz(-0.34417957) q[2];
rz(-2.2551645) q[3];
sx q[3];
rz(-0.27799806) q[3];
sx q[3];
rz(0.56604958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13720559) q[0];
sx q[0];
rz(-1.4973649) q[0];
sx q[0];
rz(1.244506) q[0];
rz(-3.0124774) q[1];
sx q[1];
rz(-1.7872417) q[1];
sx q[1];
rz(0.37277645) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0817954) q[0];
sx q[0];
rz(-0.76515388) q[0];
sx q[0];
rz(-2.9761936) q[0];
rz(1.1431085) q[2];
sx q[2];
rz(-0.83979411) q[2];
sx q[2];
rz(-1.0243624) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2671632) q[1];
sx q[1];
rz(-2.0544555) q[1];
sx q[1];
rz(2.8664385) q[1];
x q[2];
rz(-2.8126206) q[3];
sx q[3];
rz(-2.0250468) q[3];
sx q[3];
rz(-2.8153552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.63561511) q[2];
sx q[2];
rz(-1.6515235) q[2];
sx q[2];
rz(-0.90488952) q[2];
rz(0.44057009) q[3];
sx q[3];
rz(-2.6456656) q[3];
sx q[3];
rz(0.99159616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6624517) q[0];
sx q[0];
rz(-2.9859556) q[0];
sx q[0];
rz(1.2878081) q[0];
rz(1.4783391) q[1];
sx q[1];
rz(-1.1454502) q[1];
sx q[1];
rz(-0.53422654) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9984765) q[0];
sx q[0];
rz(-1.8792218) q[0];
sx q[0];
rz(-2.8240859) q[0];
rz(-pi) q[1];
x q[1];
rz(2.062837) q[2];
sx q[2];
rz(-1.8740219) q[2];
sx q[2];
rz(2.6360896) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.6492669) q[1];
sx q[1];
rz(-1.6345134) q[1];
sx q[1];
rz(1.4953514) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9910473) q[3];
sx q[3];
rz(-0.72442043) q[3];
sx q[3];
rz(2.5151593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6405032) q[2];
sx q[2];
rz(-1.2155632) q[2];
sx q[2];
rz(0.14349288) q[2];
rz(-1.3714553) q[3];
sx q[3];
rz(-1.8618795) q[3];
sx q[3];
rz(2.9218856) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4390398) q[0];
sx q[0];
rz(-0.87600231) q[0];
sx q[0];
rz(-2.904536) q[0];
rz(-1.9006231) q[1];
sx q[1];
rz(-2.3126912) q[1];
sx q[1];
rz(-1.3751078) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8538044) q[0];
sx q[0];
rz(-1.3082936) q[0];
sx q[0];
rz(0.35065035) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.077652046) q[2];
sx q[2];
rz(-0.75040557) q[2];
sx q[2];
rz(2.1466308) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.82735862) q[1];
sx q[1];
rz(-1.5463366) q[1];
sx q[1];
rz(-0.41008653) q[1];
rz(-pi) q[2];
rz(1.9938019) q[3];
sx q[3];
rz(-1.1173964) q[3];
sx q[3];
rz(2.6231433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.53283006) q[2];
sx q[2];
rz(-0.71981788) q[2];
sx q[2];
rz(2.9928845) q[2];
rz(-0.016629774) q[3];
sx q[3];
rz(-2.7745268) q[3];
sx q[3];
rz(-2.9490411) q[3];
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
rz(2.6454813) q[0];
sx q[0];
rz(-0.2922903) q[0];
sx q[0];
rz(-2.2684229) q[0];
rz(-0.9219777) q[1];
sx q[1];
rz(-2.0261804) q[1];
sx q[1];
rz(3.057664) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31939313) q[0];
sx q[0];
rz(-2.5376352) q[0];
sx q[0];
rz(1.1711867) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3533808) q[2];
sx q[2];
rz(-1.2251717) q[2];
sx q[2];
rz(-2.2923922) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0107207) q[1];
sx q[1];
rz(-2.2166703) q[1];
sx q[1];
rz(0.801553) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0725669) q[3];
sx q[3];
rz(-0.8419753) q[3];
sx q[3];
rz(2.98711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7579047) q[2];
sx q[2];
rz(-2.7010475) q[2];
sx q[2];
rz(-2.596358) q[2];
rz(0.43045726) q[3];
sx q[3];
rz(-2.0038219) q[3];
sx q[3];
rz(-2.8366413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51171821) q[0];
sx q[0];
rz(-0.015462333) q[0];
sx q[0];
rz(-1.2114725) q[0];
rz(0.97310549) q[1];
sx q[1];
rz(-0.5935697) q[1];
sx q[1];
rz(-0.94582311) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66447542) q[0];
sx q[0];
rz(-2.9070589) q[0];
sx q[0];
rz(2.7995336) q[0];
rz(-pi) q[1];
rz(2.321645) q[2];
sx q[2];
rz(-2.2689399) q[2];
sx q[2];
rz(-1.9422216) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9981873) q[1];
sx q[1];
rz(-2.7618976) q[1];
sx q[1];
rz(0.077296301) q[1];
rz(-pi) q[2];
x q[2];
rz(0.010057851) q[3];
sx q[3];
rz(-2.4754482) q[3];
sx q[3];
rz(2.2705164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8699845) q[2];
sx q[2];
rz(-1.4166069) q[2];
sx q[2];
rz(2.8611709) q[2];
rz(2.5366606) q[3];
sx q[3];
rz(-2.0623902) q[3];
sx q[3];
rz(0.98208565) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1542926) q[0];
sx q[0];
rz(-0.91006088) q[0];
sx q[0];
rz(-2.5262685) q[0];
rz(-0.92957169) q[1];
sx q[1];
rz(-0.85837448) q[1];
sx q[1];
rz(-2.5659134) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.070411365) q[0];
sx q[0];
rz(-2.1343263) q[0];
sx q[0];
rz(2.1348743) q[0];
rz(-pi) q[1];
rz(0.11328463) q[2];
sx q[2];
rz(-1.7656529) q[2];
sx q[2];
rz(1.726113) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.74734028) q[1];
sx q[1];
rz(-1.5963975) q[1];
sx q[1];
rz(-1.4634553) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.18746312) q[3];
sx q[3];
rz(-2.1037357) q[3];
sx q[3];
rz(0.77123469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.82970396) q[2];
sx q[2];
rz(-2.3788033) q[2];
sx q[2];
rz(-3.1196307) q[2];
rz(-2.9711376) q[3];
sx q[3];
rz(-1.0237834) q[3];
sx q[3];
rz(2.8767265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4158674) q[0];
sx q[0];
rz(-0.9265582) q[0];
sx q[0];
rz(-0.6341933) q[0];
rz(-3.1126853) q[1];
sx q[1];
rz(-0.78866619) q[1];
sx q[1];
rz(0.96910563) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64810753) q[0];
sx q[0];
rz(-3.0457975) q[0];
sx q[0];
rz(-0.84798725) q[0];
x q[1];
rz(1*pi/15) q[2];
sx q[2];
rz(-2.3839715) q[2];
sx q[2];
rz(2.9486738) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9358878) q[1];
sx q[1];
rz(-1.3921326) q[1];
sx q[1];
rz(0.67358576) q[1];
rz(-pi) q[2];
rz(-1.2877527) q[3];
sx q[3];
rz(-1.1453562) q[3];
sx q[3];
rz(-0.89983672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4518296) q[2];
sx q[2];
rz(-1.657106) q[2];
sx q[2];
rz(-2.7815212) q[2];
rz(1.6137971) q[3];
sx q[3];
rz(-0.66643047) q[3];
sx q[3];
rz(0.56267363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9344899) q[0];
sx q[0];
rz(-1.5705382) q[0];
sx q[0];
rz(-1.6194153) q[0];
rz(3.0974401) q[1];
sx q[1];
rz(-1.4587198) q[1];
sx q[1];
rz(-1.1062467) q[1];
rz(-1.1883696) q[2];
sx q[2];
rz(-1.2570751) q[2];
sx q[2];
rz(-3.0260069) q[2];
rz(-1.885407) q[3];
sx q[3];
rz(-1.1369858) q[3];
sx q[3];
rz(2.9997957) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
