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
rz(-2.0055327) q[1];
sx q[1];
rz(4.0772822) q[1];
sx q[1];
rz(4.7128591) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4966272) q[0];
sx q[0];
rz(-2.6500406) q[0];
sx q[0];
rz(0.77792032) q[0];
rz(0.76531305) q[2];
sx q[2];
rz(-2.1047449) q[2];
sx q[2];
rz(0.050616654) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1882602) q[1];
sx q[1];
rz(-0.42020513) q[1];
sx q[1];
rz(-0.62700595) q[1];
rz(-0.19574638) q[3];
sx q[3];
rz(-1.2139075) q[3];
sx q[3];
rz(1.9407879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2661665) q[2];
sx q[2];
rz(-0.29310075) q[2];
sx q[2];
rz(-2.0092633) q[2];
rz(-1.6752361) q[3];
sx q[3];
rz(-1.3365859) q[3];
sx q[3];
rz(2.1291389) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9448626) q[0];
sx q[0];
rz(-2.9319627) q[0];
sx q[0];
rz(-0.18584132) q[0];
rz(0.56022412) q[1];
sx q[1];
rz(-1.8461684) q[1];
sx q[1];
rz(-2.9247608) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5179948) q[0];
sx q[0];
rz(-1.2835842) q[0];
sx q[0];
rz(0.90201305) q[0];
x q[1];
rz(2.1255323) q[2];
sx q[2];
rz(-1.9811355) q[2];
sx q[2];
rz(-2.1330657) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9456957) q[1];
sx q[1];
rz(-2.2356114) q[1];
sx q[1];
rz(2.871454) q[1];
x q[2];
rz(-0.84083765) q[3];
sx q[3];
rz(-0.84078046) q[3];
sx q[3];
rz(2.9527612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.310114) q[2];
sx q[2];
rz(-2.3159413) q[2];
sx q[2];
rz(-1.2878093) q[2];
rz(2.3790322) q[3];
sx q[3];
rz(-1.1688787) q[3];
sx q[3];
rz(2.8365703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4644311) q[0];
sx q[0];
rz(-2.7966249) q[0];
sx q[0];
rz(-2.537354) q[0];
rz(1.3263946) q[1];
sx q[1];
rz(-1.7809968) q[1];
sx q[1];
rz(-0.93260971) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9661449) q[0];
sx q[0];
rz(-1.6245337) q[0];
sx q[0];
rz(-1.2530112) q[0];
rz(-0.89838018) q[2];
sx q[2];
rz(-2.8542238) q[2];
sx q[2];
rz(-2.3786366) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.052913594) q[1];
sx q[1];
rz(-2.1267849) q[1];
sx q[1];
rz(1.4312137) q[1];
x q[2];
rz(-2.2266085) q[3];
sx q[3];
rz(-0.6647771) q[3];
sx q[3];
rz(1.2802326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9937667) q[2];
sx q[2];
rz(-1.0819165) q[2];
sx q[2];
rz(-1.0926584) q[2];
rz(-2.5993733) q[3];
sx q[3];
rz(-2.0565624) q[3];
sx q[3];
rz(2.1742163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7820691) q[0];
sx q[0];
rz(-3.0451267) q[0];
sx q[0];
rz(-2.6413667) q[0];
rz(-2.3362828) q[1];
sx q[1];
rz(-1.1601245) q[1];
sx q[1];
rz(1.6436228) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.995979) q[0];
sx q[0];
rz(-2.0659475) q[0];
sx q[0];
rz(-0.33546319) q[0];
x q[1];
rz(3.0849886) q[2];
sx q[2];
rz(-1.7610234) q[2];
sx q[2];
rz(-0.20400001) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8151617) q[1];
sx q[1];
rz(-1.8306499) q[1];
sx q[1];
rz(-1.8111147) q[1];
x q[2];
rz(-2.6763776) q[3];
sx q[3];
rz(-1.9808931) q[3];
sx q[3];
rz(-2.0801534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3952289) q[2];
sx q[2];
rz(-2.5791898) q[2];
sx q[2];
rz(2.4397819) q[2];
rz(0.83135215) q[3];
sx q[3];
rz(-2.1777007) q[3];
sx q[3];
rz(-0.62197661) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9005301) q[0];
sx q[0];
rz(-0.59589544) q[0];
sx q[0];
rz(2.3262614) q[0];
rz(1.5218081) q[1];
sx q[1];
rz(-2.3074469) q[1];
sx q[1];
rz(2.0933847) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5883023) q[0];
sx q[0];
rz(-2.7748845) q[0];
sx q[0];
rz(2.328863) q[0];
rz(2.4900715) q[2];
sx q[2];
rz(-1.5610352) q[2];
sx q[2];
rz(-2.5196645) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8909) q[1];
sx q[1];
rz(-2.7437468) q[1];
sx q[1];
rz(-0.65244168) q[1];
rz(2.1861595) q[3];
sx q[3];
rz(-1.2293929) q[3];
sx q[3];
rz(-1.1036901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6158225) q[2];
sx q[2];
rz(-0.56695357) q[2];
sx q[2];
rz(-1.099951) q[2];
rz(0.82529092) q[3];
sx q[3];
rz(-1.0422948) q[3];
sx q[3];
rz(0.88551372) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0681756) q[0];
sx q[0];
rz(-0.59403479) q[0];
sx q[0];
rz(-2.2391879) q[0];
rz(-2.1249318) q[1];
sx q[1];
rz(-1.0598176) q[1];
sx q[1];
rz(3.0117603) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3467305) q[0];
sx q[0];
rz(-0.71160337) q[0];
sx q[0];
rz(0.58563389) q[0];
rz(-1.0340704) q[2];
sx q[2];
rz(-2.495129) q[2];
sx q[2];
rz(1.549364) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9784769) q[1];
sx q[1];
rz(-2.4448131) q[1];
sx q[1];
rz(0.58660581) q[1];
rz(-pi) q[2];
rz(0.31452175) q[3];
sx q[3];
rz(-2.5701227) q[3];
sx q[3];
rz(-2.275327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8292024) q[2];
sx q[2];
rz(-2.1924993) q[2];
sx q[2];
rz(2.9373346) q[2];
rz(1.2060818) q[3];
sx q[3];
rz(-1.6198502) q[3];
sx q[3];
rz(0.23541418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7234574) q[0];
sx q[0];
rz(-1.8122939) q[0];
sx q[0];
rz(-1.6947421) q[0];
rz(1.8824668) q[1];
sx q[1];
rz(-0.99021688) q[1];
sx q[1];
rz(2.4553305) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26353729) q[0];
sx q[0];
rz(-2.005308) q[0];
sx q[0];
rz(-2.6146019) q[0];
rz(-0.94481988) q[2];
sx q[2];
rz(-0.84569028) q[2];
sx q[2];
rz(-3.0996029) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2974907) q[1];
sx q[1];
rz(-2.9417848) q[1];
sx q[1];
rz(1.1872477) q[1];
x q[2];
rz(1.2495743) q[3];
sx q[3];
rz(-1.9739082) q[3];
sx q[3];
rz(-1.5954799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4454322) q[2];
sx q[2];
rz(-1.3779209) q[2];
sx q[2];
rz(-3.1398204) q[2];
rz(0.56162515) q[3];
sx q[3];
rz(-0.91149819) q[3];
sx q[3];
rz(1.6368438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5381662) q[0];
sx q[0];
rz(-0.68646938) q[0];
sx q[0];
rz(-1.6954533) q[0];
rz(0.7810477) q[1];
sx q[1];
rz(-1.8361517) q[1];
sx q[1];
rz(1.6400281) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43270375) q[0];
sx q[0];
rz(-2.6484657) q[0];
sx q[0];
rz(-1.6119484) q[0];
rz(-pi) q[1];
rz(-1.5090452) q[2];
sx q[2];
rz(-1.5973063) q[2];
sx q[2];
rz(2.3960631) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9403119) q[1];
sx q[1];
rz(-2.8021325) q[1];
sx q[1];
rz(-0.27075726) q[1];
rz(0.6286962) q[3];
sx q[3];
rz(-1.0843715) q[3];
sx q[3];
rz(1.6823671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.35187307) q[2];
sx q[2];
rz(-1.7799653) q[2];
sx q[2];
rz(-1.3191351) q[2];
rz(-1.9296648) q[3];
sx q[3];
rz(-1.2865678) q[3];
sx q[3];
rz(-2.8222728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8050352) q[0];
sx q[0];
rz(-0.55258495) q[0];
sx q[0];
rz(1.9375027) q[0];
rz(2.7583292) q[1];
sx q[1];
rz(-0.52572322) q[1];
sx q[1];
rz(-0.35167545) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6707014) q[0];
sx q[0];
rz(-2.0235217) q[0];
sx q[0];
rz(-2.7265413) q[0];
x q[1];
rz(-1.8750538) q[2];
sx q[2];
rz(-2.2396302) q[2];
sx q[2];
rz(-2.5077016) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5896776) q[1];
sx q[1];
rz(-0.81112408) q[1];
sx q[1];
rz(3.0685436) q[1];
rz(-pi) q[2];
rz(-0.10961253) q[3];
sx q[3];
rz(-1.5188367) q[3];
sx q[3];
rz(0.51652858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7982771) q[2];
sx q[2];
rz(-2.0337992) q[2];
sx q[2];
rz(-1.8593672) q[2];
rz(-1.6451689) q[3];
sx q[3];
rz(-1.6069501) q[3];
sx q[3];
rz(1.055868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6431817) q[0];
sx q[0];
rz(-1.8739941) q[0];
sx q[0];
rz(0.19432755) q[0];
rz(2.1037897) q[1];
sx q[1];
rz(-0.56832814) q[1];
sx q[1];
rz(-2.1077572) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94668418) q[0];
sx q[0];
rz(-2.2241728) q[0];
sx q[0];
rz(-0.16434591) q[0];
rz(-pi) q[1];
rz(2.6331484) q[2];
sx q[2];
rz(-0.93548453) q[2];
sx q[2];
rz(-2.9490162) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.49627134) q[1];
sx q[1];
rz(-1.4693345) q[1];
sx q[1];
rz(-2.1437777) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7887573) q[3];
sx q[3];
rz(-0.86943227) q[3];
sx q[3];
rz(0.71818128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0795435) q[2];
sx q[2];
rz(-0.94576183) q[2];
sx q[2];
rz(2.5058084) q[2];
rz(0.27030269) q[3];
sx q[3];
rz(-0.79939866) q[3];
sx q[3];
rz(-1.5283782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4476267) q[0];
sx q[0];
rz(-1.8287369) q[0];
sx q[0];
rz(1.0736314) q[0];
rz(-1.4355961) q[1];
sx q[1];
rz(-1.5626848) q[1];
sx q[1];
rz(-2.3609153) q[1];
rz(-0.96799093) q[2];
sx q[2];
rz(-1.597076) q[2];
sx q[2];
rz(1.1521641) q[2];
rz(-3.071143) q[3];
sx q[3];
rz(-1.9000713) q[3];
sx q[3];
rz(-2.6303359) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];