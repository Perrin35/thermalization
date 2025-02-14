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
rz(-1.3136366) q[0];
sx q[0];
rz(-3.0744636) q[0];
sx q[0];
rz(-0.022775291) q[0];
rz(1.0442806) q[1];
sx q[1];
rz(-2.2106946) q[1];
sx q[1];
rz(0.013962362) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9305796) q[0];
sx q[0];
rz(-1.6114116) q[0];
sx q[0];
rz(0.75675772) q[0];
x q[1];
rz(2.3910782) q[2];
sx q[2];
rz(-2.566411) q[2];
sx q[2];
rz(1.2746948) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.57052862) q[1];
sx q[1];
rz(-1.7693874) q[1];
sx q[1];
rz(1.3873439) q[1];
rz(-pi) q[2];
x q[2];
rz(0.45043378) q[3];
sx q[3];
rz(-1.8897459) q[3];
sx q[3];
rz(-1.7613086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.17249168) q[2];
sx q[2];
rz(-1.8164941) q[2];
sx q[2];
rz(-1.3410404) q[2];
rz(1.4260882) q[3];
sx q[3];
rz(-1.1172349) q[3];
sx q[3];
rz(0.30556998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7476927) q[0];
sx q[0];
rz(-1.9992398) q[0];
sx q[0];
rz(2.4161412) q[0];
rz(0.91616383) q[1];
sx q[1];
rz(-0.80892816) q[1];
sx q[1];
rz(0.61942548) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5529157) q[0];
sx q[0];
rz(-1.1631794) q[0];
sx q[0];
rz(1.0199532) q[0];
rz(-pi) q[1];
rz(0.97255637) q[2];
sx q[2];
rz(-0.91678874) q[2];
sx q[2];
rz(1.1197661) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.404285) q[1];
sx q[1];
rz(-2.5661307) q[1];
sx q[1];
rz(2.2810443) q[1];
rz(2.8579669) q[3];
sx q[3];
rz(-1.1414293) q[3];
sx q[3];
rz(0.46653433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3118423) q[2];
sx q[2];
rz(-2.1205015) q[2];
sx q[2];
rz(-2.1799083) q[2];
rz(1.1059777) q[3];
sx q[3];
rz(-1.032369) q[3];
sx q[3];
rz(0.34524125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.85213929) q[0];
sx q[0];
rz(-1.1281321) q[0];
sx q[0];
rz(-1.7732675) q[0];
rz(-0.013956919) q[1];
sx q[1];
rz(-2.0266666) q[1];
sx q[1];
rz(-1.4143292) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.208983) q[0];
sx q[0];
rz(-1.480672) q[0];
sx q[0];
rz(1.8831253) q[0];
x q[1];
rz(-1.3886202) q[2];
sx q[2];
rz(-0.37473703) q[2];
sx q[2];
rz(0.47606766) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6048879) q[1];
sx q[1];
rz(-1.9997215) q[1];
sx q[1];
rz(-2.1307441) q[1];
rz(-1.3120804) q[3];
sx q[3];
rz(-1.95041) q[3];
sx q[3];
rz(2.1624452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.76325512) q[2];
sx q[2];
rz(-1.1474643) q[2];
sx q[2];
rz(-2.9249127) q[2];
rz(2.2583151) q[3];
sx q[3];
rz(-2.7933385) q[3];
sx q[3];
rz(1.1368375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9388409) q[0];
sx q[0];
rz(-1.0067679) q[0];
sx q[0];
rz(0.83813611) q[0];
rz(0.54191598) q[1];
sx q[1];
rz(-2.7653265) q[1];
sx q[1];
rz(-3.0964877) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7823001) q[0];
sx q[0];
rz(-1.703039) q[0];
sx q[0];
rz(0.6036997) q[0];
x q[1];
rz(-2.0397268) q[2];
sx q[2];
rz(-2.6108338) q[2];
sx q[2];
rz(1.6380978) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1029658) q[1];
sx q[1];
rz(-0.47768394) q[1];
sx q[1];
rz(-1.716502) q[1];
rz(-pi) q[2];
rz(1.4724422) q[3];
sx q[3];
rz(-2.062766) q[3];
sx q[3];
rz(1.4592001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0185467) q[2];
sx q[2];
rz(-2.9728153) q[2];
sx q[2];
rz(-0.65400845) q[2];
rz(2.5096014) q[3];
sx q[3];
rz(-1.7157103) q[3];
sx q[3];
rz(-0.33647195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3376447) q[0];
sx q[0];
rz(-0.73886442) q[0];
sx q[0];
rz(1.3496572) q[0];
rz(-2.3140287) q[1];
sx q[1];
rz(-2.5065828) q[1];
sx q[1];
rz(-0.48447022) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0837528) q[0];
sx q[0];
rz(-0.50005243) q[0];
sx q[0];
rz(0.82885784) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4260277) q[2];
sx q[2];
rz(-2.5826404) q[2];
sx q[2];
rz(1.1048855) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3282991) q[1];
sx q[1];
rz(-0.44078953) q[1];
sx q[1];
rz(-2.8404124) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6788533) q[3];
sx q[3];
rz(-0.92904186) q[3];
sx q[3];
rz(-2.3433507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1822002) q[2];
sx q[2];
rz(-2.0727938) q[2];
sx q[2];
rz(0.63641316) q[2];
rz(-3.0427129) q[3];
sx q[3];
rz(-1.585588) q[3];
sx q[3];
rz(1.6695361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32475489) q[0];
sx q[0];
rz(-0.88961283) q[0];
sx q[0];
rz(-1.1725934) q[0];
rz(1.532754) q[1];
sx q[1];
rz(-1.0120665) q[1];
sx q[1];
rz(3.1408659) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4662292) q[0];
sx q[0];
rz(-0.6777938) q[0];
sx q[0];
rz(-2.2905605) q[0];
rz(-0.41623314) q[2];
sx q[2];
rz(-2.9011167) q[2];
sx q[2];
rz(1.7822303) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7757984) q[1];
sx q[1];
rz(-1.9025814) q[1];
sx q[1];
rz(-1.4514239) q[1];
rz(-pi) q[2];
rz(2.5623393) q[3];
sx q[3];
rz(-1.7369075) q[3];
sx q[3];
rz(-0.36234713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8751004) q[2];
sx q[2];
rz(-0.88783395) q[2];
sx q[2];
rz(-2.9317741) q[2];
rz(-0.80875129) q[3];
sx q[3];
rz(-2.5918312) q[3];
sx q[3];
rz(-2.2426898) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6784191) q[0];
sx q[0];
rz(-1.2043948) q[0];
sx q[0];
rz(-3.0294982) q[0];
rz(0.049292715) q[1];
sx q[1];
rz(-2.6105328) q[1];
sx q[1];
rz(1.8045527) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1111832) q[0];
sx q[0];
rz(-0.93713916) q[0];
sx q[0];
rz(-0.2588057) q[0];
rz(1.7981729) q[2];
sx q[2];
rz(-0.51926398) q[2];
sx q[2];
rz(2.8092172) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.783205) q[1];
sx q[1];
rz(-1.2239982) q[1];
sx q[1];
rz(1.5814745) q[1];
rz(-pi) q[2];
rz(-1.8932166) q[3];
sx q[3];
rz(-2.5395576) q[3];
sx q[3];
rz(-0.69794929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2992531) q[2];
sx q[2];
rz(-0.82483333) q[2];
sx q[2];
rz(0.71005589) q[2];
rz(3.1402785) q[3];
sx q[3];
rz(-0.90130663) q[3];
sx q[3];
rz(2.5406751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32427078) q[0];
sx q[0];
rz(-1.7621499) q[0];
sx q[0];
rz(-0.58260179) q[0];
rz(0.6800037) q[1];
sx q[1];
rz(-2.1425207) q[1];
sx q[1];
rz(-1.8780139) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3124426) q[0];
sx q[0];
rz(-1.1151322) q[0];
sx q[0];
rz(-1.2553041) q[0];
rz(1.3538876) q[2];
sx q[2];
rz(-1.3706651) q[2];
sx q[2];
rz(1.6871802) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.42238754) q[1];
sx q[1];
rz(-0.68049708) q[1];
sx q[1];
rz(-2.8533555) q[1];
rz(-pi) q[2];
rz(-1.6035147) q[3];
sx q[3];
rz(-1.9944291) q[3];
sx q[3];
rz(-0.7684497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5052283) q[2];
sx q[2];
rz(-1.6978846) q[2];
sx q[2];
rz(-0.12602028) q[2];
rz(-1.5382918) q[3];
sx q[3];
rz(-0.77682972) q[3];
sx q[3];
rz(2.7974424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4890323) q[0];
sx q[0];
rz(-1.6283988) q[0];
sx q[0];
rz(1.162758) q[0];
rz(1.3105505) q[1];
sx q[1];
rz(-0.73859155) q[1];
sx q[1];
rz(-1.383925) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7689631) q[0];
sx q[0];
rz(-1.2314737) q[0];
sx q[0];
rz(3.0940381) q[0];
x q[1];
rz(0.53015401) q[2];
sx q[2];
rz(-0.14214504) q[2];
sx q[2];
rz(-1.7696385) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0392929) q[1];
sx q[1];
rz(-2.1792534) q[1];
sx q[1];
rz(3.0000172) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4320763) q[3];
sx q[3];
rz(-1.528421) q[3];
sx q[3];
rz(-0.78813533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.80013529) q[2];
sx q[2];
rz(-1.8693482) q[2];
sx q[2];
rz(-2.6696491) q[2];
rz(0.74884549) q[3];
sx q[3];
rz(-2.2180836) q[3];
sx q[3];
rz(-2.3239465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6227459) q[0];
sx q[0];
rz(-2.2846344) q[0];
sx q[0];
rz(2.6527606) q[0];
rz(0.34293109) q[1];
sx q[1];
rz(-0.88645005) q[1];
sx q[1];
rz(2.0874646) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80741548) q[0];
sx q[0];
rz(-1.5332959) q[0];
sx q[0];
rz(1.6223909) q[0];
rz(2.6661886) q[2];
sx q[2];
rz(-0.66988551) q[2];
sx q[2];
rz(0.49150447) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.78275241) q[1];
sx q[1];
rz(-2.5237875) q[1];
sx q[1];
rz(1.0516404) q[1];
rz(1.098759) q[3];
sx q[3];
rz(-1.6141547) q[3];
sx q[3];
rz(1.2809917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1663345) q[2];
sx q[2];
rz(-0.20938337) q[2];
sx q[2];
rz(-3.0496822) q[2];
rz(-2.6722243) q[3];
sx q[3];
rz(-2.2769603) q[3];
sx q[3];
rz(-3.0288467) q[3];
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
rz(2.5975006) q[0];
sx q[0];
rz(-1.8202029) q[0];
sx q[0];
rz(1.2910917) q[0];
rz(0.76580936) q[1];
sx q[1];
rz(-1.3594834) q[1];
sx q[1];
rz(1.8654738) q[1];
rz(-1.2086537) q[2];
sx q[2];
rz(-1.8100783) q[2];
sx q[2];
rz(1.7892005) q[2];
rz(-2.7226483) q[3];
sx q[3];
rz(-0.57804116) q[3];
sx q[3];
rz(-3.0347435) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
