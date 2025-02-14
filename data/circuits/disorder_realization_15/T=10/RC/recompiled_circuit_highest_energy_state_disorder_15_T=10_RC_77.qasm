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
rz(3.1188174) q[0];
rz(1.0442806) q[1];
sx q[1];
rz(-2.2106946) q[1];
sx q[1];
rz(-3.1276303) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21101306) q[0];
sx q[0];
rz(-1.6114116) q[0];
sx q[0];
rz(-2.3848349) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6988552) q[2];
sx q[2];
rz(-1.95089) q[2];
sx q[2];
rz(-0.9600823) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.57052862) q[1];
sx q[1];
rz(-1.3722053) q[1];
sx q[1];
rz(-1.3873439) q[1];
rz(-pi) q[2];
rz(1.9223677) q[3];
sx q[3];
rz(-1.1446125) q[3];
sx q[3];
rz(3.1015729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.17249168) q[2];
sx q[2];
rz(-1.8164941) q[2];
sx q[2];
rz(-1.8005523) q[2];
rz(-1.7155044) q[3];
sx q[3];
rz(-1.1172349) q[3];
sx q[3];
rz(-2.8360227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39389998) q[0];
sx q[0];
rz(-1.9992398) q[0];
sx q[0];
rz(0.72545141) q[0];
rz(-2.2254288) q[1];
sx q[1];
rz(-2.3326645) q[1];
sx q[1];
rz(-0.61942548) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7325798) q[0];
sx q[0];
rz(-2.469099) q[0];
sx q[0];
rz(2.2605863) q[0];
rz(-2.393707) q[2];
sx q[2];
rz(-2.0341784) q[2];
sx q[2];
rz(-3.0836251) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.20850764) q[1];
sx q[1];
rz(-1.2080482) q[1];
sx q[1];
rz(1.1136934) q[1];
x q[2];
rz(-0.28362579) q[3];
sx q[3];
rz(-2.0001634) q[3];
sx q[3];
rz(-0.46653433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.82975036) q[2];
sx q[2];
rz(-1.0210911) q[2];
sx q[2];
rz(-2.1799083) q[2];
rz(-1.1059777) q[3];
sx q[3];
rz(-1.032369) q[3];
sx q[3];
rz(-0.34524125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.85213929) q[0];
sx q[0];
rz(-2.0134605) q[0];
sx q[0];
rz(-1.7732675) q[0];
rz(-0.013956919) q[1];
sx q[1];
rz(-1.114926) q[1];
sx q[1];
rz(-1.7272635) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9326097) q[0];
sx q[0];
rz(-1.6609207) q[0];
sx q[0];
rz(1.8831253) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0704546) q[2];
sx q[2];
rz(-1.9390328) q[2];
sx q[2];
rz(2.8609543) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5217375) q[1];
sx q[1];
rz(-0.69112366) q[1];
sx q[1];
rz(2.28165) q[1];
rz(-0.39140626) q[3];
sx q[3];
rz(-1.8107171) q[3];
sx q[3];
rz(-2.4521884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.76325512) q[2];
sx q[2];
rz(-1.1474643) q[2];
sx q[2];
rz(2.9249127) q[2];
rz(-0.8832776) q[3];
sx q[3];
rz(-2.7933385) q[3];
sx q[3];
rz(-2.0047552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20275177) q[0];
sx q[0];
rz(-2.1348248) q[0];
sx q[0];
rz(2.3034565) q[0];
rz(-2.5996767) q[1];
sx q[1];
rz(-0.37626615) q[1];
sx q[1];
rz(3.0964877) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1208266) q[0];
sx q[0];
rz(-0.97310518) q[0];
sx q[0];
rz(-1.4106013) q[0];
rz(-2.0397268) q[2];
sx q[2];
rz(-0.53075889) q[2];
sx q[2];
rz(-1.6380978) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.939203) q[1];
sx q[1];
rz(-2.0430026) q[1];
sx q[1];
rz(0.07501988) q[1];
x q[2];
rz(0.49398936) q[3];
sx q[3];
rz(-1.4841379) q[3];
sx q[3];
rz(-2.983421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1230459) q[2];
sx q[2];
rz(-2.9728153) q[2];
sx q[2];
rz(0.65400845) q[2];
rz(0.6319913) q[3];
sx q[3];
rz(-1.4258823) q[3];
sx q[3];
rz(-0.33647195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80394799) q[0];
sx q[0];
rz(-2.4027282) q[0];
sx q[0];
rz(1.7919354) q[0];
rz(0.82756394) q[1];
sx q[1];
rz(-0.63500985) q[1];
sx q[1];
rz(0.48447022) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9512099) q[0];
sx q[0];
rz(-1.9007378) q[0];
sx q[0];
rz(1.953682) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4260277) q[2];
sx q[2];
rz(-0.55895222) q[2];
sx q[2];
rz(2.0367071) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1441298) q[1];
sx q[1];
rz(-1.9904549) q[1];
sx q[1];
rz(-1.4317572) q[1];
rz(-pi) q[2];
rz(2.6788533) q[3];
sx q[3];
rz(-0.92904186) q[3];
sx q[3];
rz(0.79824191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.95939246) q[2];
sx q[2];
rz(-1.0687989) q[2];
sx q[2];
rz(2.5051795) q[2];
rz(-3.0427129) q[3];
sx q[3];
rz(-1.5560047) q[3];
sx q[3];
rz(-1.6695361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(0.32475489) q[0];
sx q[0];
rz(-2.2519798) q[0];
sx q[0];
rz(-1.1725934) q[0];
rz(-1.532754) q[1];
sx q[1];
rz(-1.0120665) q[1];
sx q[1];
rz(-3.1408659) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6470203) q[0];
sx q[0];
rz(-1.1446409) q[0];
sx q[0];
rz(2.1151353) q[0];
rz(-pi) q[1];
rz(1.6696207) q[2];
sx q[2];
rz(-1.7903869) q[2];
sx q[2];
rz(-1.7864986) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8975413) q[1];
sx q[1];
rz(-1.6836299) q[1];
sx q[1];
rz(0.33399069) q[1];
rz(-pi) q[2];
x q[2];
rz(0.57925333) q[3];
sx q[3];
rz(-1.7369075) q[3];
sx q[3];
rz(0.36234713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.26649228) q[2];
sx q[2];
rz(-0.88783395) q[2];
sx q[2];
rz(-2.9317741) q[2];
rz(0.80875129) q[3];
sx q[3];
rz(-0.54976141) q[3];
sx q[3];
rz(-2.2426898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4631735) q[0];
sx q[0];
rz(-1.9371978) q[0];
sx q[0];
rz(3.0294982) q[0];
rz(-3.0922999) q[1];
sx q[1];
rz(-2.6105328) q[1];
sx q[1];
rz(1.8045527) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4508731) q[0];
sx q[0];
rz(-2.4639122) q[0];
sx q[0];
rz(1.9059794) q[0];
x q[1];
rz(-1.7981729) q[2];
sx q[2];
rz(-2.6223287) q[2];
sx q[2];
rz(2.8092172) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2087792) q[1];
sx q[1];
rz(-1.5808388) q[1];
sx q[1];
rz(-0.34681635) q[1];
x q[2];
rz(-0.2143799) q[3];
sx q[3];
rz(-1.0037546) q[3];
sx q[3];
rz(2.0585394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2992531) q[2];
sx q[2];
rz(-2.3167593) q[2];
sx q[2];
rz(0.71005589) q[2];
rz(3.1402785) q[3];
sx q[3];
rz(-2.240286) q[3];
sx q[3];
rz(-2.5406751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8173219) q[0];
sx q[0];
rz(-1.7621499) q[0];
sx q[0];
rz(-0.58260179) q[0];
rz(-2.461589) q[1];
sx q[1];
rz(-2.1425207) q[1];
sx q[1];
rz(-1.8780139) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.82915) q[0];
sx q[0];
rz(-2.0264605) q[0];
sx q[0];
rz(1.2553041) q[0];
x q[1];
rz(-2.3266257) q[2];
sx q[2];
rz(-2.8475347) q[2];
sx q[2];
rz(-2.5240099) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7192051) q[1];
sx q[1];
rz(-2.4610956) q[1];
sx q[1];
rz(2.8533555) q[1];
rz(3.0691759) q[3];
sx q[3];
rz(-2.7167746) q[3];
sx q[3];
rz(-2.2936898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6363643) q[2];
sx q[2];
rz(-1.4437081) q[2];
sx q[2];
rz(0.12602028) q[2];
rz(-1.6033008) q[3];
sx q[3];
rz(-2.3647629) q[3];
sx q[3];
rz(2.7974424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4890323) q[0];
sx q[0];
rz(-1.5131938) q[0];
sx q[0];
rz(1.162758) q[0];
rz(1.3105505) q[1];
sx q[1];
rz(-2.4030011) q[1];
sx q[1];
rz(1.383925) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2306128) q[0];
sx q[0];
rz(-0.34251102) q[0];
sx q[0];
rz(1.4369276) q[0];
rz(-0.53015401) q[2];
sx q[2];
rz(-2.9994476) q[2];
sx q[2];
rz(1.3719541) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.38720623) q[1];
sx q[1];
rz(-1.4547567) q[1];
sx q[1];
rz(0.95761489) q[1];
rz(-0.042785809) q[3];
sx q[3];
rz(-1.4322016) q[3];
sx q[3];
rz(2.3530172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.80013529) q[2];
sx q[2];
rz(-1.8693482) q[2];
sx q[2];
rz(-2.6696491) q[2];
rz(2.3927472) q[3];
sx q[3];
rz(-2.2180836) q[3];
sx q[3];
rz(-0.81764618) q[3];
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
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6227459) q[0];
sx q[0];
rz(-2.2846344) q[0];
sx q[0];
rz(-0.48883206) q[0];
rz(2.7986616) q[1];
sx q[1];
rz(-2.2551426) q[1];
sx q[1];
rz(-1.0541281) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80741548) q[0];
sx q[0];
rz(-1.6082967) q[0];
sx q[0];
rz(1.5192017) q[0];
rz(-pi) q[1];
rz(-1.2230049) q[2];
sx q[2];
rz(-0.985983) q[2];
sx q[2];
rz(-2.0690167) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.17133896) q[1];
sx q[1];
rz(-1.043817) q[1];
sx q[1];
rz(-0.33895609) q[1];
rz(3.0929186) q[3];
sx q[3];
rz(-1.0992388) q[3];
sx q[3];
rz(2.8296628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9752581) q[2];
sx q[2];
rz(-2.9322093) q[2];
sx q[2];
rz(3.0496822) q[2];
rz(-2.6722243) q[3];
sx q[3];
rz(-2.2769603) q[3];
sx q[3];
rz(-3.0288467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
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
rz(-2.3757833) q[1];
sx q[1];
rz(-1.3594834) q[1];
sx q[1];
rz(1.8654738) q[1];
rz(-1.2086537) q[2];
sx q[2];
rz(-1.8100783) q[2];
sx q[2];
rz(1.7892005) q[2];
rz(2.7226483) q[3];
sx q[3];
rz(-2.5635515) q[3];
sx q[3];
rz(0.10684914) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
