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
rz(-0.80225575) q[0];
sx q[0];
rz(-1.7576317) q[0];
sx q[0];
rz(1.2686165) q[0];
rz(0.4624548) q[1];
sx q[1];
rz(6.3422536) q[1];
sx q[1];
rz(8.0191945) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7168306) q[0];
sx q[0];
rz(-2.4772236) q[0];
sx q[0];
rz(-0.82358255) q[0];
rz(-pi) q[1];
rz(1.8008158) q[2];
sx q[2];
rz(-0.91311306) q[2];
sx q[2];
rz(-2.7000526) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.66623233) q[1];
sx q[1];
rz(-1.1603429) q[1];
sx q[1];
rz(-0.79373856) q[1];
rz(2.2470993) q[3];
sx q[3];
rz(-1.3925806) q[3];
sx q[3];
rz(1.7591998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1543701) q[2];
sx q[2];
rz(-1.1673678) q[2];
sx q[2];
rz(-2.3517189) q[2];
rz(-0.023898276) q[3];
sx q[3];
rz(-1.8022715) q[3];
sx q[3];
rz(-2.6051615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24870482) q[0];
sx q[0];
rz(-1.4905812) q[0];
sx q[0];
rz(-1.254068) q[0];
rz(-1.6528543) q[1];
sx q[1];
rz(-2.2658927) q[1];
sx q[1];
rz(-0.48315963) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41071677) q[0];
sx q[0];
rz(-0.19598254) q[0];
sx q[0];
rz(1.3021126) q[0];
rz(-pi) q[1];
rz(0.56098361) q[2];
sx q[2];
rz(-2.3234754) q[2];
sx q[2];
rz(-2.4647692) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3846783) q[1];
sx q[1];
rz(-1.9325629) q[1];
sx q[1];
rz(1.8513308) q[1];
rz(1.3951357) q[3];
sx q[3];
rz(-1.8570199) q[3];
sx q[3];
rz(-1.3940879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2529605) q[2];
sx q[2];
rz(-1.4560207) q[2];
sx q[2];
rz(-1.6271094) q[2];
rz(-0.0557946) q[3];
sx q[3];
rz(-1.5443708) q[3];
sx q[3];
rz(1.6519215) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91419879) q[0];
sx q[0];
rz(-2.6970503) q[0];
sx q[0];
rz(3.0134873) q[0];
rz(2.5654492) q[1];
sx q[1];
rz(-0.73155254) q[1];
sx q[1];
rz(-2.3760956) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.573581) q[0];
sx q[0];
rz(-0.67493248) q[0];
sx q[0];
rz(-2.5715067) q[0];
rz(-pi) q[1];
rz(-1.9472935) q[2];
sx q[2];
rz(-0.91937477) q[2];
sx q[2];
rz(-2.5848856) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1372337) q[1];
sx q[1];
rz(-0.69789125) q[1];
sx q[1];
rz(2.3801719) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5452905) q[3];
sx q[3];
rz(-1.4273941) q[3];
sx q[3];
rz(-2.0633782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3196044) q[2];
sx q[2];
rz(-0.87279785) q[2];
sx q[2];
rz(-0.86307159) q[2];
rz(0.34267628) q[3];
sx q[3];
rz(-2.7211012) q[3];
sx q[3];
rz(-1.6531403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2494025) q[0];
sx q[0];
rz(-0.96977314) q[0];
sx q[0];
rz(-2.6595907) q[0];
rz(2.8328698) q[1];
sx q[1];
rz(-0.32544193) q[1];
sx q[1];
rz(0.6122922) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4929313) q[0];
sx q[0];
rz(-1.0928286) q[0];
sx q[0];
rz(2.1532691) q[0];
rz(-3.127549) q[2];
sx q[2];
rz(-1.4870475) q[2];
sx q[2];
rz(1.3824826) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0902113) q[1];
sx q[1];
rz(-1.5259698) q[1];
sx q[1];
rz(-0.96446891) q[1];
rz(-pi) q[2];
rz(-3.0975625) q[3];
sx q[3];
rz(-0.14942871) q[3];
sx q[3];
rz(-1.3469411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4838532) q[2];
sx q[2];
rz(-2.190399) q[2];
sx q[2];
rz(1.308002) q[2];
rz(-0.43404239) q[3];
sx q[3];
rz(-1.7900034) q[3];
sx q[3];
rz(-1.2691809) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9798715) q[0];
sx q[0];
rz(-1.5890108) q[0];
sx q[0];
rz(-0.27444926) q[0];
rz(-1.5212003) q[1];
sx q[1];
rz(-1.043964) q[1];
sx q[1];
rz(1.257198) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8306106) q[0];
sx q[0];
rz(-2.794696) q[0];
sx q[0];
rz(0.78126379) q[0];
rz(0.99165708) q[2];
sx q[2];
rz(-0.61505396) q[2];
sx q[2];
rz(-1.8293387) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2272531) q[1];
sx q[1];
rz(-1.7861331) q[1];
sx q[1];
rz(2.1244177) q[1];
rz(-2.8207079) q[3];
sx q[3];
rz(-2.2498331) q[3];
sx q[3];
rz(-1.8378375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3955128) q[2];
sx q[2];
rz(-1.6359676) q[2];
sx q[2];
rz(0.5557605) q[2];
rz(-2.3505576) q[3];
sx q[3];
rz(-1.0020703) q[3];
sx q[3];
rz(1.5033495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.898734) q[0];
sx q[0];
rz(-1.7019615) q[0];
sx q[0];
rz(-2.4295501) q[0];
rz(0.60246077) q[1];
sx q[1];
rz(-2.7280877) q[1];
sx q[1];
rz(0.079040225) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1940546) q[0];
sx q[0];
rz(-1.6350766) q[0];
sx q[0];
rz(-3.0966731) q[0];
rz(-pi) q[1];
rz(2.2964962) q[2];
sx q[2];
rz(-0.11339408) q[2];
sx q[2];
rz(1.5886943) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.81198967) q[1];
sx q[1];
rz(-1.7895849) q[1];
sx q[1];
rz(-2.3892774) q[1];
rz(1.6738123) q[3];
sx q[3];
rz(-1.8673737) q[3];
sx q[3];
rz(2.1144583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1003803) q[2];
sx q[2];
rz(-2.2104287) q[2];
sx q[2];
rz(2.1616914) q[2];
rz(1.6117217) q[3];
sx q[3];
rz(-1.2414705) q[3];
sx q[3];
rz(-0.25947586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-0.74705446) q[0];
sx q[0];
rz(-1.8600445) q[0];
sx q[0];
rz(-3.0406612) q[0];
rz(-2.586567) q[1];
sx q[1];
rz(-0.9340159) q[1];
sx q[1];
rz(1.2947882) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0880343) q[0];
sx q[0];
rz(-0.71672601) q[0];
sx q[0];
rz(1.2117366) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.567595) q[2];
sx q[2];
rz(-0.29433196) q[2];
sx q[2];
rz(-0.41121361) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.76483549) q[1];
sx q[1];
rz(-1.7143814) q[1];
sx q[1];
rz(-0.43068703) q[1];
rz(-pi) q[2];
rz(-0.28980589) q[3];
sx q[3];
rz(-2.0103526) q[3];
sx q[3];
rz(1.8686779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1780221) q[2];
sx q[2];
rz(-2.2691085) q[2];
sx q[2];
rz(-2.8998609) q[2];
rz(-2.669615) q[3];
sx q[3];
rz(-0.21502544) q[3];
sx q[3];
rz(-1.5579461) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80411512) q[0];
sx q[0];
rz(-2.1601456) q[0];
sx q[0];
rz(0.61019439) q[0];
rz(-2.923851) q[1];
sx q[1];
rz(-2.1861031) q[1];
sx q[1];
rz(-2.311923) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.037704) q[0];
sx q[0];
rz(-1.0173326) q[0];
sx q[0];
rz(-2.3174556) q[0];
x q[1];
rz(2.3102364) q[2];
sx q[2];
rz(-2.6717253) q[2];
sx q[2];
rz(0.63127764) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.024914537) q[1];
sx q[1];
rz(-0.45104154) q[1];
sx q[1];
rz(-2.7052959) q[1];
rz(-pi) q[2];
x q[2];
rz(0.85932087) q[3];
sx q[3];
rz(-1.4830657) q[3];
sx q[3];
rz(0.52857196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3024451) q[2];
sx q[2];
rz(-1.7586917) q[2];
sx q[2];
rz(-2.006532) q[2];
rz(-0.48842397) q[3];
sx q[3];
rz(-1.4819744) q[3];
sx q[3];
rz(-1.9058913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14588533) q[0];
sx q[0];
rz(-2.0304401) q[0];
sx q[0];
rz(-0.97440326) q[0];
rz(-2.3459332) q[1];
sx q[1];
rz(-2.7641422) q[1];
sx q[1];
rz(0.61765751) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5159588) q[0];
sx q[0];
rz(-0.083481073) q[0];
sx q[0];
rz(0.030103695) q[0];
rz(-pi) q[1];
rz(1.1716003) q[2];
sx q[2];
rz(-1.079487) q[2];
sx q[2];
rz(-1.0688865) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1682081) q[1];
sx q[1];
rz(-2.1323556) q[1];
sx q[1];
rz(1.2513158) q[1];
rz(-2.0418725) q[3];
sx q[3];
rz(-1.866939) q[3];
sx q[3];
rz(2.5961034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1254897) q[2];
sx q[2];
rz(-1.8718655) q[2];
sx q[2];
rz(1.3886836) q[2];
rz(-1.8359418) q[3];
sx q[3];
rz(-0.7235705) q[3];
sx q[3];
rz(-1.2828264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1284803) q[0];
sx q[0];
rz(-2.4025669) q[0];
sx q[0];
rz(-0.59984961) q[0];
rz(-0.57303095) q[1];
sx q[1];
rz(-0.60879469) q[1];
sx q[1];
rz(-1.0240239) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3442992) q[0];
sx q[0];
rz(-1.1670207) q[0];
sx q[0];
rz(2.7112194) q[0];
rz(-pi) q[1];
rz(-0.055209514) q[2];
sx q[2];
rz(-2.488236) q[2];
sx q[2];
rz(-3.0467767) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9309195) q[1];
sx q[1];
rz(-1.1751047) q[1];
sx q[1];
rz(-1.6617421) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7217595) q[3];
sx q[3];
rz(-1.1929885) q[3];
sx q[3];
rz(-1.0491766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2188501) q[2];
sx q[2];
rz(-1.3877733) q[2];
sx q[2];
rz(-2.3663523) q[2];
rz(1.5357337) q[3];
sx q[3];
rz(-1.1865059) q[3];
sx q[3];
rz(1.8839914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53104644) q[0];
sx q[0];
rz(-0.95120593) q[0];
sx q[0];
rz(-1.9134941) q[0];
rz(1.7465406) q[1];
sx q[1];
rz(-1.4735305) q[1];
sx q[1];
rz(2.7458618) q[1];
rz(-1.2634696) q[2];
sx q[2];
rz(-2.0057445) q[2];
sx q[2];
rz(2.4705171) q[2];
rz(-0.80202924) q[3];
sx q[3];
rz(-1.5939972) q[3];
sx q[3];
rz(0.93929285) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
