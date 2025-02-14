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
rz(1.8279561) q[0];
sx q[0];
rz(6.2160563) q[0];
sx q[0];
rz(9.4475533) q[0];
rz(1.0442806) q[1];
sx q[1];
rz(-2.2106946) q[1];
sx q[1];
rz(-3.1276303) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21101306) q[0];
sx q[0];
rz(-1.530181) q[0];
sx q[0];
rz(0.75675772) q[0];
rz(0.44273744) q[2];
sx q[2];
rz(-1.1907026) q[2];
sx q[2];
rz(2.1815104) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3252445) q[1];
sx q[1];
rz(-2.8720586) q[1];
sx q[1];
rz(2.4052038) q[1];
rz(-pi) q[2];
rz(0.45043378) q[3];
sx q[3];
rz(-1.2518468) q[3];
sx q[3];
rz(-1.3802841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.969101) q[2];
sx q[2];
rz(-1.3250985) q[2];
sx q[2];
rz(1.3410404) q[2];
rz(-1.4260882) q[3];
sx q[3];
rz(-1.1172349) q[3];
sx q[3];
rz(2.8360227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7476927) q[0];
sx q[0];
rz(-1.9992398) q[0];
sx q[0];
rz(-0.72545141) q[0];
rz(0.91616383) q[1];
sx q[1];
rz(-0.80892816) q[1];
sx q[1];
rz(0.61942548) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5529157) q[0];
sx q[0];
rz(-1.9784133) q[0];
sx q[0];
rz(2.1216395) q[0];
rz(-pi) q[1];
x q[1];
rz(2.393707) q[2];
sx q[2];
rz(-1.1074142) q[2];
sx q[2];
rz(0.057967535) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.933085) q[1];
sx q[1];
rz(-1.9335445) q[1];
sx q[1];
rz(2.0278992) q[1];
rz(-pi) q[2];
rz(2.8579669) q[3];
sx q[3];
rz(-2.0001634) q[3];
sx q[3];
rz(-0.46653433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.82975036) q[2];
sx q[2];
rz(-1.0210911) q[2];
sx q[2];
rz(-0.96168438) q[2];
rz(2.035615) q[3];
sx q[3];
rz(-1.032369) q[3];
sx q[3];
rz(2.7963514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2894534) q[0];
sx q[0];
rz(-2.0134605) q[0];
sx q[0];
rz(1.7732675) q[0];
rz(3.1276357) q[1];
sx q[1];
rz(-1.114926) q[1];
sx q[1];
rz(1.4143292) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0516134) q[0];
sx q[0];
rz(-0.3246626) q[0];
sx q[0];
rz(1.8568296) q[0];
rz(1.9398841) q[2];
sx q[2];
rz(-1.5044343) q[2];
sx q[2];
rz(1.8770799) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5217375) q[1];
sx q[1];
rz(-2.450469) q[1];
sx q[1];
rz(0.85994263) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5713874) q[3];
sx q[3];
rz(-0.45582882) q[3];
sx q[3];
rz(-2.7826235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.76325512) q[2];
sx q[2];
rz(-1.1474643) q[2];
sx q[2];
rz(-2.9249127) q[2];
rz(-2.2583151) q[3];
sx q[3];
rz(-2.7933385) q[3];
sx q[3];
rz(2.0047552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9388409) q[0];
sx q[0];
rz(-1.0067679) q[0];
sx q[0];
rz(-0.83813611) q[0];
rz(-0.54191598) q[1];
sx q[1];
rz(-0.37626615) q[1];
sx q[1];
rz(0.045104973) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7411557) q[0];
sx q[0];
rz(-2.525345) q[0];
sx q[0];
rz(-0.23016302) q[0];
rz(-0.25928478) q[2];
sx q[2];
rz(-1.102299) q[2];
sx q[2];
rz(1.1069654) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1029658) q[1];
sx q[1];
rz(-0.47768394) q[1];
sx q[1];
rz(1.716502) q[1];
x q[2];
rz(-1.4724422) q[3];
sx q[3];
rz(-1.0788267) q[3];
sx q[3];
rz(1.4592001) q[3];
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
rz(-2.4875842) q[2];
rz(2.5096014) q[3];
sx q[3];
rz(-1.4258823) q[3];
sx q[3];
rz(0.33647195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80394799) q[0];
sx q[0];
rz(-2.4027282) q[0];
sx q[0];
rz(-1.3496572) q[0];
rz(-2.3140287) q[1];
sx q[1];
rz(-0.63500985) q[1];
sx q[1];
rz(-2.6571224) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.890936) q[0];
sx q[0];
rz(-1.93205) q[0];
sx q[0];
rz(2.7879232) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1250317) q[2];
sx q[2];
rz(-1.6473738) q[2];
sx q[2];
rz(0.34293338) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1441298) q[1];
sx q[1];
rz(-1.9904549) q[1];
sx q[1];
rz(1.7098355) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6788533) q[3];
sx q[3];
rz(-2.2125508) q[3];
sx q[3];
rz(2.3433507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1822002) q[2];
sx q[2];
rz(-1.0687989) q[2];
sx q[2];
rz(2.5051795) q[2];
rz(3.0427129) q[3];
sx q[3];
rz(-1.585588) q[3];
sx q[3];
rz(1.4720565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.8168378) q[0];
sx q[0];
rz(-0.88961283) q[0];
sx q[0];
rz(-1.1725934) q[0];
rz(1.6088387) q[1];
sx q[1];
rz(-2.1295261) q[1];
sx q[1];
rz(-0.00072678725) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1689735) q[0];
sx q[0];
rz(-1.0797636) q[0];
sx q[0];
rz(0.48788496) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6696207) q[2];
sx q[2];
rz(-1.7903869) q[2];
sx q[2];
rz(1.7864986) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.24405133) q[1];
sx q[1];
rz(-1.6836299) q[1];
sx q[1];
rz(2.807602) q[1];
rz(-pi) q[2];
rz(2.5623393) q[3];
sx q[3];
rz(-1.4046852) q[3];
sx q[3];
rz(0.36234713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8751004) q[2];
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
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(1.6784191) q[0];
sx q[0];
rz(-1.2043948) q[0];
sx q[0];
rz(-0.11209442) q[0];
rz(3.0922999) q[1];
sx q[1];
rz(-2.6105328) q[1];
sx q[1];
rz(1.3370399) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69071952) q[0];
sx q[0];
rz(-2.4639122) q[0];
sx q[0];
rz(1.2356133) q[0];
rz(-pi) q[1];
x q[1];
rz(2.078901) q[2];
sx q[2];
rz(-1.6828949) q[2];
sx q[2];
rz(-1.0401806) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9328135) q[1];
sx q[1];
rz(-1.5808388) q[1];
sx q[1];
rz(-0.34681635) q[1];
rz(-pi) q[2];
rz(2.1483809) q[3];
sx q[3];
rz(-1.3903729) q[3];
sx q[3];
rz(-2.5374295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2992531) q[2];
sx q[2];
rz(-0.82483333) q[2];
sx q[2];
rz(-0.71005589) q[2];
rz(-0.0013141343) q[3];
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
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8173219) q[0];
sx q[0];
rz(-1.3794427) q[0];
sx q[0];
rz(0.58260179) q[0];
rz(0.6800037) q[1];
sx q[1];
rz(-2.1425207) q[1];
sx q[1];
rz(-1.8780139) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9506142) q[0];
sx q[0];
rz(-0.54784566) q[0];
sx q[0];
rz(-2.5771499) q[0];
rz(-pi) q[1];
rz(-2.3266257) q[2];
sx q[2];
rz(-0.294058) q[2];
sx q[2];
rz(-0.61758274) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.78681011) q[1];
sx q[1];
rz(-0.92325961) q[1];
sx q[1];
rz(-1.7969653) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7177592) q[3];
sx q[3];
rz(-1.5409711) q[3];
sx q[3];
rz(-0.78889293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6363643) q[2];
sx q[2];
rz(-1.6978846) q[2];
sx q[2];
rz(-0.12602028) q[2];
rz(-1.6033008) q[3];
sx q[3];
rz(-0.77682972) q[3];
sx q[3];
rz(0.34415027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-0.4890323) q[0];
sx q[0];
rz(-1.6283988) q[0];
sx q[0];
rz(-1.9788347) q[0];
rz(-1.8310422) q[1];
sx q[1];
rz(-0.73859155) q[1];
sx q[1];
rz(1.7576677) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2306128) q[0];
sx q[0];
rz(-2.7990816) q[0];
sx q[0];
rz(-1.7046651) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4985563) q[2];
sx q[2];
rz(-1.4482699) q[2];
sx q[2];
rz(-0.83736698) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1022998) q[1];
sx q[1];
rz(-0.96233923) q[1];
sx q[1];
rz(3.0000172) q[1];
rz(-1.7095164) q[3];
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
rz(-2.3414574) q[2];
sx q[2];
rz(-1.8693482) q[2];
sx q[2];
rz(0.47194353) q[2];
rz(-2.3927472) q[3];
sx q[3];
rz(-2.2180836) q[3];
sx q[3];
rz(0.81764618) q[3];
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
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6227459) q[0];
sx q[0];
rz(-0.8569583) q[0];
sx q[0];
rz(2.6527606) q[0];
rz(-2.7986616) q[1];
sx q[1];
rz(-2.2551426) q[1];
sx q[1];
rz(1.0541281) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7502405) q[0];
sx q[0];
rz(-3.0778193) q[0];
sx q[0];
rz(-0.94185726) q[0];
rz(1.9185877) q[2];
sx q[2];
rz(-0.985983) q[2];
sx q[2];
rz(1.0725759) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.78275241) q[1];
sx q[1];
rz(-0.61780518) q[1];
sx q[1];
rz(-2.0899523) q[1];
rz(-pi) q[2];
rz(-0.048674095) q[3];
sx q[3];
rz(-2.0423539) q[3];
sx q[3];
rz(0.31192985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9752581) q[2];
sx q[2];
rz(-0.20938337) q[2];
sx q[2];
rz(3.0496822) q[2];
rz(-0.46936834) q[3];
sx q[3];
rz(-0.86463237) q[3];
sx q[3];
rz(0.11274591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5975006) q[0];
sx q[0];
rz(-1.8202029) q[0];
sx q[0];
rz(1.2910917) q[0];
rz(-0.76580936) q[1];
sx q[1];
rz(-1.7821093) q[1];
sx q[1];
rz(-1.2761188) q[1];
rz(0.2551887) q[2];
sx q[2];
rz(-1.219426) q[2];
sx q[2];
rz(0.30795369) q[2];
rz(-0.41894434) q[3];
sx q[3];
rz(-2.5635515) q[3];
sx q[3];
rz(0.10684914) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
