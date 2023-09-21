OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.80823094) q[0];
sx q[0];
rz(-2.1699177) q[0];
sx q[0];
rz(1.4810286) q[0];
rz(3.0193168) q[1];
sx q[1];
rz(-3.0552157) q[1];
sx q[1];
rz(3.1187305) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6883144) q[0];
sx q[0];
rz(-2.9524321) q[0];
sx q[0];
rz(2.2384089) q[0];
rz(-pi) q[1];
rz(0.51486751) q[2];
sx q[2];
rz(-1.1866743) q[2];
sx q[2];
rz(1.9799973) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.37576807) q[1];
sx q[1];
rz(-1.4210912) q[1];
sx q[1];
rz(-1.2234115) q[1];
rz(-pi) q[2];
x q[2];
rz(1.919235) q[3];
sx q[3];
rz(-2.493353) q[3];
sx q[3];
rz(-0.1227619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7227398) q[2];
sx q[2];
rz(-0.8272233) q[2];
sx q[2];
rz(-2.375405) q[2];
rz(0.17151672) q[3];
sx q[3];
rz(-2.4092716) q[3];
sx q[3];
rz(-2.5550301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3294753) q[0];
sx q[0];
rz(-0.43918878) q[0];
sx q[0];
rz(3.0549333) q[0];
rz(0.86241972) q[1];
sx q[1];
rz(-2.5984867) q[1];
sx q[1];
rz(0.08509732) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8893196) q[0];
sx q[0];
rz(-2.5190881) q[0];
sx q[0];
rz(-0.028153367) q[0];
x q[1];
rz(-2.1727824) q[2];
sx q[2];
rz(-1.3170871) q[2];
sx q[2];
rz(1.3041376) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3955235) q[1];
sx q[1];
rz(-1.6568526) q[1];
sx q[1];
rz(-1.458191) q[1];
x q[2];
rz(-0.99382932) q[3];
sx q[3];
rz(-2.5964063) q[3];
sx q[3];
rz(2.546437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8791085) q[2];
sx q[2];
rz(-1.5510677) q[2];
sx q[2];
rz(-1.6764199) q[2];
rz(1.9654467) q[3];
sx q[3];
rz(-0.90548235) q[3];
sx q[3];
rz(0.24648497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6572606) q[0];
sx q[0];
rz(-0.13467877) q[0];
sx q[0];
rz(2.235967) q[0];
rz(-2.8088645) q[1];
sx q[1];
rz(-0.85488027) q[1];
sx q[1];
rz(1.7571626) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4577643) q[0];
sx q[0];
rz(-2.0450852) q[0];
sx q[0];
rz(1.477924) q[0];
rz(-0.39321123) q[2];
sx q[2];
rz(-1.097743) q[2];
sx q[2];
rz(-2.2676603) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.78450173) q[1];
sx q[1];
rz(-2.38584) q[1];
sx q[1];
rz(3.109039) q[1];
rz(-pi) q[2];
rz(-2.609385) q[3];
sx q[3];
rz(-1.0695219) q[3];
sx q[3];
rz(1.7621079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.24421346) q[2];
sx q[2];
rz(-2.801557) q[2];
sx q[2];
rz(1.8133694) q[2];
rz(-1.7437079) q[3];
sx q[3];
rz(-0.54026794) q[3];
sx q[3];
rz(3.0298997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.06552799) q[0];
sx q[0];
rz(-2.9632443) q[0];
sx q[0];
rz(-0.67101014) q[0];
rz(-1.3440075) q[1];
sx q[1];
rz(-1.9799415) q[1];
sx q[1];
rz(2.9342594) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5642693) q[0];
sx q[0];
rz(-0.042357001) q[0];
sx q[0];
rz(-1.7839412) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7423332) q[2];
sx q[2];
rz(-0.17645391) q[2];
sx q[2];
rz(-1.8432957) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.85992766) q[1];
sx q[1];
rz(-1.7237701) q[1];
sx q[1];
rz(1.8106736) q[1];
x q[2];
rz(0.025709318) q[3];
sx q[3];
rz(-0.51179143) q[3];
sx q[3];
rz(0.6230841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8308782) q[2];
sx q[2];
rz(-1.0881311) q[2];
sx q[2];
rz(2.3332398) q[2];
rz(0.80777848) q[3];
sx q[3];
rz(-0.72179663) q[3];
sx q[3];
rz(-2.9479153) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5508674) q[0];
sx q[0];
rz(-0.71722513) q[0];
sx q[0];
rz(-2.2461058) q[0];
rz(-2.8764309) q[1];
sx q[1];
rz(-0.83509713) q[1];
sx q[1];
rz(-0.53363824) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95414872) q[0];
sx q[0];
rz(-1.2866409) q[0];
sx q[0];
rz(2.8548293) q[0];
rz(2.3737714) q[2];
sx q[2];
rz(-1.7349093) q[2];
sx q[2];
rz(0.62190157) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.67709778) q[1];
sx q[1];
rz(-0.2650731) q[1];
sx q[1];
rz(-0.58880834) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.93612) q[3];
sx q[3];
rz(-2.1901202) q[3];
sx q[3];
rz(-1.6831786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9020033) q[2];
sx q[2];
rz(-1.5782372) q[2];
sx q[2];
rz(0.56838244) q[2];
rz(-1.2166294) q[3];
sx q[3];
rz(-2.670848) q[3];
sx q[3];
rz(0.34354982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6333273) q[0];
sx q[0];
rz(-0.87843043) q[0];
sx q[0];
rz(-1.1761965) q[0];
rz(1.5559224) q[1];
sx q[1];
rz(-2.1891179) q[1];
sx q[1];
rz(-1.0046545) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83197901) q[0];
sx q[0];
rz(-1.2011315) q[0];
sx q[0];
rz(1.2020822) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7447694) q[2];
sx q[2];
rz(-1.6453711) q[2];
sx q[2];
rz(-2.4900988) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1130044) q[1];
sx q[1];
rz(-2.5341946) q[1];
sx q[1];
rz(2.7007156) q[1];
rz(1.5435013) q[3];
sx q[3];
rz(-2.7276162) q[3];
sx q[3];
rz(-2.9881145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4040318) q[2];
sx q[2];
rz(-2.1514251) q[2];
sx q[2];
rz(0.44815865) q[2];
rz(0.29799497) q[3];
sx q[3];
rz(-2.4500676) q[3];
sx q[3];
rz(0.34860778) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3649243) q[0];
sx q[0];
rz(-1.3013327) q[0];
sx q[0];
rz(-0.65814322) q[0];
rz(-1.2843885) q[1];
sx q[1];
rz(-2.6874266) q[1];
sx q[1];
rz(-0.67214322) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2015398) q[0];
sx q[0];
rz(-1.7367898) q[0];
sx q[0];
rz(-0.0075017651) q[0];
x q[1];
rz(-2.3357046) q[2];
sx q[2];
rz(-2.561509) q[2];
sx q[2];
rz(-2.9758331) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8183648) q[1];
sx q[1];
rz(-1.0867026) q[1];
sx q[1];
rz(-1.9869884) q[1];
rz(2.4403205) q[3];
sx q[3];
rz(-1.6810732) q[3];
sx q[3];
rz(0.14227223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.60823524) q[2];
sx q[2];
rz(-0.63295263) q[2];
sx q[2];
rz(2.7427924) q[2];
rz(0.82459015) q[3];
sx q[3];
rz(-1.3922858) q[3];
sx q[3];
rz(-0.59857541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7010715) q[0];
sx q[0];
rz(-0.43061391) q[0];
sx q[0];
rz(-2.9833802) q[0];
rz(2.054706) q[1];
sx q[1];
rz(-2.465076) q[1];
sx q[1];
rz(-2.3349082) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0772484) q[0];
sx q[0];
rz(-1.4532538) q[0];
sx q[0];
rz(-2.006152) q[0];
rz(0.20571795) q[2];
sx q[2];
rz(-1.6763902) q[2];
sx q[2];
rz(-0.92668698) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.73831564) q[1];
sx q[1];
rz(-2.3936831) q[1];
sx q[1];
rz(0.43055375) q[1];
rz(1.8716378) q[3];
sx q[3];
rz(-1.5402208) q[3];
sx q[3];
rz(-1.9105063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.57758254) q[2];
sx q[2];
rz(-1.492604) q[2];
sx q[2];
rz(-0.99964833) q[2];
rz(-0.10351652) q[3];
sx q[3];
rz(-3.0045356) q[3];
sx q[3];
rz(-1.1243533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0209811) q[0];
sx q[0];
rz(-0.7779026) q[0];
sx q[0];
rz(-0.65752423) q[0];
rz(2.855037) q[1];
sx q[1];
rz(-0.92266881) q[1];
sx q[1];
rz(0.34067571) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17250241) q[0];
sx q[0];
rz(-1.889125) q[0];
sx q[0];
rz(-2.3031635) q[0];
rz(-1.0057783) q[2];
sx q[2];
rz(-1.3823595) q[2];
sx q[2];
rz(-2.4609158) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.454168) q[1];
sx q[1];
rz(-1.7523645) q[1];
sx q[1];
rz(-2.1329692) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9433603) q[3];
sx q[3];
rz(-1.7361589) q[3];
sx q[3];
rz(-2.4076622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3866773) q[2];
sx q[2];
rz(-1.9415386) q[2];
sx q[2];
rz(-2.6436451) q[2];
rz(0.30161101) q[3];
sx q[3];
rz(-2.7409654) q[3];
sx q[3];
rz(-2.6224459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72934812) q[0];
sx q[0];
rz(-3.0817741) q[0];
sx q[0];
rz(2.8630032) q[0];
rz(-2.5623698) q[1];
sx q[1];
rz(-2.2021553) q[1];
sx q[1];
rz(-3.0648807) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3007606) q[0];
sx q[0];
rz(-1.6941841) q[0];
sx q[0];
rz(-1.2542017) q[0];
rz(2.1439728) q[2];
sx q[2];
rz(-1.4263037) q[2];
sx q[2];
rz(-1.0695374) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1288651) q[1];
sx q[1];
rz(-2.7755133) q[1];
sx q[1];
rz(-0.058458316) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4823227) q[3];
sx q[3];
rz(-1.3946597) q[3];
sx q[3];
rz(0.74151553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.77438337) q[2];
sx q[2];
rz(-2.4343906) q[2];
sx q[2];
rz(-0.45551604) q[2];
rz(2.7101743) q[3];
sx q[3];
rz(-3.016267) q[3];
sx q[3];
rz(3.07807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1800304) q[0];
sx q[0];
rz(-1.7407692) q[0];
sx q[0];
rz(-2.2041007) q[0];
rz(0.58615276) q[1];
sx q[1];
rz(-1.502232) q[1];
sx q[1];
rz(-1.4486817) q[1];
rz(2.0942681) q[2];
sx q[2];
rz(-2.6101255) q[2];
sx q[2];
rz(2.5892467) q[2];
rz(-0.42701957) q[3];
sx q[3];
rz(-1.2357161) q[3];
sx q[3];
rz(-1.2824035) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];