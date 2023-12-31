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
rz(4.1132676) q[0];
sx q[0];
rz(10.905807) q[0];
rz(-0.12227585) q[1];
sx q[1];
rz(-0.08637698) q[1];
sx q[1];
rz(-3.1187305) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6883144) q[0];
sx q[0];
rz(-2.9524321) q[0];
sx q[0];
rz(0.90318371) q[0];
rz(-pi) q[1];
x q[1];
rz(0.51486751) q[2];
sx q[2];
rz(-1.1866743) q[2];
sx q[2];
rz(-1.1615953) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.37576807) q[1];
sx q[1];
rz(-1.7205015) q[1];
sx q[1];
rz(-1.9181812) q[1];
rz(-1.919235) q[3];
sx q[3];
rz(-0.64823965) q[3];
sx q[3];
rz(-0.1227619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7227398) q[2];
sx q[2];
rz(-0.8272233) q[2];
sx q[2];
rz(2.375405) q[2];
rz(0.17151672) q[3];
sx q[3];
rz(-2.4092716) q[3];
sx q[3];
rz(0.58656251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3294753) q[0];
sx q[0];
rz(-2.7024039) q[0];
sx q[0];
rz(0.086659327) q[0];
rz(2.2791729) q[1];
sx q[1];
rz(-2.5984867) q[1];
sx q[1];
rz(-0.08509732) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21762411) q[0];
sx q[0];
rz(-0.94857615) q[0];
sx q[0];
rz(-1.5909965) q[0];
x q[1];
rz(2.8367963) q[2];
sx q[2];
rz(-2.1509503) q[2];
sx q[2];
rz(3.0457029) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.47463402) q[1];
sx q[1];
rz(-2.9999795) q[1];
sx q[1];
rz(0.9160362) q[1];
rz(-2.1477633) q[3];
sx q[3];
rz(-0.54518632) q[3];
sx q[3];
rz(-0.59515566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8791085) q[2];
sx q[2];
rz(-1.5510677) q[2];
sx q[2];
rz(1.6764199) q[2];
rz(-1.1761459) q[3];
sx q[3];
rz(-0.90548235) q[3];
sx q[3];
rz(0.24648497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48433205) q[0];
sx q[0];
rz(-0.13467877) q[0];
sx q[0];
rz(0.9056257) q[0];
rz(-2.8088645) q[1];
sx q[1];
rz(-2.2867124) q[1];
sx q[1];
rz(1.3844301) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15554409) q[0];
sx q[0];
rz(-1.6533924) q[0];
sx q[0];
rz(-2.6655469) q[0];
rz(-2.076782) q[2];
sx q[2];
rz(-1.9188606) q[2];
sx q[2];
rz(2.6315174) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3789931) q[1];
sx q[1];
rz(-1.5931207) q[1];
sx q[1];
rz(0.75548817) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1372041) q[3];
sx q[3];
rz(-2.0319788) q[3];
sx q[3];
rz(0.084463488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.24421346) q[2];
sx q[2];
rz(-0.34003568) q[2];
sx q[2];
rz(-1.8133694) q[2];
rz(-1.3978847) q[3];
sx q[3];
rz(-2.6013247) q[3];
sx q[3];
rz(3.0298997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-3.0760647) q[0];
sx q[0];
rz(-2.9632443) q[0];
sx q[0];
rz(0.67101014) q[0];
rz(1.3440075) q[1];
sx q[1];
rz(-1.9799415) q[1];
sx q[1];
rz(-2.9342594) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3509388) q[0];
sx q[0];
rz(-1.5293984) q[0];
sx q[0];
rz(-3.1326276) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3992594) q[2];
sx q[2];
rz(-0.17645391) q[2];
sx q[2];
rz(-1.2982969) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9878848) q[1];
sx q[1];
rz(-2.8578836) q[1];
sx q[1];
rz(-0.99516408) q[1];
rz(-pi) q[2];
rz(-0.025709318) q[3];
sx q[3];
rz(-0.51179143) q[3];
sx q[3];
rz(-0.6230841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8308782) q[2];
sx q[2];
rz(-1.0881311) q[2];
sx q[2];
rz(-2.3332398) q[2];
rz(2.3338142) q[3];
sx q[3];
rz(-0.72179663) q[3];
sx q[3];
rz(2.9479153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5508674) q[0];
sx q[0];
rz(-2.4243675) q[0];
sx q[0];
rz(0.89548683) q[0];
rz(0.26516178) q[1];
sx q[1];
rz(-0.83509713) q[1];
sx q[1];
rz(-0.53363824) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14347178) q[0];
sx q[0];
rz(-2.7406685) q[0];
sx q[0];
rz(-2.3401767) q[0];
rz(-1.3445504) q[2];
sx q[2];
rz(-0.81586736) q[2];
sx q[2];
rz(-2.0362542) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4662019) q[1];
sx q[1];
rz(-1.7168105) q[1];
sx q[1];
rz(2.9195665) q[1];
rz(1.8495464) q[3];
sx q[3];
rz(-2.4933443) q[3];
sx q[3];
rz(1.3384782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2395893) q[2];
sx q[2];
rz(-1.5782372) q[2];
sx q[2];
rz(2.5732102) q[2];
rz(-1.9249632) q[3];
sx q[3];
rz(-2.670848) q[3];
sx q[3];
rz(-0.34354982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50826532) q[0];
sx q[0];
rz(-0.87843043) q[0];
sx q[0];
rz(-1.9653962) q[0];
rz(1.5559224) q[1];
sx q[1];
rz(-0.95247477) q[1];
sx q[1];
rz(1.0046545) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5414808) q[0];
sx q[0];
rz(-1.2280557) q[0];
sx q[0];
rz(0.39370763) q[0];
rz(-pi) q[1];
rz(2.9506358) q[2];
sx q[2];
rz(-0.40340427) q[2];
sx q[2];
rz(2.0463338) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5914454) q[1];
sx q[1];
rz(-1.0284371) q[1];
sx q[1];
rz(-1.8591327) q[1];
rz(-pi) q[2];
rz(-1.5435013) q[3];
sx q[3];
rz(-2.7276162) q[3];
sx q[3];
rz(2.9881145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7375609) q[2];
sx q[2];
rz(-2.1514251) q[2];
sx q[2];
rz(-0.44815865) q[2];
rz(2.8435977) q[3];
sx q[3];
rz(-0.69152504) q[3];
sx q[3];
rz(0.34860778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3649243) q[0];
sx q[0];
rz(-1.3013327) q[0];
sx q[0];
rz(2.4834494) q[0];
rz(1.2843885) q[1];
sx q[1];
rz(-2.6874266) q[1];
sx q[1];
rz(0.67214322) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2015398) q[0];
sx q[0];
rz(-1.7367898) q[0];
sx q[0];
rz(-3.1340909) q[0];
rz(-2.7156098) q[2];
sx q[2];
rz(-1.1642712) q[2];
sx q[2];
rz(-2.4533518) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.044683177) q[1];
sx q[1];
rz(-1.9367366) q[1];
sx q[1];
rz(0.52176042) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7147195) q[3];
sx q[3];
rz(-2.2669499) q[3];
sx q[3];
rz(-1.5211943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.60823524) q[2];
sx q[2];
rz(-0.63295263) q[2];
sx q[2];
rz(2.7427924) q[2];
rz(-2.3170025) q[3];
sx q[3];
rz(-1.7493068) q[3];
sx q[3];
rz(-2.5430172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44052112) q[0];
sx q[0];
rz(-0.43061391) q[0];
sx q[0];
rz(-0.15821247) q[0];
rz(-2.054706) q[1];
sx q[1];
rz(-2.465076) q[1];
sx q[1];
rz(2.3349082) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5935532) q[0];
sx q[0];
rz(-2.002945) q[0];
sx q[0];
rz(3.0120864) q[0];
rz(-2.9358747) q[2];
sx q[2];
rz(-1.6763902) q[2];
sx q[2];
rz(-0.92668698) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.403277) q[1];
sx q[1];
rz(-0.74790955) q[1];
sx q[1];
rz(0.43055375) q[1];
rz(1.8716378) q[3];
sx q[3];
rz(-1.6013718) q[3];
sx q[3];
rz(-1.2310864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.57758254) q[2];
sx q[2];
rz(-1.492604) q[2];
sx q[2];
rz(0.99964833) q[2];
rz(-3.0380761) q[3];
sx q[3];
rz(-3.0045356) q[3];
sx q[3];
rz(-2.0172393) q[3];
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
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0209811) q[0];
sx q[0];
rz(-0.7779026) q[0];
sx q[0];
rz(-0.65752423) q[0];
rz(-0.28655562) q[1];
sx q[1];
rz(-0.92266881) q[1];
sx q[1];
rz(-2.8009169) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0176312) q[0];
sx q[0];
rz(-2.258856) q[0];
sx q[0];
rz(0.4171564) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9129487) q[2];
sx q[2];
rz(-0.59235448) q[2];
sx q[2];
rz(0.60281384) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.77010158) q[1];
sx q[1];
rz(-2.1226468) q[1];
sx q[1];
rz(-2.9279207) q[1];
rz(-pi) q[2];
rz(0.19823234) q[3];
sx q[3];
rz(-1.7361589) q[3];
sx q[3];
rz(0.73393047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3866773) q[2];
sx q[2];
rz(-1.200054) q[2];
sx q[2];
rz(2.6436451) q[2];
rz(0.30161101) q[3];
sx q[3];
rz(-2.7409654) q[3];
sx q[3];
rz(0.51914674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4122445) q[0];
sx q[0];
rz(-0.059818581) q[0];
sx q[0];
rz(-2.8630032) q[0];
rz(-0.57922286) q[1];
sx q[1];
rz(-0.93943739) q[1];
sx q[1];
rz(-3.0648807) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3007606) q[0];
sx q[0];
rz(-1.4474086) q[0];
sx q[0];
rz(1.8873909) q[0];
rz(-pi) q[1];
rz(-1.8329343) q[2];
sx q[2];
rz(-2.5524676) q[2];
sx q[2];
rz(0.72074705) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0662688) q[1];
sx q[1];
rz(-1.2053718) q[1];
sx q[1];
rz(-1.5931904) q[1];
rz(-1.4823227) q[3];
sx q[3];
rz(-1.7469329) q[3];
sx q[3];
rz(-0.74151553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3672093) q[2];
sx q[2];
rz(-0.70720208) q[2];
sx q[2];
rz(0.45551604) q[2];
rz(-2.7101743) q[3];
sx q[3];
rz(-0.12532561) q[3];
sx q[3];
rz(-0.063522696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.9615622) q[0];
sx q[0];
rz(-1.7407692) q[0];
sx q[0];
rz(-2.2041007) q[0];
rz(-0.58615276) q[1];
sx q[1];
rz(-1.6393607) q[1];
sx q[1];
rz(1.6929109) q[1];
rz(1.0473245) q[2];
sx q[2];
rz(-0.53146711) q[2];
sx q[2];
rz(-0.55234595) q[2];
rz(0.42701957) q[3];
sx q[3];
rz(-1.9058766) q[3];
sx q[3];
rz(1.8591892) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
