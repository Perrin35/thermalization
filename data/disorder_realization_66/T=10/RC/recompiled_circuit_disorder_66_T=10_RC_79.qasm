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
rz(0.45863736) q[0];
sx q[0];
rz(-1.4541172) q[0];
sx q[0];
rz(1.7200243) q[0];
rz(-pi) q[1];
rz(0.51486751) q[2];
sx q[2];
rz(-1.1866743) q[2];
sx q[2];
rz(1.9799973) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3373128) q[1];
sx q[1];
rz(-2.7645281) q[1];
sx q[1];
rz(1.9878597) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1894737) q[3];
sx q[3];
rz(-1.7784356) q[3];
sx q[3];
rz(-1.4116956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.4188529) q[2];
sx q[2];
rz(-2.3143694) q[2];
sx q[2];
rz(2.375405) q[2];
rz(0.17151672) q[3];
sx q[3];
rz(-0.73232108) q[3];
sx q[3];
rz(2.5550301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3294753) q[0];
sx q[0];
rz(-0.43918878) q[0];
sx q[0];
rz(3.0549333) q[0];
rz(-2.2791729) q[1];
sx q[1];
rz(-0.54310596) q[1];
sx q[1];
rz(3.0564953) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9239685) q[0];
sx q[0];
rz(-2.1930165) q[0];
sx q[0];
rz(1.5505962) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.000196) q[2];
sx q[2];
rz(-0.64711249) q[2];
sx q[2];
rz(2.5246758) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.18499204) q[1];
sx q[1];
rz(-1.4586095) q[1];
sx q[1];
rz(3.0549906) q[1];
rz(-pi) q[2];
rz(-0.99382932) q[3];
sx q[3];
rz(-0.54518632) q[3];
sx q[3];
rz(0.59515566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8791085) q[2];
sx q[2];
rz(-1.5510677) q[2];
sx q[2];
rz(-1.6764199) q[2];
rz(1.1761459) q[3];
sx q[3];
rz(-2.2361103) q[3];
sx q[3];
rz(0.24648497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(2.6572606) q[0];
sx q[0];
rz(-0.13467877) q[0];
sx q[0];
rz(0.9056257) q[0];
rz(-2.8088645) q[1];
sx q[1];
rz(-0.85488027) q[1];
sx q[1];
rz(-1.3844301) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15554409) q[0];
sx q[0];
rz(-1.4882003) q[0];
sx q[0];
rz(0.47604576) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7483814) q[2];
sx q[2];
rz(-1.097743) q[2];
sx q[2];
rz(2.2676603) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3570909) q[1];
sx q[1];
rz(-0.75575268) q[1];
sx q[1];
rz(0.032553629) q[1];
rz(0.53220766) q[3];
sx q[3];
rz(-1.0695219) q[3];
sx q[3];
rz(1.7621079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.24421346) q[2];
sx q[2];
rz(-0.34003568) q[2];
sx q[2];
rz(-1.3282233) q[2];
rz(1.3978847) q[3];
sx q[3];
rz(-2.6013247) q[3];
sx q[3];
rz(0.11169294) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.06552799) q[0];
sx q[0];
rz(-2.9632443) q[0];
sx q[0];
rz(0.67101014) q[0];
rz(-1.7975851) q[1];
sx q[1];
rz(-1.1616511) q[1];
sx q[1];
rz(-0.20733325) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78051356) q[0];
sx q[0];
rz(-1.5618389) q[0];
sx q[0];
rz(-1.5293967) q[0];
rz(-pi) q[1];
rz(-1.7423332) q[2];
sx q[2];
rz(-0.17645391) q[2];
sx q[2];
rz(-1.2982969) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.281665) q[1];
sx q[1];
rz(-1.7237701) q[1];
sx q[1];
rz(1.8106736) q[1];
rz(-pi) q[2];
rz(-3.1158833) q[3];
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
rz(-0.80835289) q[2];
rz(0.80777848) q[3];
sx q[3];
rz(-0.72179663) q[3];
sx q[3];
rz(0.19367735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5508674) q[0];
sx q[0];
rz(-0.71722513) q[0];
sx q[0];
rz(-2.2461058) q[0];
rz(0.26516178) q[1];
sx q[1];
rz(-0.83509713) q[1];
sx q[1];
rz(-0.53363824) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4424608) q[0];
sx q[0];
rz(-1.8457544) q[0];
sx q[0];
rz(1.2752227) q[0];
rz(-pi) q[1];
rz(-2.3737714) q[2];
sx q[2];
rz(-1.4066833) q[2];
sx q[2];
rz(-2.5196911) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6753908) q[1];
sx q[1];
rz(-1.4247822) q[1];
sx q[1];
rz(2.9195665) q[1];
rz(-pi) q[2];
rz(1.2920462) q[3];
sx q[3];
rz(-2.4933443) q[3];
sx q[3];
rz(-1.3384782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2395893) q[2];
sx q[2];
rz(-1.5633554) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-2.6333273) q[0];
sx q[0];
rz(-2.2631622) q[0];
sx q[0];
rz(1.1761965) q[0];
rz(1.5856702) q[1];
sx q[1];
rz(-0.95247477) q[1];
sx q[1];
rz(2.1369381) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6001119) q[0];
sx q[0];
rz(-1.2280557) q[0];
sx q[0];
rz(2.747885) q[0];
rz(-pi) q[1];
rz(-0.3968233) q[2];
sx q[2];
rz(-1.6453711) q[2];
sx q[2];
rz(2.4900988) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.028588258) q[1];
sx q[1];
rz(-0.60739809) q[1];
sx q[1];
rz(0.44087704) q[1];
rz(-pi) q[2];
rz(1.9846356) q[3];
sx q[3];
rz(-1.5817747) q[3];
sx q[3];
rz(-1.699284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4040318) q[2];
sx q[2];
rz(-2.1514251) q[2];
sx q[2];
rz(-2.693434) q[2];
rz(-0.29799497) q[3];
sx q[3];
rz(-2.4500676) q[3];
sx q[3];
rz(-0.34860778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3649243) q[0];
sx q[0];
rz(-1.3013327) q[0];
sx q[0];
rz(2.4834494) q[0];
rz(-1.8572042) q[1];
sx q[1];
rz(-2.6874266) q[1];
sx q[1];
rz(-2.4694494) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2015398) q[0];
sx q[0];
rz(-1.7367898) q[0];
sx q[0];
rz(-3.1340909) q[0];
x q[1];
rz(-0.80588801) q[2];
sx q[2];
rz(-0.58008367) q[2];
sx q[2];
rz(-2.9758331) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8183648) q[1];
sx q[1];
rz(-1.0867026) q[1];
sx q[1];
rz(1.9869884) q[1];
rz(-1.4268731) q[3];
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
rz(-2.5333574) q[2];
sx q[2];
rz(-2.50864) q[2];
sx q[2];
rz(-2.7427924) q[2];
rz(-2.3170025) q[3];
sx q[3];
rz(-1.7493068) q[3];
sx q[3];
rz(-2.5430172) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7010715) q[0];
sx q[0];
rz(-2.7109787) q[0];
sx q[0];
rz(0.15821247) q[0];
rz(2.054706) q[1];
sx q[1];
rz(-2.465076) q[1];
sx q[1];
rz(0.80668443) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5935532) q[0];
sx q[0];
rz(-2.002945) q[0];
sx q[0];
rz(-0.12950626) q[0];
rz(-pi) q[1];
rz(0.20571795) q[2];
sx q[2];
rz(-1.4652025) q[2];
sx q[2];
rz(-2.2149057) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.403277) q[1];
sx q[1];
rz(-2.3936831) q[1];
sx q[1];
rz(2.7110389) q[1];
x q[2];
rz(0.032012352) q[3];
sx q[3];
rz(-1.2700998) q[3];
sx q[3];
rz(-2.7923982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5640101) q[2];
sx q[2];
rz(-1.6489886) q[2];
sx q[2];
rz(-0.99964833) q[2];
rz(0.10351652) q[3];
sx q[3];
rz(-0.13705702) q[3];
sx q[3];
rz(2.0172393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0209811) q[0];
sx q[0];
rz(-2.3636901) q[0];
sx q[0];
rz(2.4840684) q[0];
rz(2.855037) q[1];
sx q[1];
rz(-0.92266881) q[1];
sx q[1];
rz(0.34067571) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17250241) q[0];
sx q[0];
rz(-1.889125) q[0];
sx q[0];
rz(0.83842917) q[0];
rz(1.0057783) q[2];
sx q[2];
rz(-1.7592332) q[2];
sx q[2];
rz(-2.4609158) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.454168) q[1];
sx q[1];
rz(-1.7523645) q[1];
sx q[1];
rz(1.0086235) q[1];
rz(-pi) q[2];
rz(-0.19823234) q[3];
sx q[3];
rz(-1.4054338) q[3];
sx q[3];
rz(0.73393047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.75491536) q[2];
sx q[2];
rz(-1.200054) q[2];
sx q[2];
rz(0.49794751) q[2];
rz(-2.8399816) q[3];
sx q[3];
rz(-2.7409654) q[3];
sx q[3];
rz(-2.6224459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4122445) q[0];
sx q[0];
rz(-3.0817741) q[0];
sx q[0];
rz(2.8630032) q[0];
rz(-2.5623698) q[1];
sx q[1];
rz(-0.93943739) q[1];
sx q[1];
rz(3.0648807) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68966507) q[0];
sx q[0];
rz(-1.2566914) q[0];
sx q[0];
rz(-0.12977022) q[0];
x q[1];
rz(0.17148359) q[2];
sx q[2];
rz(-2.137261) q[2];
sx q[2];
rz(-2.7330074) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0662688) q[1];
sx q[1];
rz(-1.2053718) q[1];
sx q[1];
rz(-1.5931904) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6807908) q[3];
sx q[3];
rz(-2.9446903) q[3];
sx q[3];
rz(0.27289665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3672093) q[2];
sx q[2];
rz(-2.4343906) q[2];
sx q[2];
rz(-0.45551604) q[2];
rz(-0.43141836) q[3];
sx q[3];
rz(-3.016267) q[3];
sx q[3];
rz(-0.063522696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
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
rz(-0.2858328) q[2];
sx q[2];
rz(-2.0252068) q[2];
sx q[2];
rz(0.037638738) q[2];
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
