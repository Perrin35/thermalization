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
rz(0.02286214) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45863736) q[0];
sx q[0];
rz(-1.6874755) q[0];
sx q[0];
rz(1.4215683) q[0];
rz(-pi) q[1];
rz(-0.51486751) q[2];
sx q[2];
rz(-1.1866743) q[2];
sx q[2];
rz(1.1615953) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.248977) q[1];
sx q[1];
rz(-1.2274582) q[1];
sx q[1];
rz(2.9825319) q[1];
x q[2];
rz(-2.1894737) q[3];
sx q[3];
rz(-1.3631571) q[3];
sx q[3];
rz(1.729897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7227398) q[2];
sx q[2];
rz(-0.8272233) q[2];
sx q[2];
rz(2.375405) q[2];
rz(0.17151672) q[3];
sx q[3];
rz(-2.4092716) q[3];
sx q[3];
rz(-2.5550301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8121174) q[0];
sx q[0];
rz(-0.43918878) q[0];
sx q[0];
rz(-3.0549333) q[0];
rz(2.2791729) q[1];
sx q[1];
rz(-2.5984867) q[1];
sx q[1];
rz(3.0564953) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8893196) q[0];
sx q[0];
rz(-0.62250455) q[0];
sx q[0];
rz(0.028153367) q[0];
rz(-pi) q[1];
rz(-1.1413967) q[2];
sx q[2];
rz(-2.4944802) q[2];
sx q[2];
rz(2.5246758) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7460691) q[1];
sx q[1];
rz(-1.4847401) q[1];
sx q[1];
rz(-1.6834016) q[1];
x q[2];
rz(1.1005136) q[3];
sx q[3];
rz(-1.2840052) q[3];
sx q[3];
rz(-1.6581397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8791085) q[2];
sx q[2];
rz(-1.5510677) q[2];
sx q[2];
rz(-1.6764199) q[2];
rz(1.9654467) q[3];
sx q[3];
rz(-2.2361103) q[3];
sx q[3];
rz(-0.24648497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-0.48433205) q[0];
sx q[0];
rz(-3.0069139) q[0];
sx q[0];
rz(2.235967) q[0];
rz(2.8088645) q[1];
sx q[1];
rz(-0.85488027) q[1];
sx q[1];
rz(1.3844301) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4577643) q[0];
sx q[0];
rz(-2.0450852) q[0];
sx q[0];
rz(1.6636687) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.076782) q[2];
sx q[2];
rz(-1.9188606) q[2];
sx q[2];
rz(-0.51007523) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3123734) q[1];
sx q[1];
rz(-2.3260498) q[1];
sx q[1];
rz(1.6014598) q[1];
rz(-0.53220766) q[3];
sx q[3];
rz(-2.0720707) q[3];
sx q[3];
rz(-1.3794848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8973792) q[2];
sx q[2];
rz(-0.34003568) q[2];
sx q[2];
rz(-1.8133694) q[2];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.06552799) q[0];
sx q[0];
rz(-0.17834839) q[0];
sx q[0];
rz(0.67101014) q[0];
rz(1.7975851) q[1];
sx q[1];
rz(-1.9799415) q[1];
sx q[1];
rz(-0.20733325) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79065381) q[0];
sx q[0];
rz(-1.6121943) q[0];
sx q[0];
rz(-3.1326276) q[0];
rz(-pi) q[1];
rz(3.1111654) q[2];
sx q[2];
rz(-1.7446339) q[2];
sx q[2];
rz(2.0174842) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3934717) q[1];
sx q[1];
rz(-1.8078184) q[1];
sx q[1];
rz(-0.15740983) q[1];
rz(-pi) q[2];
rz(-1.585235) q[3];
sx q[3];
rz(-1.0591905) q[3];
sx q[3];
rz(-2.5479941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.31071445) q[2];
sx q[2];
rz(-1.0881311) q[2];
sx q[2];
rz(2.3332398) q[2];
rz(-2.3338142) q[3];
sx q[3];
rz(-0.72179663) q[3];
sx q[3];
rz(0.19367735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5508674) q[0];
sx q[0];
rz(-2.4243675) q[0];
sx q[0];
rz(-2.2461058) q[0];
rz(-0.26516178) q[1];
sx q[1];
rz(-0.83509713) q[1];
sx q[1];
rz(0.53363824) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4424608) q[0];
sx q[0];
rz(-1.8457544) q[0];
sx q[0];
rz(-1.2752227) q[0];
rz(-pi) q[1];
rz(2.3737714) q[2];
sx q[2];
rz(-1.7349093) q[2];
sx q[2];
rz(0.62190157) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.071760885) q[1];
sx q[1];
rz(-1.7904209) q[1];
sx q[1];
rz(-1.4211618) q[1];
rz(-pi) q[2];
x q[2];
rz(2.93612) q[3];
sx q[3];
rz(-2.1901202) q[3];
sx q[3];
rz(1.6831786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2395893) q[2];
sx q[2];
rz(-1.5633554) q[2];
sx q[2];
rz(-2.5732102) q[2];
rz(1.2166294) q[3];
sx q[3];
rz(-0.47074461) q[3];
sx q[3];
rz(0.34354982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-0.95247477) q[1];
sx q[1];
rz(1.0046545) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6509103) q[0];
sx q[0];
rz(-0.51603979) q[0];
sx q[0];
rz(2.3923621) q[0];
x q[1];
rz(0.3968233) q[2];
sx q[2];
rz(-1.6453711) q[2];
sx q[2];
rz(-2.4900988) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5914454) q[1];
sx q[1];
rz(-1.0284371) q[1];
sx q[1];
rz(1.8591327) q[1];
rz(-1.1569571) q[3];
sx q[3];
rz(-1.5817747) q[3];
sx q[3];
rz(1.4423086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4040318) q[2];
sx q[2];
rz(-2.1514251) q[2];
sx q[2];
rz(-0.44815865) q[2];
rz(2.8435977) q[3];
sx q[3];
rz(-0.69152504) q[3];
sx q[3];
rz(-2.7929849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3649243) q[0];
sx q[0];
rz(-1.8402599) q[0];
sx q[0];
rz(0.65814322) q[0];
rz(-1.8572042) q[1];
sx q[1];
rz(-2.6874266) q[1];
sx q[1];
rz(-2.4694494) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9400529) q[0];
sx q[0];
rz(-1.4048028) q[0];
sx q[0];
rz(-3.1340909) q[0];
rz(-1.1291814) q[2];
sx q[2];
rz(-1.1814983) q[2];
sx q[2];
rz(1.0600952) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0828515) q[1];
sx q[1];
rz(-0.62742678) q[1];
sx q[1];
rz(0.65545603) q[1];
rz(-pi) q[2];
x q[2];
rz(0.70127212) q[3];
sx q[3];
rz(-1.6810732) q[3];
sx q[3];
rz(2.9993204) q[3];
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
rz(2.3170025) q[3];
sx q[3];
rz(-1.3922858) q[3];
sx q[3];
rz(0.59857541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44052112) q[0];
sx q[0];
rz(-2.7109787) q[0];
sx q[0];
rz(0.15821247) q[0];
rz(1.0868866) q[1];
sx q[1];
rz(-0.67651665) q[1];
sx q[1];
rz(0.80668443) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0772484) q[0];
sx q[0];
rz(-1.4532538) q[0];
sx q[0];
rz(1.1354406) q[0];
rz(2.6629692) q[2];
sx q[2];
rz(-0.23089409) q[2];
sx q[2];
rz(-2.965197) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1572691) q[1];
sx q[1];
rz(-1.8586129) q[1];
sx q[1];
rz(-2.4411574) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.673647) q[3];
sx q[3];
rz(-2.8392483) q[3];
sx q[3];
rz(-0.24149382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5640101) q[2];
sx q[2];
rz(-1.492604) q[2];
sx q[2];
rz(-2.1419443) q[2];
rz(3.0380761) q[3];
sx q[3];
rz(-3.0045356) q[3];
sx q[3];
rz(2.0172393) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12061159) q[0];
sx q[0];
rz(-2.3636901) q[0];
sx q[0];
rz(-0.65752423) q[0];
rz(2.855037) q[1];
sx q[1];
rz(-0.92266881) q[1];
sx q[1];
rz(-2.8009169) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17250241) q[0];
sx q[0];
rz(-1.2524676) q[0];
sx q[0];
rz(-0.83842917) q[0];
rz(-1.9129487) q[2];
sx q[2];
rz(-0.59235448) q[2];
sx q[2];
rz(-2.5387788) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9790736) q[1];
sx q[1];
rz(-0.58774978) q[1];
sx q[1];
rz(-1.9025004) q[1];
x q[2];
rz(2.9433603) q[3];
sx q[3];
rz(-1.4054338) q[3];
sx q[3];
rz(-2.4076622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3866773) q[2];
sx q[2];
rz(-1.9415386) q[2];
sx q[2];
rz(-0.49794751) q[2];
rz(0.30161101) q[3];
sx q[3];
rz(-2.7409654) q[3];
sx q[3];
rz(0.51914674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72934812) q[0];
sx q[0];
rz(-3.0817741) q[0];
sx q[0];
rz(0.27858946) q[0];
rz(-0.57922286) q[1];
sx q[1];
rz(-2.2021553) q[1];
sx q[1];
rz(3.0648807) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3007606) q[0];
sx q[0];
rz(-1.6941841) q[0];
sx q[0];
rz(-1.2542017) q[0];
x q[1];
rz(-1.3086583) q[2];
sx q[2];
rz(-2.5524676) q[2];
sx q[2];
rz(-0.72074705) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1288651) q[1];
sx q[1];
rz(-2.7755133) q[1];
sx q[1];
rz(-3.0831343) q[1];
rz(-pi) q[2];
rz(0.46080188) q[3];
sx q[3];
rz(-2.9446903) q[3];
sx q[3];
rz(-2.868696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3672093) q[2];
sx q[2];
rz(-2.4343906) q[2];
sx q[2];
rz(-0.45551604) q[2];
rz(-0.43141836) q[3];
sx q[3];
rz(-0.12532561) q[3];
sx q[3];
rz(0.063522696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1800304) q[0];
sx q[0];
rz(-1.4008235) q[0];
sx q[0];
rz(0.93749198) q[0];
rz(-2.5554399) q[1];
sx q[1];
rz(-1.502232) q[1];
sx q[1];
rz(-1.4486817) q[1];
rz(1.0473245) q[2];
sx q[2];
rz(-0.53146711) q[2];
sx q[2];
rz(-0.55234595) q[2];
rz(-2.4424845) q[3];
sx q[3];
rz(-2.6052478) q[3];
sx q[3];
rz(-0.33741046) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
