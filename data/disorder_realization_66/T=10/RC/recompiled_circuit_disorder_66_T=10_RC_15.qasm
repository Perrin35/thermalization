OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3333617) q[0];
sx q[0];
rz(-0.97167492) q[0];
sx q[0];
rz(-1.4810286) q[0];
rz(-0.12227585) q[1];
sx q[1];
rz(-0.08637698) q[1];
sx q[1];
rz(-3.1187305) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4532783) q[0];
sx q[0];
rz(-0.18916057) q[0];
sx q[0];
rz(-2.2384089) q[0];
rz(-pi) q[1];
rz(-0.51486751) q[2];
sx q[2];
rz(-1.9549184) q[2];
sx q[2];
rz(1.9799973) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3373128) q[1];
sx q[1];
rz(-2.7645281) q[1];
sx q[1];
rz(1.9878597) q[1];
rz(2.888527) q[3];
sx q[3];
rz(-0.96732891) q[3];
sx q[3];
rz(2.8367708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.4188529) q[2];
sx q[2];
rz(-2.3143694) q[2];
sx q[2];
rz(-0.76618761) q[2];
rz(-0.17151672) q[3];
sx q[3];
rz(-0.73232108) q[3];
sx q[3];
rz(-2.5550301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8121174) q[0];
sx q[0];
rz(-0.43918878) q[0];
sx q[0];
rz(3.0549333) q[0];
rz(-0.86241972) q[1];
sx q[1];
rz(-2.5984867) q[1];
sx q[1];
rz(-0.08509732) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25227308) q[0];
sx q[0];
rz(-2.5190881) q[0];
sx q[0];
rz(-3.1134393) q[0];
rz(-0.9688103) q[2];
sx q[2];
rz(-1.3170871) q[2];
sx q[2];
rz(1.837455) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9566006) q[1];
sx q[1];
rz(-1.4586095) q[1];
sx q[1];
rz(3.0549906) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8220903) q[3];
sx q[3];
rz(-2.0204244) q[3];
sx q[3];
rz(-0.055469661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8791085) q[2];
sx q[2];
rz(-1.5510677) q[2];
sx q[2];
rz(-1.6764199) q[2];
rz(-1.9654467) q[3];
sx q[3];
rz(-0.90548235) q[3];
sx q[3];
rz(-0.24648497) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6572606) q[0];
sx q[0];
rz(-3.0069139) q[0];
sx q[0];
rz(-0.9056257) q[0];
rz(0.33272818) q[1];
sx q[1];
rz(-2.2867124) q[1];
sx q[1];
rz(1.3844301) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2565838) q[0];
sx q[0];
rz(-2.6589767) q[0];
sx q[0];
rz(-0.17871876) q[0];
x q[1];
rz(-0.92817523) q[2];
sx q[2];
rz(-2.536142) q[2];
sx q[2];
rz(1.5290934) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.76259957) q[1];
sx q[1];
rz(-1.5484719) q[1];
sx q[1];
rz(-0.75548817) q[1];
rz(-pi) q[2];
rz(2.1372041) q[3];
sx q[3];
rz(-2.0319788) q[3];
sx q[3];
rz(3.0571292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8973792) q[2];
sx q[2];
rz(-2.801557) q[2];
sx q[2];
rz(-1.8133694) q[2];
rz(1.3978847) q[3];
sx q[3];
rz(-2.6013247) q[3];
sx q[3];
rz(0.11169294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(3.0760647) q[0];
sx q[0];
rz(-2.9632443) q[0];
sx q[0];
rz(-0.67101014) q[0];
rz(1.7975851) q[1];
sx q[1];
rz(-1.9799415) q[1];
sx q[1];
rz(-0.20733325) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78051356) q[0];
sx q[0];
rz(-1.5618389) q[0];
sx q[0];
rz(1.5293967) q[0];
x q[1];
rz(-1.3968799) q[2];
sx q[2];
rz(-1.5408278) q[2];
sx q[2];
rz(2.7001691) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9878848) q[1];
sx q[1];
rz(-2.8578836) q[1];
sx q[1];
rz(2.1464286) q[1];
rz(1.585235) q[3];
sx q[3];
rz(-1.0591905) q[3];
sx q[3];
rz(2.5479941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8308782) q[2];
sx q[2];
rz(-1.0881311) q[2];
sx q[2];
rz(2.3332398) q[2];
rz(-0.80777848) q[3];
sx q[3];
rz(-0.72179663) q[3];
sx q[3];
rz(-0.19367735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59072524) q[0];
sx q[0];
rz(-2.4243675) q[0];
sx q[0];
rz(-2.2461058) q[0];
rz(-0.26516178) q[1];
sx q[1];
rz(-2.3064955) q[1];
sx q[1];
rz(2.6079544) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69913188) q[0];
sx q[0];
rz(-1.8457544) q[0];
sx q[0];
rz(-1.86637) q[0];
rz(-pi) q[1];
rz(-0.23405481) q[2];
sx q[2];
rz(-2.3599527) q[2];
sx q[2];
rz(2.3603338) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4662019) q[1];
sx q[1];
rz(-1.4247822) q[1];
sx q[1];
rz(-0.22202613) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2920462) q[3];
sx q[3];
rz(-2.4933443) q[3];
sx q[3];
rz(-1.8031144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2395893) q[2];
sx q[2];
rz(-1.5782372) q[2];
sx q[2];
rz(2.5732102) q[2];
rz(-1.9249632) q[3];
sx q[3];
rz(-0.47074461) q[3];
sx q[3];
rz(0.34354982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50826532) q[0];
sx q[0];
rz(-2.2631622) q[0];
sx q[0];
rz(1.1761965) q[0];
rz(-1.5856702) q[1];
sx q[1];
rz(-2.1891179) q[1];
sx q[1];
rz(2.1369381) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6509103) q[0];
sx q[0];
rz(-0.51603979) q[0];
sx q[0];
rz(-0.74923058) q[0];
rz(-pi) q[1];
rz(-1.4899646) q[2];
sx q[2];
rz(-1.1751375) q[2];
sx q[2];
rz(-0.8880907) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.028588258) q[1];
sx q[1];
rz(-0.60739809) q[1];
sx q[1];
rz(2.7007156) q[1];
rz(-1.9846356) q[3];
sx q[3];
rz(-1.5598179) q[3];
sx q[3];
rz(-1.699284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4040318) q[2];
sx q[2];
rz(-0.99016756) q[2];
sx q[2];
rz(-2.693434) q[2];
rz(-0.29799497) q[3];
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
rz(pi/2) q[3];
sx q[3];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1561688) q[0];
sx q[0];
rz(-0.16616136) q[0];
sx q[0];
rz(-1.6155433) q[0];
x q[1];
rz(-0.42598287) q[2];
sx q[2];
rz(-1.9773215) q[2];
sx q[2];
rz(0.68824088) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0587412) q[1];
sx q[1];
rz(-2.5141659) q[1];
sx q[1];
rz(2.4861366) q[1];
rz(-pi) q[2];
rz(-2.4403205) q[3];
sx q[3];
rz(-1.4605195) q[3];
sx q[3];
rz(-2.9993204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5333574) q[2];
sx q[2];
rz(-0.63295263) q[2];
sx q[2];
rz(0.39880025) q[2];
rz(-2.3170025) q[3];
sx q[3];
rz(-1.7493068) q[3];
sx q[3];
rz(-2.5430172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44052112) q[0];
sx q[0];
rz(-0.43061391) q[0];
sx q[0];
rz(0.15821247) q[0];
rz(-2.054706) q[1];
sx q[1];
rz(-0.67651665) q[1];
sx q[1];
rz(-2.3349082) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54803941) q[0];
sx q[0];
rz(-2.002945) q[0];
sx q[0];
rz(-0.12950626) q[0];
rz(-2.9358747) q[2];
sx q[2];
rz(-1.4652025) q[2];
sx q[2];
rz(0.92668698) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.17864171) q[1];
sx q[1];
rz(-2.2370403) q[1];
sx q[1];
rz(-1.2013749) q[1];
rz(-pi) q[2];
rz(0.032012352) q[3];
sx q[3];
rz(-1.2700998) q[3];
sx q[3];
rz(-2.7923982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5640101) q[2];
sx q[2];
rz(-1.492604) q[2];
sx q[2];
rz(2.1419443) q[2];
rz(-3.0380761) q[3];
sx q[3];
rz(-3.0045356) q[3];
sx q[3];
rz(-2.0172393) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12061159) q[0];
sx q[0];
rz(-2.3636901) q[0];
sx q[0];
rz(-0.65752423) q[0];
rz(-2.855037) q[1];
sx q[1];
rz(-2.2189238) q[1];
sx q[1];
rz(-2.8009169) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17250241) q[0];
sx q[0];
rz(-1.2524676) q[0];
sx q[0];
rz(2.3031635) q[0];
rz(-2.9195243) q[2];
sx q[2];
rz(-2.1246315) q[2];
sx q[2];
rz(-1.008322) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.68742467) q[1];
sx q[1];
rz(-1.7523645) q[1];
sx q[1];
rz(1.0086235) q[1];
rz(-pi) q[2];
rz(1.4021923) q[3];
sx q[3];
rz(-1.3753034) q[3];
sx q[3];
rz(-2.33778) q[3];
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
rz(-0.30161101) q[3];
sx q[3];
rz(-0.40062723) q[3];
sx q[3];
rz(-2.6224459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4122445) q[0];
sx q[0];
rz(-3.0817741) q[0];
sx q[0];
rz(-2.8630032) q[0];
rz(0.57922286) q[1];
sx q[1];
rz(-0.93943739) q[1];
sx q[1];
rz(3.0648807) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0893223) q[0];
sx q[0];
rz(-2.8025586) q[0];
sx q[0];
rz(1.9498755) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3086583) q[2];
sx q[2];
rz(-0.58912504) q[2];
sx q[2];
rz(0.72074705) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0753239) q[1];
sx q[1];
rz(-1.9362209) q[1];
sx q[1];
rz(-1.5484023) q[1];
rz(-0.17681392) q[3];
sx q[3];
rz(-1.483695) q[3];
sx q[3];
rz(-0.84482312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3672093) q[2];
sx q[2];
rz(-0.70720208) q[2];
sx q[2];
rz(2.6860766) q[2];
rz(2.7101743) q[3];
sx q[3];
rz(-0.12532561) q[3];
sx q[3];
rz(0.063522696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.099844) q[2];
sx q[2];
rz(-1.3146613) q[2];
sx q[2];
rz(1.4801499) q[2];
rz(1.9361817) q[3];
sx q[3];
rz(-1.1689417) q[3];
sx q[3];
rz(0.4369215) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
