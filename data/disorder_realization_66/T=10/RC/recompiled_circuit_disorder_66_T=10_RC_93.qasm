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
rz(-0.4532783) q[0];
sx q[0];
rz(-2.9524321) q[0];
sx q[0];
rz(-2.2384089) q[0];
rz(-pi) q[1];
rz(-2.6267251) q[2];
sx q[2];
rz(-1.9549184) q[2];
sx q[2];
rz(-1.9799973) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.37576807) q[1];
sx q[1];
rz(-1.4210912) q[1];
sx q[1];
rz(1.9181812) q[1];
x q[2];
rz(1.919235) q[3];
sx q[3];
rz(-0.64823965) q[3];
sx q[3];
rz(-3.0188308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.4188529) q[2];
sx q[2];
rz(-0.8272233) q[2];
sx q[2];
rz(-2.375405) q[2];
rz(-0.17151672) q[3];
sx q[3];
rz(-0.73232108) q[3];
sx q[3];
rz(-2.5550301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8121174) q[0];
sx q[0];
rz(-0.43918878) q[0];
sx q[0];
rz(0.086659327) q[0];
rz(-2.2791729) q[1];
sx q[1];
rz(-0.54310596) q[1];
sx q[1];
rz(3.0564953) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.800195) q[0];
sx q[0];
rz(-1.5872103) q[0];
sx q[0];
rz(-0.62231681) q[0];
rz(-pi) q[1];
rz(1.1413967) q[2];
sx q[2];
rz(-2.4944802) q[2];
sx q[2];
rz(0.61691689) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.47463402) q[1];
sx q[1];
rz(-2.9999795) q[1];
sx q[1];
rz(0.9160362) q[1];
rz(-1.1005136) q[3];
sx q[3];
rz(-1.8575875) q[3];
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
rz(-0.90548235) q[3];
sx q[3];
rz(-2.8951077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6572606) q[0];
sx q[0];
rz(-3.0069139) q[0];
sx q[0];
rz(2.235967) q[0];
rz(2.8088645) q[1];
sx q[1];
rz(-2.2867124) q[1];
sx q[1];
rz(1.7571626) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8850088) q[0];
sx q[0];
rz(-0.48261595) q[0];
sx q[0];
rz(-0.17871876) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2134174) q[2];
sx q[2];
rz(-2.536142) q[2];
sx q[2];
rz(-1.6124992) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3123734) q[1];
sx q[1];
rz(-0.81554283) q[1];
sx q[1];
rz(-1.6014598) q[1];
rz(-pi) q[2];
rz(0.53220766) q[3];
sx q[3];
rz(-2.0720707) q[3];
sx q[3];
rz(1.3794848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8973792) q[2];
sx q[2];
rz(-2.801557) q[2];
sx q[2];
rz(1.8133694) q[2];
rz(-1.3978847) q[3];
sx q[3];
rz(-0.54026794) q[3];
sx q[3];
rz(-3.0298997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
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
rz(-1.9799415) q[1];
sx q[1];
rz(0.20733325) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79065381) q[0];
sx q[0];
rz(-1.6121943) q[0];
sx q[0];
rz(3.1326276) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.030427293) q[2];
sx q[2];
rz(-1.3969587) q[2];
sx q[2];
rz(1.1241084) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.15370788) q[1];
sx q[1];
rz(-0.28370902) q[1];
sx q[1];
rz(-0.99516408) q[1];
x q[2];
rz(0.025709318) q[3];
sx q[3];
rz(-0.51179143) q[3];
sx q[3];
rz(-2.5185086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.31071445) q[2];
sx q[2];
rz(-1.0881311) q[2];
sx q[2];
rz(-0.80835289) q[2];
rz(-2.3338142) q[3];
sx q[3];
rz(-2.419796) q[3];
sx q[3];
rz(2.9479153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5508674) q[0];
sx q[0];
rz(-2.4243675) q[0];
sx q[0];
rz(-2.2461058) q[0];
rz(2.8764309) q[1];
sx q[1];
rz(-2.3064955) q[1];
sx q[1];
rz(2.6079544) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69913188) q[0];
sx q[0];
rz(-1.8457544) q[0];
sx q[0];
rz(1.2752227) q[0];
rz(-pi) q[1];
rz(1.7970423) q[2];
sx q[2];
rz(-2.3257253) q[2];
sx q[2];
rz(-1.1053384) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4662019) q[1];
sx q[1];
rz(-1.4247822) q[1];
sx q[1];
rz(-0.22202613) q[1];
rz(-0.20547262) q[3];
sx q[3];
rz(-0.95147248) q[3];
sx q[3];
rz(1.458414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2395893) q[2];
sx q[2];
rz(-1.5633554) q[2];
sx q[2];
rz(-2.5732102) q[2];
rz(-1.2166294) q[3];
sx q[3];
rz(-0.47074461) q[3];
sx q[3];
rz(-0.34354982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6333273) q[0];
sx q[0];
rz(-0.87843043) q[0];
sx q[0];
rz(-1.9653962) q[0];
rz(1.5559224) q[1];
sx q[1];
rz(-2.1891179) q[1];
sx q[1];
rz(-1.0046545) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83197901) q[0];
sx q[0];
rz(-1.9404611) q[0];
sx q[0];
rz(-1.9395104) q[0];
rz(-pi) q[1];
x q[1];
rz(0.3968233) q[2];
sx q[2];
rz(-1.4962215) q[2];
sx q[2];
rz(2.4900988) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.55014729) q[1];
sx q[1];
rz(-1.0284371) q[1];
sx q[1];
rz(-1.28246) q[1];
x q[2];
rz(-1.1569571) q[3];
sx q[3];
rz(-1.5817747) q[3];
sx q[3];
rz(1.4423086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4040318) q[2];
sx q[2];
rz(-0.99016756) q[2];
sx q[2];
rz(-2.693434) q[2];
rz(2.8435977) q[3];
sx q[3];
rz(-2.4500676) q[3];
sx q[3];
rz(-0.34860778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77666831) q[0];
sx q[0];
rz(-1.8402599) q[0];
sx q[0];
rz(-2.4834494) q[0];
rz(1.2843885) q[1];
sx q[1];
rz(-0.45416608) q[1];
sx q[1];
rz(2.4694494) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9854239) q[0];
sx q[0];
rz(-2.9754313) q[0];
sx q[0];
rz(1.5260494) q[0];
rz(-2.3357046) q[2];
sx q[2];
rz(-0.58008367) q[2];
sx q[2];
rz(-0.16575955) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3232279) q[1];
sx q[1];
rz(-1.0867026) q[1];
sx q[1];
rz(-1.9869884) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.16996202) q[3];
sx q[3];
rz(-2.4331577) q[3];
sx q[3];
rz(1.8426614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5333574) q[2];
sx q[2];
rz(-0.63295263) q[2];
sx q[2];
rz(-0.39880025) q[2];
rz(-2.3170025) q[3];
sx q[3];
rz(-1.3922858) q[3];
sx q[3];
rz(-0.59857541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.44052112) q[0];
sx q[0];
rz(-2.7109787) q[0];
sx q[0];
rz(-0.15821247) q[0];
rz(-2.054706) q[1];
sx q[1];
rz(-0.67651665) q[1];
sx q[1];
rz(-2.3349082) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54803941) q[0];
sx q[0];
rz(-1.1386477) q[0];
sx q[0];
rz(-0.12950626) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6786472) q[2];
sx q[2];
rz(-1.3662405) q[2];
sx q[2];
rz(0.66609913) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.17864171) q[1];
sx q[1];
rz(-0.9045524) q[1];
sx q[1];
rz(-1.9402177) q[1];
rz(-pi) q[2];
x q[2];
rz(0.032012352) q[3];
sx q[3];
rz(-1.2700998) q[3];
sx q[3];
rz(0.34919448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5640101) q[2];
sx q[2];
rz(-1.492604) q[2];
sx q[2];
rz(-0.99964833) q[2];
rz(0.10351652) q[3];
sx q[3];
rz(-3.0045356) q[3];
sx q[3];
rz(-2.0172393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
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
rz(-2.2189238) q[1];
sx q[1];
rz(-0.34067571) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1239615) q[0];
sx q[0];
rz(-2.258856) q[0];
sx q[0];
rz(-0.4171564) q[0];
rz(1.0057783) q[2];
sx q[2];
rz(-1.3823595) q[2];
sx q[2];
rz(2.4609158) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.454168) q[1];
sx q[1];
rz(-1.7523645) q[1];
sx q[1];
rz(1.0086235) q[1];
x q[2];
rz(1.7394003) q[3];
sx q[3];
rz(-1.7662893) q[3];
sx q[3];
rz(0.80381264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3866773) q[2];
sx q[2];
rz(-1.9415386) q[2];
sx q[2];
rz(0.49794751) q[2];
rz(-0.30161101) q[3];
sx q[3];
rz(-2.7409654) q[3];
sx q[3];
rz(-0.51914674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4122445) q[0];
sx q[0];
rz(-3.0817741) q[0];
sx q[0];
rz(-0.27858946) q[0];
rz(2.5623698) q[1];
sx q[1];
rz(-2.2021553) q[1];
sx q[1];
rz(3.0648807) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68966507) q[0];
sx q[0];
rz(-1.8849012) q[0];
sx q[0];
rz(0.12977022) q[0];
rz(2.9701091) q[2];
sx q[2];
rz(-1.0043317) q[2];
sx q[2];
rz(0.4085853) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0662688) q[1];
sx q[1];
rz(-1.9362209) q[1];
sx q[1];
rz(1.5484023) q[1];
rz(2.9647787) q[3];
sx q[3];
rz(-1.6578976) q[3];
sx q[3];
rz(0.84482312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3672093) q[2];
sx q[2];
rz(-0.70720208) q[2];
sx q[2];
rz(-2.6860766) q[2];
rz(0.43141836) q[3];
sx q[3];
rz(-3.016267) q[3];
sx q[3];
rz(-3.07807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(-1.0473245) q[2];
sx q[2];
rz(-2.6101255) q[2];
sx q[2];
rz(2.5892467) q[2];
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