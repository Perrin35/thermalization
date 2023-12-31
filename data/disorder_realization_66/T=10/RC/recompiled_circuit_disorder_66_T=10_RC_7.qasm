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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6883144) q[0];
sx q[0];
rz(-2.9524321) q[0];
sx q[0];
rz(-0.90318371) q[0];
rz(-pi) q[1];
rz(2.0055662) q[2];
sx q[2];
rz(-1.0966986) q[2];
sx q[2];
rz(-0.20027645) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.248977) q[1];
sx q[1];
rz(-1.9141344) q[1];
sx q[1];
rz(2.9825319) q[1];
rz(1.919235) q[3];
sx q[3];
rz(-2.493353) q[3];
sx q[3];
rz(3.0188308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.4188529) q[2];
sx q[2];
rz(-2.3143694) q[2];
sx q[2];
rz(-0.76618761) q[2];
rz(0.17151672) q[3];
sx q[3];
rz(-0.73232108) q[3];
sx q[3];
rz(2.5550301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3294753) q[0];
sx q[0];
rz(-0.43918878) q[0];
sx q[0];
rz(-3.0549333) q[0];
rz(2.2791729) q[1];
sx q[1];
rz(-0.54310596) q[1];
sx q[1];
rz(-3.0564953) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9239685) q[0];
sx q[0];
rz(-0.94857615) q[0];
sx q[0];
rz(-1.5505962) q[0];
rz(-pi) q[1];
x q[1];
rz(0.30479635) q[2];
sx q[2];
rz(-0.99064231) q[2];
sx q[2];
rz(3.0457029) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7460691) q[1];
sx q[1];
rz(-1.4847401) q[1];
sx q[1];
rz(-1.6834016) q[1];
rz(-0.3195023) q[3];
sx q[3];
rz(-1.1211683) q[3];
sx q[3];
rz(-0.055469661) q[3];
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
rz(-0.90548235) q[3];
sx q[3];
rz(-0.24648497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48433205) q[0];
sx q[0];
rz(-3.0069139) q[0];
sx q[0];
rz(0.9056257) q[0];
rz(2.8088645) q[1];
sx q[1];
rz(-0.85488027) q[1];
sx q[1];
rz(-1.7571626) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4577643) q[0];
sx q[0];
rz(-1.0965075) q[0];
sx q[0];
rz(1.6636687) q[0];
x q[1];
rz(0.92817523) q[2];
sx q[2];
rz(-0.60545063) q[2];
sx q[2];
rz(1.5290934) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3123734) q[1];
sx q[1];
rz(-0.81554283) q[1];
sx q[1];
rz(-1.5401328) q[1];
x q[2];
rz(2.3178187) q[3];
sx q[3];
rz(-0.71411055) q[3];
sx q[3];
rz(-0.87574524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8973792) q[2];
sx q[2];
rz(-0.34003568) q[2];
sx q[2];
rz(-1.3282233) q[2];
rz(-1.3978847) q[3];
sx q[3];
rz(-2.6013247) q[3];
sx q[3];
rz(3.0298997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0760647) q[0];
sx q[0];
rz(-2.9632443) q[0];
sx q[0];
rz(2.4705825) q[0];
rz(1.3440075) q[1];
sx q[1];
rz(-1.9799415) q[1];
sx q[1];
rz(0.20733325) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3610791) q[0];
sx q[0];
rz(-1.5797537) q[0];
sx q[0];
rz(-1.6121959) q[0];
rz(1.7423332) q[2];
sx q[2];
rz(-0.17645391) q[2];
sx q[2];
rz(1.2982969) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.15370788) q[1];
sx q[1];
rz(-0.28370902) q[1];
sx q[1];
rz(0.99516408) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.585235) q[3];
sx q[3];
rz(-1.0591905) q[3];
sx q[3];
rz(-2.5479941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59072524) q[0];
sx q[0];
rz(-2.4243675) q[0];
sx q[0];
rz(-2.2461058) q[0];
rz(2.8764309) q[1];
sx q[1];
rz(-0.83509713) q[1];
sx q[1];
rz(0.53363824) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4424608) q[0];
sx q[0];
rz(-1.8457544) q[0];
sx q[0];
rz(-1.2752227) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7970423) q[2];
sx q[2];
rz(-0.81586736) q[2];
sx q[2];
rz(-2.0362542) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.071760885) q[1];
sx q[1];
rz(-1.7904209) q[1];
sx q[1];
rz(-1.4211618) q[1];
rz(-2.93612) q[3];
sx q[3];
rz(-2.1901202) q[3];
sx q[3];
rz(1.458414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2395893) q[2];
sx q[2];
rz(-1.5782372) q[2];
sx q[2];
rz(0.56838244) q[2];
rz(-1.2166294) q[3];
sx q[3];
rz(-2.670848) q[3];
sx q[3];
rz(-2.7980428) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6333273) q[0];
sx q[0];
rz(-0.87843043) q[0];
sx q[0];
rz(-1.9653962) q[0];
rz(-1.5559224) q[1];
sx q[1];
rz(-2.1891179) q[1];
sx q[1];
rz(1.0046545) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6509103) q[0];
sx q[0];
rz(-2.6255529) q[0];
sx q[0];
rz(2.3923621) q[0];
rz(-0.19095687) q[2];
sx q[2];
rz(-2.7381884) q[2];
sx q[2];
rz(-2.0463338) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.028588258) q[1];
sx q[1];
rz(-2.5341946) q[1];
sx q[1];
rz(-2.7007156) q[1];
x q[2];
rz(-0.011990487) q[3];
sx q[3];
rz(-1.1569835) q[3];
sx q[3];
rz(0.12366611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4040318) q[2];
sx q[2];
rz(-0.99016756) q[2];
sx q[2];
rz(-0.44815865) q[2];
rz(0.29799497) q[3];
sx q[3];
rz(-2.4500676) q[3];
sx q[3];
rz(-2.7929849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77666831) q[0];
sx q[0];
rz(-1.8402599) q[0];
sx q[0];
rz(-2.4834494) q[0];
rz(1.8572042) q[1];
sx q[1];
rz(-0.45416608) q[1];
sx q[1];
rz(-2.4694494) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9854239) q[0];
sx q[0];
rz(-0.16616136) q[0];
sx q[0];
rz(-1.6155433) q[0];
rz(-pi) q[1];
rz(-2.7156098) q[2];
sx q[2];
rz(-1.1642712) q[2];
sx q[2];
rz(-2.4533518) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0828515) q[1];
sx q[1];
rz(-0.62742678) q[1];
sx q[1];
rz(2.4861366) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[1];
rz(-2.5333574) q[2];
sx q[2];
rz(-0.63295263) q[2];
sx q[2];
rz(-0.39880025) q[2];
rz(2.3170025) q[3];
sx q[3];
rz(-1.7493068) q[3];
sx q[3];
rz(2.5430172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44052112) q[0];
sx q[0];
rz(-2.7109787) q[0];
sx q[0];
rz(2.9833802) q[0];
rz(1.0868866) q[1];
sx q[1];
rz(-0.67651665) q[1];
sx q[1];
rz(-2.3349082) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8950302) q[0];
sx q[0];
rz(-2.6916305) q[0];
sx q[0];
rz(-1.8438086) q[0];
rz(-pi) q[1];
x q[1];
rz(0.20571795) q[2];
sx q[2];
rz(-1.6763902) q[2];
sx q[2];
rz(-0.92668698) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.403277) q[1];
sx q[1];
rz(-2.3936831) q[1];
sx q[1];
rz(-0.43055375) q[1];
rz(1.673647) q[3];
sx q[3];
rz(-2.8392483) q[3];
sx q[3];
rz(-2.9000988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5640101) q[2];
sx q[2];
rz(-1.492604) q[2];
sx q[2];
rz(-2.1419443) q[2];
rz(0.10351652) q[3];
sx q[3];
rz(-3.0045356) q[3];
sx q[3];
rz(1.1243533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-2.3636901) q[0];
sx q[0];
rz(0.65752423) q[0];
rz(-2.855037) q[1];
sx q[1];
rz(-2.2189238) q[1];
sx q[1];
rz(-2.8009169) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7332471) q[0];
sx q[0];
rz(-0.78662965) q[0];
sx q[0];
rz(-2.0287081) q[0];
rz(-pi) q[1];
rz(-0.22206837) q[2];
sx q[2];
rz(-1.0169612) q[2];
sx q[2];
rz(2.1332707) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.454168) q[1];
sx q[1];
rz(-1.3892281) q[1];
sx q[1];
rz(2.1329692) q[1];
rz(-pi) q[2];
rz(1.4021923) q[3];
sx q[3];
rz(-1.7662893) q[3];
sx q[3];
rz(2.33778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.75491536) q[2];
sx q[2];
rz(-1.200054) q[2];
sx q[2];
rz(0.49794751) q[2];
rz(-2.8399816) q[3];
sx q[3];
rz(-0.40062723) q[3];
sx q[3];
rz(-0.51914674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4122445) q[0];
sx q[0];
rz(-0.059818581) q[0];
sx q[0];
rz(2.8630032) q[0];
rz(-2.5623698) q[1];
sx q[1];
rz(-2.2021553) q[1];
sx q[1];
rz(0.07671193) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.840832) q[0];
sx q[0];
rz(-1.6941841) q[0];
sx q[0];
rz(-1.2542017) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1439728) q[2];
sx q[2];
rz(-1.4263037) q[2];
sx q[2];
rz(2.0720553) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0753239) q[1];
sx q[1];
rz(-1.9362209) q[1];
sx q[1];
rz(-1.5931904) q[1];
rz(-pi) q[2];
x q[2];
rz(1.65927) q[3];
sx q[3];
rz(-1.7469329) q[3];
sx q[3];
rz(2.4000771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.77438337) q[2];
sx q[2];
rz(-0.70720208) q[2];
sx q[2];
rz(2.6860766) q[2];
rz(-2.7101743) q[3];
sx q[3];
rz(-3.016267) q[3];
sx q[3];
rz(-3.07807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1800304) q[0];
sx q[0];
rz(-1.7407692) q[0];
sx q[0];
rz(-2.2041007) q[0];
rz(2.5554399) q[1];
sx q[1];
rz(-1.6393607) q[1];
sx q[1];
rz(1.6929109) q[1];
rz(2.8557599) q[2];
sx q[2];
rz(-2.0252068) q[2];
sx q[2];
rz(0.037638738) q[2];
rz(-1.9361817) q[3];
sx q[3];
rz(-1.972651) q[3];
sx q[3];
rz(-2.7046711) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
