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
rz(0.02286214) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6829553) q[0];
sx q[0];
rz(-1.6874755) q[0];
sx q[0];
rz(1.4215683) q[0];
rz(-pi) q[1];
x q[1];
rz(0.68732287) q[2];
sx q[2];
rz(-2.5097373) q[2];
sx q[2];
rz(-0.99422115) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7658246) q[1];
sx q[1];
rz(-1.7205015) q[1];
sx q[1];
rz(-1.2234115) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25306563) q[3];
sx q[3];
rz(-2.1742637) q[3];
sx q[3];
rz(2.8367708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7227398) q[2];
sx q[2];
rz(-0.8272233) q[2];
sx q[2];
rz(-2.375405) q[2];
rz(2.9700759) q[3];
sx q[3];
rz(-0.73232108) q[3];
sx q[3];
rz(0.58656251) q[3];
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
sx q[0];
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
rz(-0.08509732) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.800195) q[0];
sx q[0];
rz(-1.5872103) q[0];
sx q[0];
rz(2.5192758) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1413967) q[2];
sx q[2];
rz(-2.4944802) q[2];
sx q[2];
rz(-2.5246758) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7460691) q[1];
sx q[1];
rz(-1.6568526) q[1];
sx q[1];
rz(1.6834016) q[1];
rz(-pi) q[2];
rz(0.3195023) q[3];
sx q[3];
rz(-2.0204244) q[3];
sx q[3];
rz(-0.055469661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8791085) q[2];
sx q[2];
rz(-1.5510677) q[2];
sx q[2];
rz(-1.4651728) q[2];
rz(-1.9654467) q[3];
sx q[3];
rz(-2.2361103) q[3];
sx q[3];
rz(0.24648497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.6572606) q[0];
sx q[0];
rz(-3.0069139) q[0];
sx q[0];
rz(-0.9056257) q[0];
rz(2.8088645) q[1];
sx q[1];
rz(-2.2867124) q[1];
sx q[1];
rz(-1.3844301) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4577643) q[0];
sx q[0];
rz(-2.0450852) q[0];
sx q[0];
rz(-1.6636687) q[0];
rz(-0.39321123) q[2];
sx q[2];
rz(-2.0438497) q[2];
sx q[2];
rz(2.2676603) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3789931) q[1];
sx q[1];
rz(-1.5484719) q[1];
sx q[1];
rz(-0.75548817) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.82377394) q[3];
sx q[3];
rz(-0.71411055) q[3];
sx q[3];
rz(-0.87574524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8973792) q[2];
sx q[2];
rz(-2.801557) q[2];
sx q[2];
rz(-1.8133694) q[2];
rz(-1.3978847) q[3];
sx q[3];
rz(-0.54026794) q[3];
sx q[3];
rz(-3.0298997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
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
rz(-1.1616511) q[1];
sx q[1];
rz(0.20733325) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57732338) q[0];
sx q[0];
rz(-0.042357001) q[0];
sx q[0];
rz(1.7839412) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7423332) q[2];
sx q[2];
rz(-0.17645391) q[2];
sx q[2];
rz(1.8432957) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.15370788) q[1];
sx q[1];
rz(-2.8578836) q[1];
sx q[1];
rz(-0.99516408) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6299423) q[3];
sx q[3];
rz(-1.5582065) q[3];
sx q[3];
rz(0.97012855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.31071445) q[2];
sx q[2];
rz(-1.0881311) q[2];
sx q[2];
rz(-0.80835289) q[2];
rz(-0.80777848) q[3];
sx q[3];
rz(-2.419796) q[3];
sx q[3];
rz(0.19367735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5508674) q[0];
sx q[0];
rz(-0.71722513) q[0];
sx q[0];
rz(2.2461058) q[0];
rz(0.26516178) q[1];
sx q[1];
rz(-0.83509713) q[1];
sx q[1];
rz(-0.53363824) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9981209) q[0];
sx q[0];
rz(-0.40092418) q[0];
sx q[0];
rz(0.801416) q[0];
rz(-pi) q[1];
rz(-1.7970423) q[2];
sx q[2];
rz(-2.3257253) q[2];
sx q[2];
rz(-2.0362542) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.071760885) q[1];
sx q[1];
rz(-1.3511718) q[1];
sx q[1];
rz(1.4211618) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.94138937) q[3];
sx q[3];
rz(-1.4038868) q[3];
sx q[3];
rz(0.0080136673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9020033) q[2];
sx q[2];
rz(-1.5633554) q[2];
sx q[2];
rz(-2.5732102) q[2];
rz(1.9249632) q[3];
sx q[3];
rz(-2.670848) q[3];
sx q[3];
rz(0.34354982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(0.50826532) q[0];
sx q[0];
rz(-0.87843043) q[0];
sx q[0];
rz(-1.1761965) q[0];
rz(-1.5559224) q[1];
sx q[1];
rz(-0.95247477) q[1];
sx q[1];
rz(2.1369381) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5414808) q[0];
sx q[0];
rz(-1.2280557) q[0];
sx q[0];
rz(0.39370763) q[0];
rz(-pi) q[1];
rz(-2.9506358) q[2];
sx q[2];
rz(-2.7381884) q[2];
sx q[2];
rz(-1.0952589) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.172563) q[1];
sx q[1];
rz(-1.3247715) q[1];
sx q[1];
rz(-0.56117705) q[1];
rz(-pi) q[2];
rz(-1.1569571) q[3];
sx q[3];
rz(-1.5598179) q[3];
sx q[3];
rz(-1.4423086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4040318) q[2];
sx q[2];
rz(-2.1514251) q[2];
sx q[2];
rz(-2.693434) q[2];
rz(0.29799497) q[3];
sx q[3];
rz(-2.4500676) q[3];
sx q[3];
rz(-2.7929849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77666831) q[0];
sx q[0];
rz(-1.8402599) q[0];
sx q[0];
rz(2.4834494) q[0];
rz(-1.2843885) q[1];
sx q[1];
rz(-0.45416608) q[1];
sx q[1];
rz(-2.4694494) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1561688) q[0];
sx q[0];
rz(-2.9754313) q[0];
sx q[0];
rz(1.5260494) q[0];
x q[1];
rz(-2.0124112) q[2];
sx q[2];
rz(-1.1814983) q[2];
sx q[2];
rz(-1.0600952) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0969095) q[1];
sx q[1];
rz(-1.204856) q[1];
sx q[1];
rz(-2.6198322) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.70127212) q[3];
sx q[3];
rz(-1.4605195) q[3];
sx q[3];
rz(-0.14227223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.60823524) q[2];
sx q[2];
rz(-2.50864) q[2];
sx q[2];
rz(2.7427924) q[2];
rz(2.3170025) q[3];
sx q[3];
rz(-1.3922858) q[3];
sx q[3];
rz(-2.5430172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7010715) q[0];
sx q[0];
rz(-2.7109787) q[0];
sx q[0];
rz(0.15821247) q[0];
rz(-2.054706) q[1];
sx q[1];
rz(-0.67651665) q[1];
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
x q[1];
rz(1.4629455) q[2];
sx q[2];
rz(-1.3662405) q[2];
sx q[2];
rz(-0.66609913) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.403277) q[1];
sx q[1];
rz(-0.74790955) q[1];
sx q[1];
rz(2.7110389) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1095803) q[3];
sx q[3];
rz(-1.8714928) q[3];
sx q[3];
rz(2.7923982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.57758254) q[2];
sx q[2];
rz(-1.6489886) q[2];
sx q[2];
rz(-2.1419443) q[2];
rz(-0.10351652) q[3];
sx q[3];
rz(-3.0045356) q[3];
sx q[3];
rz(-1.1243533) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12061159) q[0];
sx q[0];
rz(-0.7779026) q[0];
sx q[0];
rz(0.65752423) q[0];
rz(-2.855037) q[1];
sx q[1];
rz(-0.92266881) q[1];
sx q[1];
rz(2.8009169) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1239615) q[0];
sx q[0];
rz(-2.258856) q[0];
sx q[0];
rz(0.4171564) q[0];
rz(-pi) q[1];
rz(-1.9129487) q[2];
sx q[2];
rz(-2.5492382) q[2];
sx q[2];
rz(-0.60281384) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9790736) q[1];
sx q[1];
rz(-0.58774978) q[1];
sx q[1];
rz(1.2390922) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4021923) q[3];
sx q[3];
rz(-1.3753034) q[3];
sx q[3];
rz(0.80381264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.75491536) q[2];
sx q[2];
rz(-1.200054) q[2];
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
rz(-pi/2) q[2];
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
rz(-0.07671193) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3007606) q[0];
sx q[0];
rz(-1.4474086) q[0];
sx q[0];
rz(-1.2542017) q[0];
rz(-pi) q[1];
rz(-0.17148359) q[2];
sx q[2];
rz(-2.137261) q[2];
sx q[2];
rz(2.7330074) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0753239) q[1];
sx q[1];
rz(-1.9362209) q[1];
sx q[1];
rz(1.5931904) q[1];
rz(-pi) q[2];
x q[2];
rz(0.46080188) q[3];
sx q[3];
rz(-0.19690234) q[3];
sx q[3];
rz(-0.27289665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.77438337) q[2];
sx q[2];
rz(-2.4343906) q[2];
sx q[2];
rz(2.6860766) q[2];
rz(0.43141836) q[3];
sx q[3];
rz(-3.016267) q[3];
sx q[3];
rz(0.063522696) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
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
rz(-2.0942681) q[2];
sx q[2];
rz(-0.53146711) q[2];
sx q[2];
rz(-0.55234595) q[2];
rz(2.4424845) q[3];
sx q[3];
rz(-0.53634488) q[3];
sx q[3];
rz(2.8041822) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
