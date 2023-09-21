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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.011933) q[0];
sx q[0];
rz(-1.4225905) q[0];
sx q[0];
rz(3.0236142) q[0];
x q[1];
rz(2.4542698) q[2];
sx q[2];
rz(-2.5097373) q[2];
sx q[2];
rz(-2.1473715) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.37576807) q[1];
sx q[1];
rz(-1.7205015) q[1];
sx q[1];
rz(1.2234115) q[1];
rz(-pi) q[2];
rz(-0.25306563) q[3];
sx q[3];
rz(-2.1742637) q[3];
sx q[3];
rz(-2.8367708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.4188529) q[2];
sx q[2];
rz(-0.8272233) q[2];
sx q[2];
rz(-0.76618761) q[2];
rz(0.17151672) q[3];
sx q[3];
rz(-2.4092716) q[3];
sx q[3];
rz(0.58656251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3294753) q[0];
sx q[0];
rz(-2.7024039) q[0];
sx q[0];
rz(-3.0549333) q[0];
rz(2.2791729) q[1];
sx q[1];
rz(-2.5984867) q[1];
sx q[1];
rz(-0.08509732) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.800195) q[0];
sx q[0];
rz(-1.5872103) q[0];
sx q[0];
rz(-0.62231681) q[0];
x q[1];
rz(2.1727824) q[2];
sx q[2];
rz(-1.8245056) q[2];
sx q[2];
rz(1.3041376) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7460691) q[1];
sx q[1];
rz(-1.4847401) q[1];
sx q[1];
rz(-1.458191) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3195023) q[3];
sx q[3];
rz(-2.0204244) q[3];
sx q[3];
rz(3.086123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.26248419) q[2];
sx q[2];
rz(-1.5905249) q[2];
sx q[2];
rz(1.4651728) q[2];
rz(1.1761459) q[3];
sx q[3];
rz(-2.2361103) q[3];
sx q[3];
rz(0.24648497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
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
rz(-0.48433205) q[0];
sx q[0];
rz(-3.0069139) q[0];
sx q[0];
rz(2.235967) q[0];
rz(-2.8088645) q[1];
sx q[1];
rz(-0.85488027) q[1];
sx q[1];
rz(-1.3844301) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6838283) q[0];
sx q[0];
rz(-2.0450852) q[0];
sx q[0];
rz(-1.6636687) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2134174) q[2];
sx q[2];
rz(-2.536142) q[2];
sx q[2];
rz(1.6124992) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.76259957) q[1];
sx q[1];
rz(-1.5484719) q[1];
sx q[1];
rz(-0.75548817) q[1];
rz(-pi) q[2];
x q[2];
rz(2.609385) q[3];
sx q[3];
rz(-2.0720707) q[3];
sx q[3];
rz(1.7621079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8973792) q[2];
sx q[2];
rz(-2.801557) q[2];
sx q[2];
rz(-1.8133694) q[2];
rz(1.3978847) q[3];
sx q[3];
rz(-0.54026794) q[3];
sx q[3];
rz(3.0298997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
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
rz(3.0760647) q[0];
sx q[0];
rz(-2.9632443) q[0];
sx q[0];
rz(-0.67101014) q[0];
rz(1.3440075) q[1];
sx q[1];
rz(-1.9799415) q[1];
sx q[1];
rz(-2.9342594) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5642693) q[0];
sx q[0];
rz(-3.0992357) q[0];
sx q[0];
rz(1.3576515) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7447128) q[2];
sx q[2];
rz(-1.5408278) q[2];
sx q[2];
rz(2.7001691) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.15370788) q[1];
sx q[1];
rz(-0.28370902) q[1];
sx q[1];
rz(2.1464286) q[1];
rz(-pi) q[2];
rz(2.6299423) q[3];
sx q[3];
rz(-1.5833862) q[3];
sx q[3];
rz(-0.97012855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8308782) q[2];
sx q[2];
rz(-2.0534616) q[2];
sx q[2];
rz(-0.80835289) q[2];
rz(0.80777848) q[3];
sx q[3];
rz(-2.419796) q[3];
sx q[3];
rz(2.9479153) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59072524) q[0];
sx q[0];
rz(-0.71722513) q[0];
sx q[0];
rz(0.89548683) q[0];
rz(2.8764309) q[1];
sx q[1];
rz(-0.83509713) q[1];
sx q[1];
rz(-2.6079544) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9981209) q[0];
sx q[0];
rz(-0.40092418) q[0];
sx q[0];
rz(-2.3401767) q[0];
rz(-pi) q[1];
rz(-2.3737714) q[2];
sx q[2];
rz(-1.7349093) q[2];
sx q[2];
rz(-0.62190157) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.071760885) q[1];
sx q[1];
rz(-1.7904209) q[1];
sx q[1];
rz(1.4211618) q[1];
rz(-pi) q[2];
x q[2];
rz(2.93612) q[3];
sx q[3];
rz(-0.95147248) q[3];
sx q[3];
rz(1.458414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9020033) q[2];
sx q[2];
rz(-1.5782372) q[2];
sx q[2];
rz(2.5732102) q[2];
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
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6333273) q[0];
sx q[0];
rz(-2.2631622) q[0];
sx q[0];
rz(1.1761965) q[0];
rz(-1.5559224) q[1];
sx q[1];
rz(-2.1891179) q[1];
sx q[1];
rz(-2.1369381) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6001119) q[0];
sx q[0];
rz(-1.913537) q[0];
sx q[0];
rz(-0.39370763) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.19095687) q[2];
sx q[2];
rz(-0.40340427) q[2];
sx q[2];
rz(2.0463338) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.028588258) q[1];
sx q[1];
rz(-2.5341946) q[1];
sx q[1];
rz(0.44087704) q[1];
x q[2];
rz(-1.5435013) q[3];
sx q[3];
rz(-0.41397646) q[3];
sx q[3];
rz(0.15347813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4040318) q[2];
sx q[2];
rz(-2.1514251) q[2];
sx q[2];
rz(2.693434) q[2];
rz(0.29799497) q[3];
sx q[3];
rz(-0.69152504) q[3];
sx q[3];
rz(-0.34860778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3649243) q[0];
sx q[0];
rz(-1.8402599) q[0];
sx q[0];
rz(-0.65814322) q[0];
rz(-1.2843885) q[1];
sx q[1];
rz(-0.45416608) q[1];
sx q[1];
rz(-2.4694494) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9400529) q[0];
sx q[0];
rz(-1.4048028) q[0];
sx q[0];
rz(-0.0075017651) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3357046) q[2];
sx q[2];
rz(-2.561509) q[2];
sx q[2];
rz(-0.16575955) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.044683177) q[1];
sx q[1];
rz(-1.9367366) q[1];
sx q[1];
rz(2.6198322) q[1];
rz(-pi) q[2];
rz(0.16996202) q[3];
sx q[3];
rz(-2.4331577) q[3];
sx q[3];
rz(1.2989312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5333574) q[2];
sx q[2];
rz(-2.50864) q[2];
sx q[2];
rz(-2.7427924) q[2];
rz(-0.82459015) q[3];
sx q[3];
rz(-1.7493068) q[3];
sx q[3];
rz(2.5430172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7010715) q[0];
sx q[0];
rz(-2.7109787) q[0];
sx q[0];
rz(-0.15821247) q[0];
rz(-1.0868866) q[1];
sx q[1];
rz(-0.67651665) q[1];
sx q[1];
rz(-0.80668443) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0772484) q[0];
sx q[0];
rz(-1.6883388) q[0];
sx q[0];
rz(2.006152) q[0];
x q[1];
rz(-2.6629692) q[2];
sx q[2];
rz(-2.9106986) q[2];
sx q[2];
rz(-2.965197) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.73831564) q[1];
sx q[1];
rz(-2.3936831) q[1];
sx q[1];
rz(0.43055375) q[1];
rz(1.673647) q[3];
sx q[3];
rz(-2.8392483) q[3];
sx q[3];
rz(-2.9000988) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0209811) q[0];
sx q[0];
rz(-0.7779026) q[0];
sx q[0];
rz(0.65752423) q[0];
rz(-2.855037) q[1];
sx q[1];
rz(-0.92266881) q[1];
sx q[1];
rz(-0.34067571) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0176312) q[0];
sx q[0];
rz(-0.88273662) q[0];
sx q[0];
rz(2.7244363) q[0];
rz(2.9195243) q[2];
sx q[2];
rz(-2.1246315) q[2];
sx q[2];
rz(-2.1332707) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9790736) q[1];
sx q[1];
rz(-2.5538429) q[1];
sx q[1];
rz(1.2390922) q[1];
rz(-pi) q[2];
rz(-2.4386028) q[3];
sx q[3];
rz(-2.8841416) q[3];
sx q[3];
rz(-1.5233745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3866773) q[2];
sx q[2];
rz(-1.9415386) q[2];
sx q[2];
rz(-2.6436451) q[2];
rz(0.30161101) q[3];
sx q[3];
rz(-0.40062723) q[3];
sx q[3];
rz(2.6224459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72934812) q[0];
sx q[0];
rz(-0.059818581) q[0];
sx q[0];
rz(-2.8630032) q[0];
rz(-2.5623698) q[1];
sx q[1];
rz(-2.2021553) q[1];
sx q[1];
rz(-3.0648807) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3007606) q[0];
sx q[0];
rz(-1.4474086) q[0];
sx q[0];
rz(1.8873909) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9701091) q[2];
sx q[2];
rz(-2.137261) q[2];
sx q[2];
rz(-0.4085853) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0662688) q[1];
sx q[1];
rz(-1.9362209) q[1];
sx q[1];
rz(1.5931904) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.65927) q[3];
sx q[3];
rz(-1.3946597) q[3];
sx q[3];
rz(2.4000771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.77438337) q[2];
sx q[2];
rz(-2.4343906) q[2];
sx q[2];
rz(-0.45551604) q[2];
rz(-2.7101743) q[3];
sx q[3];
rz(-0.12532561) q[3];
sx q[3];
rz(-0.063522696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1800304) q[0];
sx q[0];
rz(-1.4008235) q[0];
sx q[0];
rz(0.93749198) q[0];
rz(-0.58615276) q[1];
sx q[1];
rz(-1.6393607) q[1];
sx q[1];
rz(1.6929109) q[1];
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