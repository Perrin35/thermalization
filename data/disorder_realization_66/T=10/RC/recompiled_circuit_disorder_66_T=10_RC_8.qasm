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
x q[1];
rz(-1.1360264) q[2];
sx q[2];
rz(-2.044894) q[2];
sx q[2];
rz(-2.9413162) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7658246) q[1];
sx q[1];
rz(-1.7205015) q[1];
sx q[1];
rz(-1.9181812) q[1];
rz(-pi) q[2];
rz(1.2223577) q[3];
sx q[3];
rz(-0.64823965) q[3];
sx q[3];
rz(-0.1227619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.4188529) q[2];
sx q[2];
rz(-0.8272233) q[2];
sx q[2];
rz(2.375405) q[2];
rz(0.17151672) q[3];
sx q[3];
rz(-0.73232108) q[3];
sx q[3];
rz(-0.58656251) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8121174) q[0];
sx q[0];
rz(-2.7024039) q[0];
sx q[0];
rz(-3.0549333) q[0];
rz(-0.86241972) q[1];
sx q[1];
rz(-0.54310596) q[1];
sx q[1];
rz(-3.0564953) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25227308) q[0];
sx q[0];
rz(-2.5190881) q[0];
sx q[0];
rz(-3.1134393) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.000196) q[2];
sx q[2];
rz(-0.64711249) q[2];
sx q[2];
rz(-0.61691689) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7460691) q[1];
sx q[1];
rz(-1.6568526) q[1];
sx q[1];
rz(-1.458191) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3195023) q[3];
sx q[3];
rz(-2.0204244) q[3];
sx q[3];
rz(-0.055469661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.26248419) q[2];
sx q[2];
rz(-1.5510677) q[2];
sx q[2];
rz(-1.4651728) q[2];
rz(1.1761459) q[3];
sx q[3];
rz(-2.2361103) q[3];
sx q[3];
rz(-2.8951077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48433205) q[0];
sx q[0];
rz(-3.0069139) q[0];
sx q[0];
rz(2.235967) q[0];
rz(-0.33272818) q[1];
sx q[1];
rz(-2.2867124) q[1];
sx q[1];
rz(-1.3844301) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15554409) q[0];
sx q[0];
rz(-1.4882003) q[0];
sx q[0];
rz(-0.47604576) q[0];
x q[1];
rz(-2.2134174) q[2];
sx q[2];
rz(-0.60545063) q[2];
sx q[2];
rz(1.5290934) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3789931) q[1];
sx q[1];
rz(-1.5484719) q[1];
sx q[1];
rz(0.75548817) q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.8973792) q[2];
sx q[2];
rz(-2.801557) q[2];
sx q[2];
rz(-1.3282233) q[2];
rz(-1.7437079) q[3];
sx q[3];
rz(-2.6013247) q[3];
sx q[3];
rz(-3.0298997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0760647) q[0];
sx q[0];
rz(-0.17834839) q[0];
sx q[0];
rz(-2.4705825) q[0];
rz(1.7975851) q[1];
sx q[1];
rz(-1.1616511) q[1];
sx q[1];
rz(-2.9342594) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57732338) q[0];
sx q[0];
rz(-0.042357001) q[0];
sx q[0];
rz(-1.7839412) q[0];
x q[1];
rz(1.7447128) q[2];
sx q[2];
rz(-1.6007649) q[2];
sx q[2];
rz(0.44142351) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.15370788) q[1];
sx q[1];
rz(-0.28370902) q[1];
sx q[1];
rz(-2.1464286) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6299423) q[3];
sx q[3];
rz(-1.5582065) q[3];
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
rz(2.3332398) q[2];
rz(0.80777848) q[3];
sx q[3];
rz(-0.72179663) q[3];
sx q[3];
rz(0.19367735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59072524) q[0];
sx q[0];
rz(-0.71722513) q[0];
sx q[0];
rz(-0.89548683) q[0];
rz(2.8764309) q[1];
sx q[1];
rz(-2.3064955) q[1];
sx q[1];
rz(2.6079544) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14347178) q[0];
sx q[0];
rz(-0.40092418) q[0];
sx q[0];
rz(0.801416) q[0];
rz(-pi) q[1];
rz(1.3445504) q[2];
sx q[2];
rz(-2.3257253) q[2];
sx q[2];
rz(-2.0362542) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.071760885) q[1];
sx q[1];
rz(-1.7904209) q[1];
sx q[1];
rz(-1.4211618) q[1];
x q[2];
rz(2.93612) q[3];
sx q[3];
rz(-2.1901202) q[3];
sx q[3];
rz(-1.458414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9020033) q[2];
sx q[2];
rz(-1.5782372) q[2];
sx q[2];
rz(0.56838244) q[2];
rz(-1.9249632) q[3];
sx q[3];
rz(-2.670848) q[3];
sx q[3];
rz(2.7980428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50826532) q[0];
sx q[0];
rz(-2.2631622) q[0];
sx q[0];
rz(1.1761965) q[0];
rz(-1.5559224) q[1];
sx q[1];
rz(-2.1891179) q[1];
sx q[1];
rz(-2.1369381) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3096136) q[0];
sx q[0];
rz(-1.9404611) q[0];
sx q[0];
rz(1.9395104) q[0];
rz(2.7447694) q[2];
sx q[2];
rz(-1.6453711) q[2];
sx q[2];
rz(2.4900988) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9690296) q[1];
sx q[1];
rz(-1.3247715) q[1];
sx q[1];
rz(2.5804156) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9846356) q[3];
sx q[3];
rz(-1.5598179) q[3];
sx q[3];
rz(-1.699284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4040318) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3649243) q[0];
sx q[0];
rz(-1.8402599) q[0];
sx q[0];
rz(2.4834494) q[0];
rz(-1.2843885) q[1];
sx q[1];
rz(-0.45416608) q[1];
sx q[1];
rz(-2.4694494) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7710966) q[0];
sx q[0];
rz(-1.578195) q[0];
sx q[0];
rz(-1.7367944) q[0];
rz(0.42598287) q[2];
sx q[2];
rz(-1.9773215) q[2];
sx q[2];
rz(2.4533518) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0587412) q[1];
sx q[1];
rz(-0.62742678) q[1];
sx q[1];
rz(-2.4861366) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4403205) q[3];
sx q[3];
rz(-1.6810732) q[3];
sx q[3];
rz(2.9993204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.60823524) q[2];
sx q[2];
rz(-2.50864) q[2];
sx q[2];
rz(0.39880025) q[2];
rz(2.3170025) q[3];
sx q[3];
rz(-1.7493068) q[3];
sx q[3];
rz(-0.59857541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44052112) q[0];
sx q[0];
rz(-0.43061391) q[0];
sx q[0];
rz(0.15821247) q[0];
rz(1.0868866) q[1];
sx q[1];
rz(-2.465076) q[1];
sx q[1];
rz(-0.80668443) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54803941) q[0];
sx q[0];
rz(-1.1386477) q[0];
sx q[0];
rz(-0.12950626) q[0];
x q[1];
rz(-0.47862349) q[2];
sx q[2];
rz(-0.23089409) q[2];
sx q[2];
rz(0.17639562) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.17864171) q[1];
sx q[1];
rz(-2.2370403) q[1];
sx q[1];
rz(1.2013749) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2699548) q[3];
sx q[3];
rz(-1.6013718) q[3];
sx q[3];
rz(-1.9105063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5640101) q[2];
sx q[2];
rz(-1.492604) q[2];
sx q[2];
rz(-2.1419443) q[2];
rz(3.0380761) q[3];
sx q[3];
rz(-3.0045356) q[3];
sx q[3];
rz(-1.1243533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0176312) q[0];
sx q[0];
rz(-0.88273662) q[0];
sx q[0];
rz(2.7244363) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9195243) q[2];
sx q[2];
rz(-1.0169612) q[2];
sx q[2];
rz(-2.1332707) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3714911) q[1];
sx q[1];
rz(-2.1226468) q[1];
sx q[1];
rz(0.21367195) q[1];
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
rz(-0.75491536) q[2];
sx q[2];
rz(-1.200054) q[2];
sx q[2];
rz(0.49794751) q[2];
rz(-0.30161101) q[3];
sx q[3];
rz(-0.40062723) q[3];
sx q[3];
rz(0.51914674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4122445) q[0];
sx q[0];
rz(-0.059818581) q[0];
sx q[0];
rz(2.8630032) q[0];
rz(0.57922286) q[1];
sx q[1];
rz(-0.93943739) q[1];
sx q[1];
rz(-0.07671193) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0893223) q[0];
sx q[0];
rz(-0.33903402) q[0];
sx q[0];
rz(1.9498755) q[0];
rz(2.1439728) q[2];
sx q[2];
rz(-1.4263037) q[2];
sx q[2];
rz(-1.0695374) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.50347603) q[1];
sx q[1];
rz(-1.5498811) q[1];
sx q[1];
rz(-0.36550826) q[1];
x q[2];
rz(-2.6807908) q[3];
sx q[3];
rz(-0.19690234) q[3];
sx q[3];
rz(-0.27289665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3672093) q[2];
sx q[2];
rz(-0.70720208) q[2];
sx q[2];
rz(2.6860766) q[2];
rz(-2.7101743) q[3];
sx q[3];
rz(-3.016267) q[3];
sx q[3];
rz(0.063522696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9615622) q[0];
sx q[0];
rz(-1.7407692) q[0];
sx q[0];
rz(-2.2041007) q[0];
rz(2.5554399) q[1];
sx q[1];
rz(-1.6393607) q[1];
sx q[1];
rz(1.6929109) q[1];
rz(-0.2858328) q[2];
sx q[2];
rz(-2.0252068) q[2];
sx q[2];
rz(0.037638738) q[2];
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
