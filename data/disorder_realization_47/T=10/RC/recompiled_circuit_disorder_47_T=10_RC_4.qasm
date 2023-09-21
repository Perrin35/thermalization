OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7597423) q[0];
sx q[0];
rz(-2.3072825) q[0];
sx q[0];
rz(0.068339737) q[0];
rz(2.4812658) q[1];
sx q[1];
rz(3.9897444) q[1];
sx q[1];
rz(6.2453649) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28107444) q[0];
sx q[0];
rz(-1.7639419) q[0];
sx q[0];
rz(-0.079695745) q[0];
rz(-pi) q[1];
rz(0.90934609) q[2];
sx q[2];
rz(-1.6644018) q[2];
sx q[2];
rz(2.0441165) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.046207) q[1];
sx q[1];
rz(-0.74828212) q[1];
sx q[1];
rz(0.074949646) q[1];
rz(-pi) q[2];
x q[2];
rz(1.171265) q[3];
sx q[3];
rz(-0.37326187) q[3];
sx q[3];
rz(1.319862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7906856) q[2];
sx q[2];
rz(-0.51351341) q[2];
sx q[2];
rz(1.8581871) q[2];
rz(0.12617271) q[3];
sx q[3];
rz(-1.4161371) q[3];
sx q[3];
rz(0.090601966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.693817) q[0];
sx q[0];
rz(-2.2651146) q[0];
sx q[0];
rz(-2.5449261) q[0];
rz(1.5860575) q[1];
sx q[1];
rz(-1.3417473) q[1];
sx q[1];
rz(1.7780875) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1648646) q[0];
sx q[0];
rz(-2.4318027) q[0];
sx q[0];
rz(-0.68517942) q[0];
rz(-pi) q[1];
rz(0.5559276) q[2];
sx q[2];
rz(-1.5276507) q[2];
sx q[2];
rz(-2.3748929) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9047782) q[1];
sx q[1];
rz(-1.292255) q[1];
sx q[1];
rz(-2.4836471) q[1];
rz(-2.3724243) q[3];
sx q[3];
rz(-2.2393919) q[3];
sx q[3];
rz(-1.9070966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8186701) q[2];
sx q[2];
rz(-2.1408036) q[2];
sx q[2];
rz(-2.4070516) q[2];
rz(2.5143886) q[3];
sx q[3];
rz(-1.6278798) q[3];
sx q[3];
rz(-0.74497765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.7221786) q[0];
sx q[0];
rz(-0.83291554) q[0];
sx q[0];
rz(0.27221361) q[0];
rz(0.84725562) q[1];
sx q[1];
rz(-1.8224742) q[1];
sx q[1];
rz(-2.8289657) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49293672) q[0];
sx q[0];
rz(-1.3515633) q[0];
sx q[0];
rz(1.6256902) q[0];
rz(2.0414511) q[2];
sx q[2];
rz(-1.0725642) q[2];
sx q[2];
rz(-0.91980308) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4362267) q[1];
sx q[1];
rz(-1.501302) q[1];
sx q[1];
rz(1.3100998) q[1];
rz(-pi) q[2];
rz(0.39194312) q[3];
sx q[3];
rz(-0.80229811) q[3];
sx q[3];
rz(2.9585569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.12038885) q[2];
sx q[2];
rz(-2.618232) q[2];
sx q[2];
rz(-2.3266501) q[2];
rz(-2.1728544) q[3];
sx q[3];
rz(-1.0453984) q[3];
sx q[3];
rz(-2.708784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75354904) q[0];
sx q[0];
rz(-2.9409445) q[0];
sx q[0];
rz(-2.5182305) q[0];
rz(-2.3240044) q[1];
sx q[1];
rz(-0.31660429) q[1];
sx q[1];
rz(3.0923016) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.277963) q[0];
sx q[0];
rz(-2.3290714) q[0];
sx q[0];
rz(0.066594007) q[0];
rz(1.2055779) q[2];
sx q[2];
rz(-0.74539241) q[2];
sx q[2];
rz(3.0129907) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.66344391) q[1];
sx q[1];
rz(-1.899292) q[1];
sx q[1];
rz(0.25681396) q[1];
rz(-pi) q[2];
rz(-0.87818273) q[3];
sx q[3];
rz(-2.4947824) q[3];
sx q[3];
rz(-1.9424903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4108882) q[2];
sx q[2];
rz(-1.5521908) q[2];
sx q[2];
rz(-0.98497406) q[2];
rz(-2.2385521) q[3];
sx q[3];
rz(-1.1629546) q[3];
sx q[3];
rz(0.27339098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34489283) q[0];
sx q[0];
rz(-1.6051689) q[0];
sx q[0];
rz(3.0539736) q[0];
rz(0.15631974) q[1];
sx q[1];
rz(-2.5904398) q[1];
sx q[1];
rz(0.87096754) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0098457355) q[0];
sx q[0];
rz(-1.5421268) q[0];
sx q[0];
rz(2.8542551) q[0];
rz(2.9882177) q[2];
sx q[2];
rz(-1.181029) q[2];
sx q[2];
rz(0.49488059) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.761913) q[1];
sx q[1];
rz(-1.2907791) q[1];
sx q[1];
rz(-2.4680087) q[1];
rz(-pi) q[2];
rz(-0.17467588) q[3];
sx q[3];
rz(-0.81411618) q[3];
sx q[3];
rz(-1.0604309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0481723) q[2];
sx q[2];
rz(-0.8478567) q[2];
sx q[2];
rz(-2.7887153) q[2];
rz(2.2757163) q[3];
sx q[3];
rz(-2.4398949) q[3];
sx q[3];
rz(0.29233366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46309328) q[0];
sx q[0];
rz(-1.1266288) q[0];
sx q[0];
rz(3.034806) q[0];
rz(-1.1865901) q[1];
sx q[1];
rz(-1.0266961) q[1];
sx q[1];
rz(-1.0245163) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86040076) q[0];
sx q[0];
rz(-1.8587347) q[0];
sx q[0];
rz(2.6698551) q[0];
rz(0.61667911) q[2];
sx q[2];
rz(-0.70313912) q[2];
sx q[2];
rz(2.6577735) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7926327) q[1];
sx q[1];
rz(-2.0434725) q[1];
sx q[1];
rz(-0.93797586) q[1];
rz(1.2023737) q[3];
sx q[3];
rz(-2.0344067) q[3];
sx q[3];
rz(1.1045477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7496877) q[2];
sx q[2];
rz(-1.2140032) q[2];
sx q[2];
rz(-2.270703) q[2];
rz(1.332256) q[3];
sx q[3];
rz(-2.199316) q[3];
sx q[3];
rz(1.4107305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.1720599) q[0];
sx q[0];
rz(-1.4071858) q[0];
sx q[0];
rz(-2.5906738) q[0];
rz(0.46539601) q[1];
sx q[1];
rz(-2.4702563) q[1];
sx q[1];
rz(-2.904772) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2315002) q[0];
sx q[0];
rz(-1.6084558) q[0];
sx q[0];
rz(1.670027) q[0];
rz(-pi) q[1];
rz(0.15010712) q[2];
sx q[2];
rz(-1.7068212) q[2];
sx q[2];
rz(-2.0748236) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.19141087) q[1];
sx q[1];
rz(-1.6457874) q[1];
sx q[1];
rz(-1.8412776) q[1];
rz(2.3051466) q[3];
sx q[3];
rz(-1.340938) q[3];
sx q[3];
rz(-3.0695855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.44234309) q[2];
sx q[2];
rz(-2.4599059) q[2];
sx q[2];
rz(2.701475) q[2];
rz(-2.0464499) q[3];
sx q[3];
rz(-1.6710072) q[3];
sx q[3];
rz(-1.7949036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31297627) q[0];
sx q[0];
rz(-1.3422817) q[0];
sx q[0];
rz(-3.112088) q[0];
rz(2.1942031) q[1];
sx q[1];
rz(-2.3209929) q[1];
sx q[1];
rz(3.0227919) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7293852) q[0];
sx q[0];
rz(-2.8280624) q[0];
sx q[0];
rz(0.14051147) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0256836) q[2];
sx q[2];
rz(-1.9695749) q[2];
sx q[2];
rz(1.1023956) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5549705) q[1];
sx q[1];
rz(-1.277031) q[1];
sx q[1];
rz(2.657386) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0766255) q[3];
sx q[3];
rz(-1.4255187) q[3];
sx q[3];
rz(-2.6165917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7416731) q[2];
sx q[2];
rz(-0.46367773) q[2];
sx q[2];
rz(-1.5681533) q[2];
rz(1.9741156) q[3];
sx q[3];
rz(-1.7220595) q[3];
sx q[3];
rz(-1.8673816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.2385999) q[0];
sx q[0];
rz(-0.82505834) q[0];
sx q[0];
rz(0.15326823) q[0];
rz(-1.0614456) q[1];
sx q[1];
rz(-0.57466424) q[1];
sx q[1];
rz(1.9877888) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.783219) q[0];
sx q[0];
rz(-1.3010539) q[0];
sx q[0];
rz(-0.37958919) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2112054) q[2];
sx q[2];
rz(-0.79158917) q[2];
sx q[2];
rz(1.1861578) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8037387) q[1];
sx q[1];
rz(-1.6263464) q[1];
sx q[1];
rz(1.3615723) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4235751) q[3];
sx q[3];
rz(-0.69268337) q[3];
sx q[3];
rz(2.6976531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8447421) q[2];
sx q[2];
rz(-1.9694318) q[2];
sx q[2];
rz(-1.6142169) q[2];
rz(0.99689233) q[3];
sx q[3];
rz(-1.4930054) q[3];
sx q[3];
rz(-2.2629288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24755724) q[0];
sx q[0];
rz(-2.0443125) q[0];
sx q[0];
rz(0.19009185) q[0];
rz(0.67063531) q[1];
sx q[1];
rz(-1.1963444) q[1];
sx q[1];
rz(-1.6533096) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9207536) q[0];
sx q[0];
rz(-1.5399884) q[0];
sx q[0];
rz(3.0160745) q[0];
rz(-pi) q[1];
rz(1.8411631) q[2];
sx q[2];
rz(-1.9989982) q[2];
sx q[2];
rz(-2.2141475) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.30257672) q[1];
sx q[1];
rz(-0.95514983) q[1];
sx q[1];
rz(2.2018196) q[1];
x q[2];
rz(2.2404352) q[3];
sx q[3];
rz(-2.3271051) q[3];
sx q[3];
rz(-0.27064532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.0037447475) q[2];
sx q[2];
rz(-0.90128428) q[2];
sx q[2];
rz(0.89912644) q[2];
rz(-0.51504618) q[3];
sx q[3];
rz(-0.47596541) q[3];
sx q[3];
rz(-2.9639444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0777733) q[0];
sx q[0];
rz(-1.780092) q[0];
sx q[0];
rz(-0.5583981) q[0];
rz(2.8521815) q[1];
sx q[1];
rz(-2.2402973) q[1];
sx q[1];
rz(-1.4351861) q[1];
rz(0.32439705) q[2];
sx q[2];
rz(-0.73642052) q[2];
sx q[2];
rz(-3.0036075) q[2];
rz(1.8437456) q[3];
sx q[3];
rz(-1.6932586) q[3];
sx q[3];
rz(0.3441588) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
