OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3818504) q[0];
sx q[0];
rz(-0.83431017) q[0];
sx q[0];
rz(3.0732529) q[0];
rz(2.4812658) q[1];
sx q[1];
rz(-2.2934409) q[1];
sx q[1];
rz(3.1037722) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28107444) q[0];
sx q[0];
rz(-1.7639419) q[0];
sx q[0];
rz(-3.0618969) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2322466) q[2];
sx q[2];
rz(-1.6644018) q[2];
sx q[2];
rz(-2.0441165) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6112069) q[1];
sx q[1];
rz(-1.5198277) q[1];
sx q[1];
rz(-2.3947122) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1511729) q[3];
sx q[3];
rz(-1.9133948) q[3];
sx q[3];
rz(1.3959988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.35090703) q[2];
sx q[2];
rz(-0.51351341) q[2];
sx q[2];
rz(-1.2834056) q[2];
rz(-0.12617271) q[3];
sx q[3];
rz(-1.4161371) q[3];
sx q[3];
rz(3.0509907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.693817) q[0];
sx q[0];
rz(-0.87647804) q[0];
sx q[0];
rz(-0.59666657) q[0];
rz(1.5555351) q[1];
sx q[1];
rz(-1.3417473) q[1];
sx q[1];
rz(1.3635051) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1024433) q[0];
sx q[0];
rz(-1.9958695) q[0];
sx q[0];
rz(-0.58702472) q[0];
rz(-0.5559276) q[2];
sx q[2];
rz(-1.5276507) q[2];
sx q[2];
rz(-0.76669979) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.99216695) q[1];
sx q[1];
rz(-0.70632315) q[1];
sx q[1];
rz(-2.7041433) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2927106) q[3];
sx q[3];
rz(-2.1697681) q[3];
sx q[3];
rz(0.90585432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3229225) q[2];
sx q[2];
rz(-1.0007891) q[2];
sx q[2];
rz(0.73454109) q[2];
rz(0.62720403) q[3];
sx q[3];
rz(-1.5137129) q[3];
sx q[3];
rz(2.396615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7221786) q[0];
sx q[0];
rz(-0.83291554) q[0];
sx q[0];
rz(2.869379) q[0];
rz(-0.84725562) q[1];
sx q[1];
rz(-1.3191185) q[1];
sx q[1];
rz(-2.8289657) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49293672) q[0];
sx q[0];
rz(-1.7900294) q[0];
sx q[0];
rz(-1.5159025) q[0];
rz(-pi) q[1];
x q[1];
rz(0.69487822) q[2];
sx q[2];
rz(-2.4701397) q[2];
sx q[2];
rz(1.4050671) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.12007512) q[1];
sx q[1];
rz(-0.26959637) q[1];
sx q[1];
rz(-1.3070379) q[1];
rz(-pi) q[2];
x q[2];
rz(0.76287855) q[3];
sx q[3];
rz(-1.2925914) q[3];
sx q[3];
rz(-2.0335846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0212038) q[2];
sx q[2];
rz(-0.52336064) q[2];
sx q[2];
rz(-2.3266501) q[2];
rz(-2.1728544) q[3];
sx q[3];
rz(-1.0453984) q[3];
sx q[3];
rz(0.43280861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(0.75354904) q[0];
sx q[0];
rz(-2.9409445) q[0];
sx q[0];
rz(-2.5182305) q[0];
rz(-2.3240044) q[1];
sx q[1];
rz(-0.31660429) q[1];
sx q[1];
rz(-0.049291074) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8636297) q[0];
sx q[0];
rz(-0.81252126) q[0];
sx q[0];
rz(-0.066594007) q[0];
rz(-pi) q[1];
rz(0.85929112) q[2];
sx q[2];
rz(-1.8154732) q[2];
sx q[2];
rz(-1.7161075) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4781487) q[1];
sx q[1];
rz(-1.899292) q[1];
sx q[1];
rz(2.8847787) q[1];
x q[2];
rz(0.87818273) q[3];
sx q[3];
rz(-2.4947824) q[3];
sx q[3];
rz(-1.1991024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7307044) q[2];
sx q[2];
rz(-1.5521908) q[2];
sx q[2];
rz(-2.1566186) q[2];
rz(2.2385521) q[3];
sx q[3];
rz(-1.9786381) q[3];
sx q[3];
rz(0.27339098) q[3];
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
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7966998) q[0];
sx q[0];
rz(-1.5364237) q[0];
sx q[0];
rz(-0.087619089) q[0];
rz(-0.15631974) q[1];
sx q[1];
rz(-2.5904398) q[1];
sx q[1];
rz(-0.87096754) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0098457355) q[0];
sx q[0];
rz(-1.5994659) q[0];
sx q[0];
rz(2.8542551) q[0];
x q[1];
rz(-1.2147374) q[2];
sx q[2];
rz(-0.41741727) q[2];
sx q[2];
rz(-2.2603214) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4082143) q[1];
sx q[1];
rz(-0.92792643) q[1];
sx q[1];
rz(-1.9233568) q[1];
x q[2];
rz(0.17467588) q[3];
sx q[3];
rz(-0.81411618) q[3];
sx q[3];
rz(1.0604309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.093420371) q[2];
sx q[2];
rz(-2.293736) q[2];
sx q[2];
rz(2.7887153) q[2];
rz(-2.2757163) q[3];
sx q[3];
rz(-2.4398949) q[3];
sx q[3];
rz(2.849259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46309328) q[0];
sx q[0];
rz(-2.0149639) q[0];
sx q[0];
rz(0.10678664) q[0];
rz(1.1865901) q[1];
sx q[1];
rz(-2.1148966) q[1];
sx q[1];
rz(-1.0245163) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86040076) q[0];
sx q[0];
rz(-1.8587347) q[0];
sx q[0];
rz(0.47173758) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60501955) q[2];
sx q[2];
rz(-1.1875249) q[2];
sx q[2];
rz(-2.5503416) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.66631324) q[1];
sx q[1];
rz(-2.3716455) q[1];
sx q[1];
rz(-2.2837) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6496274) q[3];
sx q[3];
rz(-1.2428189) q[3];
sx q[3];
rz(-0.63719751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7496877) q[2];
sx q[2];
rz(-1.2140032) q[2];
sx q[2];
rz(-0.8708896) q[2];
rz(1.332256) q[3];
sx q[3];
rz(-2.199316) q[3];
sx q[3];
rz(-1.7308621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9695327) q[0];
sx q[0];
rz(-1.7344069) q[0];
sx q[0];
rz(0.55091888) q[0];
rz(-2.6761966) q[1];
sx q[1];
rz(-2.4702563) q[1];
sx q[1];
rz(0.23682061) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4846372) q[0];
sx q[0];
rz(-1.4716363) q[0];
sx q[0];
rz(0.037845503) q[0];
x q[1];
rz(1.4332438) q[2];
sx q[2];
rz(-1.7195065) q[2];
sx q[2];
rz(2.6580722) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3586147) q[1];
sx q[1];
rz(-1.8404984) q[1];
sx q[1];
rz(3.0637834) q[1];
x q[2];
rz(-2.3051466) q[3];
sx q[3];
rz(-1.8006547) q[3];
sx q[3];
rz(-3.0695855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.44234309) q[2];
sx q[2];
rz(-2.4599059) q[2];
sx q[2];
rz(0.44011763) q[2];
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
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(-2.8286164) q[0];
sx q[0];
rz(-1.799311) q[0];
sx q[0];
rz(-0.029504689) q[0];
rz(-0.94738952) q[1];
sx q[1];
rz(-0.82059971) q[1];
sx q[1];
rz(0.11880076) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4122075) q[0];
sx q[0];
rz(-2.8280624) q[0];
sx q[0];
rz(-3.0010812) q[0];
rz(-2.115909) q[2];
sx q[2];
rz(-1.1720177) q[2];
sx q[2];
rz(1.1023956) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1353116) q[1];
sx q[1];
rz(-1.1089919) q[1];
sx q[1];
rz(-1.2414353) q[1];
x q[2];
rz(-1.2715862) q[3];
sx q[3];
rz(-0.51338235) q[3];
sx q[3];
rz(-1.8332924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7416731) q[2];
sx q[2];
rz(-0.46367773) q[2];
sx q[2];
rz(1.5681533) q[2];
rz(-1.9741156) q[3];
sx q[3];
rz(-1.4195331) q[3];
sx q[3];
rz(-1.8673816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90299273) q[0];
sx q[0];
rz(-2.3165343) q[0];
sx q[0];
rz(2.9883244) q[0];
rz(2.080147) q[1];
sx q[1];
rz(-2.5669284) q[1];
sx q[1];
rz(-1.9877888) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.518084) q[0];
sx q[0];
rz(-0.46184807) q[0];
sx q[0];
rz(2.5005546) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2112054) q[2];
sx q[2];
rz(-2.3500035) q[2];
sx q[2];
rz(1.1861578) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.1187763) q[1];
sx q[1];
rz(-0.21636848) q[1];
sx q[1];
rz(-1.8323891) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0204569) q[3];
sx q[3];
rz(-0.88705685) q[3];
sx q[3];
rz(2.5072806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8447421) q[2];
sx q[2];
rz(-1.9694318) q[2];
sx q[2];
rz(1.5273757) q[2];
rz(0.99689233) q[3];
sx q[3];
rz(-1.6485873) q[3];
sx q[3];
rz(-0.87866384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24755724) q[0];
sx q[0];
rz(-2.0443125) q[0];
sx q[0];
rz(-2.9515008) q[0];
rz(0.67063531) q[1];
sx q[1];
rz(-1.9452483) q[1];
sx q[1];
rz(1.6533096) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9207536) q[0];
sx q[0];
rz(-1.5399884) q[0];
sx q[0];
rz(0.12551813) q[0];
rz(-pi) q[1];
rz(-0.52942099) q[2];
sx q[2];
rz(-2.6396857) q[2];
sx q[2];
rz(1.5159964) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.30257672) q[1];
sx q[1];
rz(-0.95514983) q[1];
sx q[1];
rz(-2.2018196) q[1];
rz(-2.2404352) q[3];
sx q[3];
rz(-0.81448758) q[3];
sx q[3];
rz(2.8709473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.0037447475) q[2];
sx q[2];
rz(-0.90128428) q[2];
sx q[2];
rz(2.2424662) q[2];
rz(2.6265465) q[3];
sx q[3];
rz(-0.47596541) q[3];
sx q[3];
rz(0.17764828) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0777733) q[0];
sx q[0];
rz(-1.3615006) q[0];
sx q[0];
rz(2.5831945) q[0];
rz(-0.28941119) q[1];
sx q[1];
rz(-2.2402973) q[1];
sx q[1];
rz(-1.4351861) q[1];
rz(-0.32439705) q[2];
sx q[2];
rz(-2.4051721) q[2];
sx q[2];
rz(0.13798513) q[2];
rz(3.0144721) q[3];
sx q[3];
rz(-1.8416497) q[3];
sx q[3];
rz(1.9491378) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
