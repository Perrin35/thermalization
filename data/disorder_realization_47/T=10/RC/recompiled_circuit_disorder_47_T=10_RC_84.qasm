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
rz(2.3072825) q[0];
sx q[0];
rz(6.351525) q[0];
rz(-0.66032687) q[1];
sx q[1];
rz(-0.84815174) q[1];
sx q[1];
rz(0.037820427) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2743933) q[0];
sx q[0];
rz(-1.4925856) q[0];
sx q[0];
rz(-1.377051) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1184138) q[2];
sx q[2];
rz(-0.9127494) q[2];
sx q[2];
rz(2.7409035) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6112069) q[1];
sx q[1];
rz(-1.5198277) q[1];
sx q[1];
rz(2.3947122) q[1];
x q[2];
rz(1.9170403) q[3];
sx q[3];
rz(-1.4284705) q[3];
sx q[3];
rz(-0.12366731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.35090703) q[2];
sx q[2];
rz(-0.51351341) q[2];
sx q[2];
rz(-1.2834056) q[2];
rz(3.0154199) q[3];
sx q[3];
rz(-1.7254555) q[3];
sx q[3];
rz(-3.0509907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44777563) q[0];
sx q[0];
rz(-0.87647804) q[0];
sx q[0];
rz(0.59666657) q[0];
rz(-1.5555351) q[1];
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
rz(2.4564132) q[0];
rz(0.081625799) q[2];
sx q[2];
rz(-2.584169) q[2];
sx q[2];
rz(0.73478414) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9047782) q[1];
sx q[1];
rz(-1.292255) q[1];
sx q[1];
rz(2.4836471) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3724243) q[3];
sx q[3];
rz(-2.2393919) q[3];
sx q[3];
rz(1.9070966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3229225) q[2];
sx q[2];
rz(-2.1408036) q[2];
sx q[2];
rz(2.4070516) q[2];
rz(-2.5143886) q[3];
sx q[3];
rz(-1.6278798) q[3];
sx q[3];
rz(0.74497765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7221786) q[0];
sx q[0];
rz(-2.3086771) q[0];
sx q[0];
rz(2.869379) q[0];
rz(2.294337) q[1];
sx q[1];
rz(-1.8224742) q[1];
sx q[1];
rz(-0.31262696) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4011742) q[0];
sx q[0];
rz(-2.9156988) q[0];
sx q[0];
rz(-2.9001539) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1001415) q[2];
sx q[2];
rz(-1.0725642) q[2];
sx q[2];
rz(2.2217896) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.15309139) q[1];
sx q[1];
rz(-1.3107436) q[1];
sx q[1];
rz(0.071916332) q[1];
x q[2];
rz(-2.3787141) q[3];
sx q[3];
rz(-1.2925914) q[3];
sx q[3];
rz(1.1080081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0212038) q[2];
sx q[2];
rz(-2.618232) q[2];
sx q[2];
rz(-0.81494251) q[2];
rz(-2.1728544) q[3];
sx q[3];
rz(-1.0453984) q[3];
sx q[3];
rz(-2.708784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75354904) q[0];
sx q[0];
rz(-2.9409445) q[0];
sx q[0];
rz(0.62336212) q[0];
rz(-0.81758824) q[1];
sx q[1];
rz(-0.31660429) q[1];
sx q[1];
rz(0.049291074) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.277963) q[0];
sx q[0];
rz(-2.3290714) q[0];
sx q[0];
rz(-0.066594007) q[0];
rz(-pi) q[1];
rz(-0.31844278) q[2];
sx q[2];
rz(-0.88469425) q[2];
sx q[2];
rz(2.7903914) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.019776736) q[1];
sx q[1];
rz(-0.41408086) q[1];
sx q[1];
rz(2.2112234) q[1];
rz(-pi) q[2];
rz(-0.87818273) q[3];
sx q[3];
rz(-0.64681029) q[3];
sx q[3];
rz(1.9424903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4108882) q[2];
sx q[2];
rz(-1.5521908) q[2];
sx q[2];
rz(-0.98497406) q[2];
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
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34489283) q[0];
sx q[0];
rz(-1.5364237) q[0];
sx q[0];
rz(0.087619089) q[0];
rz(-2.9852729) q[1];
sx q[1];
rz(-0.55115288) q[1];
sx q[1];
rz(2.2706251) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5694215) q[0];
sx q[0];
rz(-1.8580125) q[0];
sx q[0];
rz(1.5409018) q[0];
x q[1];
rz(0.15337495) q[2];
sx q[2];
rz(-1.9605637) q[2];
sx q[2];
rz(-2.6467121) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4082143) q[1];
sx q[1];
rz(-2.2136662) q[1];
sx q[1];
rz(1.9233568) q[1];
x q[2];
rz(-1.7528275) q[3];
sx q[3];
rz(-2.3689299) q[3];
sx q[3];
rz(2.3327737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0481723) q[2];
sx q[2];
rz(-0.8478567) q[2];
sx q[2];
rz(2.7887153) q[2];
rz(2.2757163) q[3];
sx q[3];
rz(-2.4398949) q[3];
sx q[3];
rz(0.29233366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6784994) q[0];
sx q[0];
rz(-2.0149639) q[0];
sx q[0];
rz(-3.034806) q[0];
rz(1.1865901) q[1];
sx q[1];
rz(-2.1148966) q[1];
sx q[1];
rz(2.1170763) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56652503) q[0];
sx q[0];
rz(-2.0216414) q[0];
sx q[0];
rz(-1.8917811) q[0];
x q[1];
rz(1.1149939) q[2];
sx q[2];
rz(-1.015128) q[2];
sx q[2];
rz(1.9089886) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.66631324) q[1];
sx q[1];
rz(-0.76994714) q[1];
sx q[1];
rz(2.2837) q[1];
rz(-pi) q[2];
rz(-0.6242574) q[3];
sx q[3];
rz(-0.58371021) q[3];
sx q[3];
rz(0.39238413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.391905) q[2];
sx q[2];
rz(-1.9275894) q[2];
sx q[2];
rz(2.270703) q[2];
rz(-1.8093367) q[3];
sx q[3];
rz(-0.94227666) q[3];
sx q[3];
rz(1.7308621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1720599) q[0];
sx q[0];
rz(-1.4071858) q[0];
sx q[0];
rz(2.5906738) q[0];
rz(0.46539601) q[1];
sx q[1];
rz(-2.4702563) q[1];
sx q[1];
rz(-2.904772) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.119334) q[0];
sx q[0];
rz(-3.0354781) q[0];
sx q[0];
rz(-1.9342213) q[0];
rz(-pi) q[1];
rz(-0.15010712) q[2];
sx q[2];
rz(-1.4347715) q[2];
sx q[2];
rz(-2.0748236) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6432861) q[1];
sx q[1];
rz(-0.28043881) q[1];
sx q[1];
rz(1.8449057) q[1];
x q[2];
rz(-2.3051466) q[3];
sx q[3];
rz(-1.8006547) q[3];
sx q[3];
rz(0.072007192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6992496) q[2];
sx q[2];
rz(-2.4599059) q[2];
sx q[2];
rz(-2.701475) q[2];
rz(1.0951428) q[3];
sx q[3];
rz(-1.4705855) q[3];
sx q[3];
rz(1.7949036) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8286164) q[0];
sx q[0];
rz(-1.799311) q[0];
sx q[0];
rz(3.112088) q[0];
rz(2.1942031) q[1];
sx q[1];
rz(-0.82059971) q[1];
sx q[1];
rz(-3.0227919) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4122075) q[0];
sx q[0];
rz(-0.31353024) q[0];
sx q[0];
rz(-0.14051147) q[0];
rz(-pi) q[1];
rz(-2.115909) q[2];
sx q[2];
rz(-1.9695749) q[2];
sx q[2];
rz(-1.1023956) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.48092948) q[1];
sx q[1];
rz(-2.5813563) q[1];
sx q[1];
rz(2.5653097) q[1];
x q[2];
rz(-1.0766255) q[3];
sx q[3];
rz(-1.716074) q[3];
sx q[3];
rz(0.52500099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7416731) q[2];
sx q[2];
rz(-2.6779149) q[2];
sx q[2];
rz(1.5734394) q[2];
rz(-1.167477) q[3];
sx q[3];
rz(-1.4195331) q[3];
sx q[3];
rz(-1.2742111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2385999) q[0];
sx q[0];
rz(-2.3165343) q[0];
sx q[0];
rz(2.9883244) q[0];
rz(-1.0614456) q[1];
sx q[1];
rz(-0.57466424) q[1];
sx q[1];
rz(-1.1538039) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3183347) q[0];
sx q[0];
rz(-1.9360006) q[0];
sx q[0];
rz(-1.2814786) q[0];
rz(-2.597528) q[2];
sx q[2];
rz(-2.1779034) q[2];
sx q[2];
rz(1.1407167) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8037387) q[1];
sx q[1];
rz(-1.5152463) q[1];
sx q[1];
rz(-1.3615723) q[1];
rz(2.2581402) q[3];
sx q[3];
rz(-1.4769819) q[3];
sx q[3];
rz(1.0132307) q[3];
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
rz(-1.4930054) q[3];
sx q[3];
rz(0.87866384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8940354) q[0];
sx q[0];
rz(-1.0972801) q[0];
sx q[0];
rz(-2.9515008) q[0];
rz(0.67063531) q[1];
sx q[1];
rz(-1.9452483) q[1];
sx q[1];
rz(1.6533096) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7877486) q[0];
sx q[0];
rz(-1.6962546) q[0];
sx q[0];
rz(1.6018484) q[0];
x q[1];
rz(0.44234862) q[2];
sx q[2];
rz(-1.3253691) q[2];
sx q[2];
rz(-0.52877141) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8390159) q[1];
sx q[1];
rz(-0.95514983) q[1];
sx q[1];
rz(2.2018196) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.87741239) q[3];
sx q[3];
rz(-1.1023695) q[3];
sx q[3];
rz(-1.7978158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1378479) q[2];
sx q[2];
rz(-2.2403084) q[2];
sx q[2];
rz(0.89912644) q[2];
rz(0.51504618) q[3];
sx q[3];
rz(-0.47596541) q[3];
sx q[3];
rz(2.9639444) q[3];
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
rz(0.32439705) q[2];
sx q[2];
rz(-0.73642052) q[2];
sx q[2];
rz(-3.0036075) q[2];
rz(-1.1424941) q[3];
sx q[3];
rz(-2.8430568) q[3];
sx q[3];
rz(-1.6381016) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
