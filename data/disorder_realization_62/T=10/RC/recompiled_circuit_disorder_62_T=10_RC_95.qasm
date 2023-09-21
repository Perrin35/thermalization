OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2405038) q[0];
sx q[0];
rz(-2.9641889) q[0];
sx q[0];
rz(2.0071964) q[0];
rz(1.1881243) q[1];
sx q[1];
rz(-2.1048574) q[1];
sx q[1];
rz(-0.66361767) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21270574) q[0];
sx q[0];
rz(-2.0318673) q[0];
sx q[0];
rz(1.5972208) q[0];
x q[1];
rz(2.7825836) q[2];
sx q[2];
rz(-0.81072545) q[2];
sx q[2];
rz(-2.5100978) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9259778) q[1];
sx q[1];
rz(-1.3091334) q[1];
sx q[1];
rz(1.91933) q[1];
rz(-pi) q[2];
rz(0.10856467) q[3];
sx q[3];
rz(-1.3285471) q[3];
sx q[3];
rz(0.066075174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2628281) q[2];
sx q[2];
rz(-0.44439134) q[2];
sx q[2];
rz(0.051068548) q[2];
rz(0.55705327) q[3];
sx q[3];
rz(-2.3414108) q[3];
sx q[3];
rz(-1.5548271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5550845) q[0];
sx q[0];
rz(-2.3095135) q[0];
sx q[0];
rz(-0.59659514) q[0];
rz(2.3157628) q[1];
sx q[1];
rz(-1.4412216) q[1];
sx q[1];
rz(-1.2260431) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3292424) q[0];
sx q[0];
rz(-1.7698405) q[0];
sx q[0];
rz(-2.6890486) q[0];
rz(-pi) q[1];
rz(0.82310279) q[2];
sx q[2];
rz(-0.96379333) q[2];
sx q[2];
rz(0.58128231) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5495758) q[1];
sx q[1];
rz(-1.5712275) q[1];
sx q[1];
rz(-1.7906584) q[1];
x q[2];
rz(1.5561043) q[3];
sx q[3];
rz(-0.5726632) q[3];
sx q[3];
rz(0.0088012561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.79919672) q[2];
sx q[2];
rz(-1.1725972) q[2];
sx q[2];
rz(-0.33102316) q[2];
rz(2.3349169) q[3];
sx q[3];
rz(-1.9162063) q[3];
sx q[3];
rz(1.4276918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1598635) q[0];
sx q[0];
rz(-1.230343) q[0];
sx q[0];
rz(1.249041) q[0];
rz(-0.088009134) q[1];
sx q[1];
rz(-2.0188589) q[1];
sx q[1];
rz(-1.0294611) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2217076) q[0];
sx q[0];
rz(-1.789933) q[0];
sx q[0];
rz(-2.1973781) q[0];
rz(-pi) q[1];
x q[1];
rz(0.37000334) q[2];
sx q[2];
rz(-0.68379935) q[2];
sx q[2];
rz(1.5497269) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0440812) q[1];
sx q[1];
rz(-1.4315726) q[1];
sx q[1];
rz(-0.15686762) q[1];
rz(-pi) q[2];
rz(2.401628) q[3];
sx q[3];
rz(-1.4810586) q[3];
sx q[3];
rz(1.2147853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.133698) q[2];
sx q[2];
rz(-1.7241314) q[2];
sx q[2];
rz(-0.50764817) q[2];
rz(-1.7525904) q[3];
sx q[3];
rz(-0.32998431) q[3];
sx q[3];
rz(1.9784137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65748173) q[0];
sx q[0];
rz(-2.8634475) q[0];
sx q[0];
rz(1.5456276) q[0];
rz(1.0428628) q[1];
sx q[1];
rz(-1.1735801) q[1];
sx q[1];
rz(-1.5159336) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31947485) q[0];
sx q[0];
rz(-1.7135156) q[0];
sx q[0];
rz(2.4583754) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.085725351) q[2];
sx q[2];
rz(-2.1282196) q[2];
sx q[2];
rz(2.2103708) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2559291) q[1];
sx q[1];
rz(-0.71223488) q[1];
sx q[1];
rz(0.42339143) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5528615) q[3];
sx q[3];
rz(-1.3874467) q[3];
sx q[3];
rz(-2.6453032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3778014) q[2];
sx q[2];
rz(-0.8447454) q[2];
sx q[2];
rz(1.2949004) q[2];
rz(-3.0002248) q[3];
sx q[3];
rz(-0.54261345) q[3];
sx q[3];
rz(-0.98658371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.7871053) q[0];
sx q[0];
rz(-1.961345) q[0];
sx q[0];
rz(1.5198583) q[0];
rz(-2.5095818) q[1];
sx q[1];
rz(-2.4193587) q[1];
sx q[1];
rz(2.246726) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97840727) q[0];
sx q[0];
rz(-2.1911231) q[0];
sx q[0];
rz(2.083367) q[0];
rz(-0.72563719) q[2];
sx q[2];
rz(-1.1302395) q[2];
sx q[2];
rz(-0.63647905) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.85256165) q[1];
sx q[1];
rz(-0.82419306) q[1];
sx q[1];
rz(2.186071) q[1];
x q[2];
rz(1.5630866) q[3];
sx q[3];
rz(-1.2843411) q[3];
sx q[3];
rz(-0.25659284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.616509) q[2];
sx q[2];
rz(-0.36642763) q[2];
sx q[2];
rz(-1.9990702) q[2];
rz(0.74674314) q[3];
sx q[3];
rz(-1.2544422) q[3];
sx q[3];
rz(-2.0660627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58364761) q[0];
sx q[0];
rz(-2.8201411) q[0];
sx q[0];
rz(1.7640132) q[0];
rz(-2.6691943) q[1];
sx q[1];
rz(-0.51858416) q[1];
sx q[1];
rz(-0.46498743) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0156292) q[0];
sx q[0];
rz(-1.1383801) q[0];
sx q[0];
rz(0.36884357) q[0];
rz(-pi) q[1];
rz(-1.497252) q[2];
sx q[2];
rz(-1.2300756) q[2];
sx q[2];
rz(-1.3982915) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2496693) q[1];
sx q[1];
rz(-1.6185456) q[1];
sx q[1];
rz(2.2062917) q[1];
rz(-pi) q[2];
x q[2];
rz(0.70907866) q[3];
sx q[3];
rz(-1.4104098) q[3];
sx q[3];
rz(1.3080314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.650699) q[2];
sx q[2];
rz(-1.2630254) q[2];
sx q[2];
rz(-2.6965551) q[2];
rz(-2.2079091) q[3];
sx q[3];
rz(-1.4327587) q[3];
sx q[3];
rz(2.8745108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12748195) q[0];
sx q[0];
rz(-1.5662136) q[0];
sx q[0];
rz(1.7215464) q[0];
rz(-3.1177915) q[1];
sx q[1];
rz(-0.61444608) q[1];
sx q[1];
rz(-0.15596095) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0359081) q[0];
sx q[0];
rz(-1.3996291) q[0];
sx q[0];
rz(-0.25733421) q[0];
rz(-pi) q[1];
rz(-1.9975125) q[2];
sx q[2];
rz(-1.0973755) q[2];
sx q[2];
rz(-2.0063426) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5676463) q[1];
sx q[1];
rz(-1.624755) q[1];
sx q[1];
rz(1.7458564) q[1];
x q[2];
rz(2.6050623) q[3];
sx q[3];
rz(-2.0091669) q[3];
sx q[3];
rz(2.8560672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0722787) q[2];
sx q[2];
rz(-2.5645655) q[2];
sx q[2];
rz(-1.0726661) q[2];
rz(2.8159451) q[3];
sx q[3];
rz(-1.9653392) q[3];
sx q[3];
rz(-1.6252888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9496562) q[0];
sx q[0];
rz(-2.7392116) q[0];
sx q[0];
rz(-2.8038213) q[0];
rz(-1.0900963) q[1];
sx q[1];
rz(-0.97507674) q[1];
sx q[1];
rz(2.8930194) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0524806) q[0];
sx q[0];
rz(-2.5626474) q[0];
sx q[0];
rz(-0.46778932) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0609264) q[2];
sx q[2];
rz(-1.9976127) q[2];
sx q[2];
rz(0.61444297) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.12411815) q[1];
sx q[1];
rz(-1.0867456) q[1];
sx q[1];
rz(0.89213051) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0646938) q[3];
sx q[3];
rz(-1.4014763) q[3];
sx q[3];
rz(2.3016735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8481855) q[2];
sx q[2];
rz(-0.53331393) q[2];
sx q[2];
rz(0.08671134) q[2];
rz(-2.6596206) q[3];
sx q[3];
rz(-0.50589365) q[3];
sx q[3];
rz(1.988525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(0.50865737) q[0];
sx q[0];
rz(-1.0077227) q[0];
sx q[0];
rz(0.28276643) q[0];
rz(-0.70156082) q[1];
sx q[1];
rz(-2.3200254) q[1];
sx q[1];
rz(1.823002) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2254927) q[0];
sx q[0];
rz(-2.6927104) q[0];
sx q[0];
rz(-1.7569957) q[0];
rz(-0.65638541) q[2];
sx q[2];
rz(-2.0771386) q[2];
sx q[2];
rz(-0.93751794) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2576037) q[1];
sx q[1];
rz(-2.249243) q[1];
sx q[1];
rz(0.18737327) q[1];
rz(-0.75026476) q[3];
sx q[3];
rz(-1.8599469) q[3];
sx q[3];
rz(0.17444785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41137722) q[2];
sx q[2];
rz(-2.8432379) q[2];
sx q[2];
rz(-0.39917699) q[2];
rz(-2.2579851) q[3];
sx q[3];
rz(-1.6688321) q[3];
sx q[3];
rz(-1.2214899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3783962) q[0];
sx q[0];
rz(-2.3662687) q[0];
sx q[0];
rz(-1.5989074) q[0];
rz(2.0762766) q[1];
sx q[1];
rz(-2.1677446) q[1];
sx q[1];
rz(1.261196) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6529918) q[0];
sx q[0];
rz(-1.4654667) q[0];
sx q[0];
rz(2.7306042) q[0];
rz(2.5433259) q[2];
sx q[2];
rz(-1.5891799) q[2];
sx q[2];
rz(2.604904) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.8775455) q[1];
sx q[1];
rz(-1.0730768) q[1];
sx q[1];
rz(1.4977786) q[1];
rz(-pi) q[2];
rz(-1.4114755) q[3];
sx q[3];
rz(-0.89309249) q[3];
sx q[3];
rz(0.89980984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.5579055) q[2];
sx q[2];
rz(-0.62290278) q[2];
sx q[2];
rz(-2.5718001) q[2];
rz(1.9231046) q[3];
sx q[3];
rz(-0.64703882) q[3];
sx q[3];
rz(-0.56263721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31310836) q[0];
sx q[0];
rz(-0.88043558) q[0];
sx q[0];
rz(1.2783929) q[0];
rz(-2.5333511) q[1];
sx q[1];
rz(-2.6651762) q[1];
sx q[1];
rz(-2.6574635) q[1];
rz(0.026253168) q[2];
sx q[2];
rz(-2.142341) q[2];
sx q[2];
rz(-0.061948902) q[2];
rz(2.9363587) q[3];
sx q[3];
rz(-0.60094613) q[3];
sx q[3];
rz(-0.080106674) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];