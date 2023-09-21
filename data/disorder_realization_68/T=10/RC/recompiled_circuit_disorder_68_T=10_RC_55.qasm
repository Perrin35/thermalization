OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.39419898) q[0];
sx q[0];
rz(-0.49180254) q[0];
sx q[0];
rz(-2.9536182) q[0];
rz(-1.1176874) q[1];
sx q[1];
rz(-1.517065) q[1];
sx q[1];
rz(-0.36748537) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7029019) q[0];
sx q[0];
rz(-2.5950948) q[0];
sx q[0];
rz(-1.370907) q[0];
rz(-pi) q[1];
rz(-1.3221402) q[2];
sx q[2];
rz(-0.50422943) q[2];
sx q[2];
rz(0.76489514) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.51812664) q[1];
sx q[1];
rz(-1.7382442) q[1];
sx q[1];
rz(-1.8159588) q[1];
x q[2];
rz(1.3931307) q[3];
sx q[3];
rz(-1.2347617) q[3];
sx q[3];
rz(0.48776585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1774896) q[2];
sx q[2];
rz(-0.51012817) q[2];
sx q[2];
rz(0.5509848) q[2];
rz(-1.3059113) q[3];
sx q[3];
rz(-1.6492313) q[3];
sx q[3];
rz(1.3163542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47857639) q[0];
sx q[0];
rz(-1.0000279) q[0];
sx q[0];
rz(0.4719032) q[0];
rz(-2.7117803) q[1];
sx q[1];
rz(-1.2496354) q[1];
sx q[1];
rz(2.205251) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4124356) q[0];
sx q[0];
rz(-1.658174) q[0];
sx q[0];
rz(0.23075128) q[0];
x q[1];
rz(-2.3564732) q[2];
sx q[2];
rz(-1.8765419) q[2];
sx q[2];
rz(2.608125) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1658926) q[1];
sx q[1];
rz(-0.70043889) q[1];
sx q[1];
rz(-2.9794934) q[1];
rz(-1.6132658) q[3];
sx q[3];
rz(-0.50308933) q[3];
sx q[3];
rz(-0.79723061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3669746) q[2];
sx q[2];
rz(-0.32704157) q[2];
sx q[2];
rz(-2.7152087) q[2];
rz(1.2373699) q[3];
sx q[3];
rz(-0.62785134) q[3];
sx q[3];
rz(3.1085076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8957829) q[0];
sx q[0];
rz(-1.3170467) q[0];
sx q[0];
rz(-0.93908969) q[0];
rz(-2.242873) q[1];
sx q[1];
rz(-0.4788613) q[1];
sx q[1];
rz(-2.5476707) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31375162) q[0];
sx q[0];
rz(-1.272164) q[0];
sx q[0];
rz(1.9065501) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0837768) q[2];
sx q[2];
rz(-2.2584492) q[2];
sx q[2];
rz(-0.85580326) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9582639) q[1];
sx q[1];
rz(-1.3576344) q[1];
sx q[1];
rz(-1.1413241) q[1];
x q[2];
rz(1.2265615) q[3];
sx q[3];
rz(-2.0153181) q[3];
sx q[3];
rz(1.3821186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5014191) q[2];
sx q[2];
rz(-2.4121425) q[2];
sx q[2];
rz(1.7017986) q[2];
rz(0.38763186) q[3];
sx q[3];
rz(-1.5165611) q[3];
sx q[3];
rz(-2.7534527) q[3];
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
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3751635) q[0];
sx q[0];
rz(-1.5513993) q[0];
sx q[0];
rz(2.638812) q[0];
rz(0.76820961) q[1];
sx q[1];
rz(-0.50351024) q[1];
sx q[1];
rz(0.75685135) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2198974) q[0];
sx q[0];
rz(-1.2384402) q[0];
sx q[0];
rz(-0.38145782) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1050622) q[2];
sx q[2];
rz(-1.5591991) q[2];
sx q[2];
rz(0.84601814) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6824324) q[1];
sx q[1];
rz(-1.2130514) q[1];
sx q[1];
rz(0.95174241) q[1];
rz(-pi) q[2];
rz(2.0714949) q[3];
sx q[3];
rz(-0.90161937) q[3];
sx q[3];
rz(0.6176399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.42671529) q[2];
sx q[2];
rz(-1.2299512) q[2];
sx q[2];
rz(-1.4871917) q[2];
rz(0.58250827) q[3];
sx q[3];
rz(-2.0472066) q[3];
sx q[3];
rz(2.5845161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29397598) q[0];
sx q[0];
rz(-1.0852381) q[0];
sx q[0];
rz(2.3838682) q[0];
rz(-1.2879397) q[1];
sx q[1];
rz(-0.92823354) q[1];
sx q[1];
rz(1.0505189) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8331063) q[0];
sx q[0];
rz(-1.8252488) q[0];
sx q[0];
rz(1.2480877) q[0];
x q[1];
rz(0.97425766) q[2];
sx q[2];
rz(-1.9447118) q[2];
sx q[2];
rz(-0.24335441) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.037462385) q[1];
sx q[1];
rz(-2.3298652) q[1];
sx q[1];
rz(-0.8123668) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3349808) q[3];
sx q[3];
rz(-1.1685373) q[3];
sx q[3];
rz(2.3276687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.918255) q[2];
sx q[2];
rz(-0.35327062) q[2];
sx q[2];
rz(2.5081432) q[2];
rz(1.194362) q[3];
sx q[3];
rz(-1.6703689) q[3];
sx q[3];
rz(-0.71715322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.441992) q[0];
sx q[0];
rz(-0.49015912) q[0];
sx q[0];
rz(-2.8884086) q[0];
rz(1.5340012) q[1];
sx q[1];
rz(-1.4350767) q[1];
sx q[1];
rz(1.6794499) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40186858) q[0];
sx q[0];
rz(-1.7812294) q[0];
sx q[0];
rz(-1.6519288) q[0];
rz(1.2573104) q[2];
sx q[2];
rz(-0.67337155) q[2];
sx q[2];
rz(-1.9999258) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.05789214) q[1];
sx q[1];
rz(-1.9533227) q[1];
sx q[1];
rz(1.3871357) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8263894) q[3];
sx q[3];
rz(-1.8990371) q[3];
sx q[3];
rz(-2.4344276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.2016466) q[2];
sx q[2];
rz(-0.74735171) q[2];
sx q[2];
rz(0.80491006) q[2];
rz(1.1770052) q[3];
sx q[3];
rz(-1.0369119) q[3];
sx q[3];
rz(-1.0837519) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8686304) q[0];
sx q[0];
rz(-2.0721764) q[0];
sx q[0];
rz(-0.92765635) q[0];
rz(-1.0246798) q[1];
sx q[1];
rz(-1.506348) q[1];
sx q[1];
rz(-1.0120846) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1597848) q[0];
sx q[0];
rz(-0.73043981) q[0];
sx q[0];
rz(-2.3083789) q[0];
rz(0.32832844) q[2];
sx q[2];
rz(-0.40847455) q[2];
sx q[2];
rz(0.35818737) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.95841366) q[1];
sx q[1];
rz(-1.603754) q[1];
sx q[1];
rz(0.54380137) q[1];
rz(-3.0104962) q[3];
sx q[3];
rz(-0.56250611) q[3];
sx q[3];
rz(1.2870115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.53081375) q[2];
sx q[2];
rz(-1.6616219) q[2];
sx q[2];
rz(0.42993316) q[2];
rz(-2.1271465) q[3];
sx q[3];
rz(-0.40922624) q[3];
sx q[3];
rz(2.608192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6124509) q[0];
sx q[0];
rz(-0.93944678) q[0];
sx q[0];
rz(-2.9274143) q[0];
rz(2.0902436) q[1];
sx q[1];
rz(-0.21251692) q[1];
sx q[1];
rz(2.8578551) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8026233) q[0];
sx q[0];
rz(-1.5366652) q[0];
sx q[0];
rz(-2.8404833) q[0];
x q[1];
rz(1.0546513) q[2];
sx q[2];
rz(-0.37545855) q[2];
sx q[2];
rz(-1.6106538) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.022790837) q[1];
sx q[1];
rz(-0.53794251) q[1];
sx q[1];
rz(-2.7522037) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.489336) q[3];
sx q[3];
rz(-2.1144923) q[3];
sx q[3];
rz(-0.039747681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6909137) q[2];
sx q[2];
rz(-0.51270715) q[2];
sx q[2];
rz(1.696375) q[2];
rz(-1.5971659) q[3];
sx q[3];
rz(-1.701136) q[3];
sx q[3];
rz(-0.33932313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63672367) q[0];
sx q[0];
rz(-3.0915785) q[0];
sx q[0];
rz(-0.069256393) q[0];
rz(-1.4878558) q[1];
sx q[1];
rz(-1.2830877) q[1];
sx q[1];
rz(1.5690631) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3587787) q[0];
sx q[0];
rz(-2.1920715) q[0];
sx q[0];
rz(0.57399477) q[0];
rz(-pi) q[1];
rz(-1.4885159) q[2];
sx q[2];
rz(-1.2459718) q[2];
sx q[2];
rz(0.61330739) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.931817) q[1];
sx q[1];
rz(-1.6546384) q[1];
sx q[1];
rz(-0.2340338) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8435555) q[3];
sx q[3];
rz(-2.4747304) q[3];
sx q[3];
rz(-2.42815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1853603) q[2];
sx q[2];
rz(-0.23024836) q[2];
sx q[2];
rz(-3.0017079) q[2];
rz(0.36758962) q[3];
sx q[3];
rz(-1.9544173) q[3];
sx q[3];
rz(2.1504413) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96520987) q[0];
sx q[0];
rz(-2.7503224) q[0];
sx q[0];
rz(-2.4998253) q[0];
rz(1.9104674) q[1];
sx q[1];
rz(-1.1522013) q[1];
sx q[1];
rz(-2.8737601) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56775996) q[0];
sx q[0];
rz(-1.1288252) q[0];
sx q[0];
rz(2.3964336) q[0];
x q[1];
rz(-2.1688813) q[2];
sx q[2];
rz(-2.9209666) q[2];
sx q[2];
rz(1.0704744) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1416727) q[1];
sx q[1];
rz(-1.281764) q[1];
sx q[1];
rz(-0.33860597) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0249428) q[3];
sx q[3];
rz(-0.59026679) q[3];
sx q[3];
rz(0.86436194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3315167) q[2];
sx q[2];
rz(-1.5248359) q[2];
sx q[2];
rz(-1.6646741) q[2];
rz(-0.26633513) q[3];
sx q[3];
rz(-0.24644066) q[3];
sx q[3];
rz(-0.5464856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1289566) q[0];
sx q[0];
rz(-2.2205882) q[0];
sx q[0];
rz(0.90482774) q[0];
rz(-0.77990445) q[1];
sx q[1];
rz(-0.48702469) q[1];
sx q[1];
rz(-1.3866966) q[1];
rz(-0.98942479) q[2];
sx q[2];
rz(-2.1944254) q[2];
sx q[2];
rz(2.2133322) q[2];
rz(1.7846617) q[3];
sx q[3];
rz(-2.2371117) q[3];
sx q[3];
rz(-0.37123751) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
