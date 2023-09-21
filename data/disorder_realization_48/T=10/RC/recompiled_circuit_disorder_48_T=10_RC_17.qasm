OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7735908) q[0];
sx q[0];
rz(-2.350783) q[0];
sx q[0];
rz(2.8074582) q[0];
rz(-0.45733991) q[1];
sx q[1];
rz(5.338905) q[1];
sx q[1];
rz(10.64325) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89361184) q[0];
sx q[0];
rz(-0.75151822) q[0];
sx q[0];
rz(0.052608842) q[0];
rz(-pi) q[1];
rz(-3.0604612) q[2];
sx q[2];
rz(-0.46617026) q[2];
sx q[2];
rz(2.1612338) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4301181) q[1];
sx q[1];
rz(-2.4415486) q[1];
sx q[1];
rz(-0.75573604) q[1];
rz(-pi) q[2];
rz(-0.18913194) q[3];
sx q[3];
rz(-1.7889708) q[3];
sx q[3];
rz(1.3355586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5228287) q[2];
sx q[2];
rz(-0.4814119) q[2];
sx q[2];
rz(-0.5775601) q[2];
rz(1.9918359) q[3];
sx q[3];
rz(-1.7532319) q[3];
sx q[3];
rz(-2.4770588) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98786551) q[0];
sx q[0];
rz(-0.58652121) q[0];
sx q[0];
rz(0.38744774) q[0];
rz(0.93915835) q[1];
sx q[1];
rz(-0.99717957) q[1];
sx q[1];
rz(-1.739025) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70799202) q[0];
sx q[0];
rz(-3.1118244) q[0];
sx q[0];
rz(-0.13694163) q[0];
rz(1.7895133) q[2];
sx q[2];
rz(-1.0220851) q[2];
sx q[2];
rz(2.1825841) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2991997) q[1];
sx q[1];
rz(-1.0789011) q[1];
sx q[1];
rz(-2.030034) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4757231) q[3];
sx q[3];
rz(-2.2007696) q[3];
sx q[3];
rz(2.7271987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.42276057) q[2];
sx q[2];
rz(-1.7451124) q[2];
sx q[2];
rz(-2.823901) q[2];
rz(-2.9348532) q[3];
sx q[3];
rz(-0.59967774) q[3];
sx q[3];
rz(0.81682214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3690255) q[0];
sx q[0];
rz(-1.7018397) q[0];
sx q[0];
rz(1.4136219) q[0];
rz(-0.47779045) q[1];
sx q[1];
rz(-1.7910035) q[1];
sx q[1];
rz(0.40107045) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0553592) q[0];
sx q[0];
rz(-1.8952574) q[0];
sx q[0];
rz(3.0326891) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7602073) q[2];
sx q[2];
rz(-2.2064798) q[2];
sx q[2];
rz(2.1798101) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2420826) q[1];
sx q[1];
rz(-0.84065719) q[1];
sx q[1];
rz(1.5793369) q[1];
rz(2.9455455) q[3];
sx q[3];
rz(-1.184433) q[3];
sx q[3];
rz(-0.26323174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5473189) q[2];
sx q[2];
rz(-1.5218364) q[2];
sx q[2];
rz(2.5857914) q[2];
rz(2.1650971) q[3];
sx q[3];
rz(-0.54978168) q[3];
sx q[3];
rz(0.78021375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7926086) q[0];
sx q[0];
rz(-1.5190834) q[0];
sx q[0];
rz(1.4439616) q[0];
rz(1.5199039) q[1];
sx q[1];
rz(-2.4855721) q[1];
sx q[1];
rz(-0.25340432) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69617535) q[0];
sx q[0];
rz(-1.4679969) q[0];
sx q[0];
rz(2.1131383) q[0];
x q[1];
rz(1.5379982) q[2];
sx q[2];
rz(-0.48600733) q[2];
sx q[2];
rz(-0.82644586) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6771486) q[1];
sx q[1];
rz(-2.1975937) q[1];
sx q[1];
rz(1.6897175) q[1];
x q[2];
rz(0.13417379) q[3];
sx q[3];
rz(-0.65376702) q[3];
sx q[3];
rz(1.3282446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.110934) q[2];
sx q[2];
rz(-1.3867644) q[2];
sx q[2];
rz(0.78732642) q[2];
rz(-0.91286719) q[3];
sx q[3];
rz(-2.392231) q[3];
sx q[3];
rz(-2.1319938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.426429) q[0];
sx q[0];
rz(-0.63868317) q[0];
sx q[0];
rz(-0.062967904) q[0];
rz(3.0175623) q[1];
sx q[1];
rz(-2.3359559) q[1];
sx q[1];
rz(-0.45809349) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8100909) q[0];
sx q[0];
rz(-1.1281663) q[0];
sx q[0];
rz(1.5682194) q[0];
rz(-pi) q[1];
rz(-0.45256726) q[2];
sx q[2];
rz(-2.1308225) q[2];
sx q[2];
rz(2.3134311) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0702857) q[1];
sx q[1];
rz(-2.3100393) q[1];
sx q[1];
rz(2.9830095) q[1];
x q[2];
rz(2.0714508) q[3];
sx q[3];
rz(-1.6541012) q[3];
sx q[3];
rz(1.3798151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1725585) q[2];
sx q[2];
rz(-2.2183552) q[2];
sx q[2];
rz(-2.9120973) q[2];
rz(0.0028006639) q[3];
sx q[3];
rz(-2.2715748) q[3];
sx q[3];
rz(1.3389448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5620419) q[0];
sx q[0];
rz(-0.27844772) q[0];
sx q[0];
rz(0.91947412) q[0];
rz(-3.0793076) q[1];
sx q[1];
rz(-1.0039763) q[1];
sx q[1];
rz(-1.2671635) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4290659) q[0];
sx q[0];
rz(-2.292284) q[0];
sx q[0];
rz(0.52961403) q[0];
rz(2.9847449) q[2];
sx q[2];
rz(-2.3914797) q[2];
sx q[2];
rz(-1.7240806) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2515182) q[1];
sx q[1];
rz(-1.8971844) q[1];
sx q[1];
rz(-2.9337487) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1250161) q[3];
sx q[3];
rz(-2.6187754) q[3];
sx q[3];
rz(2.0857874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.9617812) q[2];
sx q[2];
rz(-0.87819019) q[2];
sx q[2];
rz(-0.58376694) q[2];
rz(2.4328655) q[3];
sx q[3];
rz(-1.3112336) q[3];
sx q[3];
rz(3.1183929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0631183) q[0];
sx q[0];
rz(-0.45409504) q[0];
sx q[0];
rz(-2.069058) q[0];
rz(-2.5947) q[1];
sx q[1];
rz(-1.2439589) q[1];
sx q[1];
rz(1.1118719) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3355013) q[0];
sx q[0];
rz(-1.6643235) q[0];
sx q[0];
rz(1.085824) q[0];
rz(-pi) q[1];
rz(-0.56356168) q[2];
sx q[2];
rz(-2.5346018) q[2];
sx q[2];
rz(0.8286455) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0791196) q[1];
sx q[1];
rz(-1.8603431) q[1];
sx q[1];
rz(1.6505961) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0734378) q[3];
sx q[3];
rz(-1.916269) q[3];
sx q[3];
rz(0.84514602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.303858) q[2];
sx q[2];
rz(-1.4759109) q[2];
sx q[2];
rz(1.973935) q[2];
rz(-1.5363103) q[3];
sx q[3];
rz(-1.4669908) q[3];
sx q[3];
rz(1.4060098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7325571) q[0];
sx q[0];
rz(-1.5696101) q[0];
sx q[0];
rz(0.73079601) q[0];
rz(-2.2413975) q[1];
sx q[1];
rz(-0.80454818) q[1];
sx q[1];
rz(-2.3866167) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2728111) q[0];
sx q[0];
rz(-2.4065354) q[0];
sx q[0];
rz(-1.3520157) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5078012) q[2];
sx q[2];
rz(-0.83665028) q[2];
sx q[2];
rz(1.447669) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4432115) q[1];
sx q[1];
rz(-1.4581919) q[1];
sx q[1];
rz(0.57971445) q[1];
rz(-2.2488392) q[3];
sx q[3];
rz(-2.4340981) q[3];
sx q[3];
rz(-2.1219818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3770611) q[2];
sx q[2];
rz(-1.3724644) q[2];
sx q[2];
rz(-0.63684741) q[2];
rz(2.8751255) q[3];
sx q[3];
rz(-1.0576495) q[3];
sx q[3];
rz(-1.586097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4144142) q[0];
sx q[0];
rz(-2.0120912) q[0];
sx q[0];
rz(-1.138858) q[0];
rz(0.75421929) q[1];
sx q[1];
rz(-2.8051839) q[1];
sx q[1];
rz(3.1220904) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50197983) q[0];
sx q[0];
rz(-1.6507971) q[0];
sx q[0];
rz(-1.7777068) q[0];
rz(2.2150061) q[2];
sx q[2];
rz(-2.8928061) q[2];
sx q[2];
rz(0.34005806) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.86917415) q[1];
sx q[1];
rz(-0.65009102) q[1];
sx q[1];
rz(-1.8723349) q[1];
rz(-pi) q[2];
rz(-1.8832302) q[3];
sx q[3];
rz(-1.4133246) q[3];
sx q[3];
rz(1.0851932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.0043682178) q[2];
sx q[2];
rz(-1.7251816) q[2];
sx q[2];
rz(0.84890378) q[2];
rz(-0.38765872) q[3];
sx q[3];
rz(-2.013423) q[3];
sx q[3];
rz(1.5415812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3432817) q[0];
sx q[0];
rz(-0.16769519) q[0];
sx q[0];
rz(-2.6570901) q[0];
rz(-1.3867406) q[1];
sx q[1];
rz(-1.7157028) q[1];
sx q[1];
rz(-1.1482931) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5507817) q[0];
sx q[0];
rz(-1.592144) q[0];
sx q[0];
rz(3.000893) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8414367) q[2];
sx q[2];
rz(-0.65320063) q[2];
sx q[2];
rz(2.4094827) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.53459586) q[1];
sx q[1];
rz(-1.3998704) q[1];
sx q[1];
rz(-0.89764885) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.34928068) q[3];
sx q[3];
rz(-1.5216856) q[3];
sx q[3];
rz(1.9625361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2293573) q[2];
sx q[2];
rz(-1.2939913) q[2];
sx q[2];
rz(2.7764376) q[2];
rz(0.12864104) q[3];
sx q[3];
rz(-1.235685) q[3];
sx q[3];
rz(-2.685759) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0317595) q[0];
sx q[0];
rz(-0.84072996) q[0];
sx q[0];
rz(1.6054556) q[0];
rz(2.1784492) q[1];
sx q[1];
rz(-1.2711202) q[1];
sx q[1];
rz(-1.0585379) q[1];
rz(1.9671494) q[2];
sx q[2];
rz(-0.94559961) q[2];
sx q[2];
rz(-0.0018975817) q[2];
rz(-1.5118602) q[3];
sx q[3];
rz(-2.5806576) q[3];
sx q[3];
rz(2.9604119) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
