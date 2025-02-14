OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0951618) q[0];
sx q[0];
rz(-2.996063) q[0];
sx q[0];
rz(-0.50330436) q[0];
rz(1.1777999) q[1];
sx q[1];
rz(-1.5947394) q[1];
sx q[1];
rz(1.4454747) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4325949) q[0];
sx q[0];
rz(-1.1721969) q[0];
sx q[0];
rz(1.523827) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2022871) q[2];
sx q[2];
rz(-2.21647) q[2];
sx q[2];
rz(-1.1935357) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7912476) q[1];
sx q[1];
rz(-2.4295632) q[1];
sx q[1];
rz(-0.51598926) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.30685356) q[3];
sx q[3];
rz(-3.0899991) q[3];
sx q[3];
rz(-1.2449698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4989) q[2];
sx q[2];
rz(-1.91012) q[2];
sx q[2];
rz(-1.6373681) q[2];
rz(-0.73689342) q[3];
sx q[3];
rz(-1.6156018) q[3];
sx q[3];
rz(-2.1604497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-2.7351643) q[0];
sx q[0];
rz(-0.46877113) q[0];
sx q[0];
rz(2.6306187) q[0];
rz(1.3276395) q[1];
sx q[1];
rz(-1.4013441) q[1];
sx q[1];
rz(2.8151292) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5559306) q[0];
sx q[0];
rz(-1.8343975) q[0];
sx q[0];
rz(-1.8469514) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9316177) q[2];
sx q[2];
rz(-2.2433503) q[2];
sx q[2];
rz(2.5598516) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.49681155) q[1];
sx q[1];
rz(-2.4678951) q[1];
sx q[1];
rz(-2.069887) q[1];
x q[2];
rz(-1.1835008) q[3];
sx q[3];
rz(-1.3069469) q[3];
sx q[3];
rz(1.4572168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.1702801) q[2];
sx q[2];
rz(-0.22394094) q[2];
sx q[2];
rz(1.8453321) q[2];
rz(0.41804677) q[3];
sx q[3];
rz(-2.1635735) q[3];
sx q[3];
rz(2.4372098) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.068785) q[0];
sx q[0];
rz(-2.7588221) q[0];
sx q[0];
rz(-2.5991154) q[0];
rz(-1.3756649) q[1];
sx q[1];
rz(-2.723697) q[1];
sx q[1];
rz(-0.03104041) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.870793) q[0];
sx q[0];
rz(-0.36211553) q[0];
sx q[0];
rz(-2.5848081) q[0];
x q[1];
rz(-1.8607742) q[2];
sx q[2];
rz(-2.463468) q[2];
sx q[2];
rz(-1.4416308) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.52306226) q[1];
sx q[1];
rz(-1.3253322) q[1];
sx q[1];
rz(0.95644585) q[1];
x q[2];
rz(-3.1150041) q[3];
sx q[3];
rz(-0.55517653) q[3];
sx q[3];
rz(2.7257277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5151908) q[2];
sx q[2];
rz(-0.73651892) q[2];
sx q[2];
rz(2.314563) q[2];
rz(0.51860297) q[3];
sx q[3];
rz(-1.1122455) q[3];
sx q[3];
rz(-1.5145068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-2.9943635) q[0];
sx q[0];
rz(-1.3059068) q[0];
sx q[0];
rz(-1.2116785) q[0];
rz(1.0481102) q[1];
sx q[1];
rz(-2.4217889) q[1];
sx q[1];
rz(1.8009708) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3723243) q[0];
sx q[0];
rz(-2.4256673) q[0];
sx q[0];
rz(0.47874079) q[0];
x q[1];
rz(-0.2366613) q[2];
sx q[2];
rz(-1.4259031) q[2];
sx q[2];
rz(1.3587111) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0405492) q[1];
sx q[1];
rz(-1.2770997) q[1];
sx q[1];
rz(0.59567506) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9817438) q[3];
sx q[3];
rz(-1.2839497) q[3];
sx q[3];
rz(2.8887859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2010605) q[2];
sx q[2];
rz(-0.17720711) q[2];
sx q[2];
rz(2.000467) q[2];
rz(1.5308258) q[3];
sx q[3];
rz(-2.174236) q[3];
sx q[3];
rz(-2.358986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.57946) q[0];
sx q[0];
rz(-1.9714404) q[0];
sx q[0];
rz(0.60254565) q[0];
rz(2.5613979) q[1];
sx q[1];
rz(-1.6289026) q[1];
sx q[1];
rz(0.7811195) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21602042) q[0];
sx q[0];
rz(-1.4013774) q[0];
sx q[0];
rz(-1.2792619) q[0];
x q[1];
rz(-0.76483043) q[2];
sx q[2];
rz(-2.0553608) q[2];
sx q[2];
rz(2.7861905) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.32051099) q[1];
sx q[1];
rz(-2.2901504) q[1];
sx q[1];
rz(3.0266552) q[1];
x q[2];
rz(-1.7479728) q[3];
sx q[3];
rz(-2.2886638) q[3];
sx q[3];
rz(-2.126463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.96131229) q[2];
sx q[2];
rz(-1.140241) q[2];
sx q[2];
rz(0.22892924) q[2];
rz(1.3198352) q[3];
sx q[3];
rz(-1.6596551) q[3];
sx q[3];
rz(2.2975217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22885403) q[0];
sx q[0];
rz(-1.5874533) q[0];
sx q[0];
rz(1.9629021) q[0];
rz(-1.2294058) q[1];
sx q[1];
rz(-0.39437672) q[1];
sx q[1];
rz(-1.2961402) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73979171) q[0];
sx q[0];
rz(-1.8515046) q[0];
sx q[0];
rz(0.14174353) q[0];
rz(2.2887694) q[2];
sx q[2];
rz(-2.1666424) q[2];
sx q[2];
rz(-2.9862474) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.069433) q[1];
sx q[1];
rz(-1.0899223) q[1];
sx q[1];
rz(-2.062501) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9518106) q[3];
sx q[3];
rz(-0.49473195) q[3];
sx q[3];
rz(-3.1056946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.487315) q[2];
sx q[2];
rz(-2.2389905) q[2];
sx q[2];
rz(-2.8793092) q[2];
rz(-2.8413963) q[3];
sx q[3];
rz(-1.9191977) q[3];
sx q[3];
rz(-2.9330971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7009785) q[0];
sx q[0];
rz(-0.3955667) q[0];
sx q[0];
rz(1.8390919) q[0];
rz(2.473096) q[1];
sx q[1];
rz(-1.786247) q[1];
sx q[1];
rz(1.0350593) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.221091) q[0];
sx q[0];
rz(-1.971444) q[0];
sx q[0];
rz(-0.23200881) q[0];
x q[1];
rz(1.1167489) q[2];
sx q[2];
rz(-0.40698689) q[2];
sx q[2];
rz(1.5429516) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0332843) q[1];
sx q[1];
rz(-2.7893745) q[1];
sx q[1];
rz(-1.7914821) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1627102) q[3];
sx q[3];
rz(-2.22284) q[3];
sx q[3];
rz(2.121821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.11789007) q[2];
sx q[2];
rz(-2.1390476) q[2];
sx q[2];
rz(-1.4855851) q[2];
rz(2.8591136) q[3];
sx q[3];
rz(-1.3586724) q[3];
sx q[3];
rz(-0.43454596) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78541237) q[0];
sx q[0];
rz(-1.0419351) q[0];
sx q[0];
rz(2.8670512) q[0];
rz(-1.2046332) q[1];
sx q[1];
rz(-0.58012539) q[1];
sx q[1];
rz(0.56142941) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.954079) q[0];
sx q[0];
rz(-0.8773548) q[0];
sx q[0];
rz(2.2827143) q[0];
rz(-0.37158575) q[2];
sx q[2];
rz(-2.4306261) q[2];
sx q[2];
rz(-1.9869876) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.16374396) q[1];
sx q[1];
rz(-1.9789961) q[1];
sx q[1];
rz(-0.331649) q[1];
rz(-pi) q[2];
rz(-0.92989489) q[3];
sx q[3];
rz(-0.42503438) q[3];
sx q[3];
rz(2.2512521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.56208912) q[2];
sx q[2];
rz(-1.1056489) q[2];
sx q[2];
rz(-1.4037464) q[2];
rz(-1.7956519) q[3];
sx q[3];
rz(-0.74508777) q[3];
sx q[3];
rz(-2.6534206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3812934) q[0];
sx q[0];
rz(-2.5385222) q[0];
sx q[0];
rz(2.2316933) q[0];
rz(-1.9445885) q[1];
sx q[1];
rz(-2.0884114) q[1];
sx q[1];
rz(-2.2534456) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5407046) q[0];
sx q[0];
rz(-1.6645169) q[0];
sx q[0];
rz(-1.4613495) q[0];
rz(2.4817564) q[2];
sx q[2];
rz(-1.5232695) q[2];
sx q[2];
rz(-2.8328676) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3067813) q[1];
sx q[1];
rz(-1.1843006) q[1];
sx q[1];
rz(0.67910925) q[1];
rz(-pi) q[2];
rz(1.8664594) q[3];
sx q[3];
rz(-1.5348892) q[3];
sx q[3];
rz(-1.1334238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1278648) q[2];
sx q[2];
rz(-1.700054) q[2];
sx q[2];
rz(2.0590651) q[2];
rz(-0.23076375) q[3];
sx q[3];
rz(-1.328238) q[3];
sx q[3];
rz(0.48399353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8312663) q[0];
sx q[0];
rz(-0.75208298) q[0];
sx q[0];
rz(-2.5600774) q[0];
rz(0.61559081) q[1];
sx q[1];
rz(-2.1938727) q[1];
sx q[1];
rz(-1.8448578) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0727507) q[0];
sx q[0];
rz(-1.5889935) q[0];
sx q[0];
rz(-0.010404603) q[0];
rz(-0.026785568) q[2];
sx q[2];
rz(-2.7973632) q[2];
sx q[2];
rz(2.8452498) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.77380005) q[1];
sx q[1];
rz(-1.8270711) q[1];
sx q[1];
rz(2.0092177) q[1];
rz(-0.6720242) q[3];
sx q[3];
rz(-1.8906396) q[3];
sx q[3];
rz(-0.32957382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7221308) q[2];
sx q[2];
rz(-0.1408793) q[2];
sx q[2];
rz(-2.6922743) q[2];
rz(0.32014534) q[3];
sx q[3];
rz(-1.3195427) q[3];
sx q[3];
rz(1.8951353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37353361) q[0];
sx q[0];
rz(-0.67714416) q[0];
sx q[0];
rz(-1.9376391) q[0];
rz(-1.3846579) q[1];
sx q[1];
rz(-2.90381) q[1];
sx q[1];
rz(-0.79868383) q[1];
rz(0.95296994) q[2];
sx q[2];
rz(-1.7838799) q[2];
sx q[2];
rz(-1.5418216) q[2];
rz(-2.2438335) q[3];
sx q[3];
rz(-2.5644292) q[3];
sx q[3];
rz(0.226365) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
