OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.49139872) q[0];
sx q[0];
rz(-0.2645275) q[0];
sx q[0];
rz(2.7471623) q[0];
rz(-3.1354304) q[1];
sx q[1];
rz(-2.8013464) q[1];
sx q[1];
rz(-1.9415829) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8631247) q[0];
sx q[0];
rz(-1.5764563) q[0];
sx q[0];
rz(1.6186884) q[0];
rz(-pi) q[1];
rz(-0.33485246) q[2];
sx q[2];
rz(-1.9960253) q[2];
sx q[2];
rz(-1.877117) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.29772273) q[1];
sx q[1];
rz(-0.6134609) q[1];
sx q[1];
rz(0.86617275) q[1];
x q[2];
rz(1.5267738) q[3];
sx q[3];
rz(-1.1929973) q[3];
sx q[3];
rz(-1.6954741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3068984) q[2];
sx q[2];
rz(-1.4935741) q[2];
sx q[2];
rz(-2.8519894) q[2];
rz(-0.87537193) q[3];
sx q[3];
rz(-1.0062199) q[3];
sx q[3];
rz(3.0818821) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5263379) q[0];
sx q[0];
rz(-0.86741388) q[0];
sx q[0];
rz(0.34399024) q[0];
rz(3.0572609) q[1];
sx q[1];
rz(-2.4721959) q[1];
sx q[1];
rz(1.3551691) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2750759) q[0];
sx q[0];
rz(-1.6356902) q[0];
sx q[0];
rz(-2.1823723) q[0];
rz(0.64135735) q[2];
sx q[2];
rz(-0.20198447) q[2];
sx q[2];
rz(-2.5833677) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1234839) q[1];
sx q[1];
rz(-1.5068441) q[1];
sx q[1];
rz(-2.0204087) q[1];
x q[2];
rz(-1.387272) q[3];
sx q[3];
rz(-1.4944544) q[3];
sx q[3];
rz(2.9294088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8460059) q[2];
sx q[2];
rz(-0.83335525) q[2];
sx q[2];
rz(0.53768349) q[2];
rz(-0.50283557) q[3];
sx q[3];
rz(-0.42481315) q[3];
sx q[3];
rz(1.0104377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4780592) q[0];
sx q[0];
rz(-0.53776598) q[0];
sx q[0];
rz(-0.64087254) q[0];
rz(-0.74869853) q[1];
sx q[1];
rz(-1.0083895) q[1];
sx q[1];
rz(-1.0650939) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0707866) q[0];
sx q[0];
rz(-2.7822128) q[0];
sx q[0];
rz(-1.5887567) q[0];
x q[1];
rz(-0.89945729) q[2];
sx q[2];
rz(-2.3231635) q[2];
sx q[2];
rz(-2.1607272) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6667337) q[1];
sx q[1];
rz(-1.534728) q[1];
sx q[1];
rz(-2.1875847) q[1];
rz(0.25032708) q[3];
sx q[3];
rz(-2.2776789) q[3];
sx q[3];
rz(1.3644497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.37725267) q[2];
sx q[2];
rz(-1.9621907) q[2];
sx q[2];
rz(-0.71737814) q[2];
rz(2.453089) q[3];
sx q[3];
rz(-2.513956) q[3];
sx q[3];
rz(-0.29754105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82729572) q[0];
sx q[0];
rz(-1.2001487) q[0];
sx q[0];
rz(-2.8919019) q[0];
rz(1.0149792) q[1];
sx q[1];
rz(-2.8453638) q[1];
sx q[1];
rz(0.011118523) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4269037) q[0];
sx q[0];
rz(-0.62781292) q[0];
sx q[0];
rz(0.83616242) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3035994) q[2];
sx q[2];
rz(-2.0760771) q[2];
sx q[2];
rz(-2.879564) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.74777714) q[1];
sx q[1];
rz(-1.2443466) q[1];
sx q[1];
rz(2.087489) q[1];
rz(0.31828493) q[3];
sx q[3];
rz(-1.3878229) q[3];
sx q[3];
rz(-2.0103679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.49542385) q[2];
sx q[2];
rz(-1.8858706) q[2];
sx q[2];
rz(0.28309506) q[2];
rz(-0.66343534) q[3];
sx q[3];
rz(-2.5484271) q[3];
sx q[3];
rz(0.88808131) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11113142) q[0];
sx q[0];
rz(-0.24065329) q[0];
sx q[0];
rz(-0.33183137) q[0];
rz(0.49452531) q[1];
sx q[1];
rz(-1.2978413) q[1];
sx q[1];
rz(-1.8146851) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0622612) q[0];
sx q[0];
rz(-2.3752897) q[0];
sx q[0];
rz(2.9192231) q[0];
rz(-0.58381501) q[2];
sx q[2];
rz(-0.8875672) q[2];
sx q[2];
rz(-0.77606397) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2697619) q[1];
sx q[1];
rz(-2.7360536) q[1];
sx q[1];
rz(-2.521442) q[1];
x q[2];
rz(-2.1760686) q[3];
sx q[3];
rz(-11*pi/13) q[3];
sx q[3];
rz(1.0799288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1422687) q[2];
sx q[2];
rz(-0.48196718) q[2];
sx q[2];
rz(-1.2456606) q[2];
rz(1.2549531) q[3];
sx q[3];
rz(-0.84147036) q[3];
sx q[3];
rz(0.83827034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6923043) q[0];
sx q[0];
rz(-0.0031539991) q[0];
sx q[0];
rz(-2.4601049) q[0];
rz(0.20755126) q[1];
sx q[1];
rz(-0.46336585) q[1];
sx q[1];
rz(-2.025827) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.241248) q[0];
sx q[0];
rz(-1.1660518) q[0];
sx q[0];
rz(-0.83165283) q[0];
x q[1];
rz(0.90066465) q[2];
sx q[2];
rz(-1.9158944) q[2];
sx q[2];
rz(2.5010441) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.15397729) q[1];
sx q[1];
rz(-2.1591641) q[1];
sx q[1];
rz(-1.5920559) q[1];
x q[2];
rz(0.22478215) q[3];
sx q[3];
rz(-2.4596679) q[3];
sx q[3];
rz(-0.080760591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.21268022) q[2];
sx q[2];
rz(-1.2968411) q[2];
sx q[2];
rz(0.35432717) q[2];
rz(-0.3195233) q[3];
sx q[3];
rz(-1.1525681) q[3];
sx q[3];
rz(0.37187809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-1.0091244) q[0];
sx q[0];
rz(-2.9102944) q[0];
sx q[0];
rz(-0.67434597) q[0];
rz(1.1122423) q[1];
sx q[1];
rz(-2.4770885) q[1];
sx q[1];
rz(-2.5792714) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6816662) q[0];
sx q[0];
rz(-2.473218) q[0];
sx q[0];
rz(1.5947123) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.855905) q[2];
sx q[2];
rz(-2.7966768) q[2];
sx q[2];
rz(-0.069318511) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6866236) q[1];
sx q[1];
rz(-1.7640055) q[1];
sx q[1];
rz(-1.7186233) q[1];
rz(-pi) q[2];
rz(0.9741707) q[3];
sx q[3];
rz(-1.6430292) q[3];
sx q[3];
rz(2.7330287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0397296) q[2];
sx q[2];
rz(-0.9938643) q[2];
sx q[2];
rz(-2.8015461) q[2];
rz(-0.21751054) q[3];
sx q[3];
rz(-0.9483996) q[3];
sx q[3];
rz(0.31869179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44889221) q[0];
sx q[0];
rz(-3.0954439) q[0];
sx q[0];
rz(-2.7451519) q[0];
rz(-3.0026644) q[1];
sx q[1];
rz(-2.6815092) q[1];
sx q[1];
rz(1.6202392) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28946653) q[0];
sx q[0];
rz(-0.20848256) q[0];
sx q[0];
rz(-2.8058488) q[0];
x q[1];
rz(-2.2659726) q[2];
sx q[2];
rz(-0.5860098) q[2];
sx q[2];
rz(-2.3919174) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7594787) q[1];
sx q[1];
rz(-1.5112875) q[1];
sx q[1];
rz(-2.3630523) q[1];
x q[2];
rz(2.2500854) q[3];
sx q[3];
rz(-1.4818947) q[3];
sx q[3];
rz(0.31739435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5788995) q[2];
sx q[2];
rz(-1.0572628) q[2];
sx q[2];
rz(2.8472624) q[2];
rz(-1.1307905) q[3];
sx q[3];
rz(-1.3785988) q[3];
sx q[3];
rz(-2.14595) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49333736) q[0];
sx q[0];
rz(-0.91518891) q[0];
sx q[0];
rz(-2.7822568) q[0];
rz(2.1954779) q[1];
sx q[1];
rz(-2.7455536) q[1];
sx q[1];
rz(-0.27063453) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6217057) q[0];
sx q[0];
rz(-1.6019551) q[0];
sx q[0];
rz(-3.1057538) q[0];
rz(-pi) q[1];
rz(2.6793849) q[2];
sx q[2];
rz(-0.86420977) q[2];
sx q[2];
rz(2.8414937) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.088085727) q[1];
sx q[1];
rz(-3.0778031) q[1];
sx q[1];
rz(-2.3699058) q[1];
x q[2];
rz(-0.6587894) q[3];
sx q[3];
rz(-2.3156392) q[3];
sx q[3];
rz(-0.39289075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.82751194) q[2];
sx q[2];
rz(-2.9120047) q[2];
sx q[2];
rz(-2.6861526) q[2];
rz(0.81196249) q[3];
sx q[3];
rz(-2.2596695) q[3];
sx q[3];
rz(-2.5922095) q[3];
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
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51046002) q[0];
sx q[0];
rz(-1.6332508) q[0];
sx q[0];
rz(-2.4023138) q[0];
rz(-0.23070681) q[1];
sx q[1];
rz(-0.46805996) q[1];
sx q[1];
rz(2.646692) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5737168) q[0];
sx q[0];
rz(-2.0567237) q[0];
sx q[0];
rz(-0.042851187) q[0];
rz(-1.5636744) q[2];
sx q[2];
rz(-1.2429503) q[2];
sx q[2];
rz(0.42051007) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.010023821) q[1];
sx q[1];
rz(-1.6820757) q[1];
sx q[1];
rz(2.2113423) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.095330843) q[3];
sx q[3];
rz(-2.2653326) q[3];
sx q[3];
rz(-1.8180083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0080002) q[2];
sx q[2];
rz(-2.0501037) q[2];
sx q[2];
rz(0.32785329) q[2];
rz(2.949529) q[3];
sx q[3];
rz(-2.8938507) q[3];
sx q[3];
rz(1.0333992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0192169) q[0];
sx q[0];
rz(-1.5126956) q[0];
sx q[0];
rz(1.4737286) q[0];
rz(-2.3090251) q[1];
sx q[1];
rz(-1.7201798) q[1];
sx q[1];
rz(-1.7725772) q[1];
rz(2.1140425) q[2];
sx q[2];
rz(-1.7474568) q[2];
sx q[2];
rz(-2.2251868) q[2];
rz(-2.8752747) q[3];
sx q[3];
rz(-1.8177633) q[3];
sx q[3];
rz(0.37533356) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
