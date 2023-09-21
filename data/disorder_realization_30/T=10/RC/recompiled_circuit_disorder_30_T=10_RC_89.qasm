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
rz(2.8770652) q[0];
sx q[0];
rz(9.8192083) q[0];
rz(0.0061622942) q[1];
sx q[1];
rz(-0.34024629) q[1];
sx q[1];
rz(1.9415829) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.278468) q[0];
sx q[0];
rz(-1.5764563) q[0];
sx q[0];
rz(-1.5229043) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0179022) q[2];
sx q[2];
rz(-1.8748218) q[2];
sx q[2];
rz(0.44888874) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.29772273) q[1];
sx q[1];
rz(-0.6134609) q[1];
sx q[1];
rz(2.2754199) q[1];
rz(-2.7634611) q[3];
sx q[3];
rz(-1.5298801) q[3];
sx q[3];
rz(-0.10842987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8346943) q[2];
sx q[2];
rz(-1.4935741) q[2];
sx q[2];
rz(-2.8519894) q[2];
rz(-0.87537193) q[3];
sx q[3];
rz(-2.1353728) q[3];
sx q[3];
rz(0.059710596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61525476) q[0];
sx q[0];
rz(-2.2741788) q[0];
sx q[0];
rz(-2.7976024) q[0];
rz(3.0572609) q[1];
sx q[1];
rz(-2.4721959) q[1];
sx q[1];
rz(1.3551691) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8913169) q[0];
sx q[0];
rz(-2.180897) q[0];
sx q[0];
rz(0.079205714) q[0];
x q[1];
rz(2.5002353) q[2];
sx q[2];
rz(-0.20198447) q[2];
sx q[2];
rz(2.5833677) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.018108798) q[1];
sx q[1];
rz(-1.6347486) q[1];
sx q[1];
rz(2.0204087) q[1];
rz(1.1739028) q[3];
sx q[3];
rz(-2.9429884) q[3];
sx q[3];
rz(-1.393115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.29558674) q[2];
sx q[2];
rz(-0.83335525) q[2];
sx q[2];
rz(0.53768349) q[2];
rz(2.6387571) q[3];
sx q[3];
rz(-2.7167795) q[3];
sx q[3];
rz(-1.0104377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66353345) q[0];
sx q[0];
rz(-2.6038267) q[0];
sx q[0];
rz(0.64087254) q[0];
rz(-0.74869853) q[1];
sx q[1];
rz(-2.1332032) q[1];
sx q[1];
rz(1.0650939) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.070806064) q[0];
sx q[0];
rz(-0.35937989) q[0];
sx q[0];
rz(-1.552836) q[0];
rz(-pi) q[1];
rz(0.89945729) q[2];
sx q[2];
rz(-0.81842917) q[2];
sx q[2];
rz(-2.1607272) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9948431) q[1];
sx q[1];
rz(-2.5238876) q[1];
sx q[1];
rz(1.6330994) q[1];
rz(-0.84824003) q[3];
sx q[3];
rz(-1.7602929) q[3];
sx q[3];
rz(0.041785985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.37725267) q[2];
sx q[2];
rz(-1.179402) q[2];
sx q[2];
rz(-0.71737814) q[2];
rz(-2.453089) q[3];
sx q[3];
rz(-0.62763667) q[3];
sx q[3];
rz(2.8440516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3142969) q[0];
sx q[0];
rz(-1.941444) q[0];
sx q[0];
rz(-0.24969077) q[0];
rz(2.1266134) q[1];
sx q[1];
rz(-2.8453638) q[1];
sx q[1];
rz(-0.011118523) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.714689) q[0];
sx q[0];
rz(-0.62781292) q[0];
sx q[0];
rz(-2.3054302) q[0];
rz(-pi) q[1];
rz(2.3035994) q[2];
sx q[2];
rz(-1.0655155) q[2];
sx q[2];
rz(0.2620286) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3938155) q[1];
sx q[1];
rz(-1.897246) q[1];
sx q[1];
rz(1.0541037) q[1];
x q[2];
rz(-0.53400455) q[3];
sx q[3];
rz(-0.36557331) q[3];
sx q[3];
rz(-0.064985736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.49542385) q[2];
sx q[2];
rz(-1.8858706) q[2];
sx q[2];
rz(-2.8584976) q[2];
rz(0.66343534) q[3];
sx q[3];
rz(-2.5484271) q[3];
sx q[3];
rz(2.2535113) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11113142) q[0];
sx q[0];
rz(-0.24065329) q[0];
sx q[0];
rz(-2.8097613) q[0];
rz(-0.49452531) q[1];
sx q[1];
rz(-1.8437513) q[1];
sx q[1];
rz(1.3269075) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0793314) q[0];
sx q[0];
rz(-2.3752897) q[0];
sx q[0];
rz(-2.9192231) q[0];
x q[1];
rz(-2.1660216) q[2];
sx q[2];
rz(-0.86704463) q[2];
sx q[2];
rz(3.1095568) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2111645) q[1];
sx q[1];
rz(-1.8976364) q[1];
sx q[1];
rz(-1.8153166) q[1];
rz(-2.1760686) q[3];
sx q[3];
rz(-11*pi/13) q[3];
sx q[3];
rz(1.0799288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1422687) q[2];
sx q[2];
rz(-2.6596255) q[2];
sx q[2];
rz(1.2456606) q[2];
rz(1.2549531) q[3];
sx q[3];
rz(-2.3001223) q[3];
sx q[3];
rz(2.3033223) q[3];
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
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44928837) q[0];
sx q[0];
rz(-3.1384387) q[0];
sx q[0];
rz(0.6814878) q[0];
rz(-0.20755126) q[1];
sx q[1];
rz(-2.6782268) q[1];
sx q[1];
rz(1.1157657) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67416699) q[0];
sx q[0];
rz(-0.90303991) q[0];
sx q[0];
rz(-0.52533124) q[0];
rz(1.0461147) q[2];
sx q[2];
rz(-0.74138734) q[2];
sx q[2];
rz(-1.3337097) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7129732) q[1];
sx q[1];
rz(-1.553112) q[1];
sx q[1];
rz(2.5531205) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.22478215) q[3];
sx q[3];
rz(-0.68192476) q[3];
sx q[3];
rz(-0.080760591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9289124) q[2];
sx q[2];
rz(-1.2968411) q[2];
sx q[2];
rz(-0.35432717) q[2];
rz(2.8220693) q[3];
sx q[3];
rz(-1.1525681) q[3];
sx q[3];
rz(0.37187809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0091244) q[0];
sx q[0];
rz(-0.23129825) q[0];
sx q[0];
rz(-2.4672467) q[0];
rz(-2.0293503) q[1];
sx q[1];
rz(-2.4770885) q[1];
sx q[1];
rz(0.56232125) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6816662) q[0];
sx q[0];
rz(-2.473218) q[0];
sx q[0];
rz(1.5468803) q[0];
rz(-pi) q[1];
rz(-1.4698896) q[2];
sx q[2];
rz(-1.9011874) q[2];
sx q[2];
rz(2.908387) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3462785) q[1];
sx q[1];
rz(-0.242713) q[1];
sx q[1];
rz(-0.64530356) q[1];
rz(-pi) q[2];
rz(3.0543442) q[3];
sx q[3];
rz(-2.1656519) q[3];
sx q[3];
rz(2.0283386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0397296) q[2];
sx q[2];
rz(-2.1477284) q[2];
sx q[2];
rz(-0.34004655) q[2];
rz(-2.9240821) q[3];
sx q[3];
rz(-2.1931931) q[3];
sx q[3];
rz(0.31869179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6927004) q[0];
sx q[0];
rz(-0.046148766) q[0];
sx q[0];
rz(-2.7451519) q[0];
rz(3.0026644) q[1];
sx q[1];
rz(-0.46008343) q[1];
sx q[1];
rz(-1.5213535) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5312708) q[0];
sx q[0];
rz(-1.5025508) q[0];
sx q[0];
rz(-2.9444429) q[0];
rz(-pi) q[1];
rz(-2.0422158) q[2];
sx q[2];
rz(-1.9328914) q[2];
sx q[2];
rz(-2.9277756) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8926881) q[1];
sx q[1];
rz(-0.78033328) q[1];
sx q[1];
rz(-3.0569539) q[1];
x q[2];
rz(3.0275214) q[3];
sx q[3];
rz(-2.2469006) q[3];
sx q[3];
rz(1.3249719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5788995) q[2];
sx q[2];
rz(-2.0843299) q[2];
sx q[2];
rz(0.29433027) q[2];
rz(1.1307905) q[3];
sx q[3];
rz(-1.3785988) q[3];
sx q[3];
rz(2.14595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
rz(0.49333736) q[0];
sx q[0];
rz(-0.91518891) q[0];
sx q[0];
rz(-0.35933581) q[0];
rz(-0.94611478) q[1];
sx q[1];
rz(-0.39603907) q[1];
sx q[1];
rz(-2.8709581) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5198869) q[0];
sx q[0];
rz(-1.5396376) q[0];
sx q[0];
rz(3.1057538) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.052202) q[2];
sx q[2];
rz(-0.82197661) q[2];
sx q[2];
rz(-0.35441986) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.86079019) q[1];
sx q[1];
rz(-1.5250912) q[1];
sx q[1];
rz(-1.6153107) q[1];
x q[2];
rz(2.4828033) q[3];
sx q[3];
rz(-2.3156392) q[3];
sx q[3];
rz(-0.39289075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3140807) q[2];
sx q[2];
rz(-2.9120047) q[2];
sx q[2];
rz(-0.45544004) q[2];
rz(0.81196249) q[3];
sx q[3];
rz(-2.2596695) q[3];
sx q[3];
rz(0.5493831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51046002) q[0];
sx q[0];
rz(-1.5083418) q[0];
sx q[0];
rz(2.4023138) q[0];
rz(-0.23070681) q[1];
sx q[1];
rz(-0.46805996) q[1];
sx q[1];
rz(2.646692) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56787581) q[0];
sx q[0];
rz(-2.0567237) q[0];
sx q[0];
rz(0.042851187) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5636744) q[2];
sx q[2];
rz(-1.2429503) q[2];
sx q[2];
rz(2.7210826) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1315688) q[1];
sx q[1];
rz(-1.6820757) q[1];
sx q[1];
rz(2.2113423) q[1];
x q[2];
rz(3.0462618) q[3];
sx q[3];
rz(-2.2653326) q[3];
sx q[3];
rz(-1.8180083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.13359244) q[2];
sx q[2];
rz(-2.0501037) q[2];
sx q[2];
rz(2.8137394) q[2];
rz(2.949529) q[3];
sx q[3];
rz(-2.8938507) q[3];
sx q[3];
rz(-2.1081934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1223758) q[0];
sx q[0];
rz(-1.5126956) q[0];
sx q[0];
rz(1.4737286) q[0];
rz(-2.3090251) q[1];
sx q[1];
rz(-1.7201798) q[1];
sx q[1];
rz(-1.7725772) q[1];
rz(-0.20559786) q[2];
sx q[2];
rz(-2.1046706) q[2];
sx q[2];
rz(-0.54866366) q[2];
rz(-1.8264063) q[3];
sx q[3];
rz(-1.8288463) q[3];
sx q[3];
rz(2.0127206) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];