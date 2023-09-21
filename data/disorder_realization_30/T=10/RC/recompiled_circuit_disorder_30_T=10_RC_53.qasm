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
rz(-0.39443031) q[0];
rz(0.0061622942) q[1];
sx q[1];
rz(-0.34024629) q[1];
sx q[1];
rz(1.9415829) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8495356) q[0];
sx q[0];
rz(-1.5229051) q[0];
sx q[0];
rz(-0.0056664771) q[0];
x q[1];
rz(-1.1236905) q[2];
sx q[2];
rz(-1.8748218) q[2];
sx q[2];
rz(2.6927039) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0389509) q[1];
sx q[1];
rz(-1.1167553) q[1];
sx q[1];
rz(-2.7137043) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6148189) q[3];
sx q[3];
rz(-1.9485954) q[3];
sx q[3];
rz(1.6954741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8346943) q[2];
sx q[2];
rz(-1.6480185) q[2];
sx q[2];
rz(-2.8519894) q[2];
rz(0.87537193) q[3];
sx q[3];
rz(-1.0062199) q[3];
sx q[3];
rz(0.059710596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61525476) q[0];
sx q[0];
rz(-0.86741388) q[0];
sx q[0];
rz(-0.34399024) q[0];
rz(-0.084331766) q[1];
sx q[1];
rz(-0.66939676) q[1];
sx q[1];
rz(-1.3551691) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2750759) q[0];
sx q[0];
rz(-1.5059024) q[0];
sx q[0];
rz(-0.95922031) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6927035) q[2];
sx q[2];
rz(-1.4093471) q[2];
sx q[2];
rz(0.09300692) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.018108798) q[1];
sx q[1];
rz(-1.5068441) q[1];
sx q[1];
rz(-1.121184) q[1];
rz(1.7543206) q[3];
sx q[3];
rz(-1.6471383) q[3];
sx q[3];
rz(-2.9294088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8460059) q[2];
sx q[2];
rz(-0.83335525) q[2];
sx q[2];
rz(2.6039092) q[2];
rz(0.50283557) q[3];
sx q[3];
rz(-0.42481315) q[3];
sx q[3];
rz(2.1311549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4780592) q[0];
sx q[0];
rz(-2.6038267) q[0];
sx q[0];
rz(-0.64087254) q[0];
rz(-2.3928941) q[1];
sx q[1];
rz(-2.1332032) q[1];
sx q[1];
rz(-1.0650939) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6247892) q[0];
sx q[0];
rz(-1.5644801) q[0];
sx q[0];
rz(-1.2114695) q[0];
rz(-2.2674019) q[2];
sx q[2];
rz(-1.0993996) q[2];
sx q[2];
rz(-0.092560571) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6667337) q[1];
sx q[1];
rz(-1.6068646) q[1];
sx q[1];
rz(0.95400793) q[1];
rz(1.8530811) q[3];
sx q[3];
rz(-0.74263393) q[3];
sx q[3];
rz(-1.4020855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.37725267) q[2];
sx q[2];
rz(-1.179402) q[2];
sx q[2];
rz(0.71737814) q[2];
rz(-0.68850368) q[3];
sx q[3];
rz(-0.62763667) q[3];
sx q[3];
rz(-2.8440516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82729572) q[0];
sx q[0];
rz(-1.941444) q[0];
sx q[0];
rz(0.24969077) q[0];
rz(-1.0149792) q[1];
sx q[1];
rz(-2.8453638) q[1];
sx q[1];
rz(3.1304741) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7751559) q[0];
sx q[0];
rz(-1.9754793) q[0];
sx q[0];
rz(1.0767656) q[0];
rz(-pi) q[1];
rz(-2.501802) q[2];
sx q[2];
rz(-0.94546972) q[2];
sx q[2];
rz(1.7196136) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.74777714) q[1];
sx q[1];
rz(-1.897246) q[1];
sx q[1];
rz(2.087489) q[1];
x q[2];
rz(1.3783781) q[3];
sx q[3];
rz(-1.2580066) q[3];
sx q[3];
rz(2.6421412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.49542385) q[2];
sx q[2];
rz(-1.255722) q[2];
sx q[2];
rz(-0.28309506) q[2];
rz(-0.66343534) q[3];
sx q[3];
rz(-0.59316558) q[3];
sx q[3];
rz(2.2535113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11113142) q[0];
sx q[0];
rz(-2.9009394) q[0];
sx q[0];
rz(-2.8097613) q[0];
rz(2.6470673) q[1];
sx q[1];
rz(-1.2978413) q[1];
sx q[1];
rz(-1.3269075) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0793314) q[0];
sx q[0];
rz(-2.3752897) q[0];
sx q[0];
rz(-0.22236951) q[0];
x q[1];
rz(-0.97557108) q[2];
sx q[2];
rz(-0.86704463) q[2];
sx q[2];
rz(0.032035839) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2697619) q[1];
sx q[1];
rz(-0.40553906) q[1];
sx q[1];
rz(-0.62015066) q[1];
rz(-pi) q[2];
rz(2.1760686) q[3];
sx q[3];
rz(-11*pi/13) q[3];
sx q[3];
rz(-1.0799288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1422687) q[2];
sx q[2];
rz(-2.6596255) q[2];
sx q[2];
rz(-1.8959321) q[2];
rz(1.8866395) q[3];
sx q[3];
rz(-0.84147036) q[3];
sx q[3];
rz(2.3033223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6923043) q[0];
sx q[0];
rz(-3.1384387) q[0];
sx q[0];
rz(2.4601049) q[0];
rz(0.20755126) q[1];
sx q[1];
rz(-2.6782268) q[1];
sx q[1];
rz(-1.1157657) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.078243144) q[0];
sx q[0];
rz(-2.3176498) q[0];
sx q[0];
rz(1.0043762) q[0];
rz(-pi) q[1];
rz(1.0461147) q[2];
sx q[2];
rz(-0.74138734) q[2];
sx q[2];
rz(-1.3337097) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4286194) q[1];
sx q[1];
rz(-1.553112) q[1];
sx q[1];
rz(-2.5531205) q[1];
rz(-pi) q[2];
rz(0.66949087) q[3];
sx q[3];
rz(-1.4298425) q[3];
sx q[3];
rz(1.475875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9289124) q[2];
sx q[2];
rz(-1.2968411) q[2];
sx q[2];
rz(-0.35432717) q[2];
rz(0.3195233) q[3];
sx q[3];
rz(-1.1525681) q[3];
sx q[3];
rz(2.7697146) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0091244) q[0];
sx q[0];
rz(-2.9102944) q[0];
sx q[0];
rz(0.67434597) q[0];
rz(1.1122423) q[1];
sx q[1];
rz(-0.66450417) q[1];
sx q[1];
rz(-0.56232125) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0494941) q[0];
sx q[0];
rz(-1.5856165) q[0];
sx q[0];
rz(-0.90256079) q[0];
rz(-1.4698896) q[2];
sx q[2];
rz(-1.9011874) q[2];
sx q[2];
rz(2.908387) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3462785) q[1];
sx q[1];
rz(-0.242713) q[1];
sx q[1];
rz(-2.4962891) q[1];
x q[2];
rz(-0.9741707) q[3];
sx q[3];
rz(-1.4985634) q[3];
sx q[3];
rz(2.7330287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.101863) q[2];
sx q[2];
rz(-0.9938643) q[2];
sx q[2];
rz(0.34004655) q[2];
rz(-0.21751054) q[3];
sx q[3];
rz(-0.9483996) q[3];
sx q[3];
rz(-2.8229009) q[3];
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
rz(-2.6927004) q[0];
sx q[0];
rz(-0.046148766) q[0];
sx q[0];
rz(-2.7451519) q[0];
rz(-0.13892826) q[1];
sx q[1];
rz(-2.6815092) q[1];
sx q[1];
rz(-1.6202392) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0884468) q[0];
sx q[0];
rz(-1.3741115) q[0];
sx q[0];
rz(-1.501207) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0993768) q[2];
sx q[2];
rz(-1.2087012) q[2];
sx q[2];
rz(0.21381703) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.24890451) q[1];
sx q[1];
rz(-0.78033328) q[1];
sx q[1];
rz(0.08463879) q[1];
x q[2];
rz(0.89150724) q[3];
sx q[3];
rz(-1.4818947) q[3];
sx q[3];
rz(2.8241983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.56269318) q[2];
sx q[2];
rz(-2.0843299) q[2];
sx q[2];
rz(0.29433027) q[2];
rz(-2.0108022) q[3];
sx q[3];
rz(-1.3785988) q[3];
sx q[3];
rz(2.14595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6482553) q[0];
sx q[0];
rz(-2.2264037) q[0];
sx q[0];
rz(-0.35933581) q[0];
rz(-0.94611478) q[1];
sx q[1];
rz(-2.7455536) q[1];
sx q[1];
rz(-0.27063453) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5198869) q[0];
sx q[0];
rz(-1.5396376) q[0];
sx q[0];
rz(-3.1057538) q[0];
x q[1];
rz(1.0893906) q[2];
sx q[2];
rz(-0.82197661) q[2];
sx q[2];
rz(2.7871728) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2808025) q[1];
sx q[1];
rz(-1.6165015) q[1];
sx q[1];
rz(1.526282) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.7089013) q[3];
sx q[3];
rz(-1.1039675) q[3];
sx q[3];
rz(2.4469568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.82751194) q[2];
sx q[2];
rz(-0.22958799) q[2];
sx q[2];
rz(-0.45544004) q[2];
rz(-0.81196249) q[3];
sx q[3];
rz(-2.2596695) q[3];
sx q[3];
rz(-0.5493831) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51046002) q[0];
sx q[0];
rz(-1.6332508) q[0];
sx q[0];
rz(0.73927885) q[0];
rz(2.9108858) q[1];
sx q[1];
rz(-0.46805996) q[1];
sx q[1];
rz(-0.49490067) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98289821) q[0];
sx q[0];
rz(-1.532908) q[0];
sx q[0];
rz(1.0844896) q[0];
rz(1.5636744) q[2];
sx q[2];
rz(-1.8986423) q[2];
sx q[2];
rz(-2.7210826) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.478185) q[1];
sx q[1];
rz(-0.93485281) q[1];
sx q[1];
rz(-3.0031167) q[1];
rz(-pi) q[2];
rz(0.87402113) q[3];
sx q[3];
rz(-1.4975944) q[3];
sx q[3];
rz(-0.30833581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1223758) q[0];
sx q[0];
rz(-1.5126956) q[0];
sx q[0];
rz(1.4737286) q[0];
rz(2.3090251) q[1];
sx q[1];
rz(-1.4214129) q[1];
sx q[1];
rz(1.3690154) q[1];
rz(1.238263) q[2];
sx q[2];
rz(-0.56849545) q[2];
sx q[2];
rz(-0.93760437) q[2];
rz(-0.26631793) q[3];
sx q[3];
rz(-1.3238293) q[3];
sx q[3];
rz(-2.7662591) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];