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
rz(-1.2000097) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8495356) q[0];
sx q[0];
rz(-1.5229051) q[0];
sx q[0];
rz(-0.0056664771) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33485246) q[2];
sx q[2];
rz(-1.9960253) q[2];
sx q[2];
rz(-1.877117) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8438699) q[1];
sx q[1];
rz(-2.5281318) q[1];
sx q[1];
rz(2.2754199) q[1];
rz(-pi) q[2];
rz(-2.7634611) q[3];
sx q[3];
rz(-1.6117125) q[3];
sx q[3];
rz(-3.0331628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8346943) q[2];
sx q[2];
rz(-1.4935741) q[2];
sx q[2];
rz(-2.8519894) q[2];
rz(0.87537193) q[3];
sx q[3];
rz(-2.1353728) q[3];
sx q[3];
rz(-0.059710596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-0.61525476) q[0];
sx q[0];
rz(-0.86741388) q[0];
sx q[0];
rz(0.34399024) q[0];
rz(-3.0572609) q[1];
sx q[1];
rz(-2.4721959) q[1];
sx q[1];
rz(-1.3551691) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25027572) q[0];
sx q[0];
rz(-2.180897) q[0];
sx q[0];
rz(3.0623869) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4488892) q[2];
sx q[2];
rz(-1.7322455) q[2];
sx q[2];
rz(3.0485857) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4210216) q[1];
sx q[1];
rz(-0.45383006) q[1];
sx q[1];
rz(1.4245016) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.077640688) q[3];
sx q[3];
rz(-1.3878126) q[3];
sx q[3];
rz(-1.7971355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8460059) q[2];
sx q[2];
rz(-2.3082374) q[2];
sx q[2];
rz(2.6039092) q[2];
rz(-0.50283557) q[3];
sx q[3];
rz(-0.42481315) q[3];
sx q[3];
rz(-2.1311549) q[3];
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
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
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
rz(-1.0083895) q[1];
sx q[1];
rz(1.0650939) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6247892) q[0];
sx q[0];
rz(-1.5771126) q[0];
sx q[0];
rz(1.9301231) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2421354) q[2];
sx q[2];
rz(-0.81842917) q[2];
sx q[2];
rz(-2.1607272) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0703735) q[1];
sx q[1];
rz(-0.95446903) q[1];
sx q[1];
rz(3.0973869) q[1];
x q[2];
rz(0.25032708) q[3];
sx q[3];
rz(-0.86391376) q[3];
sx q[3];
rz(1.7771429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.76434) q[2];
sx q[2];
rz(-1.179402) q[2];
sx q[2];
rz(0.71737814) q[2];
rz(-2.453089) q[3];
sx q[3];
rz(-2.513956) q[3];
sx q[3];
rz(0.29754105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82729572) q[0];
sx q[0];
rz(-1.2001487) q[0];
sx q[0];
rz(2.8919019) q[0];
rz(2.1266134) q[1];
sx q[1];
rz(-2.8453638) q[1];
sx q[1];
rz(-0.011118523) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4269037) q[0];
sx q[0];
rz(-2.5137797) q[0];
sx q[0];
rz(2.3054302) q[0];
rz(-pi) q[1];
rz(-0.63979062) q[2];
sx q[2];
rz(-0.94546972) q[2];
sx q[2];
rz(1.421979) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1383458) q[1];
sx q[1];
rz(-2.0577288) q[1];
sx q[1];
rz(0.37133118) q[1];
rz(-pi) q[2];
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
rz(-2.6461688) q[2];
sx q[2];
rz(-1.255722) q[2];
sx q[2];
rz(-0.28309506) q[2];
rz(0.66343534) q[3];
sx q[3];
rz(-2.5484271) q[3];
sx q[3];
rz(2.2535113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11113142) q[0];
sx q[0];
rz(-0.24065329) q[0];
sx q[0];
rz(0.33183137) q[0];
rz(-2.6470673) q[1];
sx q[1];
rz(-1.8437513) q[1];
sx q[1];
rz(-1.3269075) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4715695) q[0];
sx q[0];
rz(-1.7243392) q[0];
sx q[0];
rz(-0.75385401) q[0];
rz(-pi) q[1];
rz(-0.7977428) q[2];
sx q[2];
rz(-2.0125055) q[2];
sx q[2];
rz(-1.1898578) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9304282) q[1];
sx q[1];
rz(-1.8976364) q[1];
sx q[1];
rz(-1.8153166) q[1];
rz(-1.9782449) q[3];
sx q[3];
rz(-1.8383998) q[3];
sx q[3];
rz(2.101055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.999324) q[2];
sx q[2];
rz(-0.48196718) q[2];
sx q[2];
rz(-1.2456606) q[2];
rz(-1.8866395) q[3];
sx q[3];
rz(-2.3001223) q[3];
sx q[3];
rz(2.3033223) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6923043) q[0];
sx q[0];
rz(-3.1384387) q[0];
sx q[0];
rz(-0.6814878) q[0];
rz(0.20755126) q[1];
sx q[1];
rz(-0.46336585) q[1];
sx q[1];
rz(-2.025827) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.241248) q[0];
sx q[0];
rz(-1.1660518) q[0];
sx q[0];
rz(-2.3099398) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.240928) q[2];
sx q[2];
rz(-1.9158944) q[2];
sx q[2];
rz(-0.64054856) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.11568497) q[1];
sx q[1];
rz(-0.58870643) q[1];
sx q[1];
rz(-0.031850978) q[1];
rz(-pi) q[2];
x q[2];
rz(0.66949087) q[3];
sx q[3];
rz(-1.7117501) q[3];
sx q[3];
rz(1.6657176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
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
rz(-2.1324683) q[0];
sx q[0];
rz(-0.23129825) q[0];
sx q[0];
rz(-0.67434597) q[0];
rz(1.1122423) q[1];
sx q[1];
rz(-2.4770885) q[1];
sx q[1];
rz(-2.5792714) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4599265) q[0];
sx q[0];
rz(-2.473218) q[0];
sx q[0];
rz(-1.5468803) q[0];
x q[1];
rz(2.855905) q[2];
sx q[2];
rz(-2.7966768) q[2];
sx q[2];
rz(0.069318511) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9971804) q[1];
sx q[1];
rz(-1.4257396) q[1];
sx q[1];
rz(-2.9463065) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4427156) q[3];
sx q[3];
rz(-2.5411385) q[3];
sx q[3];
rz(1.8734991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.101863) q[2];
sx q[2];
rz(-2.1477284) q[2];
sx q[2];
rz(0.34004655) q[2];
rz(-2.9240821) q[3];
sx q[3];
rz(-0.9483996) q[3];
sx q[3];
rz(2.8229009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6927004) q[0];
sx q[0];
rz(-3.0954439) q[0];
sx q[0];
rz(-0.39644077) q[0];
rz(-3.0026644) q[1];
sx q[1];
rz(-0.46008343) q[1];
sx q[1];
rz(-1.6202392) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5312708) q[0];
sx q[0];
rz(-1.6390419) q[0];
sx q[0];
rz(2.9444429) q[0];
rz(-pi) q[1];
rz(-2.2659726) q[2];
sx q[2];
rz(-2.5555829) q[2];
sx q[2];
rz(2.3919174) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.382114) q[1];
sx q[1];
rz(-1.6303051) q[1];
sx q[1];
rz(0.77854034) q[1];
rz(-1.429854) q[3];
sx q[3];
rz(-2.4574276) q[3];
sx q[3];
rz(-1.9977026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.56269318) q[2];
sx q[2];
rz(-2.0843299) q[2];
sx q[2];
rz(2.8472624) q[2];
rz(-2.0108022) q[3];
sx q[3];
rz(-1.7629938) q[3];
sx q[3];
rz(0.99564266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6482553) q[0];
sx q[0];
rz(-0.91518891) q[0];
sx q[0];
rz(2.7822568) q[0];
rz(0.94611478) q[1];
sx q[1];
rz(-0.39603907) q[1];
sx q[1];
rz(-0.27063453) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3753189) q[0];
sx q[0];
rz(-3.0941071) q[0];
sx q[0];
rz(2.4256698) q[0];
rz(-pi) q[1];
rz(-0.46220772) q[2];
sx q[2];
rz(-2.2773829) q[2];
sx q[2];
rz(0.30009899) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4295514) q[1];
sx q[1];
rz(-1.5263285) q[1];
sx q[1];
rz(0.045750381) q[1];
rz(-pi) q[2];
rz(-2.1569096) q[3];
sx q[3];
rz(-2.1911748) q[3];
sx q[3];
rz(1.2445205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.82751194) q[2];
sx q[2];
rz(-0.22958799) q[2];
sx q[2];
rz(-2.6861526) q[2];
rz(-0.81196249) q[3];
sx q[3];
rz(-0.88192314) q[3];
sx q[3];
rz(0.5493831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51046002) q[0];
sx q[0];
rz(-1.6332508) q[0];
sx q[0];
rz(2.4023138) q[0];
rz(0.23070681) q[1];
sx q[1];
rz(-2.6735327) q[1];
sx q[1];
rz(-0.49490067) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5737168) q[0];
sx q[0];
rz(-2.0567237) q[0];
sx q[0];
rz(-3.0987415) q[0];
rz(-pi) q[1];
rz(-1.5779183) q[2];
sx q[2];
rz(-1.2429503) q[2];
sx q[2];
rz(-0.42051007) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.478185) q[1];
sx q[1];
rz(-2.2067398) q[1];
sx q[1];
rz(0.138476) q[1];
rz(0.095330843) q[3];
sx q[3];
rz(-2.2653326) q[3];
sx q[3];
rz(1.8180083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.13359244) q[2];
sx q[2];
rz(-2.0501037) q[2];
sx q[2];
rz(-0.32785329) q[2];
rz(-2.949529) q[3];
sx q[3];
rz(-2.8938507) q[3];
sx q[3];
rz(2.1081934) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0192169) q[0];
sx q[0];
rz(-1.6288971) q[0];
sx q[0];
rz(-1.6678641) q[0];
rz(-2.3090251) q[1];
sx q[1];
rz(-1.7201798) q[1];
sx q[1];
rz(-1.7725772) q[1];
rz(1.238263) q[2];
sx q[2];
rz(-0.56849545) q[2];
sx q[2];
rz(-0.93760437) q[2];
rz(-0.76392382) q[3];
sx q[3];
rz(-0.36119701) q[3];
sx q[3];
rz(1.2154538) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
