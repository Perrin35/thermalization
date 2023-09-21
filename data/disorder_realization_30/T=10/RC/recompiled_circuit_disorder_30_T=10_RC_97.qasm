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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8495356) q[0];
sx q[0];
rz(-1.5229051) q[0];
sx q[0];
rz(0.0056664771) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33485246) q[2];
sx q[2];
rz(-1.1455673) q[2];
sx q[2];
rz(-1.2644757) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.29772273) q[1];
sx q[1];
rz(-0.6134609) q[1];
sx q[1];
rz(-2.2754199) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6148189) q[3];
sx q[3];
rz(-1.1929973) q[3];
sx q[3];
rz(-1.4461185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3068984) q[2];
sx q[2];
rz(-1.4935741) q[2];
sx q[2];
rz(-2.8519894) q[2];
rz(-2.2662207) q[3];
sx q[3];
rz(-1.0062199) q[3];
sx q[3];
rz(0.059710596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61525476) q[0];
sx q[0];
rz(-2.2741788) q[0];
sx q[0];
rz(2.7976024) q[0];
rz(-3.0572609) q[1];
sx q[1];
rz(-2.4721959) q[1];
sx q[1];
rz(-1.3551691) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25027572) q[0];
sx q[0];
rz(-2.180897) q[0];
sx q[0];
rz(-3.0623869) q[0];
rz(-pi) q[1];
x q[1];
rz(0.16263527) q[2];
sx q[2];
rz(-1.4504823) q[2];
sx q[2];
rz(-1.6441117) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1234839) q[1];
sx q[1];
rz(-1.6347486) q[1];
sx q[1];
rz(2.0204087) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.387272) q[3];
sx q[3];
rz(-1.6471383) q[3];
sx q[3];
rz(0.21218382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8460059) q[2];
sx q[2];
rz(-0.83335525) q[2];
sx q[2];
rz(2.6039092) q[2];
rz(0.50283557) q[3];
sx q[3];
rz(-2.7167795) q[3];
sx q[3];
rz(-2.1311549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66353345) q[0];
sx q[0];
rz(-2.6038267) q[0];
sx q[0];
rz(-2.5007201) q[0];
rz(-2.3928941) q[1];
sx q[1];
rz(-2.1332032) q[1];
sx q[1];
rz(2.0764988) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-2.2674019) q[2];
sx q[2];
rz(-1.0993996) q[2];
sx q[2];
rz(3.0490321) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9948431) q[1];
sx q[1];
rz(-0.61770505) q[1];
sx q[1];
rz(1.5084933) q[1];
rz(-pi) q[2];
rz(0.84824003) q[3];
sx q[3];
rz(-1.7602929) q[3];
sx q[3];
rz(3.0998067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.76434) q[2];
sx q[2];
rz(-1.9621907) q[2];
sx q[2];
rz(-0.71737814) q[2];
rz(-0.68850368) q[3];
sx q[3];
rz(-2.513956) q[3];
sx q[3];
rz(2.8440516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3142969) q[0];
sx q[0];
rz(-1.2001487) q[0];
sx q[0];
rz(-0.24969077) q[0];
rz(2.1266134) q[1];
sx q[1];
rz(-0.29622886) q[1];
sx q[1];
rz(0.011118523) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.714689) q[0];
sx q[0];
rz(-0.62781292) q[0];
sx q[0];
rz(2.3054302) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2617433) q[2];
sx q[2];
rz(-2.2789311) q[2];
sx q[2];
rz(0.81530064) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8324082) q[1];
sx q[1];
rz(-0.60317457) q[1];
sx q[1];
rz(-2.1716154) q[1];
rz(-pi) q[2];
rz(-2.8233077) q[3];
sx q[3];
rz(-1.7537698) q[3];
sx q[3];
rz(2.0103679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6461688) q[2];
sx q[2];
rz(-1.8858706) q[2];
sx q[2];
rz(0.28309506) q[2];
rz(-2.4781573) q[3];
sx q[3];
rz(-0.59316558) q[3];
sx q[3];
rz(-2.2535113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(3.0304612) q[0];
sx q[0];
rz(-2.9009394) q[0];
sx q[0];
rz(2.8097613) q[0];
rz(2.6470673) q[1];
sx q[1];
rz(-1.2978413) q[1];
sx q[1];
rz(1.8146851) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0793314) q[0];
sx q[0];
rz(-2.3752897) q[0];
sx q[0];
rz(0.22236951) q[0];
rz(-0.7977428) q[2];
sx q[2];
rz(-1.1290871) q[2];
sx q[2];
rz(-1.9517348) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8618968) q[1];
sx q[1];
rz(-1.8021291) q[1];
sx q[1];
rz(2.8054603) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1760686) q[3];
sx q[3];
rz(-11*pi/13) q[3];
sx q[3];
rz(-1.0799288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.999324) q[2];
sx q[2];
rz(-2.6596255) q[2];
sx q[2];
rz(-1.2456606) q[2];
rz(1.8866395) q[3];
sx q[3];
rz(-0.84147036) q[3];
sx q[3];
rz(2.3033223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-0.46336585) q[1];
sx q[1];
rz(-1.1157657) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0633495) q[0];
sx q[0];
rz(-2.3176498) q[0];
sx q[0];
rz(-1.0043762) q[0];
rz(-pi) q[1];
rz(0.43004604) q[2];
sx q[2];
rz(-2.1950245) q[2];
sx q[2];
rz(-2.4732694) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.15397729) q[1];
sx q[1];
rz(-0.9824285) q[1];
sx q[1];
rz(1.5920559) q[1];
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
sx q[1];
rz(-pi/2) q[1];
rz(-0.21268022) q[2];
sx q[2];
rz(-1.8447515) q[2];
sx q[2];
rz(-0.35432717) q[2];
rz(-0.3195233) q[3];
sx q[3];
rz(-1.9890246) q[3];
sx q[3];
rz(2.7697146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1324683) q[0];
sx q[0];
rz(-2.9102944) q[0];
sx q[0];
rz(2.4672467) q[0];
rz(2.0293503) q[1];
sx q[1];
rz(-2.4770885) q[1];
sx q[1];
rz(-0.56232125) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4903957) q[0];
sx q[0];
rz(-2.2389452) q[0];
sx q[0];
rz(0.018880318) q[0];
x q[1];
rz(-1.6717031) q[2];
sx q[2];
rz(-1.9011874) q[2];
sx q[2];
rz(-2.908387) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.14441227) q[1];
sx q[1];
rz(-1.715853) q[1];
sx q[1];
rz(-0.1952862) q[1];
rz(-pi) q[2];
rz(-1.4427156) q[3];
sx q[3];
rz(-2.5411385) q[3];
sx q[3];
rz(-1.2680935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.101863) q[2];
sx q[2];
rz(-0.9938643) q[2];
sx q[2];
rz(2.8015461) q[2];
rz(0.21751054) q[3];
sx q[3];
rz(-2.1931931) q[3];
sx q[3];
rz(-2.8229009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6927004) q[0];
sx q[0];
rz(-3.0954439) q[0];
sx q[0];
rz(-0.39644077) q[0];
rz(-0.13892826) q[1];
sx q[1];
rz(-2.6815092) q[1];
sx q[1];
rz(1.5213535) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6103219) q[0];
sx q[0];
rz(-1.5025508) q[0];
sx q[0];
rz(-2.9444429) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0422158) q[2];
sx q[2];
rz(-1.2087012) q[2];
sx q[2];
rz(-2.9277756) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8926881) q[1];
sx q[1];
rz(-2.3612594) q[1];
sx q[1];
rz(-0.08463879) q[1];
rz(-pi) q[2];
rz(-0.89150724) q[3];
sx q[3];
rz(-1.4818947) q[3];
sx q[3];
rz(0.31739435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5788995) q[2];
sx q[2];
rz(-1.0572628) q[2];
sx q[2];
rz(2.8472624) q[2];
rz(2.0108022) q[3];
sx q[3];
rz(-1.7629938) q[3];
sx q[3];
rz(-0.99564266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6482553) q[0];
sx q[0];
rz(-0.91518891) q[0];
sx q[0];
rz(-2.7822568) q[0];
rz(-0.94611478) q[1];
sx q[1];
rz(-0.39603907) q[1];
sx q[1];
rz(-2.8709581) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7662738) q[0];
sx q[0];
rz(-3.0941071) q[0];
sx q[0];
rz(2.4256698) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.46220772) q[2];
sx q[2];
rz(-2.2773829) q[2];
sx q[2];
rz(-2.8414937) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2808025) q[1];
sx q[1];
rz(-1.6165015) q[1];
sx q[1];
rz(1.6153107) q[1];
rz(-pi) q[2];
rz(0.7089013) q[3];
sx q[3];
rz(-2.0376251) q[3];
sx q[3];
rz(-0.69463581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.82751194) q[2];
sx q[2];
rz(-0.22958799) q[2];
sx q[2];
rz(0.45544004) q[2];
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
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.6332508) q[0];
sx q[0];
rz(0.73927885) q[0];
rz(0.23070681) q[1];
sx q[1];
rz(-2.6735327) q[1];
sx q[1];
rz(-0.49490067) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5737168) q[0];
sx q[0];
rz(-1.084869) q[0];
sx q[0];
rz(-3.0987415) q[0];
rz(1.5779183) q[2];
sx q[2];
rz(-1.2429503) q[2];
sx q[2];
rz(0.42051007) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1315688) q[1];
sx q[1];
rz(-1.6820757) q[1];
sx q[1];
rz(2.2113423) q[1];
rz(-pi) q[2];
x q[2];
rz(0.87402113) q[3];
sx q[3];
rz(-1.6439983) q[3];
sx q[3];
rz(0.30833581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13359244) q[2];
sx q[2];
rz(-1.091489) q[2];
sx q[2];
rz(0.32785329) q[2];
rz(-2.949529) q[3];
sx q[3];
rz(-2.8938507) q[3];
sx q[3];
rz(-1.0333992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1223758) q[0];
sx q[0];
rz(-1.6288971) q[0];
sx q[0];
rz(-1.6678641) q[0];
rz(0.83256759) q[1];
sx q[1];
rz(-1.7201798) q[1];
sx q[1];
rz(-1.7725772) q[1];
rz(-1.9033296) q[2];
sx q[2];
rz(-0.56849545) q[2];
sx q[2];
rz(-0.93760437) q[2];
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
