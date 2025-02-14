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
rz(-2.3186853) q[0];
sx q[0];
rz(3.5003852) q[0];
sx q[0];
rz(8.556463) q[0];
rz(-1.4153642) q[1];
sx q[1];
rz(-1.1338898) q[1];
sx q[1];
rz(2.147832) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.225634) q[0];
sx q[0];
rz(-1.4580112) q[0];
sx q[0];
rz(0.40333545) q[0];
rz(-pi) q[1];
rz(-2.0103942) q[2];
sx q[2];
rz(-2.0313259) q[2];
sx q[2];
rz(-1.6701513) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6958625) q[1];
sx q[1];
rz(-0.63803405) q[1];
sx q[1];
rz(0.7971493) q[1];
rz(-pi) q[2];
rz(0.24561974) q[3];
sx q[3];
rz(-1.9960476) q[3];
sx q[3];
rz(-2.5568145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.913784) q[2];
sx q[2];
rz(-1.8081534) q[2];
sx q[2];
rz(-0.58161962) q[2];
rz(0.73451129) q[3];
sx q[3];
rz(-1.651265) q[3];
sx q[3];
rz(2.7233126) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2468579) q[0];
sx q[0];
rz(-1.3792091) q[0];
sx q[0];
rz(-2.6439164) q[0];
rz(1.0408164) q[1];
sx q[1];
rz(-0.40319315) q[1];
sx q[1];
rz(-1.2373479) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.055641551) q[0];
sx q[0];
rz(-0.72362542) q[0];
sx q[0];
rz(-3.0776204) q[0];
rz(-pi) q[1];
rz(-1.2096268) q[2];
sx q[2];
rz(-1.3945082) q[2];
sx q[2];
rz(0.90106264) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9927382) q[1];
sx q[1];
rz(-1.9991367) q[1];
sx q[1];
rz(-0.15495877) q[1];
x q[2];
rz(2.4604843) q[3];
sx q[3];
rz(-0.33675413) q[3];
sx q[3];
rz(0.82267534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.441338) q[2];
sx q[2];
rz(-2.4105218) q[2];
sx q[2];
rz(-2.8311484) q[2];
rz(-1.1566409) q[3];
sx q[3];
rz(-2.0168596) q[3];
sx q[3];
rz(1.9748851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.1336841) q[0];
sx q[0];
rz(-0.5846566) q[0];
sx q[0];
rz(0.34580082) q[0];
rz(2.8969104) q[1];
sx q[1];
rz(-2.209765) q[1];
sx q[1];
rz(1.7049047) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0999683) q[0];
sx q[0];
rz(-1.1505373) q[0];
sx q[0];
rz(-2.7970275) q[0];
x q[1];
rz(-1.1780924) q[2];
sx q[2];
rz(-2.1059193) q[2];
sx q[2];
rz(1.8491883) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6356698) q[1];
sx q[1];
rz(-1.782592) q[1];
sx q[1];
rz(-0.79896547) q[1];
x q[2];
rz(-2.6639943) q[3];
sx q[3];
rz(-1.2354038) q[3];
sx q[3];
rz(0.23303495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1338542) q[2];
sx q[2];
rz(-2.014092) q[2];
sx q[2];
rz(0.63068843) q[2];
rz(2.808029) q[3];
sx q[3];
rz(-1.0015229) q[3];
sx q[3];
rz(1.001531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0665322) q[0];
sx q[0];
rz(-0.94828951) q[0];
sx q[0];
rz(1.7370268) q[0];
rz(-2.7724077) q[1];
sx q[1];
rz(-1.4063947) q[1];
sx q[1];
rz(-0.085748347) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2959901) q[0];
sx q[0];
rz(-1.1628502) q[0];
sx q[0];
rz(1.8970117) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6811738) q[2];
sx q[2];
rz(-1.9398089) q[2];
sx q[2];
rz(-2.2031914) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4032674) q[1];
sx q[1];
rz(-1.8083296) q[1];
sx q[1];
rz(2.0473785) q[1];
rz(0.39601456) q[3];
sx q[3];
rz(-1.2761444) q[3];
sx q[3];
rz(2.9673607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.839445) q[2];
sx q[2];
rz(-1.2371233) q[2];
sx q[2];
rz(-0.5298003) q[2];
rz(-0.038330404) q[3];
sx q[3];
rz(-0.72961346) q[3];
sx q[3];
rz(-1.5420325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2235276) q[0];
sx q[0];
rz(-2.6730838) q[0];
sx q[0];
rz(1.9388306) q[0];
rz(-2.862152) q[1];
sx q[1];
rz(-1.025082) q[1];
sx q[1];
rz(2.3675945) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27085486) q[0];
sx q[0];
rz(-1.8841198) q[0];
sx q[0];
rz(-0.82672755) q[0];
rz(-pi) q[1];
rz(2.8061637) q[2];
sx q[2];
rz(-2.895439) q[2];
sx q[2];
rz(-0.85573643) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2629649) q[1];
sx q[1];
rz(-1.557918) q[1];
sx q[1];
rz(1.7098544) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0353885) q[3];
sx q[3];
rz(-1.5683953) q[3];
sx q[3];
rz(2.6794499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0224132) q[2];
sx q[2];
rz(-3.000562) q[2];
sx q[2];
rz(-0.5640344) q[2];
rz(-0.65308475) q[3];
sx q[3];
rz(-1.1354732) q[3];
sx q[3];
rz(-0.45026067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3611203) q[0];
sx q[0];
rz(-2.5887964) q[0];
sx q[0];
rz(2.4984388) q[0];
rz(1.9505352) q[1];
sx q[1];
rz(-1.4570313) q[1];
sx q[1];
rz(0.65470421) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7793286) q[0];
sx q[0];
rz(-1.6342499) q[0];
sx q[0];
rz(-0.16880798) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.68393647) q[2];
sx q[2];
rz(-1.6003055) q[2];
sx q[2];
rz(-2.3060407) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.0400078) q[1];
sx q[1];
rz(-1.5052649) q[1];
sx q[1];
rz(1.0270018) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1804462) q[3];
sx q[3];
rz(-1.6359513) q[3];
sx q[3];
rz(-2.3303243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6341298) q[2];
sx q[2];
rz(-2.7644988) q[2];
sx q[2];
rz(-1.6667574) q[2];
rz(-2.6569488) q[3];
sx q[3];
rz(-2.1438997) q[3];
sx q[3];
rz(-1.2887597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6396879) q[0];
sx q[0];
rz(-0.26308331) q[0];
sx q[0];
rz(-2.4581773) q[0];
rz(0.13038334) q[1];
sx q[1];
rz(-1.5510635) q[1];
sx q[1];
rz(0.032616671) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6382173) q[0];
sx q[0];
rz(-2.1299358) q[0];
sx q[0];
rz(2.9342164) q[0];
x q[1];
rz(-0.68565418) q[2];
sx q[2];
rz(-1.7259806) q[2];
sx q[2];
rz(2.6557166) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8396388) q[1];
sx q[1];
rz(-2.8263013) q[1];
sx q[1];
rz(-1.006542) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86883162) q[3];
sx q[3];
rz(-1.371939) q[3];
sx q[3];
rz(1.6384038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.42328295) q[2];
sx q[2];
rz(-2.0330567) q[2];
sx q[2];
rz(-0.56524593) q[2];
rz(-1.7806753) q[3];
sx q[3];
rz(-0.2667242) q[3];
sx q[3];
rz(2.3274073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3703506) q[0];
sx q[0];
rz(-0.057436198) q[0];
sx q[0];
rz(2.8420319) q[0];
rz(1.6948505) q[1];
sx q[1];
rz(-0.87211496) q[1];
sx q[1];
rz(2.1902693) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0977391) q[0];
sx q[0];
rz(-1.365775) q[0];
sx q[0];
rz(-0.033173843) q[0];
rz(1.3592654) q[2];
sx q[2];
rz(-1.5617687) q[2];
sx q[2];
rz(0.8789076) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.48288667) q[1];
sx q[1];
rz(-1.1466007) q[1];
sx q[1];
rz(-1.5484323) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0537576) q[3];
sx q[3];
rz(-1.3099652) q[3];
sx q[3];
rz(1.5276599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1424554) q[2];
sx q[2];
rz(-0.74644011) q[2];
sx q[2];
rz(2.8738521) q[2];
rz(-2.1382051) q[3];
sx q[3];
rz(-0.95971862) q[3];
sx q[3];
rz(-0.80823922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7830911) q[0];
sx q[0];
rz(-2.6434904) q[0];
sx q[0];
rz(-0.95712334) q[0];
rz(-0.3793017) q[1];
sx q[1];
rz(-0.4117659) q[1];
sx q[1];
rz(-2.9764825) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3527213) q[0];
sx q[0];
rz(-1.9130978) q[0];
sx q[0];
rz(-0.71705937) q[0];
rz(-0.050641955) q[2];
sx q[2];
rz(-1.0989185) q[2];
sx q[2];
rz(0.34220055) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.86805008) q[1];
sx q[1];
rz(-2.1434577) q[1];
sx q[1];
rz(0.21577253) q[1];
x q[2];
rz(1.2312789) q[3];
sx q[3];
rz(-2.5105308) q[3];
sx q[3];
rz(1.5724044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.86753201) q[2];
sx q[2];
rz(-1.910285) q[2];
sx q[2];
rz(0.77862281) q[2];
rz(0.33729956) q[3];
sx q[3];
rz(-1.0236579) q[3];
sx q[3];
rz(-1.5571099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.857665) q[0];
sx q[0];
rz(-2.0511257) q[0];
sx q[0];
rz(1.0887867) q[0];
rz(-0.42824832) q[1];
sx q[1];
rz(-1.0232404) q[1];
sx q[1];
rz(-2.4931989) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26864949) q[0];
sx q[0];
rz(-0.50794426) q[0];
sx q[0];
rz(-2.5769039) q[0];
rz(-0.044610046) q[2];
sx q[2];
rz(-1.9084435) q[2];
sx q[2];
rz(0.57429796) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9251483) q[1];
sx q[1];
rz(-2.2183462) q[1];
sx q[1];
rz(-1.5179894) q[1];
rz(2.9620578) q[3];
sx q[3];
rz(-0.57099062) q[3];
sx q[3];
rz(1.4654311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.26583656) q[2];
sx q[2];
rz(-2.0529604) q[2];
sx q[2];
rz(-0.17975532) q[2];
rz(1.9836551) q[3];
sx q[3];
rz(-1.6854743) q[3];
sx q[3];
rz(1.2402844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4614048) q[0];
sx q[0];
rz(-2.9869933) q[0];
sx q[0];
rz(-2.3416478) q[0];
rz(-1.1663306) q[1];
sx q[1];
rz(-0.98465289) q[1];
sx q[1];
rz(-2.226895) q[1];
rz(-2.8958252) q[2];
sx q[2];
rz(-1.2817597) q[2];
sx q[2];
rz(1.3195932) q[2];
rz(2.8020482) q[3];
sx q[3];
rz(-1.3673269) q[3];
sx q[3];
rz(0.23585503) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
