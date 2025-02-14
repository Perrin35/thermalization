OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.71275467) q[0];
sx q[0];
rz(4.6729597) q[0];
sx q[0];
rz(10.529411) q[0];
rz(1.334335) q[1];
sx q[1];
rz(-2.4185138) q[1];
sx q[1];
rz(1.877797) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9845147) q[0];
sx q[0];
rz(-2.4369168) q[0];
sx q[0];
rz(-1.9542171) q[0];
rz(-pi) q[1];
x q[1];
rz(0.99743263) q[2];
sx q[2];
rz(-1.2094398) q[2];
sx q[2];
rz(-0.45623764) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7967148) q[1];
sx q[1];
rz(-1.895322) q[1];
sx q[1];
rz(-3.0561563) q[1];
x q[2];
rz(-2.8975211) q[3];
sx q[3];
rz(-1.9979949) q[3];
sx q[3];
rz(1.106316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9870712) q[2];
sx q[2];
rz(-1.0769341) q[2];
sx q[2];
rz(-1.3686352) q[2];
rz(-2.9030419) q[3];
sx q[3];
rz(-2.7261901) q[3];
sx q[3];
rz(-0.11428782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(-2.7411165) q[0];
sx q[0];
rz(-2.0023161) q[0];
sx q[0];
rz(-1.1503295) q[0];
rz(0.087609619) q[1];
sx q[1];
rz(-1.8116415) q[1];
sx q[1];
rz(-1.5709343) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51975361) q[0];
sx q[0];
rz(-1.6770419) q[0];
sx q[0];
rz(-0.48522093) q[0];
rz(0.17878647) q[2];
sx q[2];
rz(-1.6331722) q[2];
sx q[2];
rz(1.5449926) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0391866) q[1];
sx q[1];
rz(-2.3433609) q[1];
sx q[1];
rz(-0.95603966) q[1];
rz(0.47111311) q[3];
sx q[3];
rz(-2.2918252) q[3];
sx q[3];
rz(1.5136994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7972083) q[2];
sx q[2];
rz(-0.71343652) q[2];
sx q[2];
rz(2.9551282) q[2];
rz(2.3421085) q[3];
sx q[3];
rz(-1.5954285) q[3];
sx q[3];
rz(-2.1605087) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8568273) q[0];
sx q[0];
rz(-1.3815877) q[0];
sx q[0];
rz(-3.0322266) q[0];
rz(2.5378387) q[1];
sx q[1];
rz(-1.8021288) q[1];
sx q[1];
rz(-1.0341136) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7138192) q[0];
sx q[0];
rz(-1.7805011) q[0];
sx q[0];
rz(-0.13596491) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.500152) q[2];
sx q[2];
rz(-1.6598668) q[2];
sx q[2];
rz(-1.3180817) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2455403) q[1];
sx q[1];
rz(-2.8421092) q[1];
sx q[1];
rz(0.36548945) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9963861) q[3];
sx q[3];
rz(-2.3646486) q[3];
sx q[3];
rz(0.29175419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5321396) q[2];
sx q[2];
rz(-0.9708465) q[2];
sx q[2];
rz(-2.5965221) q[2];
rz(-2.3405781) q[3];
sx q[3];
rz(-0.89371926) q[3];
sx q[3];
rz(0.092122294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6545749) q[0];
sx q[0];
rz(-1.4734522) q[0];
sx q[0];
rz(-0.45904485) q[0];
rz(-2.0252939) q[1];
sx q[1];
rz(-2.3479159) q[1];
sx q[1];
rz(0.8265411) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7918043) q[0];
sx q[0];
rz(-1.9902181) q[0];
sx q[0];
rz(0.38695199) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.097657875) q[2];
sx q[2];
rz(-1.4726761) q[2];
sx q[2];
rz(-0.94933921) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.53163995) q[1];
sx q[1];
rz(-2.5594257) q[1];
sx q[1];
rz(1.9332631) q[1];
rz(2.8070368) q[3];
sx q[3];
rz(-1.4589196) q[3];
sx q[3];
rz(2.3528683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7894342) q[2];
sx q[2];
rz(-0.23098478) q[2];
sx q[2];
rz(-2.3918772) q[2];
rz(-2.0189144) q[3];
sx q[3];
rz(-1.2238945) q[3];
sx q[3];
rz(2.1813755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4753251) q[0];
sx q[0];
rz(-0.98639494) q[0];
sx q[0];
rz(-3.0295897) q[0];
rz(2.9169967) q[1];
sx q[1];
rz(-1.9738395) q[1];
sx q[1];
rz(-0.78132838) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8664426) q[0];
sx q[0];
rz(-0.40497447) q[0];
sx q[0];
rz(-1.5030131) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8728169) q[2];
sx q[2];
rz(-0.34160638) q[2];
sx q[2];
rz(-0.20526055) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0722149) q[1];
sx q[1];
rz(-0.88000127) q[1];
sx q[1];
rz(1.0051954) q[1];
rz(-pi) q[2];
rz(-1.8086241) q[3];
sx q[3];
rz(-2.5356511) q[3];
sx q[3];
rz(1.1693418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3096699) q[2];
sx q[2];
rz(-0.90193844) q[2];
sx q[2];
rz(0.93264467) q[2];
rz(2.0617088) q[3];
sx q[3];
rz(-0.44411689) q[3];
sx q[3];
rz(-0.035695765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90750736) q[0];
sx q[0];
rz(-1.1312753) q[0];
sx q[0];
rz(-1.750741) q[0];
rz(2.8098409) q[1];
sx q[1];
rz(-1.2650047) q[1];
sx q[1];
rz(1.5914241) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6734579) q[0];
sx q[0];
rz(-2.4979257) q[0];
sx q[0];
rz(1.629384) q[0];
rz(3.0671547) q[2];
sx q[2];
rz(-2.6808969) q[2];
sx q[2];
rz(-1.4975394) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0968714) q[1];
sx q[1];
rz(-1.7270178) q[1];
sx q[1];
rz(0.12466431) q[1];
x q[2];
rz(1.9523296) q[3];
sx q[3];
rz(-1.9100827) q[3];
sx q[3];
rz(-1.8041704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8288237) q[2];
sx q[2];
rz(-2.2569816) q[2];
sx q[2];
rz(1.620232) q[2];
rz(-0.96794266) q[3];
sx q[3];
rz(-1.5588375) q[3];
sx q[3];
rz(2.2255285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5172326) q[0];
sx q[0];
rz(-0.42651287) q[0];
sx q[0];
rz(-0.17159167) q[0];
rz(1.0393556) q[1];
sx q[1];
rz(-1.8828705) q[1];
sx q[1];
rz(1.7003869) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65214163) q[0];
sx q[0];
rz(-1.2496557) q[0];
sx q[0];
rz(-1.4399066) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4588548) q[2];
sx q[2];
rz(-2.0840692) q[2];
sx q[2];
rz(1.7626732) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.68908697) q[1];
sx q[1];
rz(-1.2156665) q[1];
sx q[1];
rz(-2.0562812) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.14136858) q[3];
sx q[3];
rz(-1.0118359) q[3];
sx q[3];
rz(-1.5593004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.27199304) q[2];
sx q[2];
rz(-1.7337493) q[2];
sx q[2];
rz(0.54455152) q[2];
rz(-0.49992391) q[3];
sx q[3];
rz(-1.0285503) q[3];
sx q[3];
rz(0.47206363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2917824) q[0];
sx q[0];
rz(-0.53684679) q[0];
sx q[0];
rz(-2.612402) q[0];
rz(-0.57580194) q[1];
sx q[1];
rz(-1.2628097) q[1];
sx q[1];
rz(-2.0226488) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42903439) q[0];
sx q[0];
rz(-0.055744113) q[0];
sx q[0];
rz(-0.38614614) q[0];
rz(-pi) q[1];
rz(-2.7706657) q[2];
sx q[2];
rz(-1.7417522) q[2];
sx q[2];
rz(1.8425187) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.24209472) q[1];
sx q[1];
rz(-0.78259727) q[1];
sx q[1];
rz(-0.67983869) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0980074) q[3];
sx q[3];
rz(-1.3889179) q[3];
sx q[3];
rz(1.393569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.240856) q[2];
sx q[2];
rz(-0.46515981) q[2];
sx q[2];
rz(-0.71211234) q[2];
rz(1.5322878) q[3];
sx q[3];
rz(-2.1963547) q[3];
sx q[3];
rz(-1.3473264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79494548) q[0];
sx q[0];
rz(-1.1279339) q[0];
sx q[0];
rz(-2.0704863) q[0];
rz(-1.2062997) q[1];
sx q[1];
rz(-1.0164398) q[1];
sx q[1];
rz(-1.6546904) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3275571) q[0];
sx q[0];
rz(-0.90891713) q[0];
sx q[0];
rz(1.3332086) q[0];
x q[1];
rz(0.98574443) q[2];
sx q[2];
rz(-0.92454708) q[2];
sx q[2];
rz(-2.1304634) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2675954) q[1];
sx q[1];
rz(-1.2800299) q[1];
sx q[1];
rz(-1.7520755) q[1];
x q[2];
rz(2.1110299) q[3];
sx q[3];
rz(-2.6298454) q[3];
sx q[3];
rz(-0.47357163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.508076) q[2];
sx q[2];
rz(-2.241892) q[2];
sx q[2];
rz(0.077795204) q[2];
rz(0.22732321) q[3];
sx q[3];
rz(-1.1319755) q[3];
sx q[3];
rz(1.0556861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28854293) q[0];
sx q[0];
rz(-0.74497861) q[0];
sx q[0];
rz(2.7689834) q[0];
rz(2.695072) q[1];
sx q[1];
rz(-1.6852854) q[1];
sx q[1];
rz(2.2850697) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83607446) q[0];
sx q[0];
rz(-1.5139607) q[0];
sx q[0];
rz(-1.6147805) q[0];
rz(3.0356016) q[2];
sx q[2];
rz(-1.4376831) q[2];
sx q[2];
rz(0.0078545257) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0622934) q[1];
sx q[1];
rz(-2.7253299) q[1];
sx q[1];
rz(1.8658616) q[1];
rz(-pi) q[2];
rz(0.89337279) q[3];
sx q[3];
rz(-2.2489298) q[3];
sx q[3];
rz(-3.107389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7591758) q[2];
sx q[2];
rz(-2.0903812) q[2];
sx q[2];
rz(-2.5860533) q[2];
rz(2.7555452) q[3];
sx q[3];
rz(-1.4945533) q[3];
sx q[3];
rz(-2.4773795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0059218) q[0];
sx q[0];
rz(-1.4501403) q[0];
sx q[0];
rz(-1.8474664) q[0];
rz(0.80264965) q[1];
sx q[1];
rz(-1.4603271) q[1];
sx q[1];
rz(-2.7535798) q[1];
rz(0.31617185) q[2];
sx q[2];
rz(-0.36973047) q[2];
sx q[2];
rz(0.75512259) q[2];
rz(-0.070128154) q[3];
sx q[3];
rz(-0.84134103) q[3];
sx q[3];
rz(0.76920912) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
