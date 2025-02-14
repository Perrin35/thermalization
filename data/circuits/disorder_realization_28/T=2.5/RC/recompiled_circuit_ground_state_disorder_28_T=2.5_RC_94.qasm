OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.82792038) q[0];
sx q[0];
rz(-1.6452687) q[0];
sx q[0];
rz(-0.21284719) q[0];
rz(-1.4053474) q[1];
sx q[1];
rz(-1.5264629) q[1];
sx q[1];
rz(2.5636173) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5151857) q[0];
sx q[0];
rz(-0.66045633) q[0];
sx q[0];
rz(-0.68645262) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.87066381) q[2];
sx q[2];
rz(-1.7144702) q[2];
sx q[2];
rz(-1.1149017) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2702944) q[1];
sx q[1];
rz(-1.9828856) q[1];
sx q[1];
rz(-0.006448556) q[1];
rz(-0.36823456) q[3];
sx q[3];
rz(-1.7126178) q[3];
sx q[3];
rz(-2.5960245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2026244) q[2];
sx q[2];
rz(-1.6309996) q[2];
sx q[2];
rz(2.5228339) q[2];
rz(-0.42801157) q[3];
sx q[3];
rz(-1.908554) q[3];
sx q[3];
rz(-1.4324043) q[3];
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
rz(1.0851058) q[0];
sx q[0];
rz(-3.1054057) q[0];
sx q[0];
rz(-0.65271839) q[0];
rz(-2.808049) q[1];
sx q[1];
rz(-2.3338552) q[1];
sx q[1];
rz(-0.52887598) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6655114) q[0];
sx q[0];
rz(-1.9383069) q[0];
sx q[0];
rz(-1.0239081) q[0];
rz(-0.11153259) q[2];
sx q[2];
rz(-2.1443488) q[2];
sx q[2];
rz(-0.83441624) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3470001) q[1];
sx q[1];
rz(-2.9459822) q[1];
sx q[1];
rz(0.83949714) q[1];
rz(-pi) q[2];
rz(0.06270285) q[3];
sx q[3];
rz(-2.3093105) q[3];
sx q[3];
rz(-0.073918561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.35398206) q[2];
sx q[2];
rz(-1.3926316) q[2];
sx q[2];
rz(-2.7667513) q[2];
rz(-1.2104872) q[3];
sx q[3];
rz(-0.46606246) q[3];
sx q[3];
rz(-1.2208285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1061363) q[0];
sx q[0];
rz(-1.9566256) q[0];
sx q[0];
rz(0.55738604) q[0];
rz(-2.4222971) q[1];
sx q[1];
rz(-0.44984111) q[1];
sx q[1];
rz(2.5296899) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5267964) q[0];
sx q[0];
rz(-0.4810534) q[0];
sx q[0];
rz(0.78562268) q[0];
x q[1];
rz(0.8863897) q[2];
sx q[2];
rz(-1.7409313) q[2];
sx q[2];
rz(-1.1332133) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.53781834) q[1];
sx q[1];
rz(-0.99939954) q[1];
sx q[1];
rz(0.091544108) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7541998) q[3];
sx q[3];
rz(-1.9678183) q[3];
sx q[3];
rz(-2.3698636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8039852) q[2];
sx q[2];
rz(-1.3472202) q[2];
sx q[2];
rz(-0.95532974) q[2];
rz(2.5721926) q[3];
sx q[3];
rz(-2.0208713) q[3];
sx q[3];
rz(-1.7795405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78205458) q[0];
sx q[0];
rz(-0.29490092) q[0];
sx q[0];
rz(-0.094757946) q[0];
rz(1.5358198) q[1];
sx q[1];
rz(-1.1631807) q[1];
sx q[1];
rz(-0.43704978) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9511918) q[0];
sx q[0];
rz(-3.0110208) q[0];
sx q[0];
rz(-2.269362) q[0];
rz(1.3895459) q[2];
sx q[2];
rz(-1.1409352) q[2];
sx q[2];
rz(1.4409122) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6833675) q[1];
sx q[1];
rz(-1.5705854) q[1];
sx q[1];
rz(1.3092759) q[1];
rz(-pi) q[2];
rz(-2.0548204) q[3];
sx q[3];
rz(-2.4216702) q[3];
sx q[3];
rz(2.7663224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7729418) q[2];
sx q[2];
rz(-2.3029885) q[2];
sx q[2];
rz(1.3402026) q[2];
rz(2.7784427) q[3];
sx q[3];
rz(-2.4996417) q[3];
sx q[3];
rz(-2.6974881) q[3];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7813613) q[0];
sx q[0];
rz(-2.4960127) q[0];
sx q[0];
rz(-0.050882291) q[0];
rz(0.52967349) q[1];
sx q[1];
rz(-2.5310204) q[1];
sx q[1];
rz(-0.018459056) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0063361703) q[0];
sx q[0];
rz(-1.0632035) q[0];
sx q[0];
rz(-0.35180845) q[0];
x q[1];
rz(1.7713806) q[2];
sx q[2];
rz(-0.86435741) q[2];
sx q[2];
rz(-2.4141224) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.94973719) q[1];
sx q[1];
rz(-2.4337447) q[1];
sx q[1];
rz(1.1337277) q[1];
x q[2];
rz(-0.2996188) q[3];
sx q[3];
rz(-0.54237759) q[3];
sx q[3];
rz(0.49132916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9336443) q[2];
sx q[2];
rz(-1.4241781) q[2];
sx q[2];
rz(-1.8682182) q[2];
rz(0.4452855) q[3];
sx q[3];
rz(-2.4280426) q[3];
sx q[3];
rz(-0.85550296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2263366) q[0];
sx q[0];
rz(-2.7610918) q[0];
sx q[0];
rz(-0.29391995) q[0];
rz(1.2433012) q[1];
sx q[1];
rz(-1.2970122) q[1];
sx q[1];
rz(-1.6506724) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4478233) q[0];
sx q[0];
rz(-2.4943879) q[0];
sx q[0];
rz(-1.8497516) q[0];
rz(0.32604172) q[2];
sx q[2];
rz(-1.4801822) q[2];
sx q[2];
rz(2.2590274) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6676644) q[1];
sx q[1];
rz(-1.3753624) q[1];
sx q[1];
rz(0.21904314) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0206679) q[3];
sx q[3];
rz(-0.48512019) q[3];
sx q[3];
rz(-2.7435722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2095498) q[2];
sx q[2];
rz(-1.8609293) q[2];
sx q[2];
rz(0.10188421) q[2];
rz(0.79002506) q[3];
sx q[3];
rz(-0.70370379) q[3];
sx q[3];
rz(1.2557282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1387966) q[0];
sx q[0];
rz(-3.0468472) q[0];
sx q[0];
rz(-0.62762678) q[0];
rz(-1.7812642) q[1];
sx q[1];
rz(-1.5831213) q[1];
sx q[1];
rz(0.21242177) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4110574) q[0];
sx q[0];
rz(-0.90170398) q[0];
sx q[0];
rz(-2.6691666) q[0];
rz(0.27842267) q[2];
sx q[2];
rz(-2.1489193) q[2];
sx q[2];
rz(2.9591564) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2523035) q[1];
sx q[1];
rz(-2.6037137) q[1];
sx q[1];
rz(-1.1982615) q[1];
rz(-pi) q[2];
rz(0.48177367) q[3];
sx q[3];
rz(-1.9984233) q[3];
sx q[3];
rz(-2.4419344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3820496) q[2];
sx q[2];
rz(-0.86882773) q[2];
sx q[2];
rz(2.7810435) q[2];
rz(-0.77066317) q[3];
sx q[3];
rz(-1.9375075) q[3];
sx q[3];
rz(0.26309052) q[3];
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
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78779546) q[0];
sx q[0];
rz(-2.0833092) q[0];
sx q[0];
rz(1.298792) q[0];
rz(-2.6020488) q[1];
sx q[1];
rz(-1.6863457) q[1];
sx q[1];
rz(1.65666) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7138707) q[0];
sx q[0];
rz(-2.1146647) q[0];
sx q[0];
rz(1.9440249) q[0];
x q[1];
rz(0.46631821) q[2];
sx q[2];
rz(-2.6728485) q[2];
sx q[2];
rz(-3.0006204) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.114552) q[1];
sx q[1];
rz(-1.5590347) q[1];
sx q[1];
rz(1.8742754) q[1];
rz(-pi) q[2];
rz(1.1670145) q[3];
sx q[3];
rz(-2.2444668) q[3];
sx q[3];
rz(0.63577494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3761042) q[2];
sx q[2];
rz(-1.8643943) q[2];
sx q[2];
rz(-0.45846024) q[2];
rz(-1.1485398) q[3];
sx q[3];
rz(-0.12980041) q[3];
sx q[3];
rz(2.5965918) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5371567) q[0];
sx q[0];
rz(-0.24188365) q[0];
sx q[0];
rz(2.3271374) q[0];
rz(-0.7704598) q[1];
sx q[1];
rz(-1.5372246) q[1];
sx q[1];
rz(1.8155712) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6961937) q[0];
sx q[0];
rz(-0.32055453) q[0];
sx q[0];
rz(0.93648367) q[0];
rz(-pi) q[1];
rz(3.0203462) q[2];
sx q[2];
rz(-0.7077982) q[2];
sx q[2];
rz(-1.540465) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0596366) q[1];
sx q[1];
rz(-0.15742144) q[1];
sx q[1];
rz(-1.6379959) q[1];
x q[2];
rz(2.2758851) q[3];
sx q[3];
rz(-2.9031347) q[3];
sx q[3];
rz(-1.2442279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.068437964) q[2];
sx q[2];
rz(-2.1775553) q[2];
sx q[2];
rz(2.5835719) q[2];
rz(-0.26992118) q[3];
sx q[3];
rz(-0.25756613) q[3];
sx q[3];
rz(0.71815467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59543264) q[0];
sx q[0];
rz(-1.3256925) q[0];
sx q[0];
rz(-0.85920715) q[0];
rz(-0.72256207) q[1];
sx q[1];
rz(-0.95046202) q[1];
sx q[1];
rz(1.439646) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42992698) q[0];
sx q[0];
rz(-1.5646805) q[0];
sx q[0];
rz(-1.5625009) q[0];
rz(-pi) q[1];
rz(0.53368081) q[2];
sx q[2];
rz(-0.42842487) q[2];
sx q[2];
rz(2.3103586) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5240092) q[1];
sx q[1];
rz(-1.7477711) q[1];
sx q[1];
rz(2.3664775) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3297746) q[3];
sx q[3];
rz(-1.9621578) q[3];
sx q[3];
rz(-2.4956839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.667111) q[2];
sx q[2];
rz(-1.2789395) q[2];
sx q[2];
rz(-1.4943592) q[2];
rz(0.4161559) q[3];
sx q[3];
rz(-1.4978724) q[3];
sx q[3];
rz(-0.070629899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1304929) q[0];
sx q[0];
rz(-0.53642219) q[0];
sx q[0];
rz(-0.82905967) q[0];
rz(-2.2566055) q[1];
sx q[1];
rz(-2.5288455) q[1];
sx q[1];
rz(1.7987342) q[1];
rz(3.0961453) q[2];
sx q[2];
rz(-0.54879594) q[2];
sx q[2];
rz(-1.9898228) q[2];
rz(0.78444278) q[3];
sx q[3];
rz(-1.7517437) q[3];
sx q[3];
rz(-2.4117398) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
