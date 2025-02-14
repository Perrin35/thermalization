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
rz(2.2344196) q[0];
sx q[0];
rz(-0.89651674) q[0];
sx q[0];
rz(3.1007015) q[0];
rz(-1.2603124) q[1];
sx q[1];
rz(-2.6949096) q[1];
sx q[1];
rz(0.094951542) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4584158) q[0];
sx q[0];
rz(-1.9703034) q[0];
sx q[0];
rz(1.3241121) q[0];
rz(-pi) q[1];
rz(-1.3407732) q[2];
sx q[2];
rz(-1.025295) q[2];
sx q[2];
rz(-0.70822424) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5168001) q[1];
sx q[1];
rz(-1.2007371) q[1];
sx q[1];
rz(2.3462252) q[1];
x q[2];
rz(2.0435752) q[3];
sx q[3];
rz(-1.1679139) q[3];
sx q[3];
rz(-1.1927746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8591259) q[2];
sx q[2];
rz(-2.5390415) q[2];
sx q[2];
rz(-1.9762543) q[2];
rz(2.2230542) q[3];
sx q[3];
rz(-0.10401741) q[3];
sx q[3];
rz(1.9281841) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7391881) q[0];
sx q[0];
rz(-0.63888752) q[0];
sx q[0];
rz(0.30493394) q[0];
rz(-1.0324427) q[1];
sx q[1];
rz(-2.86125) q[1];
sx q[1];
rz(-0.84411821) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2518678) q[0];
sx q[0];
rz(-1.8093523) q[0];
sx q[0];
rz(-2.9846334) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3621037) q[2];
sx q[2];
rz(-1.450453) q[2];
sx q[2];
rz(-0.99882767) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1010437) q[1];
sx q[1];
rz(-1.1660468) q[1];
sx q[1];
rz(-1.3084662) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5145766) q[3];
sx q[3];
rz(-1.5582784) q[3];
sx q[3];
rz(2.1345958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1673163) q[2];
sx q[2];
rz(-0.68334371) q[2];
sx q[2];
rz(2.5565476) q[2];
rz(-0.85917464) q[3];
sx q[3];
rz(-1.1480568) q[3];
sx q[3];
rz(-2.008321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1187196) q[0];
sx q[0];
rz(-0.82094231) q[0];
sx q[0];
rz(-0.098544772) q[0];
rz(-2.5755644) q[1];
sx q[1];
rz(-0.84605828) q[1];
sx q[1];
rz(2.6354094) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9509256) q[0];
sx q[0];
rz(-1.902421) q[0];
sx q[0];
rz(2.5965461) q[0];
rz(0.82019851) q[2];
sx q[2];
rz(-1.2389822) q[2];
sx q[2];
rz(-0.92727509) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0876956) q[1];
sx q[1];
rz(-1.3511041) q[1];
sx q[1];
rz(3.0151574) q[1];
x q[2];
rz(0.39301707) q[3];
sx q[3];
rz(-1.0342433) q[3];
sx q[3];
rz(-0.69325939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.026406) q[2];
sx q[2];
rz(-1.9280484) q[2];
sx q[2];
rz(0.085065993) q[2];
rz(-2.3885942) q[3];
sx q[3];
rz(-2.4716061) q[3];
sx q[3];
rz(-1.5299214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-1.6104777) q[0];
sx q[0];
rz(-1.5014638) q[0];
sx q[0];
rz(0.53989545) q[0];
rz(-3.1119697) q[1];
sx q[1];
rz(-0.74631515) q[1];
sx q[1];
rz(1.7866887) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82096174) q[0];
sx q[0];
rz(-2.3665049) q[0];
sx q[0];
rz(1.8863999) q[0];
rz(-pi) q[1];
rz(0.37432713) q[2];
sx q[2];
rz(-0.74639635) q[2];
sx q[2];
rz(2.057892) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.78614178) q[1];
sx q[1];
rz(-1.0075354) q[1];
sx q[1];
rz(-1.2652629) q[1];
x q[2];
rz(-2.3514495) q[3];
sx q[3];
rz(-1.6340268) q[3];
sx q[3];
rz(2.9823398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.70298355) q[2];
sx q[2];
rz(-0.97038022) q[2];
sx q[2];
rz(-2.2461829) q[2];
rz(-1.4488975) q[3];
sx q[3];
rz(-1.1619032) q[3];
sx q[3];
rz(0.40941516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92622906) q[0];
sx q[0];
rz(-2.138593) q[0];
sx q[0];
rz(-0.0005501752) q[0];
rz(-2.6161361) q[1];
sx q[1];
rz(-0.84392396) q[1];
sx q[1];
rz(0.97132436) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4068953) q[0];
sx q[0];
rz(-0.79888035) q[0];
sx q[0];
rz(-2.6968234) q[0];
x q[1];
rz(-2.2699666) q[2];
sx q[2];
rz(-0.52025929) q[2];
sx q[2];
rz(1.546788) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3974053) q[1];
sx q[1];
rz(-2.1402855) q[1];
sx q[1];
rz(2.2389212) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8841142) q[3];
sx q[3];
rz(-1.7755812) q[3];
sx q[3];
rz(1.9792894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9690669) q[2];
sx q[2];
rz(-2.064164) q[2];
sx q[2];
rz(1.7871008) q[2];
rz(-0.6012249) q[3];
sx q[3];
rz(-2.0982592) q[3];
sx q[3];
rz(-1.575298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71156597) q[0];
sx q[0];
rz(-0.30421782) q[0];
sx q[0];
rz(-2.8416204) q[0];
rz(-0.9777588) q[1];
sx q[1];
rz(-2.561196) q[1];
sx q[1];
rz(0.93394867) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47563206) q[0];
sx q[0];
rz(-2.7802659) q[0];
sx q[0];
rz(0.40478171) q[0];
x q[1];
rz(2.5061692) q[2];
sx q[2];
rz(-2.4197856) q[2];
sx q[2];
rz(1.3363163) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.63864795) q[1];
sx q[1];
rz(-2.3343122) q[1];
sx q[1];
rz(2.817383) q[1];
rz(-1.9225695) q[3];
sx q[3];
rz(-1.9493503) q[3];
sx q[3];
rz(1.564725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0705491) q[2];
sx q[2];
rz(-2.1328378) q[2];
sx q[2];
rz(-2.0096931) q[2];
rz(-1.2567358) q[3];
sx q[3];
rz(-2.3102424) q[3];
sx q[3];
rz(-1.7251714) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32135949) q[0];
sx q[0];
rz(-1.2437404) q[0];
sx q[0];
rz(-2.5309122) q[0];
rz(-2.3848379) q[1];
sx q[1];
rz(-2.7532531) q[1];
sx q[1];
rz(-1.6441708) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36726481) q[0];
sx q[0];
rz(-1.5455565) q[0];
sx q[0];
rz(-0.12406524) q[0];
rz(-pi) q[1];
x q[1];
rz(1.146422) q[2];
sx q[2];
rz(-1.904663) q[2];
sx q[2];
rz(0.41347143) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.037379384) q[1];
sx q[1];
rz(-0.92475212) q[1];
sx q[1];
rz(1.437633) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.11388643) q[3];
sx q[3];
rz(-2.4166757) q[3];
sx q[3];
rz(2.8020086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.18621592) q[2];
sx q[2];
rz(-2.4775041) q[2];
sx q[2];
rz(-0.60626283) q[2];
rz(-2.9210505) q[3];
sx q[3];
rz(-1.8044148) q[3];
sx q[3];
rz(-2.7748599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1425988) q[0];
sx q[0];
rz(-1.6479011) q[0];
sx q[0];
rz(2.2338474) q[0];
rz(-1.6084464) q[1];
sx q[1];
rz(-0.95450675) q[1];
sx q[1];
rz(-0.018891637) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9241727) q[0];
sx q[0];
rz(-2.0402363) q[0];
sx q[0];
rz(-1.4408235) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.18757815) q[2];
sx q[2];
rz(-2.8370428) q[2];
sx q[2];
rz(-1.9001324) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.42378473) q[1];
sx q[1];
rz(-0.46324364) q[1];
sx q[1];
rz(-0.06991212) q[1];
rz(-2.3109679) q[3];
sx q[3];
rz(-2.0181832) q[3];
sx q[3];
rz(0.4448286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3733526) q[2];
sx q[2];
rz(-1.588593) q[2];
sx q[2];
rz(0.24660435) q[2];
rz(1.8433579) q[3];
sx q[3];
rz(-1.4539098) q[3];
sx q[3];
rz(1.2687782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75803718) q[0];
sx q[0];
rz(-2.0691431) q[0];
sx q[0];
rz(-0.86395907) q[0];
rz(-1.2147238) q[1];
sx q[1];
rz(-2.0236969) q[1];
sx q[1];
rz(2.5606959) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6462517) q[0];
sx q[0];
rz(-2.2682307) q[0];
sx q[0];
rz(2.4861885) q[0];
rz(-1.6132108) q[2];
sx q[2];
rz(-2.4047432) q[2];
sx q[2];
rz(2.9880817) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.249596) q[1];
sx q[1];
rz(-0.70887762) q[1];
sx q[1];
rz(2.6061771) q[1];
x q[2];
rz(-2.8618699) q[3];
sx q[3];
rz(-2.3424853) q[3];
sx q[3];
rz(-0.38995686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2555799) q[2];
sx q[2];
rz(-0.46147999) q[2];
sx q[2];
rz(-0.18730051) q[2];
rz(2.1646132) q[3];
sx q[3];
rz(-1.6410442) q[3];
sx q[3];
rz(-1.2788844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63000694) q[0];
sx q[0];
rz(-2.4749909) q[0];
sx q[0];
rz(-1.2641719) q[0];
rz(2.1695747) q[1];
sx q[1];
rz(-2.3897901) q[1];
sx q[1];
rz(-1.1690296) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80315932) q[0];
sx q[0];
rz(-1.2375323) q[0];
sx q[0];
rz(0.14305556) q[0];
x q[1];
rz(0.75836597) q[2];
sx q[2];
rz(-2.6347199) q[2];
sx q[2];
rz(-3.1257983) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0824521) q[1];
sx q[1];
rz(-2.2603717) q[1];
sx q[1];
rz(-3.0060958) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8253978) q[3];
sx q[3];
rz(-1.4619816) q[3];
sx q[3];
rz(1.8121882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.4740037) q[2];
sx q[2];
rz(-1.7550125) q[2];
sx q[2];
rz(0.37707314) q[2];
rz(-0.22370473) q[3];
sx q[3];
rz(-2.5535899) q[3];
sx q[3];
rz(1.912775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8857464) q[0];
sx q[0];
rz(-1.7485913) q[0];
sx q[0];
rz(2.8843256) q[0];
rz(-2.5480351) q[1];
sx q[1];
rz(-2.3496353) q[1];
sx q[1];
rz(-1.4812462) q[1];
rz(-2.42057) q[2];
sx q[2];
rz(-2.1446054) q[2];
sx q[2];
rz(-2.9531324) q[2];
rz(-2.878876) q[3];
sx q[3];
rz(-1.1843984) q[3];
sx q[3];
rz(-0.57412) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
