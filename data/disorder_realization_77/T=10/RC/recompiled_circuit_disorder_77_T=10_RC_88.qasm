OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.66184008) q[0];
sx q[0];
rz(-0.84364426) q[0];
sx q[0];
rz(0.16790976) q[0];
rz(1.1711988) q[1];
sx q[1];
rz(3.436915) q[1];
sx q[1];
rz(9.480939) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4958772) q[0];
sx q[0];
rz(-1.0951395) q[0];
sx q[0];
rz(-0.23947421) q[0];
rz(2.216823) q[2];
sx q[2];
rz(-1.3999108) q[2];
sx q[2];
rz(-0.31121635) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9819298) q[1];
sx q[1];
rz(-0.61239457) q[1];
sx q[1];
rz(2.0484522) q[1];
rz(2.9763016) q[3];
sx q[3];
rz(-0.2632907) q[3];
sx q[3];
rz(-2.1804682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7636259) q[2];
sx q[2];
rz(-0.28181919) q[2];
sx q[2];
rz(0.4326694) q[2];
rz(1.1928605) q[3];
sx q[3];
rz(-1.2377219) q[3];
sx q[3];
rz(-2.7584934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52779657) q[0];
sx q[0];
rz(-2.6531117) q[0];
sx q[0];
rz(1.312785) q[0];
rz(-2.9361172) q[1];
sx q[1];
rz(-2.165129) q[1];
sx q[1];
rz(1.1516494) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8173556) q[0];
sx q[0];
rz(-1.2756057) q[0];
sx q[0];
rz(-2.2833707) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0552911) q[2];
sx q[2];
rz(-2.0995579) q[2];
sx q[2];
rz(2.8690086) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.18560219) q[1];
sx q[1];
rz(-0.62528175) q[1];
sx q[1];
rz(-2.37466) q[1];
rz(-pi) q[2];
x q[2];
rz(0.63668164) q[3];
sx q[3];
rz(-2.5425306) q[3];
sx q[3];
rz(2.3707795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.0097222086) q[2];
sx q[2];
rz(-1.4902318) q[2];
sx q[2];
rz(-2.9197664) q[2];
rz(-2.7644073) q[3];
sx q[3];
rz(-0.42755869) q[3];
sx q[3];
rz(0.77243531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8310228) q[0];
sx q[0];
rz(-0.092386827) q[0];
sx q[0];
rz(0.036852766) q[0];
rz(2.316078) q[1];
sx q[1];
rz(-1.3157536) q[1];
sx q[1];
rz(0.056578606) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7650334) q[0];
sx q[0];
rz(-2.3925836) q[0];
sx q[0];
rz(0.88699938) q[0];
rz(-pi) q[1];
rz(0.47358863) q[2];
sx q[2];
rz(-2.8550672) q[2];
sx q[2];
rz(-0.63602704) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.11675662) q[1];
sx q[1];
rz(-0.36326888) q[1];
sx q[1];
rz(0.58961745) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9583086) q[3];
sx q[3];
rz(-0.96252493) q[3];
sx q[3];
rz(-1.1063948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.27292192) q[2];
sx q[2];
rz(-1.4977095) q[2];
sx q[2];
rz(0.92612129) q[2];
rz(-2.5849294) q[3];
sx q[3];
rz(-2.848048) q[3];
sx q[3];
rz(-2.0986957) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91519231) q[0];
sx q[0];
rz(-2.4214348) q[0];
sx q[0];
rz(2.3994989) q[0];
rz(1.1391976) q[1];
sx q[1];
rz(-0.4793872) q[1];
sx q[1];
rz(0.46359584) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6620561) q[0];
sx q[0];
rz(-2.2502796) q[0];
sx q[0];
rz(0.38328538) q[0];
rz(-pi) q[1];
rz(2.8088403) q[2];
sx q[2];
rz(-1.6289662) q[2];
sx q[2];
rz(-2.1259049) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.61958085) q[1];
sx q[1];
rz(-1.3995692) q[1];
sx q[1];
rz(-0.79230688) q[1];
rz(-pi) q[2];
rz(2.1257524) q[3];
sx q[3];
rz(-1.98588) q[3];
sx q[3];
rz(-0.89158981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7745557) q[2];
sx q[2];
rz(-0.15731263) q[2];
sx q[2];
rz(3.0920933) q[2];
rz(0.1285304) q[3];
sx q[3];
rz(-1.5934207) q[3];
sx q[3];
rz(0.11894225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13609919) q[0];
sx q[0];
rz(-2.7042784) q[0];
sx q[0];
rz(2.8438925) q[0];
rz(2.659335) q[1];
sx q[1];
rz(-0.75459701) q[1];
sx q[1];
rz(2.1972426) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34236318) q[0];
sx q[0];
rz(-0.21393299) q[0];
sx q[0];
rz(1.7844723) q[0];
rz(1.1733426) q[2];
sx q[2];
rz(-0.83565088) q[2];
sx q[2];
rz(0.022692516) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7625092) q[1];
sx q[1];
rz(-0.06113872) q[1];
sx q[1];
rz(-1.6054543) q[1];
rz(-pi) q[2];
x q[2];
rz(0.43290187) q[3];
sx q[3];
rz(-1.8042943) q[3];
sx q[3];
rz(-1.3732861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.258761) q[2];
sx q[2];
rz(-1.0927039) q[2];
sx q[2];
rz(-3.0333701) q[2];
rz(0.0023068874) q[3];
sx q[3];
rz(-1.5294411) q[3];
sx q[3];
rz(-0.32430696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43679431) q[0];
sx q[0];
rz(-0.44610766) q[0];
sx q[0];
rz(2.5571402) q[0];
rz(-0.8862409) q[1];
sx q[1];
rz(-0.61683547) q[1];
sx q[1];
rz(-3.086673) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7391101) q[0];
sx q[0];
rz(-1.4882898) q[0];
sx q[0];
rz(-1.8399747) q[0];
x q[1];
rz(-1.5149649) q[2];
sx q[2];
rz(-2.8179114) q[2];
sx q[2];
rz(-1.2602381) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1850486) q[1];
sx q[1];
rz(-1.946432) q[1];
sx q[1];
rz(2.5513785) q[1];
x q[2];
rz(0.46338007) q[3];
sx q[3];
rz(-1.0372835) q[3];
sx q[3];
rz(2.2802071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3871258) q[2];
sx q[2];
rz(-3.0209164) q[2];
sx q[2];
rz(-1.0167271) q[2];
rz(-2.5975442) q[3];
sx q[3];
rz(-0.35990158) q[3];
sx q[3];
rz(-1.3325161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.6761557) q[0];
sx q[0];
rz(-2.1570719) q[0];
sx q[0];
rz(0.28453919) q[0];
rz(0.94447213) q[1];
sx q[1];
rz(-1.1962793) q[1];
sx q[1];
rz(2.231266) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.401424) q[0];
sx q[0];
rz(-1.6356042) q[0];
sx q[0];
rz(-0.054697371) q[0];
rz(2.2482713) q[2];
sx q[2];
rz(-0.51572323) q[2];
sx q[2];
rz(0.90781462) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.42156223) q[1];
sx q[1];
rz(-1.7517462) q[1];
sx q[1];
rz(-2.6471495) q[1];
rz(-pi) q[2];
rz(0.99366412) q[3];
sx q[3];
rz(-2.2054407) q[3];
sx q[3];
rz(-1.5821379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3900782) q[2];
sx q[2];
rz(-0.063515924) q[2];
sx q[2];
rz(-0.92203036) q[2];
rz(2.5743124) q[3];
sx q[3];
rz(-1.4138979) q[3];
sx q[3];
rz(1.0197619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2475125) q[0];
sx q[0];
rz(-2.4649354) q[0];
sx q[0];
rz(-0.12938736) q[0];
rz(-2.5091876) q[1];
sx q[1];
rz(-2.1148732) q[1];
sx q[1];
rz(-2.8410889) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8757272) q[0];
sx q[0];
rz(-0.89571307) q[0];
sx q[0];
rz(-2.6595594) q[0];
x q[1];
rz(2.5007162) q[2];
sx q[2];
rz(-2.2705728) q[2];
sx q[2];
rz(1.7388294) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.955519) q[1];
sx q[1];
rz(-1.3961853) q[1];
sx q[1];
rz(-2.7956635) q[1];
x q[2];
rz(2.825533) q[3];
sx q[3];
rz(-0.83871597) q[3];
sx q[3];
rz(0.72533208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.58632103) q[2];
sx q[2];
rz(-2.1466612) q[2];
sx q[2];
rz(-0.78197455) q[2];
rz(-2.590495) q[3];
sx q[3];
rz(-1.7588153) q[3];
sx q[3];
rz(-0.15792318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5732116) q[0];
sx q[0];
rz(-1.9848354) q[0];
sx q[0];
rz(3.0138299) q[0];
rz(0.54221517) q[1];
sx q[1];
rz(-0.95710373) q[1];
sx q[1];
rz(2.382747) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0046783) q[0];
sx q[0];
rz(-1.5537964) q[0];
sx q[0];
rz(-1.9949811) q[0];
rz(2.5326469) q[2];
sx q[2];
rz(-1.7296089) q[2];
sx q[2];
rz(1.2283404) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.89275937) q[1];
sx q[1];
rz(-1.4151238) q[1];
sx q[1];
rz(-2.4397736) q[1];
x q[2];
rz(0.90518732) q[3];
sx q[3];
rz(-1.5822516) q[3];
sx q[3];
rz(0.36597914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0163429) q[2];
sx q[2];
rz(-1.7813851) q[2];
sx q[2];
rz(-0.49003595) q[2];
rz(-1.4222493) q[3];
sx q[3];
rz(-1.9336721) q[3];
sx q[3];
rz(-1.0629883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35995099) q[0];
sx q[0];
rz(-2.521823) q[0];
sx q[0];
rz(-3.066257) q[0];
rz(2.244859) q[1];
sx q[1];
rz(-1.1743841) q[1];
sx q[1];
rz(0.60992253) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9655351) q[0];
sx q[0];
rz(-2.3241204) q[0];
sx q[0];
rz(-0.39359351) q[0];
rz(2.7650325) q[2];
sx q[2];
rz(-1.3137523) q[2];
sx q[2];
rz(0.10274796) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6590609) q[1];
sx q[1];
rz(-2.3606803) q[1];
sx q[1];
rz(0.34300967) q[1];
rz(-pi) q[2];
rz(2.6033953) q[3];
sx q[3];
rz(-0.81848577) q[3];
sx q[3];
rz(-1.324211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.909409) q[2];
sx q[2];
rz(-0.82513088) q[2];
sx q[2];
rz(2.4278736) q[2];
rz(2.7632726) q[3];
sx q[3];
rz(-2.648073) q[3];
sx q[3];
rz(-0.87987125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80355766) q[0];
sx q[0];
rz(-1.150158) q[0];
sx q[0];
rz(-1.5858571) q[0];
rz(-0.65080416) q[1];
sx q[1];
rz(-1.4918068) q[1];
sx q[1];
rz(3.0202958) q[1];
rz(1.4177633) q[2];
sx q[2];
rz(-2.9433708) q[2];
sx q[2];
rz(2.3797258) q[2];
rz(-2.8888632) q[3];
sx q[3];
rz(-1.1128294) q[3];
sx q[3];
rz(2.2313234) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
