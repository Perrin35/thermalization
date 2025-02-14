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
rz(-0.77684075) q[0];
sx q[0];
rz(-0.87640327) q[0];
sx q[0];
rz(2.8887698) q[0];
rz(1.1480992) q[1];
sx q[1];
rz(4.0568772) q[1];
sx q[1];
rz(8.6203909) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0774121) q[0];
sx q[0];
rz(-1.8837116) q[0];
sx q[0];
rz(1.1699647) q[0];
rz(-pi) q[1];
rz(1.8041274) q[2];
sx q[2];
rz(-1.1510282) q[2];
sx q[2];
rz(0.45316089) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.48015311) q[1];
sx q[1];
rz(-2.6240908) q[1];
sx q[1];
rz(2.3913942) q[1];
x q[2];
rz(0.31320235) q[3];
sx q[3];
rz(-2.7814908) q[3];
sx q[3];
rz(1.2822156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.048451) q[2];
sx q[2];
rz(-0.96258771) q[2];
sx q[2];
rz(-2.2691881) q[2];
rz(2.8060272) q[3];
sx q[3];
rz(-1.7458956) q[3];
sx q[3];
rz(3.057835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5035079) q[0];
sx q[0];
rz(-0.26569772) q[0];
sx q[0];
rz(0.22915325) q[0];
rz(2.9511662) q[1];
sx q[1];
rz(-1.8016022) q[1];
sx q[1];
rz(0.71358877) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3772954) q[0];
sx q[0];
rz(-2.4141001) q[0];
sx q[0];
rz(1.9801674) q[0];
rz(-pi) q[1];
x q[1];
rz(0.53411463) q[2];
sx q[2];
rz(-0.91098173) q[2];
sx q[2];
rz(2.5995863) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.59871626) q[1];
sx q[1];
rz(-2.2440254) q[1];
sx q[1];
rz(1.190541) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.76477625) q[3];
sx q[3];
rz(-1.4524609) q[3];
sx q[3];
rz(-1.8272994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4552292) q[2];
sx q[2];
rz(-0.31193048) q[2];
sx q[2];
rz(2.6658106) q[2];
rz(1.0698211) q[3];
sx q[3];
rz(-1.0082303) q[3];
sx q[3];
rz(2.3750335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0723202) q[0];
sx q[0];
rz(-1.9444436) q[0];
sx q[0];
rz(-1.9092165) q[0];
rz(0.45066372) q[1];
sx q[1];
rz(-1.181299) q[1];
sx q[1];
rz(-2.252069) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4159195) q[0];
sx q[0];
rz(-1.9751722) q[0];
sx q[0];
rz(-2.5017991) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0126152) q[2];
sx q[2];
rz(-2.0348437) q[2];
sx q[2];
rz(-1.2166263) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0625713) q[1];
sx q[1];
rz(-1.3600742) q[1];
sx q[1];
rz(-0.91075588) q[1];
rz(-pi) q[2];
x q[2];
rz(0.74340762) q[3];
sx q[3];
rz(-0.49170676) q[3];
sx q[3];
rz(1.3431026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.687261) q[2];
sx q[2];
rz(-0.4250409) q[2];
sx q[2];
rz(-1.8827776) q[2];
rz(1.2339633) q[3];
sx q[3];
rz(-2.0089269) q[3];
sx q[3];
rz(-3.1335355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(2.4312129) q[0];
sx q[0];
rz(-2.2358535) q[0];
sx q[0];
rz(1.6718965) q[0];
rz(1.1471033) q[1];
sx q[1];
rz(-1.0018145) q[1];
sx q[1];
rz(-0.9009487) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6916136) q[0];
sx q[0];
rz(-1.2581345) q[0];
sx q[0];
rz(-1.481196) q[0];
rz(-pi) q[1];
rz(-2.35975) q[2];
sx q[2];
rz(-2.3580209) q[2];
sx q[2];
rz(-1.9279566) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.21877737) q[1];
sx q[1];
rz(-2.4115857) q[1];
sx q[1];
rz(-0.14132146) q[1];
rz(-pi) q[2];
rz(0.73293682) q[3];
sx q[3];
rz(-1.4519435) q[3];
sx q[3];
rz(-0.87866966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.49742302) q[2];
sx q[2];
rz(-1.3966509) q[2];
sx q[2];
rz(-2.4141342) q[2];
rz(0.71508956) q[3];
sx q[3];
rz(-2.6810724) q[3];
sx q[3];
rz(0.49811825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6775976) q[0];
sx q[0];
rz(-1.3901187) q[0];
sx q[0];
rz(-0.736262) q[0];
rz(-2.5212506) q[1];
sx q[1];
rz(-0.54519975) q[1];
sx q[1];
rz(1.5249407) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.072983532) q[0];
sx q[0];
rz(-1.6260379) q[0];
sx q[0];
rz(1.5806356) q[0];
rz(-2.9102566) q[2];
sx q[2];
rz(-1.8978737) q[2];
sx q[2];
rz(1.2755204) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9774919) q[1];
sx q[1];
rz(-0.79704801) q[1];
sx q[1];
rz(1.0352831) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3518538) q[3];
sx q[3];
rz(-0.92806714) q[3];
sx q[3];
rz(-1.5830884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6058495) q[2];
sx q[2];
rz(-2.3652786) q[2];
sx q[2];
rz(0.95174754) q[2];
rz(-0.4176628) q[3];
sx q[3];
rz(-0.2499191) q[3];
sx q[3];
rz(-1.9865659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8010537) q[0];
sx q[0];
rz(-1.8908353) q[0];
sx q[0];
rz(2.387555) q[0];
rz(0.89318371) q[1];
sx q[1];
rz(-1.040753) q[1];
sx q[1];
rz(2.975614) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6245859) q[0];
sx q[0];
rz(-2.4501194) q[0];
sx q[0];
rz(-0.16384478) q[0];
rz(-pi) q[1];
rz(0.31678172) q[2];
sx q[2];
rz(-2.6133279) q[2];
sx q[2];
rz(-2.1342017) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.41478911) q[1];
sx q[1];
rz(-2.1926452) q[1];
sx q[1];
rz(2.1305103) q[1];
x q[2];
rz(-0.13940553) q[3];
sx q[3];
rz(-1.2981112) q[3];
sx q[3];
rz(-2.8806339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7665427) q[2];
sx q[2];
rz(-2.0039717) q[2];
sx q[2];
rz(3.1362015) q[2];
rz(0.088717669) q[3];
sx q[3];
rz(-1.5327449) q[3];
sx q[3];
rz(0.79571342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2583112) q[0];
sx q[0];
rz(-1.9323876) q[0];
sx q[0];
rz(2.9364371) q[0];
rz(2.2329277) q[1];
sx q[1];
rz(-0.87037194) q[1];
sx q[1];
rz(0.044513449) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9712898) q[0];
sx q[0];
rz(-1.4409541) q[0];
sx q[0];
rz(-1.4128039) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8327863) q[2];
sx q[2];
rz(-0.14523187) q[2];
sx q[2];
rz(-2.0042208) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8537112) q[1];
sx q[1];
rz(-1.3559623) q[1];
sx q[1];
rz(-0.98835215) q[1];
x q[2];
rz(-0.022734108) q[3];
sx q[3];
rz(-1.6904545) q[3];
sx q[3];
rz(-0.59198635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.32555106) q[2];
sx q[2];
rz(-0.52246919) q[2];
sx q[2];
rz(-2.6204056) q[2];
rz(-2.9442287) q[3];
sx q[3];
rz(-1.3857931) q[3];
sx q[3];
rz(0.50281966) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7407783) q[0];
sx q[0];
rz(-2.9476808) q[0];
sx q[0];
rz(2.8863696) q[0];
rz(0.1420282) q[1];
sx q[1];
rz(-1.4165001) q[1];
sx q[1];
rz(-0.35370383) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.557705) q[0];
sx q[0];
rz(-1.9206128) q[0];
sx q[0];
rz(-0.34829287) q[0];
rz(-0.50019294) q[2];
sx q[2];
rz(-1.0076773) q[2];
sx q[2];
rz(1.4266071) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.6981909) q[1];
sx q[1];
rz(-0.56487067) q[1];
sx q[1];
rz(-2.2400212) q[1];
rz(-pi) q[2];
rz(-0.59732893) q[3];
sx q[3];
rz(-0.84395614) q[3];
sx q[3];
rz(2.8387808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3069309) q[2];
sx q[2];
rz(-0.31535172) q[2];
sx q[2];
rz(1.6174512) q[2];
rz(-2.2751685) q[3];
sx q[3];
rz(-1.6774991) q[3];
sx q[3];
rz(2.5518937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50188142) q[0];
sx q[0];
rz(-0.15637936) q[0];
sx q[0];
rz(1.7891275) q[0];
rz(-1.2347219) q[1];
sx q[1];
rz(-2.9094978) q[1];
sx q[1];
rz(0.51188767) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0842132) q[0];
sx q[0];
rz(-1.6239927) q[0];
sx q[0];
rz(-1.4131111) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4134679) q[2];
sx q[2];
rz(-2.3300588) q[2];
sx q[2];
rz(-3.1026156) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3583957) q[1];
sx q[1];
rz(-1.8640717) q[1];
sx q[1];
rz(-0.44289987) q[1];
rz(-2.1549822) q[3];
sx q[3];
rz(-1.8711149) q[3];
sx q[3];
rz(-0.14726135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0074244) q[2];
sx q[2];
rz(-1.3205426) q[2];
sx q[2];
rz(0.39719886) q[2];
rz(2.126179) q[3];
sx q[3];
rz(-2.1404603) q[3];
sx q[3];
rz(0.5526244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39887244) q[0];
sx q[0];
rz(-1.9341368) q[0];
sx q[0];
rz(-0.43750986) q[0];
rz(0.73879009) q[1];
sx q[1];
rz(-0.80016017) q[1];
sx q[1];
rz(-1.4791666) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0453607) q[0];
sx q[0];
rz(-0.49043617) q[0];
sx q[0];
rz(-3.0231721) q[0];
rz(-pi) q[1];
rz(1.9301038) q[2];
sx q[2];
rz(-1.4622765) q[2];
sx q[2];
rz(1.6720275) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.91875118) q[1];
sx q[1];
rz(-1.6770308) q[1];
sx q[1];
rz(1.5083471) q[1];
x q[2];
rz(2.4303834) q[3];
sx q[3];
rz(-0.91806245) q[3];
sx q[3];
rz(-2.0026596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2207569) q[2];
sx q[2];
rz(-1.3584542) q[2];
sx q[2];
rz(2.9662568) q[2];
rz(1.0026503) q[3];
sx q[3];
rz(-0.23313871) q[3];
sx q[3];
rz(-1.233915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47068448) q[0];
sx q[0];
rz(-1.4574454) q[0];
sx q[0];
rz(-1.1048143) q[0];
rz(-2.9910174) q[1];
sx q[1];
rz(-0.87599788) q[1];
sx q[1];
rz(-1.8465975) q[1];
rz(-2.9433123) q[2];
sx q[2];
rz(-1.56326) q[2];
sx q[2];
rz(2.9247912) q[2];
rz(2.8908875) q[3];
sx q[3];
rz(-1.7964994) q[3];
sx q[3];
rz(-2.9356706) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
