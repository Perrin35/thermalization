OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.82436615) q[0];
sx q[0];
rz(5.1685652) q[0];
sx q[0];
rz(9.4246372) q[0];
rz(1.3340985) q[1];
sx q[1];
rz(-2.1773832) q[1];
sx q[1];
rz(1.1934086) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50085917) q[0];
sx q[0];
rz(-1.2841604) q[0];
sx q[0];
rz(-2.9564234) q[0];
rz(-pi) q[1];
rz(-2.5932181) q[2];
sx q[2];
rz(-1.8273241) q[2];
sx q[2];
rz(-0.89474364) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5692917) q[1];
sx q[1];
rz(-2.3094059) q[1];
sx q[1];
rz(2.4898847) q[1];
x q[2];
rz(3.0317806) q[3];
sx q[3];
rz(-1.7870652) q[3];
sx q[3];
rz(-3.0307378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6821735) q[2];
sx q[2];
rz(-3.1176304) q[2];
sx q[2];
rz(-1.2288644) q[2];
rz(-1.4131644) q[3];
sx q[3];
rz(-2.0404405) q[3];
sx q[3];
rz(-1.6536973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6035778) q[0];
sx q[0];
rz(-1.6390272) q[0];
sx q[0];
rz(2.1287825) q[0];
rz(-3.1139328) q[1];
sx q[1];
rz(-0.67359567) q[1];
sx q[1];
rz(-1.123463) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2259953) q[0];
sx q[0];
rz(-1.5132656) q[0];
sx q[0];
rz(1.1094088) q[0];
rz(2.3617619) q[2];
sx q[2];
rz(-1.5260328) q[2];
sx q[2];
rz(-1.1080527) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9065735) q[1];
sx q[1];
rz(-2.0779013) q[1];
sx q[1];
rz(-0.97096918) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4853737) q[3];
sx q[3];
rz(-1.1074293) q[3];
sx q[3];
rz(-0.046063395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3479487) q[2];
sx q[2];
rz(-1.0898033) q[2];
sx q[2];
rz(-2.222555) q[2];
rz(0.67409003) q[3];
sx q[3];
rz(-2.489311) q[3];
sx q[3];
rz(1.526171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8640901) q[0];
sx q[0];
rz(-0.16177495) q[0];
sx q[0];
rz(1.8664237) q[0];
rz(-2.4480942) q[1];
sx q[1];
rz(-1.8854515) q[1];
sx q[1];
rz(2.0085874) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4988574) q[0];
sx q[0];
rz(-2.9992636) q[0];
sx q[0];
rz(-1.5475153) q[0];
rz(-pi) q[1];
rz(2.3511231) q[2];
sx q[2];
rz(-1.1970453) q[2];
sx q[2];
rz(-0.30465301) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9085711) q[1];
sx q[1];
rz(-1.034522) q[1];
sx q[1];
rz(-0.33645333) q[1];
x q[2];
rz(2.3660925) q[3];
sx q[3];
rz(-2.3470078) q[3];
sx q[3];
rz(2.9426108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8901849) q[2];
sx q[2];
rz(-2.3501985) q[2];
sx q[2];
rz(1.8481002) q[2];
rz(-3.1022762) q[3];
sx q[3];
rz(-1.9226363) q[3];
sx q[3];
rz(-1.2600651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2599729) q[0];
sx q[0];
rz(-3.0631174) q[0];
sx q[0];
rz(1.9807293) q[0];
rz(2.2456031) q[1];
sx q[1];
rz(-1.4410102) q[1];
sx q[1];
rz(-0.13555759) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7318152) q[0];
sx q[0];
rz(-1.8827794) q[0];
sx q[0];
rz(0.72809763) q[0];
x q[1];
rz(0.074684871) q[2];
sx q[2];
rz(-1.1876145) q[2];
sx q[2];
rz(0.89786868) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.63771026) q[1];
sx q[1];
rz(-1.4117068) q[1];
sx q[1];
rz(-1.3824944) q[1];
rz(1.2961779) q[3];
sx q[3];
rz(-1.7227968) q[3];
sx q[3];
rz(-0.27384588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9049412) q[2];
sx q[2];
rz(-2.1950978) q[2];
sx q[2];
rz(2.2616852) q[2];
rz(-3.0974292) q[3];
sx q[3];
rz(-1.6396089) q[3];
sx q[3];
rz(2.8529609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
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
rz(-1.0376461) q[0];
sx q[0];
rz(-0.3750616) q[0];
sx q[0];
rz(-2.1283545) q[0];
rz(-3.0918616) q[1];
sx q[1];
rz(-0.91369349) q[1];
sx q[1];
rz(1.0838881) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0439107) q[0];
sx q[0];
rz(-1.327276) q[0];
sx q[0];
rz(2.9527412) q[0];
x q[1];
rz(-0.82917825) q[2];
sx q[2];
rz(-0.36311705) q[2];
sx q[2];
rz(-1.5311637) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0005972) q[1];
sx q[1];
rz(-0.53783572) q[1];
sx q[1];
rz(-1.3586033) q[1];
rz(-pi) q[2];
rz(-2.6277222) q[3];
sx q[3];
rz(-1.6146982) q[3];
sx q[3];
rz(1.1843137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9087387) q[2];
sx q[2];
rz(-0.32662699) q[2];
sx q[2];
rz(2.8971635) q[2];
rz(2.7092253) q[3];
sx q[3];
rz(-1.3997388) q[3];
sx q[3];
rz(-2.6385245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8571092) q[0];
sx q[0];
rz(-1.720022) q[0];
sx q[0];
rz(3.0474512) q[0];
rz(-2.969818) q[1];
sx q[1];
rz(-2.005902) q[1];
sx q[1];
rz(-0.89541268) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58127922) q[0];
sx q[0];
rz(-2.7972097) q[0];
sx q[0];
rz(3.0292105) q[0];
x q[1];
rz(-1.8108098) q[2];
sx q[2];
rz(-0.95108205) q[2];
sx q[2];
rz(2.0680708) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4737282) q[1];
sx q[1];
rz(-1.6147991) q[1];
sx q[1];
rz(-1.7394789) q[1];
x q[2];
rz(3.1392519) q[3];
sx q[3];
rz(-1.6454835) q[3];
sx q[3];
rz(2.5403288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.133693) q[2];
sx q[2];
rz(-0.40863016) q[2];
sx q[2];
rz(-0.80319476) q[2];
rz(-1.9512272) q[3];
sx q[3];
rz(-1.2322216) q[3];
sx q[3];
rz(2.7289594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.068709277) q[0];
sx q[0];
rz(-0.16462737) q[0];
sx q[0];
rz(2.6224526) q[0];
rz(0.58147645) q[1];
sx q[1];
rz(-2.0362208) q[1];
sx q[1];
rz(1.2566459) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0518258) q[0];
sx q[0];
rz(-1.8221812) q[0];
sx q[0];
rz(-3.0999523) q[0];
rz(-pi) q[1];
rz(1.9728327) q[2];
sx q[2];
rz(-0.49905825) q[2];
sx q[2];
rz(-1.4025584) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.17649594) q[1];
sx q[1];
rz(-0.28897044) q[1];
sx q[1];
rz(-0.22775905) q[1];
rz(-1.890896) q[3];
sx q[3];
rz(-0.88722908) q[3];
sx q[3];
rz(-2.2724831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1533623) q[2];
sx q[2];
rz(-2.1116657) q[2];
sx q[2];
rz(1.3640277) q[2];
rz(2.2310232) q[3];
sx q[3];
rz(-1.986859) q[3];
sx q[3];
rz(1.6114657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0751188) q[0];
sx q[0];
rz(-2.5771038) q[0];
sx q[0];
rz(-0.30817729) q[0];
rz(-0.072487436) q[1];
sx q[1];
rz(-1.0132353) q[1];
sx q[1];
rz(-2.7546308) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9403801) q[0];
sx q[0];
rz(-1.7965172) q[0];
sx q[0];
rz(2.1696027) q[0];
x q[1];
rz(0.87256356) q[2];
sx q[2];
rz(-1.1485032) q[2];
sx q[2];
rz(0.97908212) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8291694) q[1];
sx q[1];
rz(-0.90065354) q[1];
sx q[1];
rz(-1.6664684) q[1];
rz(-pi) q[2];
rz(0.39051315) q[3];
sx q[3];
rz(-2.6115341) q[3];
sx q[3];
rz(-0.9583677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.62362921) q[2];
sx q[2];
rz(-1.7771746) q[2];
sx q[2];
rz(2.7015838) q[2];
rz(-0.7157588) q[3];
sx q[3];
rz(-1.4322759) q[3];
sx q[3];
rz(-1.0796775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(3.0004262) q[0];
sx q[0];
rz(-0.74581242) q[0];
sx q[0];
rz(-1.0986885) q[0];
rz(2.4138342) q[1];
sx q[1];
rz(-0.37574238) q[1];
sx q[1];
rz(-3.0922906) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44137529) q[0];
sx q[0];
rz(-2.2451631) q[0];
sx q[0];
rz(-0.81654878) q[0];
rz(-pi) q[1];
rz(2.6463037) q[2];
sx q[2];
rz(-2.1914748) q[2];
sx q[2];
rz(-2.4923313) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.8763435) q[1];
sx q[1];
rz(-1.1117522) q[1];
sx q[1];
rz(-2.6190119) q[1];
rz(-pi) q[2];
rz(-0.52600577) q[3];
sx q[3];
rz(-0.81323871) q[3];
sx q[3];
rz(-0.40532743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0631642) q[2];
sx q[2];
rz(-2.5421263) q[2];
sx q[2];
rz(-0.72193974) q[2];
rz(-0.94349629) q[3];
sx q[3];
rz(-2.3908581) q[3];
sx q[3];
rz(-2.8872484) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9938875) q[0];
sx q[0];
rz(-1.1575971) q[0];
sx q[0];
rz(-2.06185) q[0];
rz(1.059277) q[1];
sx q[1];
rz(-0.22288999) q[1];
sx q[1];
rz(-1.4019029) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83090529) q[0];
sx q[0];
rz(-1.3584104) q[0];
sx q[0];
rz(3.0124245) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.18677588) q[2];
sx q[2];
rz(-1.5439856) q[2];
sx q[2];
rz(-1.1550127) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2682174) q[1];
sx q[1];
rz(-2.2443319) q[1];
sx q[1];
rz(2.7482277) q[1];
rz(-pi) q[2];
x q[2];
rz(2.510342) q[3];
sx q[3];
rz(-1.4025941) q[3];
sx q[3];
rz(-1.2390618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.4828651) q[2];
sx q[2];
rz(-1.7752703) q[2];
sx q[2];
rz(1.520291) q[2];
rz(0.55082095) q[3];
sx q[3];
rz(-2.3362624) q[3];
sx q[3];
rz(0.6974535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.993492) q[0];
sx q[0];
rz(-1.8363331) q[0];
sx q[0];
rz(1.6114417) q[0];
rz(2.2254754) q[1];
sx q[1];
rz(-2.5506908) q[1];
sx q[1];
rz(2.5509902) q[1];
rz(0.82952164) q[2];
sx q[2];
rz(-2.1567232) q[2];
sx q[2];
rz(0.17270252) q[2];
rz(1.5394474) q[3];
sx q[3];
rz(-0.51732224) q[3];
sx q[3];
rz(2.8696685) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
