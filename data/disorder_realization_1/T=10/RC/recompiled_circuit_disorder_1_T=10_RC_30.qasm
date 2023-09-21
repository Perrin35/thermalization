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
rz(-1.1146201) q[0];
sx q[0];
rz(-0.00014076509) q[0];
rz(1.3340985) q[1];
sx q[1];
rz(4.1058022) q[1];
sx q[1];
rz(10.618187) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1228468) q[0];
sx q[0];
rz(-1.7483286) q[0];
sx q[0];
rz(-1.2794504) q[0];
x q[1];
rz(1.2725856) q[2];
sx q[2];
rz(-1.0422921) q[2];
sx q[2];
rz(2.3117711) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.57230091) q[1];
sx q[1];
rz(-2.3094059) q[1];
sx q[1];
rz(2.4898847) q[1];
rz(-pi) q[2];
rz(-3.0317806) q[3];
sx q[3];
rz(-1.7870652) q[3];
sx q[3];
rz(3.0307378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6821735) q[2];
sx q[2];
rz(-3.1176304) q[2];
sx q[2];
rz(-1.9127282) q[2];
rz(-1.7284283) q[3];
sx q[3];
rz(-2.0404405) q[3];
sx q[3];
rz(1.6536973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5380149) q[0];
sx q[0];
rz(-1.6390272) q[0];
sx q[0];
rz(2.1287825) q[0];
rz(3.1139328) q[1];
sx q[1];
rz(-0.67359567) q[1];
sx q[1];
rz(1.123463) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2259953) q[0];
sx q[0];
rz(-1.628327) q[0];
sx q[0];
rz(1.1094088) q[0];
x q[1];
rz(-0.77983071) q[2];
sx q[2];
rz(-1.5260328) q[2];
sx q[2];
rz(2.0335399) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.23501913) q[1];
sx q[1];
rz(-2.0779013) q[1];
sx q[1];
rz(2.1706235) q[1];
rz(-pi) q[2];
rz(0.65621891) q[3];
sx q[3];
rz(-2.0341633) q[3];
sx q[3];
rz(3.0955293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3479487) q[2];
sx q[2];
rz(-2.0517893) q[2];
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
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27750257) q[0];
sx q[0];
rz(-0.16177495) q[0];
sx q[0];
rz(1.8664237) q[0];
rz(-2.4480942) q[1];
sx q[1];
rz(-1.2561412) q[1];
sx q[1];
rz(-2.0085874) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61921652) q[0];
sx q[0];
rz(-1.4285061) q[0];
sx q[0];
rz(-0.0033358047) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6373367) q[2];
sx q[2];
rz(-2.2849053) q[2];
sx q[2];
rz(-2.2222663) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9085711) q[1];
sx q[1];
rz(-1.034522) q[1];
sx q[1];
rz(-0.33645333) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1902309) q[3];
sx q[3];
rz(-1.0361443) q[3];
sx q[3];
rz(1.149328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.2514078) q[2];
sx q[2];
rz(-2.3501985) q[2];
sx q[2];
rz(-1.8481002) q[2];
rz(3.1022762) q[3];
sx q[3];
rz(-1.2189564) q[3];
sx q[3];
rz(-1.2600651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2599729) q[0];
sx q[0];
rz(-3.0631174) q[0];
sx q[0];
rz(1.9807293) q[0];
rz(0.89598957) q[1];
sx q[1];
rz(-1.7005824) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-3.0669078) q[2];
sx q[2];
rz(-1.1876145) q[2];
sx q[2];
rz(0.89786868) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.23952661) q[1];
sx q[1];
rz(-2.8956928) q[1];
sx q[1];
rz(0.86218254) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0566063) q[3];
sx q[3];
rz(-2.8286472) q[3];
sx q[3];
rz(1.7900975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.23665145) q[2];
sx q[2];
rz(-0.94649482) q[2];
sx q[2];
rz(2.2616852) q[2];
rz(0.044163477) q[3];
sx q[3];
rz(-1.6396089) q[3];
sx q[3];
rz(2.8529609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1039466) q[0];
sx q[0];
rz(-0.3750616) q[0];
sx q[0];
rz(-1.0132382) q[0];
rz(-0.049731072) q[1];
sx q[1];
rz(-0.91369349) q[1];
sx q[1];
rz(2.0577046) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0439107) q[0];
sx q[0];
rz(-1.327276) q[0];
sx q[0];
rz(2.9527412) q[0];
rz(-pi) q[1];
rz(-0.82917825) q[2];
sx q[2];
rz(-2.7784756) q[2];
sx q[2];
rz(-1.6104289) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.10478445) q[1];
sx q[1];
rz(-1.0462865) q[1];
sx q[1];
rz(0.12496897) q[1];
rz(-pi) q[2];
x q[2];
rz(0.089133457) q[3];
sx q[3];
rz(-2.6260178) q[3];
sx q[3];
rz(2.6775132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.23285398) q[2];
sx q[2];
rz(-2.8149657) q[2];
sx q[2];
rz(0.24442913) q[2];
rz(2.7092253) q[3];
sx q[3];
rz(-1.7418539) q[3];
sx q[3];
rz(-0.50306815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2844834) q[0];
sx q[0];
rz(-1.4215707) q[0];
sx q[0];
rz(3.0474512) q[0];
rz(2.969818) q[1];
sx q[1];
rz(-2.005902) q[1];
sx q[1];
rz(0.89541268) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58127922) q[0];
sx q[0];
rz(-2.7972097) q[0];
sx q[0];
rz(-3.0292105) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63352872) q[2];
sx q[2];
rz(-1.765536) q[2];
sx q[2];
rz(-0.63846904) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.15570607) q[1];
sx q[1];
rz(-2.967318) q[1];
sx q[1];
rz(1.3143015) q[1];
rz(-pi) q[2];
rz(1.5395245) q[3];
sx q[3];
rz(-3.0668689) q[3];
sx q[3];
rz(0.63262313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.133693) q[2];
sx q[2];
rz(-0.40863016) q[2];
sx q[2];
rz(0.80319476) q[2];
rz(-1.1903654) q[3];
sx q[3];
rz(-1.2322216) q[3];
sx q[3];
rz(-2.7289594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(3.0728834) q[0];
sx q[0];
rz(-2.9769653) q[0];
sx q[0];
rz(2.6224526) q[0];
rz(-2.5601162) q[1];
sx q[1];
rz(-2.0362208) q[1];
sx q[1];
rz(-1.8849467) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089766895) q[0];
sx q[0];
rz(-1.8221812) q[0];
sx q[0];
rz(0.041640394) q[0];
rz(-pi) q[1];
rz(2.9314552) q[2];
sx q[2];
rz(-1.1147095) q[2];
sx q[2];
rz(-1.8535341) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9650967) q[1];
sx q[1];
rz(-2.8526222) q[1];
sx q[1];
rz(-2.9138336) q[1];
x q[2];
rz(-1.890896) q[3];
sx q[3];
rz(-2.2543636) q[3];
sx q[3];
rz(2.2724831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.98823035) q[2];
sx q[2];
rz(-2.1116657) q[2];
sx q[2];
rz(-1.777565) q[2];
rz(0.91056943) q[3];
sx q[3];
rz(-1.1547337) q[3];
sx q[3];
rz(1.6114657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0751188) q[0];
sx q[0];
rz(-0.56448889) q[0];
sx q[0];
rz(-0.30817729) q[0];
rz(0.072487436) q[1];
sx q[1];
rz(-2.1283573) q[1];
sx q[1];
rz(0.38696188) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9235619) q[0];
sx q[0];
rz(-2.1523928) q[0];
sx q[0];
rz(-0.27115718) q[0];
rz(-pi) q[1];
rz(-2.1808124) q[2];
sx q[2];
rz(-2.3443601) q[2];
sx q[2];
rz(1.0459935) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3124233) q[1];
sx q[1];
rz(-2.2409391) q[1];
sx q[1];
rz(1.6664684) q[1];
x q[2];
rz(2.6450063) q[3];
sx q[3];
rz(-1.3771309) q[3];
sx q[3];
rz(0.95369875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5179634) q[2];
sx q[2];
rz(-1.7771746) q[2];
sx q[2];
rz(-0.44000885) q[2];
rz(-0.7157588) q[3];
sx q[3];
rz(-1.4322759) q[3];
sx q[3];
rz(2.0619152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0004262) q[0];
sx q[0];
rz(-0.74581242) q[0];
sx q[0];
rz(2.0429042) q[0];
rz(2.4138342) q[1];
sx q[1];
rz(-2.7658503) q[1];
sx q[1];
rz(-0.049302014) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44137529) q[0];
sx q[0];
rz(-2.2451631) q[0];
sx q[0];
rz(0.81654878) q[0];
rz(-pi) q[1];
x q[1];
rz(0.88843139) q[2];
sx q[2];
rz(-1.1738136) q[2];
sx q[2];
rz(1.9156485) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.8763435) q[1];
sx q[1];
rz(-1.1117522) q[1];
sx q[1];
rz(0.52258073) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6155869) q[3];
sx q[3];
rz(-0.81323871) q[3];
sx q[3];
rz(0.40532743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0631642) q[2];
sx q[2];
rz(-0.59946632) q[2];
sx q[2];
rz(0.72193974) q[2];
rz(0.94349629) q[3];
sx q[3];
rz(-0.75073457) q[3];
sx q[3];
rz(0.25434428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9938875) q[0];
sx q[0];
rz(-1.9839956) q[0];
sx q[0];
rz(1.0797427) q[0];
rz(2.0823157) q[1];
sx q[1];
rz(-0.22288999) q[1];
sx q[1];
rz(1.4019029) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3831543) q[0];
sx q[0];
rz(-0.24807319) q[0];
sx q[0];
rz(1.0323348) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5980814) q[2];
sx q[2];
rz(-1.7575043) q[2];
sx q[2];
rz(2.7207431) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8733752) q[1];
sx q[1];
rz(-2.2443319) q[1];
sx q[1];
rz(2.7482277) q[1];
rz(2.8614282) q[3];
sx q[3];
rz(-2.4912842) q[3];
sx q[3];
rz(-0.10661099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6587276) q[2];
sx q[2];
rz(-1.3663224) q[2];
sx q[2];
rz(1.520291) q[2];
rz(-0.55082095) q[3];
sx q[3];
rz(-0.80533022) q[3];
sx q[3];
rz(-2.4441392) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.993492) q[0];
sx q[0];
rz(-1.3052595) q[0];
sx q[0];
rz(-1.530151) q[0];
rz(2.2254754) q[1];
sx q[1];
rz(-2.5506908) q[1];
sx q[1];
rz(2.5509902) q[1];
rz(-0.79402906) q[2];
sx q[2];
rz(-0.90894884) q[2];
sx q[2];
rz(-0.85469645) q[2];
rz(-1.5394474) q[3];
sx q[3];
rz(-2.6242704) q[3];
sx q[3];
rz(-0.27192413) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];