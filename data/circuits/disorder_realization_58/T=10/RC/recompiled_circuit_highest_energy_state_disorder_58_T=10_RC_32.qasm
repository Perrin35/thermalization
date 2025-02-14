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
rz(-0.11197055) q[0];
sx q[0];
rz(-2.114871) q[0];
sx q[0];
rz(1.8486899) q[0];
rz(-1.0678043) q[1];
sx q[1];
rz(-0.81967241) q[1];
sx q[1];
rz(-1.9424865) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1333251) q[0];
sx q[0];
rz(-0.12524167) q[0];
sx q[0];
rz(-0.52405806) q[0];
rz(-pi) q[1];
x q[1];
rz(0.22817518) q[2];
sx q[2];
rz(-0.68636218) q[2];
sx q[2];
rz(2.82968) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1242697) q[1];
sx q[1];
rz(-0.60397479) q[1];
sx q[1];
rz(1.9159957) q[1];
x q[2];
rz(-2.5254123) q[3];
sx q[3];
rz(-0.95762223) q[3];
sx q[3];
rz(-2.4626436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6370411) q[2];
sx q[2];
rz(-1.4013638) q[2];
sx q[2];
rz(-1.441628) q[2];
rz(1.0124892) q[3];
sx q[3];
rz(-1.1463405) q[3];
sx q[3];
rz(0.47154021) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83990324) q[0];
sx q[0];
rz(-0.56115264) q[0];
sx q[0];
rz(-2.1121693) q[0];
rz(1.7816211) q[1];
sx q[1];
rz(-1.1861035) q[1];
sx q[1];
rz(-2.634868) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27924867) q[0];
sx q[0];
rz(-1.5900668) q[0];
sx q[0];
rz(-1.6750245) q[0];
x q[1];
rz(0.36594351) q[2];
sx q[2];
rz(-2.2323713) q[2];
sx q[2];
rz(-0.55398527) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4480865) q[1];
sx q[1];
rz(-0.58219456) q[1];
sx q[1];
rz(1.4513998) q[1];
rz(-pi) q[2];
rz(-0.84552879) q[3];
sx q[3];
rz(-1.3432028) q[3];
sx q[3];
rz(-1.9967772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8257251) q[2];
sx q[2];
rz(-2.5434912) q[2];
sx q[2];
rz(0.21710795) q[2];
rz(-1.0202967) q[3];
sx q[3];
rz(-1.812457) q[3];
sx q[3];
rz(-2.1609658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8751136) q[0];
sx q[0];
rz(-1.3631692) q[0];
sx q[0];
rz(-0.76817051) q[0];
rz(0.29113302) q[1];
sx q[1];
rz(-1.7765287) q[1];
sx q[1];
rz(1.1870144) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6757557) q[0];
sx q[0];
rz(-1.4872941) q[0];
sx q[0];
rz(-2.9299987) q[0];
rz(-1.6441543) q[2];
sx q[2];
rz(-1.6147476) q[2];
sx q[2];
rz(0.32273656) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.62115034) q[1];
sx q[1];
rz(-1.598473) q[1];
sx q[1];
rz(2.7872183) q[1];
rz(-0.94908917) q[3];
sx q[3];
rz(-1.9636077) q[3];
sx q[3];
rz(-1.2248685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8698296) q[2];
sx q[2];
rz(-0.84438476) q[2];
sx q[2];
rz(-2.1738906) q[2];
rz(-2.9076231) q[3];
sx q[3];
rz(-0.7015737) q[3];
sx q[3];
rz(-2.7847737) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1135603) q[0];
sx q[0];
rz(-1.5505294) q[0];
sx q[0];
rz(2.0618942) q[0];
rz(-0.30403852) q[1];
sx q[1];
rz(-1.1540776) q[1];
sx q[1];
rz(-1.8255723) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9213778) q[0];
sx q[0];
rz(-2.1002033) q[0];
sx q[0];
rz(1.7142943) q[0];
x q[1];
rz(1.448668) q[2];
sx q[2];
rz(-2.5425362) q[2];
sx q[2];
rz(2.0640822) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6027863) q[1];
sx q[1];
rz(-2.1486308) q[1];
sx q[1];
rz(-2.2972107) q[1];
x q[2];
rz(1.1445266) q[3];
sx q[3];
rz(-1.5139587) q[3];
sx q[3];
rz(2.6714966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6198081) q[2];
sx q[2];
rz(-1.2336171) q[2];
sx q[2];
rz(-0.141315) q[2];
rz(-0.58051336) q[3];
sx q[3];
rz(-1.6655191) q[3];
sx q[3];
rz(-0.22835246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2572131) q[0];
sx q[0];
rz(-2.3452106) q[0];
sx q[0];
rz(-1.4759395) q[0];
rz(0.93302226) q[1];
sx q[1];
rz(-2.4304183) q[1];
sx q[1];
rz(-1.3915541) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6835502) q[0];
sx q[0];
rz(-1.453896) q[0];
sx q[0];
rz(0.83302541) q[0];
x q[1];
rz(2.8093908) q[2];
sx q[2];
rz(-2.6966288) q[2];
sx q[2];
rz(1.6138168) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7793624) q[1];
sx q[1];
rz(-1.0500776) q[1];
sx q[1];
rz(-3.1073276) q[1];
rz(-pi) q[2];
rz(0.33008195) q[3];
sx q[3];
rz(-2.0795855) q[3];
sx q[3];
rz(-1.1710081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.89620227) q[2];
sx q[2];
rz(-2.7707477) q[2];
sx q[2];
rz(-2.89213) q[2];
rz(-1.4977411) q[3];
sx q[3];
rz(-1.7353053) q[3];
sx q[3];
rz(-2.7264285) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45667085) q[0];
sx q[0];
rz(-2.6455854) q[0];
sx q[0];
rz(-1.3859092) q[0];
rz(-1.2617525) q[1];
sx q[1];
rz(-1.933814) q[1];
sx q[1];
rz(-2.6127846) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.945707) q[0];
sx q[0];
rz(-1.56371) q[0];
sx q[0];
rz(-1.12557) q[0];
rz(2.4086921) q[2];
sx q[2];
rz(-2.356988) q[2];
sx q[2];
rz(0.32151383) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0046407) q[1];
sx q[1];
rz(-1.6378504) q[1];
sx q[1];
rz(-2.2100558) q[1];
x q[2];
rz(-0.15915079) q[3];
sx q[3];
rz(-2.4503631) q[3];
sx q[3];
rz(-2.7214512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.95278198) q[2];
sx q[2];
rz(-0.47241259) q[2];
sx q[2];
rz(-1.3713651) q[2];
rz(-2.93907) q[3];
sx q[3];
rz(-1.4684418) q[3];
sx q[3];
rz(1.1892345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8093569) q[0];
sx q[0];
rz(-1.7928596) q[0];
sx q[0];
rz(-0.15383823) q[0];
rz(2.4065252) q[1];
sx q[1];
rz(-2.0497597) q[1];
sx q[1];
rz(0.14911266) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47220889) q[0];
sx q[0];
rz(-2.4824004) q[0];
sx q[0];
rz(1.7071595) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9682069) q[2];
sx q[2];
rz(-1.2116222) q[2];
sx q[2];
rz(-2.2196291) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9184409) q[1];
sx q[1];
rz(-2.0479124) q[1];
sx q[1];
rz(-2.3972307) q[1];
x q[2];
rz(1.7343765) q[3];
sx q[3];
rz(-0.50967607) q[3];
sx q[3];
rz(-1.0496751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4601124) q[2];
sx q[2];
rz(-1.511938) q[2];
sx q[2];
rz(-0.45664772) q[2];
rz(-1.418142) q[3];
sx q[3];
rz(-2.2599615) q[3];
sx q[3];
rz(-0.68896967) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5099455) q[0];
sx q[0];
rz(-0.25196415) q[0];
sx q[0];
rz(3.0889567) q[0];
rz(-1.3098199) q[1];
sx q[1];
rz(-2.8927264) q[1];
sx q[1];
rz(1.2059258) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.238097) q[0];
sx q[0];
rz(-0.66199979) q[0];
sx q[0];
rz(-2.2308626) q[0];
x q[1];
rz(-0.56059273) q[2];
sx q[2];
rz(-1.3729323) q[2];
sx q[2];
rz(2.3677111) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5284023) q[1];
sx q[1];
rz(-2.2397016) q[1];
sx q[1];
rz(1.7714798) q[1];
x q[2];
rz(-1.8926335) q[3];
sx q[3];
rz(-0.6475237) q[3];
sx q[3];
rz(0.81784883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.84411088) q[2];
sx q[2];
rz(-1.7606807) q[2];
sx q[2];
rz(-2.5816176) q[2];
rz(1.0567793) q[3];
sx q[3];
rz(-0.70928514) q[3];
sx q[3];
rz(2.3287676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5071204) q[0];
sx q[0];
rz(-1.7422603) q[0];
sx q[0];
rz(-1.5768453) q[0];
rz(1.5722081) q[1];
sx q[1];
rz(-1.1610616) q[1];
sx q[1];
rz(-2.3326468) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2574111) q[0];
sx q[0];
rz(-2.0657259) q[0];
sx q[0];
rz(-1.5280426) q[0];
rz(0.81124146) q[2];
sx q[2];
rz(-1.277207) q[2];
sx q[2];
rz(-2.9305803) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8508243) q[1];
sx q[1];
rz(-1.5838669) q[1];
sx q[1];
rz(2.1816861) q[1];
rz(1.21502) q[3];
sx q[3];
rz(-2.9529533) q[3];
sx q[3];
rz(-1.1805746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.481015) q[2];
sx q[2];
rz(-0.95260859) q[2];
sx q[2];
rz(3.0926404) q[2];
rz(1.8598716) q[3];
sx q[3];
rz(-0.47544026) q[3];
sx q[3];
rz(-1.4648645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82185164) q[0];
sx q[0];
rz(-2.8806683) q[0];
sx q[0];
rz(1.9191746) q[0];
rz(-1.6659196) q[1];
sx q[1];
rz(-1.9404575) q[1];
sx q[1];
rz(-2.6840721) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2777785) q[0];
sx q[0];
rz(-0.84054986) q[0];
sx q[0];
rz(-1.449683) q[0];
x q[1];
rz(-2.7399212) q[2];
sx q[2];
rz(-2.0435512) q[2];
sx q[2];
rz(-0.061701802) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.88120115) q[1];
sx q[1];
rz(-0.40801755) q[1];
sx q[1];
rz(0.53149207) q[1];
rz(-0.68660929) q[3];
sx q[3];
rz(-1.2401738) q[3];
sx q[3];
rz(-2.7329684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7484625) q[2];
sx q[2];
rz(-0.25105432) q[2];
sx q[2];
rz(-1.040323) q[2];
rz(2.6535502) q[3];
sx q[3];
rz(-1.7010242) q[3];
sx q[3];
rz(-2.603781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3671065) q[0];
sx q[0];
rz(-2.0699061) q[0];
sx q[0];
rz(0.62181428) q[0];
rz(1.8472916) q[1];
sx q[1];
rz(-0.58302561) q[1];
sx q[1];
rz(2.2517712) q[1];
rz(-1.1414083) q[2];
sx q[2];
rz(-0.63776897) q[2];
sx q[2];
rz(-2.706131) q[2];
rz(-0.46066416) q[3];
sx q[3];
rz(-2.3973) q[3];
sx q[3];
rz(2.2826791) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
