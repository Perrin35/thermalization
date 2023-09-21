OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4797526) q[0];
sx q[0];
rz(-2.2979484) q[0];
sx q[0];
rz(2.9736829) q[0];
rz(1.1711988) q[1];
sx q[1];
rz(3.436915) q[1];
sx q[1];
rz(9.480939) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6457155) q[0];
sx q[0];
rz(-2.0464532) q[0];
sx q[0];
rz(-0.23947421) q[0];
x q[1];
rz(0.92476966) q[2];
sx q[2];
rz(-1.7416818) q[2];
sx q[2];
rz(2.8303763) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.72370428) q[1];
sx q[1];
rz(-1.0350409) q[1];
sx q[1];
rz(0.31236155) q[1];
x q[2];
rz(-2.9763016) q[3];
sx q[3];
rz(-0.2632907) q[3];
sx q[3];
rz(-0.96112448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.37796676) q[2];
sx q[2];
rz(-0.28181919) q[2];
sx q[2];
rz(0.4326694) q[2];
rz(-1.1928605) q[3];
sx q[3];
rz(-1.2377219) q[3];
sx q[3];
rz(-0.38309923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52779657) q[0];
sx q[0];
rz(-0.48848099) q[0];
sx q[0];
rz(-1.8288076) q[0];
rz(-0.20547543) q[1];
sx q[1];
rz(-0.97646362) q[1];
sx q[1];
rz(1.1516494) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8173556) q[0];
sx q[0];
rz(-1.2756057) q[0];
sx q[0];
rz(0.8582219) q[0];
rz(-2.4685681) q[2];
sx q[2];
rz(-2.4403799) q[2];
sx q[2];
rz(-0.53403026) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.18560219) q[1];
sx q[1];
rz(-0.62528175) q[1];
sx q[1];
rz(-2.37466) q[1];
rz(-pi) q[2];
rz(-1.9563975) q[3];
sx q[3];
rz(-1.10023) q[3];
sx q[3];
rz(0.040599559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1318704) q[2];
sx q[2];
rz(-1.6513609) q[2];
sx q[2];
rz(-2.9197664) q[2];
rz(2.7644073) q[3];
sx q[3];
rz(-0.42755869) q[3];
sx q[3];
rz(2.3691573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8310228) q[0];
sx q[0];
rz(-0.092386827) q[0];
sx q[0];
rz(0.036852766) q[0];
rz(-0.82551461) q[1];
sx q[1];
rz(-1.8258391) q[1];
sx q[1];
rz(-0.056578606) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3765592) q[0];
sx q[0];
rz(-2.3925836) q[0];
sx q[0];
rz(-0.88699938) q[0];
rz(-pi) q[1];
rz(-0.47358863) q[2];
sx q[2];
rz(-0.28652546) q[2];
sx q[2];
rz(2.5055656) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6371582) q[1];
sx q[1];
rz(-1.8706026) q[1];
sx q[1];
rz(1.3624886) q[1];
rz(-0.49719663) q[3];
sx q[3];
rz(-0.70780863) q[3];
sx q[3];
rz(-1.7266112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8686707) q[2];
sx q[2];
rz(-1.6438831) q[2];
sx q[2];
rz(-2.2154714) q[2];
rz(0.55666322) q[3];
sx q[3];
rz(-0.29354468) q[3];
sx q[3];
rz(-1.042897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2264003) q[0];
sx q[0];
rz(-2.4214348) q[0];
sx q[0];
rz(-0.74209374) q[0];
rz(-2.0023951) q[1];
sx q[1];
rz(-2.6622055) q[1];
sx q[1];
rz(-0.46359584) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84309834) q[0];
sx q[0];
rz(-1.8659741) q[0];
sx q[0];
rz(-0.85423268) q[0];
rz(-pi) q[1];
x q[1];
rz(0.17642994) q[2];
sx q[2];
rz(-2.8039805) q[2];
sx q[2];
rz(-2.7531429) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5220118) q[1];
sx q[1];
rz(-1.7420235) q[1];
sx q[1];
rz(2.3492858) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6632518) q[3];
sx q[3];
rz(-2.0739177) q[3];
sx q[3];
rz(-0.92418811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7745557) q[2];
sx q[2];
rz(-0.15731263) q[2];
sx q[2];
rz(-3.0920933) q[2];
rz(-0.1285304) q[3];
sx q[3];
rz(-1.5481719) q[3];
sx q[3];
rz(-3.0226504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0054935) q[0];
sx q[0];
rz(-0.43731421) q[0];
sx q[0];
rz(0.29770011) q[0];
rz(2.659335) q[1];
sx q[1];
rz(-0.75459701) q[1];
sx q[1];
rz(-0.94435) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7992295) q[0];
sx q[0];
rz(-2.9276597) q[0];
sx q[0];
rz(-1.3571204) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3659533) q[2];
sx q[2];
rz(-1.2795942) q[2];
sx q[2];
rz(-1.8679384) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7625092) q[1];
sx q[1];
rz(-3.0804539) q[1];
sx q[1];
rz(-1.5361384) q[1];
rz(-pi) q[2];
rz(-0.51575757) q[3];
sx q[3];
rz(-2.6532647) q[3];
sx q[3];
rz(-2.8749136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8828316) q[2];
sx q[2];
rz(-1.0927039) q[2];
sx q[2];
rz(3.0333701) q[2];
rz(-0.0023068874) q[3];
sx q[3];
rz(-1.5294411) q[3];
sx q[3];
rz(0.32430696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43679431) q[0];
sx q[0];
rz(-2.695485) q[0];
sx q[0];
rz(-0.58445245) q[0];
rz(2.2553518) q[1];
sx q[1];
rz(-2.5247572) q[1];
sx q[1];
rz(-0.054919682) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6829837) q[0];
sx q[0];
rz(-2.8603472) q[0];
sx q[0];
rz(-1.8722697) q[0];
rz(-pi) q[1];
rz(-3.1228742) q[2];
sx q[2];
rz(-1.2476377) q[2];
sx q[2];
rz(-1.2013555) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1850486) q[1];
sx q[1];
rz(-1.946432) q[1];
sx q[1];
rz(-0.59021414) q[1];
rz(-pi) q[2];
rz(-0.92298569) q[3];
sx q[3];
rz(-0.69159782) q[3];
sx q[3];
rz(-1.5036316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3871258) q[2];
sx q[2];
rz(-0.12067623) q[2];
sx q[2];
rz(1.0167271) q[2];
rz(-0.54404849) q[3];
sx q[3];
rz(-2.7816911) q[3];
sx q[3];
rz(-1.3325161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6761557) q[0];
sx q[0];
rz(-2.1570719) q[0];
sx q[0];
rz(2.8570535) q[0];
rz(-0.94447213) q[1];
sx q[1];
rz(-1.1962793) q[1];
sx q[1];
rz(0.91032666) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0383354) q[0];
sx q[0];
rz(-3.0568125) q[0];
sx q[0];
rz(-0.8707365) q[0];
x q[1];
rz(-1.9867284) q[2];
sx q[2];
rz(-1.2565194) q[2];
sx q[2];
rz(-1.8679801) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4714204) q[1];
sx q[1];
rz(-0.52392611) q[1];
sx q[1];
rz(-2.7736204) q[1];
x q[2];
rz(-0.99366412) q[3];
sx q[3];
rz(-0.93615195) q[3];
sx q[3];
rz(-1.5821379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3900782) q[2];
sx q[2];
rz(-3.0780767) q[2];
sx q[2];
rz(0.92203036) q[2];
rz(2.5743124) q[3];
sx q[3];
rz(-1.4138979) q[3];
sx q[3];
rz(1.0197619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2475125) q[0];
sx q[0];
rz(-2.4649354) q[0];
sx q[0];
rz(0.12938736) q[0];
rz(-2.5091876) q[1];
sx q[1];
rz(-1.0267195) q[1];
sx q[1];
rz(2.8410889) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26586543) q[0];
sx q[0];
rz(-0.89571307) q[0];
sx q[0];
rz(-2.6595594) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1883165) q[2];
sx q[2];
rz(-2.2308908) q[2];
sx q[2];
rz(-0.88027871) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.955519) q[1];
sx q[1];
rz(-1.3961853) q[1];
sx q[1];
rz(0.34592918) q[1];
rz(-pi) q[2];
rz(-0.31605966) q[3];
sx q[3];
rz(-0.83871597) q[3];
sx q[3];
rz(0.72533208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5552716) q[2];
sx q[2];
rz(-2.1466612) q[2];
sx q[2];
rz(-0.78197455) q[2];
rz(0.55109763) q[3];
sx q[3];
rz(-1.7588153) q[3];
sx q[3];
rz(-0.15792318) q[3];
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
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5683811) q[0];
sx q[0];
rz(-1.1567572) q[0];
sx q[0];
rz(-3.0138299) q[0];
rz(-2.5993775) q[1];
sx q[1];
rz(-0.95710373) q[1];
sx q[1];
rz(2.382747) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6700867) q[0];
sx q[0];
rz(-2.717088) q[0];
sx q[0];
rz(1.5295117) q[0];
rz(-pi) q[1];
x q[1];
rz(0.27300948) q[2];
sx q[2];
rz(-2.5148279) q[2];
sx q[2];
rz(-3.0221456) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2488333) q[1];
sx q[1];
rz(-1.4151238) q[1];
sx q[1];
rz(-2.4397736) q[1];
rz(3.127029) q[3];
sx q[3];
rz(-2.2363538) q[3];
sx q[3];
rz(-1.927782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0163429) q[2];
sx q[2];
rz(-1.3602076) q[2];
sx q[2];
rz(2.6515567) q[2];
rz(1.4222493) q[3];
sx q[3];
rz(-1.2079206) q[3];
sx q[3];
rz(2.0786044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35995099) q[0];
sx q[0];
rz(-2.521823) q[0];
sx q[0];
rz(-0.075335659) q[0];
rz(2.244859) q[1];
sx q[1];
rz(-1.1743841) q[1];
sx q[1];
rz(-2.5316701) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88203726) q[0];
sx q[0];
rz(-1.2872739) q[0];
sx q[0];
rz(-0.77772227) q[0];
x q[1];
rz(-1.2953193) q[2];
sx q[2];
rz(-1.9343978) q[2];
sx q[2];
rz(1.3678577) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.016644195) q[1];
sx q[1];
rz(-0.84608191) q[1];
sx q[1];
rz(1.8925341) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0697332) q[3];
sx q[3];
rz(-0.89322972) q[3];
sx q[3];
rz(-2.5354405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.23218368) q[2];
sx q[2];
rz(-0.82513088) q[2];
sx q[2];
rz(-2.4278736) q[2];
rz(-2.7632726) q[3];
sx q[3];
rz(-2.648073) q[3];
sx q[3];
rz(0.87987125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80355766) q[0];
sx q[0];
rz(-1.9914347) q[0];
sx q[0];
rz(1.5557355) q[0];
rz(0.65080416) q[1];
sx q[1];
rz(-1.6497859) q[1];
sx q[1];
rz(-0.12129687) q[1];
rz(-1.4177633) q[2];
sx q[2];
rz(-0.19822181) q[2];
sx q[2];
rz(-0.76186686) q[2];
rz(-0.25272947) q[3];
sx q[3];
rz(-2.0287632) q[3];
sx q[3];
rz(-0.91026929) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];