OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0754492) q[0];
sx q[0];
rz(-1.0439405) q[0];
sx q[0];
rz(-3.1312842) q[0];
rz(-0.97996867) q[1];
sx q[1];
rz(-1.6966532) q[1];
sx q[1];
rz(0.57656062) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4701472) q[0];
sx q[0];
rz(-1.2630442) q[0];
sx q[0];
rz(-1.6977915) q[0];
x q[1];
rz(-2.5288183) q[2];
sx q[2];
rz(-0.98721993) q[2];
sx q[2];
rz(2.539361) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.51839169) q[1];
sx q[1];
rz(-1.3956337) q[1];
sx q[1];
rz(0.53304146) q[1];
rz(-pi) q[2];
x q[2];
rz(0.16657942) q[3];
sx q[3];
rz(-2.2969679) q[3];
sx q[3];
rz(2.8960901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.52446857) q[2];
sx q[2];
rz(-1.6813797) q[2];
sx q[2];
rz(-1.8560393) q[2];
rz(1.5517976) q[3];
sx q[3];
rz(-1.0714622) q[3];
sx q[3];
rz(-2.0901399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1262421) q[0];
sx q[0];
rz(-1.9168357) q[0];
sx q[0];
rz(-2.8958564) q[0];
rz(1.0579146) q[1];
sx q[1];
rz(-1.3163047) q[1];
sx q[1];
rz(-0.33448514) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34821389) q[0];
sx q[0];
rz(-2.0248027) q[0];
sx q[0];
rz(0.55310849) q[0];
x q[1];
rz(-2.7203015) q[2];
sx q[2];
rz(-2.3347046) q[2];
sx q[2];
rz(2.5786006) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2776371) q[1];
sx q[1];
rz(-0.57687981) q[1];
sx q[1];
rz(0.83630348) q[1];
rz(0.40615079) q[3];
sx q[3];
rz(-0.63780071) q[3];
sx q[3];
rz(-2.8075308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.9574531) q[2];
sx q[2];
rz(-2.2440971) q[2];
sx q[2];
rz(1.755836) q[2];
rz(-2.1615248) q[3];
sx q[3];
rz(-2.6762784) q[3];
sx q[3];
rz(3.0906299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7053213) q[0];
sx q[0];
rz(-2.1846117) q[0];
sx q[0];
rz(1.3901688) q[0];
rz(0.35274371) q[1];
sx q[1];
rz(-0.91068641) q[1];
sx q[1];
rz(2.5211451) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1245956) q[0];
sx q[0];
rz(-0.11717883) q[0];
sx q[0];
rz(2.4524053) q[0];
x q[1];
rz(-2.1587055) q[2];
sx q[2];
rz(-1.6574142) q[2];
sx q[2];
rz(2.4482705) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5689743) q[1];
sx q[1];
rz(-0.44562045) q[1];
sx q[1];
rz(-2.3072412) q[1];
rz(0.87832344) q[3];
sx q[3];
rz(-1.1583503) q[3];
sx q[3];
rz(-2.6869122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0188401) q[2];
sx q[2];
rz(-2.7285125) q[2];
sx q[2];
rz(1.2472461) q[2];
rz(0.16417575) q[3];
sx q[3];
rz(-1.434606) q[3];
sx q[3];
rz(0.629614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71488798) q[0];
sx q[0];
rz(-2.1737104) q[0];
sx q[0];
rz(1.8545275) q[0];
rz(0.60028589) q[1];
sx q[1];
rz(-1.7900107) q[1];
sx q[1];
rz(-2.2379025) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5716887) q[0];
sx q[0];
rz(-1.6780417) q[0];
sx q[0];
rz(-1.6374169) q[0];
rz(-pi) q[1];
rz(0.37995423) q[2];
sx q[2];
rz(-1.6340874) q[2];
sx q[2];
rz(1.5587057) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.43634448) q[1];
sx q[1];
rz(-1.3911934) q[1];
sx q[1];
rz(-0.257538) q[1];
rz(-pi) q[2];
x q[2];
rz(0.65151188) q[3];
sx q[3];
rz(-2.6389696) q[3];
sx q[3];
rz(-2.9732957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.51957447) q[2];
sx q[2];
rz(-2.3081686) q[2];
sx q[2];
rz(2.3681417) q[2];
rz(1.8526239) q[3];
sx q[3];
rz(-0.52923146) q[3];
sx q[3];
rz(-0.73498631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2499823) q[0];
sx q[0];
rz(-1.0050499) q[0];
sx q[0];
rz(0.49215677) q[0];
rz(-0.53897578) q[1];
sx q[1];
rz(-1.4424126) q[1];
sx q[1];
rz(-2.0416562) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66456074) q[0];
sx q[0];
rz(-0.43433055) q[0];
sx q[0];
rz(0.6566027) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5861584) q[2];
sx q[2];
rz(-2.5298497) q[2];
sx q[2];
rz(-1.3604915) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7485986) q[1];
sx q[1];
rz(-1.6637319) q[1];
sx q[1];
rz(-0.8108906) q[1];
rz(-pi) q[2];
rz(-0.48305579) q[3];
sx q[3];
rz(-1.8641346) q[3];
sx q[3];
rz(1.4774069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.34053549) q[2];
sx q[2];
rz(-0.33581442) q[2];
sx q[2];
rz(0.56279969) q[2];
rz(-0.69989145) q[3];
sx q[3];
rz(-1.8507345) q[3];
sx q[3];
rz(0.25115299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3444779) q[0];
sx q[0];
rz(-0.7239224) q[0];
sx q[0];
rz(-2.7864454) q[0];
rz(-1.2190602) q[1];
sx q[1];
rz(-1.2390169) q[1];
sx q[1];
rz(-0.452279) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1817976) q[0];
sx q[0];
rz(-2.1451412) q[0];
sx q[0];
rz(2.6622165) q[0];
rz(-pi) q[1];
rz(-1.1638648) q[2];
sx q[2];
rz(-1.2135047) q[2];
sx q[2];
rz(-0.99321625) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4411583) q[1];
sx q[1];
rz(-1.5467073) q[1];
sx q[1];
rz(-1.8932996) q[1];
rz(-pi) q[2];
rz(-0.62884738) q[3];
sx q[3];
rz(-1.6862685) q[3];
sx q[3];
rz(-0.94146282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5111115) q[2];
sx q[2];
rz(-2.6177572) q[2];
sx q[2];
rz(0.13709489) q[2];
rz(1.5824205) q[3];
sx q[3];
rz(-1.987792) q[3];
sx q[3];
rz(0.77663511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95098507) q[0];
sx q[0];
rz(-0.76222104) q[0];
sx q[0];
rz(-1.5451587) q[0];
rz(-0.51721382) q[1];
sx q[1];
rz(-1.9904174) q[1];
sx q[1];
rz(0.98701611) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68970976) q[0];
sx q[0];
rz(-1.5496043) q[0];
sx q[0];
rz(-3.0107493) q[0];
rz(-pi) q[1];
rz(-2.0275063) q[2];
sx q[2];
rz(-2.4825077) q[2];
sx q[2];
rz(-0.50625077) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.034711866) q[1];
sx q[1];
rz(-1.8000523) q[1];
sx q[1];
rz(0.27193141) q[1];
rz(-pi) q[2];
rz(-0.084264755) q[3];
sx q[3];
rz(-1.2532506) q[3];
sx q[3];
rz(1.0061044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.78642693) q[2];
sx q[2];
rz(-1.0309018) q[2];
sx q[2];
rz(-3.1332341) q[2];
rz(2.9110294) q[3];
sx q[3];
rz(-1.2059261) q[3];
sx q[3];
rz(-1.4768538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1281328) q[0];
sx q[0];
rz(-3.1210493) q[0];
sx q[0];
rz(1.531456) q[0];
rz(1.5229335) q[1];
sx q[1];
rz(-1.5459272) q[1];
sx q[1];
rz(1.5628372) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7195936) q[0];
sx q[0];
rz(-0.61924705) q[0];
sx q[0];
rz(-0.25281711) q[0];
rz(-0.37603907) q[2];
sx q[2];
rz(-1.4362772) q[2];
sx q[2];
rz(-1.0512811) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7111499) q[1];
sx q[1];
rz(-0.64666574) q[1];
sx q[1];
rz(2.4025687) q[1];
x q[2];
rz(-2.4736797) q[3];
sx q[3];
rz(-1.3156943) q[3];
sx q[3];
rz(2.9817443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.83850399) q[2];
sx q[2];
rz(-1.9946626) q[2];
sx q[2];
rz(-3.1040891) q[2];
rz(-0.23669067) q[3];
sx q[3];
rz(-2.762837) q[3];
sx q[3];
rz(1.2934575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26838747) q[0];
sx q[0];
rz(-0.78998843) q[0];
sx q[0];
rz(1.0733806) q[0];
rz(0.34183303) q[1];
sx q[1];
rz(-0.57917246) q[1];
sx q[1];
rz(-0.62172186) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9699696) q[0];
sx q[0];
rz(-1.8299654) q[0];
sx q[0];
rz(0.39879946) q[0];
rz(-1.694181) q[2];
sx q[2];
rz(-0.91431352) q[2];
sx q[2];
rz(0.8379762) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.67111193) q[1];
sx q[1];
rz(-0.59146229) q[1];
sx q[1];
rz(-2.5745113) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.55523606) q[3];
sx q[3];
rz(-2.976417) q[3];
sx q[3];
rz(0.092008807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6217893) q[2];
sx q[2];
rz(-2.9534464) q[2];
sx q[2];
rz(0.56358799) q[2];
rz(1.8300736) q[3];
sx q[3];
rz(-1.9026285) q[3];
sx q[3];
rz(-0.81609503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1821063) q[0];
sx q[0];
rz(-2.2725548) q[0];
sx q[0];
rz(2.2667789) q[0];
rz(-1.4080217) q[1];
sx q[1];
rz(-2.4751016) q[1];
sx q[1];
rz(2.5061238) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4548774) q[0];
sx q[0];
rz(-1.2736763) q[0];
sx q[0];
rz(-2.3700506) q[0];
rz(-pi) q[1];
rz(0.27443011) q[2];
sx q[2];
rz(-1.9727025) q[2];
sx q[2];
rz(2.2257471) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8567896) q[1];
sx q[1];
rz(-1.9724047) q[1];
sx q[1];
rz(-2.1339586) q[1];
rz(-2.1346774) q[3];
sx q[3];
rz(-2.0479408) q[3];
sx q[3];
rz(-1.3729707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1982939) q[2];
sx q[2];
rz(-1.0903) q[2];
sx q[2];
rz(-2.3401006) q[2];
rz(-1.9258026) q[3];
sx q[3];
rz(-2.2299168) q[3];
sx q[3];
rz(-0.041778684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6112678) q[0];
sx q[0];
rz(-1.4401191) q[0];
sx q[0];
rz(-2.8339207) q[0];
rz(-1.2904185) q[1];
sx q[1];
rz(-2.1967874) q[1];
sx q[1];
rz(-2.8312942) q[1];
rz(-1.4353095) q[2];
sx q[2];
rz(-1.3904962) q[2];
sx q[2];
rz(2.4615859) q[2];
rz(-1.4361868) q[3];
sx q[3];
rz(-2.2231839) q[3];
sx q[3];
rz(-0.25683944) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
