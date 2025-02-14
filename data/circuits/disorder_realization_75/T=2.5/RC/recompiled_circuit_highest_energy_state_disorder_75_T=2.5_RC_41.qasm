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
rz(1.6504352) q[0];
sx q[0];
rz(2.7661134) q[0];
sx q[0];
rz(9.8162415) q[0];
rz(0.23671167) q[1];
sx q[1];
rz(-2.1032636) q[1];
sx q[1];
rz(-2.2957323) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4028735) q[0];
sx q[0];
rz(-2.4493626) q[0];
sx q[0];
rz(-2.8701001) q[0];
rz(-2.3416145) q[2];
sx q[2];
rz(-2.5350179) q[2];
sx q[2];
rz(0.457636) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1717901) q[1];
sx q[1];
rz(-2.7493021) q[1];
sx q[1];
rz(-0.7902625) q[1];
x q[2];
rz(-2.1448946) q[3];
sx q[3];
rz(-0.59060589) q[3];
sx q[3];
rz(2.8656538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5819431) q[2];
sx q[2];
rz(-2.1008284) q[2];
sx q[2];
rz(-2.9259658) q[2];
rz(2.1590479) q[3];
sx q[3];
rz(-1.6590786) q[3];
sx q[3];
rz(-1.2561579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4172281) q[0];
sx q[0];
rz(-3.0475782) q[0];
sx q[0];
rz(2.7650058) q[0];
rz(1.4812034) q[1];
sx q[1];
rz(-1.9638289) q[1];
sx q[1];
rz(-2.1362163) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1621012) q[0];
sx q[0];
rz(-2.8164356) q[0];
sx q[0];
rz(1.5473015) q[0];
x q[1];
rz(2.9769276) q[2];
sx q[2];
rz(-0.90064183) q[2];
sx q[2];
rz(-2.7014794) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5498916) q[1];
sx q[1];
rz(-1.6049859) q[1];
sx q[1];
rz(2.4728165) q[1];
x q[2];
rz(2.5672417) q[3];
sx q[3];
rz(-0.6308517) q[3];
sx q[3];
rz(1.7973981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5542095) q[2];
sx q[2];
rz(-0.46450928) q[2];
sx q[2];
rz(1.0760388) q[2];
rz(-0.46766034) q[3];
sx q[3];
rz(-2.3105919) q[3];
sx q[3];
rz(-2.9369798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6429546) q[0];
sx q[0];
rz(-2.0641646) q[0];
sx q[0];
rz(2.0306008) q[0];
rz(-0.30771646) q[1];
sx q[1];
rz(-1.6650763) q[1];
sx q[1];
rz(1.841338) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0205524) q[0];
sx q[0];
rz(-2.5835134) q[0];
sx q[0];
rz(2.9738725) q[0];
rz(1.6069534) q[2];
sx q[2];
rz(-0.99315182) q[2];
sx q[2];
rz(-2.4473178) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8474147) q[1];
sx q[1];
rz(-1.9929033) q[1];
sx q[1];
rz(2.2314784) q[1];
x q[2];
rz(-1.9228705) q[3];
sx q[3];
rz(-1.0961367) q[3];
sx q[3];
rz(-2.3941306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.51851455) q[2];
sx q[2];
rz(-1.0479505) q[2];
sx q[2];
rz(2.9497362) q[2];
rz(1.1040556) q[3];
sx q[3];
rz(-2.7147229) q[3];
sx q[3];
rz(-1.2087315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0693102) q[0];
sx q[0];
rz(-0.59232124) q[0];
sx q[0];
rz(0.88822547) q[0];
rz(-2.8690673) q[1];
sx q[1];
rz(-1.8452294) q[1];
sx q[1];
rz(-1.9857508) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4711535) q[0];
sx q[0];
rz(-1.5616172) q[0];
sx q[0];
rz(-1.8646452) q[0];
x q[1];
rz(-1.5290909) q[2];
sx q[2];
rz(-1.0805849) q[2];
sx q[2];
rz(-2.797567) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.19570505) q[1];
sx q[1];
rz(-2.3263756) q[1];
sx q[1];
rz(2.7570711) q[1];
rz(-0.0083752092) q[3];
sx q[3];
rz(-0.44784689) q[3];
sx q[3];
rz(-0.26219246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3027975) q[2];
sx q[2];
rz(-0.28442997) q[2];
sx q[2];
rz(0.75308853) q[2];
rz(0.49790844) q[3];
sx q[3];
rz(-1.4830736) q[3];
sx q[3];
rz(-2.8411617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4860151) q[0];
sx q[0];
rz(-0.97633728) q[0];
sx q[0];
rz(1.2855541) q[0];
rz(-1.9517508) q[1];
sx q[1];
rz(-0.68584502) q[1];
sx q[1];
rz(1.5600342) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.536895) q[0];
sx q[0];
rz(-1.6486537) q[0];
sx q[0];
rz(1.6998864) q[0];
rz(-pi) q[1];
rz(-2.6597775) q[2];
sx q[2];
rz(-0.31793943) q[2];
sx q[2];
rz(0.16127333) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.469857) q[1];
sx q[1];
rz(-1.3570045) q[1];
sx q[1];
rz(2.695347) q[1];
rz(-1.3693344) q[3];
sx q[3];
rz(-1.3115495) q[3];
sx q[3];
rz(-0.5539971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2994284) q[2];
sx q[2];
rz(-2.2101768) q[2];
sx q[2];
rz(-0.13775873) q[2];
rz(-0.44451851) q[3];
sx q[3];
rz(-1.3714182) q[3];
sx q[3];
rz(2.3270512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.21823068) q[0];
sx q[0];
rz(-0.87319279) q[0];
sx q[0];
rz(0.55106226) q[0];
rz(-0.29523826) q[1];
sx q[1];
rz(-0.72160882) q[1];
sx q[1];
rz(-2.3982184) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0449042) q[0];
sx q[0];
rz(-1.4375234) q[0];
sx q[0];
rz(-0.11479423) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3347537) q[2];
sx q[2];
rz(-1.9034263) q[2];
sx q[2];
rz(1.1548551) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2908926) q[1];
sx q[1];
rz(-0.85324436) q[1];
sx q[1];
rz(3.0102171) q[1];
x q[2];
rz(-2.4287534) q[3];
sx q[3];
rz(-2.3037315) q[3];
sx q[3];
rz(0.43213613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.729852) q[2];
sx q[2];
rz(-0.19062947) q[2];
sx q[2];
rz(-2.8594678) q[2];
rz(0.43635803) q[3];
sx q[3];
rz(-1.3151508) q[3];
sx q[3];
rz(0.25562975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9333164) q[0];
sx q[0];
rz(-1.0641119) q[0];
sx q[0];
rz(0.59260416) q[0];
rz(1.1366049) q[1];
sx q[1];
rz(-1.8698144) q[1];
sx q[1];
rz(2.073854) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16780218) q[0];
sx q[0];
rz(-0.67989698) q[0];
sx q[0];
rz(0.25088422) q[0];
rz(-pi) q[1];
rz(2.7726658) q[2];
sx q[2];
rz(-1.4539945) q[2];
sx q[2];
rz(-0.78250767) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0760666) q[1];
sx q[1];
rz(-2.5812831) q[1];
sx q[1];
rz(1.0275155) q[1];
rz(-pi) q[2];
rz(-0.41530825) q[3];
sx q[3];
rz(-1.7303932) q[3];
sx q[3];
rz(-2.2006765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9390255) q[2];
sx q[2];
rz(-1.86684) q[2];
sx q[2];
rz(-2.1620046) q[2];
rz(-3.0136287) q[3];
sx q[3];
rz(-0.94240677) q[3];
sx q[3];
rz(1.6894107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99280438) q[0];
sx q[0];
rz(-2.1908741) q[0];
sx q[0];
rz(3.0606781) q[0];
rz(-0.12123904) q[1];
sx q[1];
rz(-2.464005) q[1];
sx q[1];
rz(-2.0176719) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8351562) q[0];
sx q[0];
rz(-2.298931) q[0];
sx q[0];
rz(-0.50424285) q[0];
rz(-1.6155682) q[2];
sx q[2];
rz(-1.977619) q[2];
sx q[2];
rz(-1.3699832) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.96059166) q[1];
sx q[1];
rz(-2.7762526) q[1];
sx q[1];
rz(-0.2802556) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8550265) q[3];
sx q[3];
rz(-2.5919302) q[3];
sx q[3];
rz(0.46297234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.27440444) q[2];
sx q[2];
rz(-0.97325456) q[2];
sx q[2];
rz(1.6483866) q[2];
rz(0.46227208) q[3];
sx q[3];
rz(-2.3267764) q[3];
sx q[3];
rz(-1.7260684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82861519) q[0];
sx q[0];
rz(-1.3884437) q[0];
sx q[0];
rz(-1.1653362) q[0];
rz(-0.041821592) q[1];
sx q[1];
rz(-2.1526497) q[1];
sx q[1];
rz(1.3335386) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3305369) q[0];
sx q[0];
rz(-1.0687318) q[0];
sx q[0];
rz(-2.3105225) q[0];
rz(0.012537738) q[2];
sx q[2];
rz(-0.31883966) q[2];
sx q[2];
rz(2.7818305) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0153253) q[1];
sx q[1];
rz(-0.42540144) q[1];
sx q[1];
rz(0.22169561) q[1];
x q[2];
rz(-1.6707375) q[3];
sx q[3];
rz(-1.458711) q[3];
sx q[3];
rz(-2.374927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.24451438) q[2];
sx q[2];
rz(-1.3099193) q[2];
sx q[2];
rz(-0.46373996) q[2];
rz(2.5044299) q[3];
sx q[3];
rz(-1.623268) q[3];
sx q[3];
rz(-2.7630828) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31507444) q[0];
sx q[0];
rz(-1.0239064) q[0];
sx q[0];
rz(-1.6465323) q[0];
rz(-0.82978326) q[1];
sx q[1];
rz(-1.8779495) q[1];
sx q[1];
rz(1.8216546) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49176106) q[0];
sx q[0];
rz(-1.4433089) q[0];
sx q[0];
rz(0.37212917) q[0];
rz(-1.7989203) q[2];
sx q[2];
rz(-1.14883) q[2];
sx q[2];
rz(0.97944469) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4530427) q[1];
sx q[1];
rz(-0.31089766) q[1];
sx q[1];
rz(-2.4255358) q[1];
rz(-pi) q[2];
rz(1.0089031) q[3];
sx q[3];
rz(-1.3504538) q[3];
sx q[3];
rz(1.6957078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7497471) q[2];
sx q[2];
rz(-1.4186991) q[2];
sx q[2];
rz(0.5113655) q[2];
rz(2.5857207) q[3];
sx q[3];
rz(-1.116773) q[3];
sx q[3];
rz(-2.6218124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2340672) q[0];
sx q[0];
rz(-3.0341442) q[0];
sx q[0];
rz(-2.7035614) q[0];
rz(3.0829433) q[1];
sx q[1];
rz(-1.6390683) q[1];
sx q[1];
rz(-1.0674089) q[1];
rz(1.628834) q[2];
sx q[2];
rz(-0.87067247) q[2];
sx q[2];
rz(-1.6280328) q[2];
rz(-2.7300446) q[3];
sx q[3];
rz(-0.66933142) q[3];
sx q[3];
rz(0.29585024) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
