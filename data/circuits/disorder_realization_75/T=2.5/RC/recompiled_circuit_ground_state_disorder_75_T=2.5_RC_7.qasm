OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.92575443) q[0];
sx q[0];
rz(-0.32045066) q[0];
sx q[0];
rz(-0.12838636) q[0];
rz(4.0965962) q[1];
sx q[1];
rz(8.7648865) q[1];
sx q[1];
rz(13.153491) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.059537236) q[0];
sx q[0];
rz(-1.1063516) q[0];
sx q[0];
rz(1.4091253) q[0];
rz(-2.0287747) q[2];
sx q[2];
rz(-2.0484148) q[2];
sx q[2];
rz(-2.7736563) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.93361357) q[1];
sx q[1];
rz(-2.1449403) q[1];
sx q[1];
rz(1.0453392) q[1];
rz(-pi) q[2];
rz(0.0026851459) q[3];
sx q[3];
rz(-1.955282) q[3];
sx q[3];
rz(-0.17385305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.2653653) q[2];
sx q[2];
rz(-2.9076125) q[2];
sx q[2];
rz(2.2205676) q[2];
rz(1.6169351) q[3];
sx q[3];
rz(-1.2911258) q[3];
sx q[3];
rz(-2.9558712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8601473) q[0];
sx q[0];
rz(-2.7466819) q[0];
sx q[0];
rz(-1.349378) q[0];
rz(0.34740627) q[1];
sx q[1];
rz(-1.1888622) q[1];
sx q[1];
rz(-1.4030392) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8257997) q[0];
sx q[0];
rz(-0.99195882) q[0];
sx q[0];
rz(-1.3114503) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3187693) q[2];
sx q[2];
rz(-2.0206656) q[2];
sx q[2];
rz(2.8577386) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7561556) q[1];
sx q[1];
rz(-1.7240579) q[1];
sx q[1];
rz(2.4244244) q[1];
rz(0.88256695) q[3];
sx q[3];
rz(-1.360713) q[3];
sx q[3];
rz(-0.20245598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8189524) q[2];
sx q[2];
rz(-1.3078657) q[2];
sx q[2];
rz(-0.15882203) q[2];
rz(-2.9239376) q[3];
sx q[3];
rz(-1.7526151) q[3];
sx q[3];
rz(-1.3505107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36667103) q[0];
sx q[0];
rz(-1.8407624) q[0];
sx q[0];
rz(-1.8721254) q[0];
rz(-0.41995755) q[1];
sx q[1];
rz(-2.1407514) q[1];
sx q[1];
rz(-1.6023844) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7363971) q[0];
sx q[0];
rz(-1.1292845) q[0];
sx q[0];
rz(0.34607743) q[0];
x q[1];
rz(1.3871269) q[2];
sx q[2];
rz(-1.3546126) q[2];
sx q[2];
rz(-0.64858299) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.91130891) q[1];
sx q[1];
rz(-0.3989987) q[1];
sx q[1];
rz(1.4330652) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5747803) q[3];
sx q[3];
rz(-1.2897964) q[3];
sx q[3];
rz(-2.5917201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.082108214) q[2];
sx q[2];
rz(-2.2993645) q[2];
sx q[2];
rz(-0.19169894) q[2];
rz(-2.5890403) q[3];
sx q[3];
rz(-1.6165761) q[3];
sx q[3];
rz(1.0244757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(1.7315652) q[0];
sx q[0];
rz(-2.4931694) q[0];
sx q[0];
rz(-2.5372274) q[0];
rz(-1.2454237) q[1];
sx q[1];
rz(-2.3912997) q[1];
sx q[1];
rz(2.2818458) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3533136) q[0];
sx q[0];
rz(-1.9082359) q[0];
sx q[0];
rz(-0.91633015) q[0];
x q[1];
rz(3.1273975) q[2];
sx q[2];
rz(-1.4687339) q[2];
sx q[2];
rz(2.101574) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5936232) q[1];
sx q[1];
rz(-1.7523578) q[1];
sx q[1];
rz(-1.6214002) q[1];
x q[2];
rz(-0.28713496) q[3];
sx q[3];
rz(-0.74112149) q[3];
sx q[3];
rz(2.5484619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2802281) q[2];
sx q[2];
rz(-2.1275529) q[2];
sx q[2];
rz(-1.7029943) q[2];
rz(1.8493308) q[3];
sx q[3];
rz(-2.3137213) q[3];
sx q[3];
rz(-1.0546225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(1.0629405) q[0];
sx q[0];
rz(-0.27314726) q[0];
sx q[0];
rz(-0.31785059) q[0];
rz(1.3392797) q[1];
sx q[1];
rz(-0.41792089) q[1];
sx q[1];
rz(0.97871614) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62417108) q[0];
sx q[0];
rz(-1.540543) q[0];
sx q[0];
rz(-0.4010648) q[0];
x q[1];
rz(2.7182104) q[2];
sx q[2];
rz(-0.21373323) q[2];
sx q[2];
rz(-0.20704421) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3826344) q[1];
sx q[1];
rz(-2.3178604) q[1];
sx q[1];
rz(-1.7428223) q[1];
rz(-pi) q[2];
x q[2];
rz(0.78345244) q[3];
sx q[3];
rz(-2.8718567) q[3];
sx q[3];
rz(1.3470284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9432482) q[2];
sx q[2];
rz(-0.94503108) q[2];
sx q[2];
rz(-1.1172969) q[2];
rz(0.18103389) q[3];
sx q[3];
rz(-1.9046611) q[3];
sx q[3];
rz(1.1972637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4365874) q[0];
sx q[0];
rz(-2.8745108) q[0];
sx q[0];
rz(-1.3096814) q[0];
rz(-1.0218703) q[1];
sx q[1];
rz(-2.3049057) q[1];
sx q[1];
rz(-1.4885363) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2383136) q[0];
sx q[0];
rz(-2.2198703) q[0];
sx q[0];
rz(-2.7191976) q[0];
x q[1];
rz(2.4632881) q[2];
sx q[2];
rz(-1.3875752) q[2];
sx q[2];
rz(-1.2228633) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6147731) q[1];
sx q[1];
rz(-2.062791) q[1];
sx q[1];
rz(1.4089877) q[1];
x q[2];
rz(1.7506787) q[3];
sx q[3];
rz(-2.5733786) q[3];
sx q[3];
rz(-3.1364322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.66270193) q[2];
sx q[2];
rz(-1.7544489) q[2];
sx q[2];
rz(2.877511) q[2];
rz(-1.4763907) q[3];
sx q[3];
rz(-2.4614406) q[3];
sx q[3];
rz(2.7217854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5298117) q[0];
sx q[0];
rz(-1.6407069) q[0];
sx q[0];
rz(2.0773326) q[0];
rz(0.09919676) q[1];
sx q[1];
rz(-1.7925037) q[1];
sx q[1];
rz(1.681021) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49584606) q[0];
sx q[0];
rz(-2.6529713) q[0];
sx q[0];
rz(3.0054128) q[0];
rz(-pi) q[1];
rz(1.6978092) q[2];
sx q[2];
rz(-0.30108967) q[2];
sx q[2];
rz(-1.8179229) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.34014186) q[1];
sx q[1];
rz(-2.6281118) q[1];
sx q[1];
rz(-2.5380847) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4574473) q[3];
sx q[3];
rz(-0.6421291) q[3];
sx q[3];
rz(0.76434842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7907685) q[2];
sx q[2];
rz(-2.0227573) q[2];
sx q[2];
rz(-0.59949818) q[2];
rz(-2.2439469) q[3];
sx q[3];
rz(-1.7220327) q[3];
sx q[3];
rz(-3.1035778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81383234) q[0];
sx q[0];
rz(-1.0603511) q[0];
sx q[0];
rz(-1.4014442) q[0];
rz(0.071488149) q[1];
sx q[1];
rz(-1.6782327) q[1];
sx q[1];
rz(1.00114) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9545549) q[0];
sx q[0];
rz(-1.1251377) q[0];
sx q[0];
rz(1.2902687) q[0];
rz(-pi) q[1];
rz(-0.66364786) q[2];
sx q[2];
rz(-1.8816684) q[2];
sx q[2];
rz(3.1093017) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3197131) q[1];
sx q[1];
rz(-2.1997169) q[1];
sx q[1];
rz(-1.7340388) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4372479) q[3];
sx q[3];
rz(-2.337269) q[3];
sx q[3];
rz(-0.75596228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9643758) q[2];
sx q[2];
rz(-1.588725) q[2];
sx q[2];
rz(0.71844086) q[2];
rz(0.046772379) q[3];
sx q[3];
rz(-2.434157) q[3];
sx q[3];
rz(-2.5717946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30963323) q[0];
sx q[0];
rz(-2.9871812) q[0];
sx q[0];
rz(2.4511448) q[0];
rz(2.9613103) q[1];
sx q[1];
rz(-2.0803662) q[1];
sx q[1];
rz(-2.8188425) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6647341) q[0];
sx q[0];
rz(-2.9551934) q[0];
sx q[0];
rz(2.5405681) q[0];
rz(-pi) q[1];
rz(-0.006470764) q[2];
sx q[2];
rz(-1.0186884) q[2];
sx q[2];
rz(-3.0478549) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.018498713) q[1];
sx q[1];
rz(-2.1040175) q[1];
sx q[1];
rz(-2.4936704) q[1];
rz(-pi) q[2];
rz(2.460367) q[3];
sx q[3];
rz(-1.0325357) q[3];
sx q[3];
rz(-2.1665239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.47101578) q[2];
sx q[2];
rz(-1.3479503) q[2];
sx q[2];
rz(2.1237109) q[2];
rz(0.0035303591) q[3];
sx q[3];
rz(-2.1719833) q[3];
sx q[3];
rz(-1.9297011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3280535) q[0];
sx q[0];
rz(-2.4040451) q[0];
sx q[0];
rz(0.92593431) q[0];
rz(3.0042341) q[1];
sx q[1];
rz(-0.67247144) q[1];
sx q[1];
rz(-0.59991178) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1178499) q[0];
sx q[0];
rz(-1.568092) q[0];
sx q[0];
rz(-2.5976624) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1258862) q[2];
sx q[2];
rz(-2.1246464) q[2];
sx q[2];
rz(2.8987543) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9194738) q[1];
sx q[1];
rz(-1.5403693) q[1];
sx q[1];
rz(0.98201507) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0447015) q[3];
sx q[3];
rz(-1.8623433) q[3];
sx q[3];
rz(2.4959559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7291193) q[2];
sx q[2];
rz(-2.1226661) q[2];
sx q[2];
rz(2.4998383) q[2];
rz(-1.1563835) q[3];
sx q[3];
rz(-0.76120794) q[3];
sx q[3];
rz(1.6183287) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.079035096) q[0];
sx q[0];
rz(-2.2618444) q[0];
sx q[0];
rz(-2.3656144) q[0];
rz(1.7017801) q[1];
sx q[1];
rz(-1.6701313) q[1];
sx q[1];
rz(-3.0727542) q[1];
rz(-0.25664888) q[2];
sx q[2];
rz(-1.7636074) q[2];
sx q[2];
rz(-3.1008813) q[2];
rz(-0.49700789) q[3];
sx q[3];
rz(-1.0773226) q[3];
sx q[3];
rz(1.1763431) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
