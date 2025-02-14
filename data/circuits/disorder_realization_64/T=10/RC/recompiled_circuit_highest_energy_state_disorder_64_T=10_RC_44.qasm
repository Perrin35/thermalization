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
rz(2.0848701) q[0];
sx q[0];
rz(4.4216006) q[0];
sx q[0];
rz(13.292424) q[0];
rz(-2.5441406) q[1];
sx q[1];
rz(-0.51550454) q[1];
sx q[1];
rz(-1.0322303) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7643578) q[0];
sx q[0];
rz(-2.4774158) q[0];
sx q[0];
rz(2.9773877) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0529673) q[2];
sx q[2];
rz(-1.5368688) q[2];
sx q[2];
rz(-2.9133493) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.5121442) q[1];
sx q[1];
rz(-1.8820861) q[1];
sx q[1];
rz(1.6988847) q[1];
x q[2];
rz(-1.0606403) q[3];
sx q[3];
rz(-2.3054696) q[3];
sx q[3];
rz(1.5165129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.44077474) q[2];
sx q[2];
rz(-2.5371607) q[2];
sx q[2];
rz(-3.0536998) q[2];
rz(1.4676189) q[3];
sx q[3];
rz(-1.2940977) q[3];
sx q[3];
rz(-0.99883336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0198332) q[0];
sx q[0];
rz(-1.3061981) q[0];
sx q[0];
rz(-1.4922967) q[0];
rz(0.78798931) q[1];
sx q[1];
rz(-1.9580656) q[1];
sx q[1];
rz(0.96510395) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0168646) q[0];
sx q[0];
rz(-1.4788606) q[0];
sx q[0];
rz(-2.4107736) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4048058) q[2];
sx q[2];
rz(-1.6766785) q[2];
sx q[2];
rz(-0.69241619) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.31307855) q[1];
sx q[1];
rz(-1.113849) q[1];
sx q[1];
rz(-0.20707737) q[1];
rz(-pi) q[2];
rz(0.4697462) q[3];
sx q[3];
rz(-1.6516017) q[3];
sx q[3];
rz(1.684379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1284156) q[2];
sx q[2];
rz(-0.33856496) q[2];
sx q[2];
rz(-1.3236807) q[2];
rz(-0.88661083) q[3];
sx q[3];
rz(-1.9804136) q[3];
sx q[3];
rz(0.19865856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24432467) q[0];
sx q[0];
rz(-2.1710945) q[0];
sx q[0];
rz(-2.4460728) q[0];
rz(2.3059402) q[1];
sx q[1];
rz(-1.4202159) q[1];
sx q[1];
rz(0.64170352) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4280689) q[0];
sx q[0];
rz(-1.4260573) q[0];
sx q[0];
rz(-0.83065303) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4303248) q[2];
sx q[2];
rz(-1.4431653) q[2];
sx q[2];
rz(2.7421599) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1759626) q[1];
sx q[1];
rz(-1.2568736) q[1];
sx q[1];
rz(2.0759275) q[1];
x q[2];
rz(-0.68797994) q[3];
sx q[3];
rz(-2.0350581) q[3];
sx q[3];
rz(-1.2076499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.16331638) q[2];
sx q[2];
rz(-2.0127608) q[2];
sx q[2];
rz(1.8070492) q[2];
rz(-1.4000019) q[3];
sx q[3];
rz(-1.497523) q[3];
sx q[3];
rz(1.1388206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9982933) q[0];
sx q[0];
rz(-2.2704953) q[0];
sx q[0];
rz(-2.368108) q[0];
rz(-0.086291226) q[1];
sx q[1];
rz(-2.6373865) q[1];
sx q[1];
rz(-0.48826826) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21742709) q[0];
sx q[0];
rz(-2.0323951) q[0];
sx q[0];
rz(1.4946926) q[0];
x q[1];
rz(2.5163842) q[2];
sx q[2];
rz(-2.3397988) q[2];
sx q[2];
rz(-1.2259353) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.5544071) q[1];
sx q[1];
rz(-0.76495612) q[1];
sx q[1];
rz(-1.2423963) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7375324) q[3];
sx q[3];
rz(-1.9480138) q[3];
sx q[3];
rz(-0.67507832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6733072) q[2];
sx q[2];
rz(-1.1136592) q[2];
sx q[2];
rz(-0.16806531) q[2];
rz(-2.7189861) q[3];
sx q[3];
rz(-2.1674619) q[3];
sx q[3];
rz(0.15338038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76822686) q[0];
sx q[0];
rz(-0.90665561) q[0];
sx q[0];
rz(-1.1905097) q[0];
rz(2.6137784) q[1];
sx q[1];
rz(-1.0820791) q[1];
sx q[1];
rz(3.1324918) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6884228) q[0];
sx q[0];
rz(-2.1707782) q[0];
sx q[0];
rz(-2.1196516) q[0];
rz(0.78872719) q[2];
sx q[2];
rz(-3.0556779) q[2];
sx q[2];
rz(-0.14363657) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.76705383) q[1];
sx q[1];
rz(-2.3696345) q[1];
sx q[1];
rz(-2.1943387) q[1];
rz(-pi) q[2];
rz(2.6713773) q[3];
sx q[3];
rz(-2.145916) q[3];
sx q[3];
rz(1.0660389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3147543) q[2];
sx q[2];
rz(-0.49586168) q[2];
sx q[2];
rz(2.4763988) q[2];
rz(-2.1151755) q[3];
sx q[3];
rz(-1.0746936) q[3];
sx q[3];
rz(0.78072602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.6219532) q[0];
sx q[0];
rz(-0.53549796) q[0];
sx q[0];
rz(0.39805472) q[0];
rz(-1.6710501) q[1];
sx q[1];
rz(-2.2649951) q[1];
sx q[1];
rz(0.7801396) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1127284) q[0];
sx q[0];
rz(-2.5327589) q[0];
sx q[0];
rz(-0.91797773) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.69177912) q[2];
sx q[2];
rz(-2.8304812) q[2];
sx q[2];
rz(1.9115314) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2237423) q[1];
sx q[1];
rz(-0.7484352) q[1];
sx q[1];
rz(0.42389943) q[1];
rz(-pi) q[2];
rz(2.8242802) q[3];
sx q[3];
rz(-2.2307591) q[3];
sx q[3];
rz(2.8803006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.857343) q[2];
sx q[2];
rz(-1.1025068) q[2];
sx q[2];
rz(-2.9798689) q[2];
rz(-3.0297847) q[3];
sx q[3];
rz(-0.72607741) q[3];
sx q[3];
rz(2.4017754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1424471) q[0];
sx q[0];
rz(-1.7113926) q[0];
sx q[0];
rz(0.69851843) q[0];
rz(-2.0493719) q[1];
sx q[1];
rz(-2.3869546) q[1];
sx q[1];
rz(0.05365595) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2630477) q[0];
sx q[0];
rz(-1.8302655) q[0];
sx q[0];
rz(2.3058913) q[0];
rz(-0.8648134) q[2];
sx q[2];
rz(-2.5891782) q[2];
sx q[2];
rz(-0.68589003) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0567017) q[1];
sx q[1];
rz(-0.47353256) q[1];
sx q[1];
rz(-3.0242306) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7683623) q[3];
sx q[3];
rz(-2.4854492) q[3];
sx q[3];
rz(-3.0246317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1412389) q[2];
sx q[2];
rz(-0.69811368) q[2];
sx q[2];
rz(2.9290283) q[2];
rz(2.8138748) q[3];
sx q[3];
rz(-2.1543584) q[3];
sx q[3];
rz(1.7879965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89606365) q[0];
sx q[0];
rz(-1.117417) q[0];
sx q[0];
rz(0.32972202) q[0];
rz(-2.3700736) q[1];
sx q[1];
rz(-0.78116575) q[1];
sx q[1];
rz(0.82569295) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5485848) q[0];
sx q[0];
rz(-2.4720044) q[0];
sx q[0];
rz(-0.14051814) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4554899) q[2];
sx q[2];
rz(-0.59944154) q[2];
sx q[2];
rz(1.2931461) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9553884) q[1];
sx q[1];
rz(-2.3729747) q[1];
sx q[1];
rz(2.9849418) q[1];
rz(-pi) q[2];
rz(0.65123043) q[3];
sx q[3];
rz(-2.8222601) q[3];
sx q[3];
rz(-0.2624707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9416435) q[2];
sx q[2];
rz(-1.7295094) q[2];
sx q[2];
rz(-1.4818209) q[2];
rz(-3.0485349) q[3];
sx q[3];
rz(-1.2991644) q[3];
sx q[3];
rz(2.602128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2411497) q[0];
sx q[0];
rz(-0.41541442) q[0];
sx q[0];
rz(-2.7442617) q[0];
rz(2.1516402) q[1];
sx q[1];
rz(-1.3497458) q[1];
sx q[1];
rz(2.1053402) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91401228) q[0];
sx q[0];
rz(-1.5591994) q[0];
sx q[0];
rz(-2.7368403) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7781899) q[2];
sx q[2];
rz(-2.7663603) q[2];
sx q[2];
rz(0.64363942) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3237988) q[1];
sx q[1];
rz(-1.3518055) q[1];
sx q[1];
rz(-1.6621672) q[1];
rz(0.37452368) q[3];
sx q[3];
rz(-1.0820884) q[3];
sx q[3];
rz(-0.54218369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.022126023) q[2];
sx q[2];
rz(-2.4194722) q[2];
sx q[2];
rz(-1.8252581) q[2];
rz(0.25257603) q[3];
sx q[3];
rz(-1.794869) q[3];
sx q[3];
rz(-0.36309567) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13755688) q[0];
sx q[0];
rz(-2.3513849) q[0];
sx q[0];
rz(0.23605119) q[0];
rz(3.0421742) q[1];
sx q[1];
rz(-0.75629083) q[1];
sx q[1];
rz(-1.258446) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44348587) q[0];
sx q[0];
rz(-2.5763566) q[0];
sx q[0];
rz(0.14552571) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.070932) q[2];
sx q[2];
rz(-2.1720083) q[2];
sx q[2];
rz(-0.76155969) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7459877) q[1];
sx q[1];
rz(-1.4561936) q[1];
sx q[1];
rz(0.68273165) q[1];
rz(-pi) q[2];
rz(-2.3317565) q[3];
sx q[3];
rz(-2.5812529) q[3];
sx q[3];
rz(2.5446825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1666169) q[2];
sx q[2];
rz(-0.48305837) q[2];
sx q[2];
rz(1.0413337) q[2];
rz(0.58639041) q[3];
sx q[3];
rz(-0.47724637) q[3];
sx q[3];
rz(0.81048107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36622421) q[0];
sx q[0];
rz(-2.0851705) q[0];
sx q[0];
rz(2.0023517) q[0];
rz(0.99209039) q[1];
sx q[1];
rz(-1.6012259) q[1];
sx q[1];
rz(1.59902) q[1];
rz(-0.31789657) q[2];
sx q[2];
rz(-1.7902725) q[2];
sx q[2];
rz(-0.27223311) q[2];
rz(-2.548389) q[3];
sx q[3];
rz(-1.2013669) q[3];
sx q[3];
rz(-0.47040924) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
