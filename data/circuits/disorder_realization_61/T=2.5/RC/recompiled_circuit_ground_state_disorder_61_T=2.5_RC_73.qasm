OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9899848) q[0];
sx q[0];
rz(-2.2014872) q[0];
sx q[0];
rz(-0.54036933) q[0];
rz(-2.6481533) q[1];
sx q[1];
rz(-0.72055888) q[1];
sx q[1];
rz(-2.5277353) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9073528) q[0];
sx q[0];
rz(-2.1932903) q[0];
sx q[0];
rz(-1.3433775) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1356642) q[2];
sx q[2];
rz(-1.2790171) q[2];
sx q[2];
rz(0.31729441) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7481374) q[1];
sx q[1];
rz(-2.2917213) q[1];
sx q[1];
rz(1.9490521) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0560232) q[3];
sx q[3];
rz(-0.24805476) q[3];
sx q[3];
rz(0.55094592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2962239) q[2];
sx q[2];
rz(-1.2426528) q[2];
sx q[2];
rz(0.63429147) q[2];
rz(2.2058709) q[3];
sx q[3];
rz(-2.8959385) q[3];
sx q[3];
rz(-0.74265695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.3212386) q[0];
sx q[0];
rz(-1.5411493) q[0];
sx q[0];
rz(0.4441922) q[0];
rz(-0.76002899) q[1];
sx q[1];
rz(-1.1385695) q[1];
sx q[1];
rz(0.98145032) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4430005) q[0];
sx q[0];
rz(-2.5410497) q[0];
sx q[0];
rz(0.53437676) q[0];
x q[1];
rz(-0.75656105) q[2];
sx q[2];
rz(-2.8943099) q[2];
sx q[2];
rz(-0.73723388) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1552298) q[1];
sx q[1];
rz(-1.6587509) q[1];
sx q[1];
rz(-2.7025239) q[1];
rz(-2.7090453) q[3];
sx q[3];
rz(-1.8125696) q[3];
sx q[3];
rz(1.6530619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0280219) q[2];
sx q[2];
rz(-0.61442033) q[2];
sx q[2];
rz(1.2051955) q[2];
rz(2.6049854) q[3];
sx q[3];
rz(-1.9034932) q[3];
sx q[3];
rz(-1.4046148) q[3];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2236915) q[0];
sx q[0];
rz(-1.8603928) q[0];
sx q[0];
rz(-2.3028288) q[0];
rz(2.5054848) q[1];
sx q[1];
rz(-1.6517703) q[1];
sx q[1];
rz(0.020523358) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3394525) q[0];
sx q[0];
rz(-2.1198556) q[0];
sx q[0];
rz(2.2744176) q[0];
rz(-pi) q[1];
rz(-1.1032365) q[2];
sx q[2];
rz(-0.9976495) q[2];
sx q[2];
rz(-1.3064885) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.518259) q[1];
sx q[1];
rz(-1.9837399) q[1];
sx q[1];
rz(1.5625801) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0666615) q[3];
sx q[3];
rz(-0.32204667) q[3];
sx q[3];
rz(-2.5888458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3402349) q[2];
sx q[2];
rz(-1.6131718) q[2];
sx q[2];
rz(-0.94432962) q[2];
rz(-1.0229735) q[3];
sx q[3];
rz(-1.2228271) q[3];
sx q[3];
rz(-1.1378707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.2706547) q[0];
sx q[0];
rz(-1.174467) q[0];
sx q[0];
rz(-1.9842072) q[0];
rz(1.4986787) q[1];
sx q[1];
rz(-1.4721556) q[1];
sx q[1];
rz(-1.1882943) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.37876) q[0];
sx q[0];
rz(-0.68183696) q[0];
sx q[0];
rz(-2.3781611) q[0];
rz(-3.068014) q[2];
sx q[2];
rz(-0.54983222) q[2];
sx q[2];
rz(-2.4940707) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.77515652) q[1];
sx q[1];
rz(-0.85298733) q[1];
sx q[1];
rz(-2.3597673) q[1];
x q[2];
rz(-0.16309901) q[3];
sx q[3];
rz(-1.2360611) q[3];
sx q[3];
rz(-1.7645477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0231861) q[2];
sx q[2];
rz(-1.964317) q[2];
sx q[2];
rz(2.4510621) q[2];
rz(-2.0265419) q[3];
sx q[3];
rz(-2.3670022) q[3];
sx q[3];
rz(1.5007277) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.046415064) q[0];
sx q[0];
rz(-1.7061808) q[0];
sx q[0];
rz(2.114356) q[0];
rz(1.3014303) q[1];
sx q[1];
rz(-0.65168989) q[1];
sx q[1];
rz(-2.8177736) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0831858) q[0];
sx q[0];
rz(-0.91256419) q[0];
sx q[0];
rz(-0.51256521) q[0];
rz(-pi) q[1];
rz(1.0179881) q[2];
sx q[2];
rz(-1.5828642) q[2];
sx q[2];
rz(2.7928305) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.71621543) q[1];
sx q[1];
rz(-1.0310804) q[1];
sx q[1];
rz(1.9464689) q[1];
rz(-1.845076) q[3];
sx q[3];
rz(-2.6135332) q[3];
sx q[3];
rz(1.5825281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.41574898) q[2];
sx q[2];
rz(-0.21012935) q[2];
sx q[2];
rz(2.6591163) q[2];
rz(-1.1680565) q[3];
sx q[3];
rz(-1.2169633) q[3];
sx q[3];
rz(-0.67679685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93550682) q[0];
sx q[0];
rz(-0.85010234) q[0];
sx q[0];
rz(0.8859984) q[0];
rz(-2.50792) q[1];
sx q[1];
rz(-2.251667) q[1];
sx q[1];
rz(2.450313) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1544428) q[0];
sx q[0];
rz(-1.6073174) q[0];
sx q[0];
rz(0.64572389) q[0];
x q[1];
rz(-0.69181594) q[2];
sx q[2];
rz(-1.6199281) q[2];
sx q[2];
rz(3.0494339) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9554841) q[1];
sx q[1];
rz(-2.3209474) q[1];
sx q[1];
rz(3.0226743) q[1];
x q[2];
rz(-1.8227303) q[3];
sx q[3];
rz(-0.79753424) q[3];
sx q[3];
rz(-3.0424527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1401691) q[2];
sx q[2];
rz(-0.836335) q[2];
sx q[2];
rz(2.8720065) q[2];
rz(-2.1932898) q[3];
sx q[3];
rz(-1.5323261) q[3];
sx q[3];
rz(-1.3986826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7798994) q[0];
sx q[0];
rz(-2.7044856) q[0];
sx q[0];
rz(-2.1345188) q[0];
rz(-2.7359447) q[1];
sx q[1];
rz(-2.5463153) q[1];
sx q[1];
rz(2.5880623) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6104975) q[0];
sx q[0];
rz(-0.26108867) q[0];
sx q[0];
rz(1.7946662) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0231452) q[2];
sx q[2];
rz(-2.0549462) q[2];
sx q[2];
rz(0.99768703) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.86846734) q[1];
sx q[1];
rz(-2.1842064) q[1];
sx q[1];
rz(-2.1347743) q[1];
x q[2];
rz(-1.5958435) q[3];
sx q[3];
rz(-1.950693) q[3];
sx q[3];
rz(-2.8632426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8290528) q[2];
sx q[2];
rz(-1.092814) q[2];
sx q[2];
rz(1.9276169) q[2];
rz(-2.35516) q[3];
sx q[3];
rz(-1.7460456) q[3];
sx q[3];
rz(0.023155183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19602747) q[0];
sx q[0];
rz(-2.8493311) q[0];
sx q[0];
rz(-1.3319525) q[0];
rz(-2.5281483) q[1];
sx q[1];
rz(-2.1211801) q[1];
sx q[1];
rz(-2.6920998) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7004823) q[0];
sx q[0];
rz(-1.4876517) q[0];
sx q[0];
rz(1.0821728) q[0];
rz(-1.9067326) q[2];
sx q[2];
rz(-1.6523696) q[2];
sx q[2];
rz(2.5334266) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.66287884) q[1];
sx q[1];
rz(-2.3509074) q[1];
sx q[1];
rz(-1.3776758) q[1];
rz(-0.21829101) q[3];
sx q[3];
rz(-2.8779753) q[3];
sx q[3];
rz(1.5052312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.31676644) q[2];
sx q[2];
rz(-2.4591441) q[2];
sx q[2];
rz(-0.60834926) q[2];
rz(1.1634722) q[3];
sx q[3];
rz(-1.7721662) q[3];
sx q[3];
rz(2.5782862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1083531) q[0];
sx q[0];
rz(-2.2305363) q[0];
sx q[0];
rz(2.5392927) q[0];
rz(0.96744084) q[1];
sx q[1];
rz(-2.2106705) q[1];
sx q[1];
rz(1.1036576) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8062144) q[0];
sx q[0];
rz(-2.4455482) q[0];
sx q[0];
rz(-2.2048414) q[0];
rz(-pi) q[1];
x q[1];
rz(3.121762) q[2];
sx q[2];
rz(-2.5049372) q[2];
sx q[2];
rz(3.0002468) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5491268) q[1];
sx q[1];
rz(-1.5696206) q[1];
sx q[1];
rz(-1.1196613) q[1];
rz(2.1488701) q[3];
sx q[3];
rz(-1.6274656) q[3];
sx q[3];
rz(2.4538159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2063107) q[2];
sx q[2];
rz(-1.8359567) q[2];
sx q[2];
rz(0.3375816) q[2];
rz(0.28389367) q[3];
sx q[3];
rz(-0.68113911) q[3];
sx q[3];
rz(0.18686992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5481446) q[0];
sx q[0];
rz(-1.633506) q[0];
sx q[0];
rz(3.0737851) q[0];
rz(-1.9920805) q[1];
sx q[1];
rz(-1.5510473) q[1];
sx q[1];
rz(0.4745208) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10013469) q[0];
sx q[0];
rz(-0.23915072) q[0];
sx q[0];
rz(1.6275703) q[0];
rz(-pi) q[1];
rz(-0.24093012) q[2];
sx q[2];
rz(-1.4388771) q[2];
sx q[2];
rz(1.9384188) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.671512) q[1];
sx q[1];
rz(-1.5639515) q[1];
sx q[1];
rz(1.6575302) q[1];
x q[2];
rz(0.81409295) q[3];
sx q[3];
rz(-2.11907) q[3];
sx q[3];
rz(0.89422885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.48156753) q[2];
sx q[2];
rz(-2.2109172) q[2];
sx q[2];
rz(1.0732667) q[2];
rz(0.16452161) q[3];
sx q[3];
rz(-1.8142895) q[3];
sx q[3];
rz(-2.8209414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5545223) q[0];
sx q[0];
rz(-1.5656492) q[0];
sx q[0];
rz(1.5026305) q[0];
rz(1.5837689) q[1];
sx q[1];
rz(-1.0777892) q[1];
sx q[1];
rz(-0.52660175) q[1];
rz(2.6804994) q[2];
sx q[2];
rz(-1.245022) q[2];
sx q[2];
rz(1.0427708) q[2];
rz(2.6216636) q[3];
sx q[3];
rz(-0.16362301) q[3];
sx q[3];
rz(2.5273821) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
