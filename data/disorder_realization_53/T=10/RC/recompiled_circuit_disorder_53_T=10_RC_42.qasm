OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.18908137) q[0];
sx q[0];
rz(-0.12726769) q[0];
sx q[0];
rz(0.98841086) q[0];
rz(-0.53108162) q[1];
sx q[1];
rz(-2.792882) q[1];
sx q[1];
rz(1.7928064) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.759401) q[0];
sx q[0];
rz(-2.9392951) q[0];
sx q[0];
rz(-2.5281744) q[0];
rz(-pi) q[1];
rz(-0.30795745) q[2];
sx q[2];
rz(-0.27543682) q[2];
sx q[2];
rz(-2.3759885) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3218282) q[1];
sx q[1];
rz(-2.1121703) q[1];
sx q[1];
rz(2.2775047) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3710255) q[3];
sx q[3];
rz(-1.7598745) q[3];
sx q[3];
rz(0.13669361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7621883) q[2];
sx q[2];
rz(-2.4193802) q[2];
sx q[2];
rz(2.7963426) q[2];
rz(-0.18940645) q[3];
sx q[3];
rz(-2.6510986) q[3];
sx q[3];
rz(2.1729443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9703366) q[0];
sx q[0];
rz(-0.67983627) q[0];
sx q[0];
rz(-2.7217857) q[0];
rz(2.7833815) q[1];
sx q[1];
rz(-2.5102291) q[1];
sx q[1];
rz(1.0158687) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3623747) q[0];
sx q[0];
rz(-1.172116) q[0];
sx q[0];
rz(2.6108512) q[0];
rz(2.935479) q[2];
sx q[2];
rz(-1.3473251) q[2];
sx q[2];
rz(-0.4414562) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7125883) q[1];
sx q[1];
rz(-1.1947558) q[1];
sx q[1];
rz(-1.4274548) q[1];
rz(-pi) q[2];
rz(-2.0738828) q[3];
sx q[3];
rz(-1.9670608) q[3];
sx q[3];
rz(0.5644507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.35280716) q[2];
sx q[2];
rz(-2.8217227) q[2];
sx q[2];
rz(-2.0324198) q[2];
rz(-0.43168133) q[3];
sx q[3];
rz(-0.3698529) q[3];
sx q[3];
rz(-0.060401827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9839086) q[0];
sx q[0];
rz(-1.8439872) q[0];
sx q[0];
rz(2.2614959) q[0];
rz(1.8799211) q[1];
sx q[1];
rz(-2.547956) q[1];
sx q[1];
rz(-0.18149158) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49849579) q[0];
sx q[0];
rz(-2.4884014) q[0];
sx q[0];
rz(-0.80723395) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.7736189) q[2];
sx q[2];
rz(-0.81589375) q[2];
sx q[2];
rz(-1.9739763) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5955334) q[1];
sx q[1];
rz(-1.5946348) q[1];
sx q[1];
rz(-1.8438575) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6604026) q[3];
sx q[3];
rz(-1.7627197) q[3];
sx q[3];
rz(2.6709675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1168388) q[2];
sx q[2];
rz(-2.0334058) q[2];
sx q[2];
rz(1.3428358) q[2];
rz(-1.8163619) q[3];
sx q[3];
rz(-0.67458761) q[3];
sx q[3];
rz(-1.0935812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6670068) q[0];
sx q[0];
rz(-2.6285567) q[0];
sx q[0];
rz(2.4530607) q[0];
rz(-0.21903285) q[1];
sx q[1];
rz(-2.7929247) q[1];
sx q[1];
rz(0.15904388) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6200136) q[0];
sx q[0];
rz(-1.0431246) q[0];
sx q[0];
rz(-0.85711993) q[0];
rz(-pi) q[1];
rz(2.7646388) q[2];
sx q[2];
rz(-1.2231584) q[2];
sx q[2];
rz(-0.24573791) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.70119584) q[1];
sx q[1];
rz(-0.21322589) q[1];
sx q[1];
rz(1.6885944) q[1];
x q[2];
rz(-0.36985107) q[3];
sx q[3];
rz(-2.5858013) q[3];
sx q[3];
rz(0.37214798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.29545414) q[2];
sx q[2];
rz(-1.4825772) q[2];
sx q[2];
rz(-2.7159178) q[2];
rz(-3.1221636) q[3];
sx q[3];
rz(-0.21605505) q[3];
sx q[3];
rz(-0.69452906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6762125) q[0];
sx q[0];
rz(-3.1349482) q[0];
sx q[0];
rz(3.0131969) q[0];
rz(2.45576) q[1];
sx q[1];
rz(-0.99629712) q[1];
sx q[1];
rz(-2.593186) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5760561) q[0];
sx q[0];
rz(-0.57665529) q[0];
sx q[0];
rz(1.499349) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3850645) q[2];
sx q[2];
rz(-2.3302632) q[2];
sx q[2];
rz(-0.059543691) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9644331) q[1];
sx q[1];
rz(-1.1017283) q[1];
sx q[1];
rz(-0.068084929) q[1];
rz(-pi) q[2];
rz(-1.8091651) q[3];
sx q[3];
rz(-0.41394627) q[3];
sx q[3];
rz(-0.64869374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2942865) q[2];
sx q[2];
rz(-2.8319478) q[2];
sx q[2];
rz(-1.4667286) q[2];
rz(-1.0855801) q[3];
sx q[3];
rz(-1.836136) q[3];
sx q[3];
rz(-0.87695688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10228957) q[0];
sx q[0];
rz(-1.1941432) q[0];
sx q[0];
rz(-2.0423245) q[0];
rz(2.8052203) q[1];
sx q[1];
rz(-2.4772494) q[1];
sx q[1];
rz(-2.1469595) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.97904) q[0];
sx q[0];
rz(-2.131922) q[0];
sx q[0];
rz(-2.8812376) q[0];
rz(-2.6380831) q[2];
sx q[2];
rz(-1.5739872) q[2];
sx q[2];
rz(1.9838331) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5632538) q[1];
sx q[1];
rz(-1.564338) q[1];
sx q[1];
rz(1.9087113) q[1];
rz(-0.46835132) q[3];
sx q[3];
rz(-2.3845209) q[3];
sx q[3];
rz(-3.0401331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8424592) q[2];
sx q[2];
rz(-2.2392539) q[2];
sx q[2];
rz(-2.6966406) q[2];
rz(-2.3343202) q[3];
sx q[3];
rz(-3.0326796) q[3];
sx q[3];
rz(-2.3256433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.826236) q[0];
sx q[0];
rz(-0.5744136) q[0];
sx q[0];
rz(0.89609599) q[0];
rz(-0.90944666) q[1];
sx q[1];
rz(-2.8912631) q[1];
sx q[1];
rz(-3.0665841) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79338193) q[0];
sx q[0];
rz(-0.36455867) q[0];
sx q[0];
rz(1.7478463) q[0];
rz(-2.3415065) q[2];
sx q[2];
rz(-1.1598088) q[2];
sx q[2];
rz(1.5223741) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9662465) q[1];
sx q[1];
rz(-0.65755492) q[1];
sx q[1];
rz(1.6772126) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.92618561) q[3];
sx q[3];
rz(-1.965596) q[3];
sx q[3];
rz(-1.3046164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9748777) q[2];
sx q[2];
rz(-2.0090943) q[2];
sx q[2];
rz(-1.9141076) q[2];
rz(-0.04315367) q[3];
sx q[3];
rz(-1.6601325) q[3];
sx q[3];
rz(0.22668049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.9880923) q[0];
sx q[0];
rz(-0.38689125) q[0];
sx q[0];
rz(0.41326997) q[0];
rz(2.5832672) q[1];
sx q[1];
rz(-0.84086001) q[1];
sx q[1];
rz(2.4024898) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.055187125) q[0];
sx q[0];
rz(-1.0079181) q[0];
sx q[0];
rz(1.6291717) q[0];
rz(-pi) q[1];
rz(-3.1355255) q[2];
sx q[2];
rz(-2.0268315) q[2];
sx q[2];
rz(2.9059682) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6369789) q[1];
sx q[1];
rz(-0.74090545) q[1];
sx q[1];
rz(-1.8085338) q[1];
x q[2];
rz(-1.7750752) q[3];
sx q[3];
rz(-0.26976997) q[3];
sx q[3];
rz(-1.1360053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.83453137) q[2];
sx q[2];
rz(-0.31690159) q[2];
sx q[2];
rz(-1.1018264) q[2];
rz(-0.39673355) q[3];
sx q[3];
rz(-1.4361897) q[3];
sx q[3];
rz(-2.326899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48801625) q[0];
sx q[0];
rz(-0.62960136) q[0];
sx q[0];
rz(-0.6189515) q[0];
rz(2.1221819) q[1];
sx q[1];
rz(-1.5827725) q[1];
sx q[1];
rz(2.6143262) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70452481) q[0];
sx q[0];
rz(-1.2906115) q[0];
sx q[0];
rz(1.2824476) q[0];
rz(-pi) q[1];
rz(0.68248827) q[2];
sx q[2];
rz(-2.3253369) q[2];
sx q[2];
rz(1.3860821) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3026915) q[1];
sx q[1];
rz(-2.5040002) q[1];
sx q[1];
rz(-0.80877766) q[1];
rz(-1.8760975) q[3];
sx q[3];
rz(-0.68317181) q[3];
sx q[3];
rz(0.20753577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8573389) q[2];
sx q[2];
rz(-1.3353835) q[2];
sx q[2];
rz(2.2100892) q[2];
rz(-2.3305317) q[3];
sx q[3];
rz(-2.5274726) q[3];
sx q[3];
rz(1.5739937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5321524) q[0];
sx q[0];
rz(-2.092822) q[0];
sx q[0];
rz(-2.6623181) q[0];
rz(-2.2562064) q[1];
sx q[1];
rz(-0.65941864) q[1];
sx q[1];
rz(0.17818174) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70885926) q[0];
sx q[0];
rz(-1.0646001) q[0];
sx q[0];
rz(2.3082255) q[0];
rz(-pi) q[1];
x q[1];
rz(1.222625) q[2];
sx q[2];
rz(-0.91297075) q[2];
sx q[2];
rz(2.1994176) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.61160117) q[1];
sx q[1];
rz(-2.0833263) q[1];
sx q[1];
rz(2.8793053) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7111669) q[3];
sx q[3];
rz(-1.4595375) q[3];
sx q[3];
rz(-0.15669565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.46897727) q[2];
sx q[2];
rz(-0.72300935) q[2];
sx q[2];
rz(2.318312) q[2];
rz(0.12100425) q[3];
sx q[3];
rz(-2.3779317) q[3];
sx q[3];
rz(-0.96414375) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2330033) q[0];
sx q[0];
rz(-2.067726) q[0];
sx q[0];
rz(2.4698972) q[0];
rz(1.2188777) q[1];
sx q[1];
rz(-1.3296483) q[1];
sx q[1];
rz(1.7760361) q[1];
rz(-1.6505966) q[2];
sx q[2];
rz(-2.7814293) q[2];
sx q[2];
rz(-0.89520892) q[2];
rz(-0.34006313) q[3];
sx q[3];
rz(-2.1652514) q[3];
sx q[3];
rz(1.3596331) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
