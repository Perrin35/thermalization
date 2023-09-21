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
rz(-2.1531818) q[0];
rz(2.610511) q[1];
sx q[1];
rz(-0.34871066) q[1];
sx q[1];
rz(-1.7928064) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.759401) q[0];
sx q[0];
rz(-2.9392951) q[0];
sx q[0];
rz(-0.6134183) q[0];
rz(-pi) q[1];
rz(0.30795745) q[2];
sx q[2];
rz(-0.27543682) q[2];
sx q[2];
rz(-0.76560417) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.29404902) q[1];
sx q[1];
rz(-2.2806892) q[1];
sx q[1];
rz(-0.82378806) q[1];
x q[2];
rz(0.19282135) q[3];
sx q[3];
rz(-1.76696) q[3];
sx q[3];
rz(1.3960658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3794043) q[2];
sx q[2];
rz(-0.72221243) q[2];
sx q[2];
rz(0.34525004) q[2];
rz(-2.9521862) q[3];
sx q[3];
rz(-2.6510986) q[3];
sx q[3];
rz(-2.1729443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9703366) q[0];
sx q[0];
rz(-2.4617564) q[0];
sx q[0];
rz(2.7217857) q[0];
rz(-2.7833815) q[1];
sx q[1];
rz(-2.5102291) q[1];
sx q[1];
rz(2.125724) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43243877) q[0];
sx q[0];
rz(-1.0854939) q[0];
sx q[0];
rz(-1.1164467) q[0];
x q[1];
rz(0.83769669) q[2];
sx q[2];
rz(-0.3028377) q[2];
sx q[2];
rz(-0.31485117) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.8034755) q[1];
sx q[1];
rz(-2.7403767) q[1];
sx q[1];
rz(2.7944399) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2855371) q[3];
sx q[3];
rz(-0.62969172) q[3];
sx q[3];
rz(-2.7468908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7887855) q[2];
sx q[2];
rz(-2.8217227) q[2];
sx q[2];
rz(-1.1091728) q[2];
rz(0.43168133) q[3];
sx q[3];
rz(-2.7717398) q[3];
sx q[3];
rz(-0.060401827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1576841) q[0];
sx q[0];
rz(-1.2976054) q[0];
sx q[0];
rz(-0.88009673) q[0];
rz(-1.2616715) q[1];
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
rz(2.3343587) q[0];
rz(-pi) q[1];
rz(-2.3679738) q[2];
sx q[2];
rz(-0.81589375) q[2];
sx q[2];
rz(1.9739763) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5460593) q[1];
sx q[1];
rz(-1.5469578) q[1];
sx q[1];
rz(-1.2977352) q[1];
rz(-pi) q[2];
rz(0.39748945) q[3];
sx q[3];
rz(-2.6263413) q[3];
sx q[3];
rz(1.6911563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.024753831) q[2];
sx q[2];
rz(-1.1081868) q[2];
sx q[2];
rz(-1.7987569) q[2];
rz(1.8163619) q[3];
sx q[3];
rz(-2.467005) q[3];
sx q[3];
rz(2.0480115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4745859) q[0];
sx q[0];
rz(-0.51303595) q[0];
sx q[0];
rz(2.4530607) q[0];
rz(-0.21903285) q[1];
sx q[1];
rz(-2.7929247) q[1];
sx q[1];
rz(-2.9825488) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5656887) q[0];
sx q[0];
rz(-2.2825147) q[0];
sx q[0];
rz(-2.2982161) q[0];
rz(0.37695388) q[2];
sx q[2];
rz(-1.9184343) q[2];
sx q[2];
rz(2.8958547) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5608983) q[1];
sx q[1];
rz(-1.782522) q[1];
sx q[1];
rz(-3.1161518) q[1];
rz(-pi) q[2];
rz(-2.6166797) q[3];
sx q[3];
rz(-1.7626926) q[3];
sx q[3];
rz(-1.6247941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8461385) q[2];
sx q[2];
rz(-1.6590154) q[2];
sx q[2];
rz(-0.42567483) q[2];
rz(3.1221636) q[3];
sx q[3];
rz(-2.9255376) q[3];
sx q[3];
rz(-0.69452906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46538019) q[0];
sx q[0];
rz(-0.0066444962) q[0];
sx q[0];
rz(3.0131969) q[0];
rz(0.68583268) q[1];
sx q[1];
rz(-0.99629712) q[1];
sx q[1];
rz(-0.54840666) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0869285) q[0];
sx q[0];
rz(-1.6097277) q[0];
sx q[0];
rz(0.995308) q[0];
rz(-pi) q[1];
rz(0.76813278) q[2];
sx q[2];
rz(-1.7051201) q[2];
sx q[2];
rz(1.3825934) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.027443176) q[1];
sx q[1];
rz(-2.6679734) q[1];
sx q[1];
rz(1.4373535) q[1];
rz(-3.0382285) q[3];
sx q[3];
rz(-1.9723537) q[3];
sx q[3];
rz(0.38927024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2942865) q[2];
sx q[2];
rz(-2.8319478) q[2];
sx q[2];
rz(-1.6748641) q[2];
rz(1.0855801) q[3];
sx q[3];
rz(-1.3054566) q[3];
sx q[3];
rz(-0.87695688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10228957) q[0];
sx q[0];
rz(-1.1941432) q[0];
sx q[0];
rz(-2.0423245) q[0];
rz(-2.8052203) q[1];
sx q[1];
rz(-2.4772494) q[1];
sx q[1];
rz(2.1469595) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5925245) q[0];
sx q[0];
rz(-1.3510834) q[0];
sx q[0];
rz(0.99411221) q[0];
rz(0.50350952) q[2];
sx q[2];
rz(-1.5676055) q[2];
sx q[2];
rz(1.1577595) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.99018807) q[1];
sx q[1];
rz(-1.908704) q[1];
sx q[1];
rz(-3.1347472) q[1];
rz(-1.9739705) q[3];
sx q[3];
rz(-2.2304428) q[3];
sx q[3];
rz(-2.6350104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8424592) q[2];
sx q[2];
rz(-2.2392539) q[2];
sx q[2];
rz(0.4449521) q[2];
rz(2.3343202) q[3];
sx q[3];
rz(-0.10891309) q[3];
sx q[3];
rz(0.81594938) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31535661) q[0];
sx q[0];
rz(-2.5671791) q[0];
sx q[0];
rz(-2.2454967) q[0];
rz(0.90944666) q[1];
sx q[1];
rz(-2.8912631) q[1];
sx q[1];
rz(3.0665841) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5298115) q[0];
sx q[0];
rz(-1.6336332) q[0];
sx q[0];
rz(-1.2114552) q[0];
rz(-1.0117815) q[2];
sx q[2];
rz(-0.85306877) q[2];
sx q[2];
rz(-2.7028524) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9662465) q[1];
sx q[1];
rz(-2.4840377) q[1];
sx q[1];
rz(1.46438) q[1];
rz(-pi) q[2];
rz(-2.6610664) q[3];
sx q[3];
rz(-0.98283813) q[3];
sx q[3];
rz(-2.5939536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9748777) q[2];
sx q[2];
rz(-1.1324984) q[2];
sx q[2];
rz(-1.2274851) q[2];
rz(-0.04315367) q[3];
sx q[3];
rz(-1.4814601) q[3];
sx q[3];
rz(-0.22668049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1535004) q[0];
sx q[0];
rz(-2.7547014) q[0];
sx q[0];
rz(-0.41326997) q[0];
rz(-0.55832541) q[1];
sx q[1];
rz(-0.84086001) q[1];
sx q[1];
rz(2.4024898) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.055187125) q[0];
sx q[0];
rz(-2.1336745) q[0];
sx q[0];
rz(1.6291717) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5584281) q[2];
sx q[2];
rz(-0.45607273) q[2];
sx q[2];
rz(-0.2218483) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6369789) q[1];
sx q[1];
rz(-2.4006872) q[1];
sx q[1];
rz(1.3330589) q[1];
rz(1.8352175) q[3];
sx q[3];
rz(-1.5167055) q[3];
sx q[3];
rz(-0.63187481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3070613) q[2];
sx q[2];
rz(-0.31690159) q[2];
sx q[2];
rz(2.0397662) q[2];
rz(-2.7448591) q[3];
sx q[3];
rz(-1.4361897) q[3];
sx q[3];
rz(-0.81469369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6535764) q[0];
sx q[0];
rz(-0.62960136) q[0];
sx q[0];
rz(-0.6189515) q[0];
rz(-1.0194107) q[1];
sx q[1];
rz(-1.5827725) q[1];
sx q[1];
rz(2.6143262) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70452481) q[0];
sx q[0];
rz(-1.2906115) q[0];
sx q[0];
rz(1.2824476) q[0];
x q[1];
rz(0.69006069) q[2];
sx q[2];
rz(-2.0482716) q[2];
sx q[2];
rz(0.32327393) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4317324) q[1];
sx q[1];
rz(-1.1255956) q[1];
sx q[1];
rz(2.6688337) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23993581) q[3];
sx q[3];
rz(-0.92471189) q[3];
sx q[3];
rz(-2.5480888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.28425372) q[2];
sx q[2];
rz(-1.8062091) q[2];
sx q[2];
rz(-0.93150345) q[2];
rz(-2.3305317) q[3];
sx q[3];
rz(-0.61412007) q[3];
sx q[3];
rz(1.567599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5321524) q[0];
sx q[0];
rz(-1.0487707) q[0];
sx q[0];
rz(0.47927454) q[0];
rz(-0.88538623) q[1];
sx q[1];
rz(-0.65941864) q[1];
sx q[1];
rz(-0.17818174) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7893716) q[0];
sx q[0];
rz(-0.86666115) q[0];
sx q[0];
rz(-2.2602918) q[0];
x q[1];
rz(-2.7257596) q[2];
sx q[2];
rz(-0.73199474) q[2];
sx q[2];
rz(-0.40643613) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0900967) q[1];
sx q[1];
rz(-1.3428543) q[1];
sx q[1];
rz(-1.0432613) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7111669) q[3];
sx q[3];
rz(-1.4595375) q[3];
sx q[3];
rz(0.15669565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.46897727) q[2];
sx q[2];
rz(-2.4185833) q[2];
sx q[2];
rz(2.318312) q[2];
rz(3.0205884) q[3];
sx q[3];
rz(-0.76366097) q[3];
sx q[3];
rz(2.1774489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9085893) q[0];
sx q[0];
rz(-2.067726) q[0];
sx q[0];
rz(2.4698972) q[0];
rz(-1.9227149) q[1];
sx q[1];
rz(-1.3296483) q[1];
sx q[1];
rz(1.7760361) q[1];
rz(1.6505966) q[2];
sx q[2];
rz(-0.36016338) q[2];
sx q[2];
rz(2.2463837) q[2];
rz(2.1929019) q[3];
sx q[3];
rz(-1.8507675) q[3];
sx q[3];
rz(-0.40678195) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];