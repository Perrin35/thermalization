OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.9525113) q[0];
sx q[0];
rz(-3.014325) q[0];
sx q[0];
rz(-0.98841086) q[0];
rz(-0.53108162) q[1];
sx q[1];
rz(-2.792882) q[1];
sx q[1];
rz(1.7928064) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.792359) q[0];
sx q[0];
rz(-1.6867189) q[0];
sx q[0];
rz(2.9754292) q[0];
rz(-pi) q[1];
rz(2.8336352) q[2];
sx q[2];
rz(-0.27543682) q[2];
sx q[2];
rz(-2.3759885) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4781293) q[1];
sx q[1];
rz(-0.98065286) q[1];
sx q[1];
rz(0.66901916) q[1];
x q[2];
rz(-0.80356055) q[3];
sx q[3];
rz(-2.8674012) q[3];
sx q[3];
rz(-2.1823332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7621883) q[2];
sx q[2];
rz(-0.72221243) q[2];
sx q[2];
rz(0.34525004) q[2];
rz(-0.18940645) q[3];
sx q[3];
rz(-0.49049401) q[3];
sx q[3];
rz(-2.1729443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(-2.9703366) q[0];
sx q[0];
rz(-0.67983627) q[0];
sx q[0];
rz(0.41980699) q[0];
rz(2.7833815) q[1];
sx q[1];
rz(-0.63136357) q[1];
sx q[1];
rz(-1.0158687) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.779218) q[0];
sx q[0];
rz(-1.9694766) q[0];
sx q[0];
rz(2.6108512) q[0];
rz(-1.7989356) q[2];
sx q[2];
rz(-1.7717138) q[2];
sx q[2];
rz(-1.9659496) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7125883) q[1];
sx q[1];
rz(-1.1947558) q[1];
sx q[1];
rz(-1.7141378) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0677098) q[3];
sx q[3];
rz(-1.1745319) q[3];
sx q[3];
rz(0.5644507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.35280716) q[2];
sx q[2];
rz(-2.8217227) q[2];
sx q[2];
rz(-1.1091728) q[2];
rz(-0.43168133) q[3];
sx q[3];
rz(-0.3698529) q[3];
sx q[3];
rz(-0.060401827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9839086) q[0];
sx q[0];
rz(-1.8439872) q[0];
sx q[0];
rz(-0.88009673) q[0];
rz(-1.8799211) q[1];
sx q[1];
rz(-2.547956) q[1];
sx q[1];
rz(-2.9601011) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3767552) q[0];
sx q[0];
rz(-1.1163045) q[0];
sx q[0];
rz(2.6548813) q[0];
rz(-pi) q[1];
rz(-0.9319839) q[2];
sx q[2];
rz(-1.0227232) q[2];
sx q[2];
rz(-1.0149479) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.03141244) q[1];
sx q[1];
rz(-1.8437779) q[1];
sx q[1];
rz(-3.1168373) q[1];
rz(-pi) q[2];
rz(-0.39748945) q[3];
sx q[3];
rz(-0.51525138) q[3];
sx q[3];
rz(-1.4504364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1168388) q[2];
sx q[2];
rz(-2.0334058) q[2];
sx q[2];
rz(1.3428358) q[2];
rz(-1.3252307) q[3];
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
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4745859) q[0];
sx q[0];
rz(-0.51303595) q[0];
sx q[0];
rz(2.4530607) q[0];
rz(-2.9225598) q[1];
sx q[1];
rz(-0.34866798) q[1];
sx q[1];
rz(-2.9825488) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.521579) q[0];
sx q[0];
rz(-1.0431246) q[0];
sx q[0];
rz(-0.85711993) q[0];
x q[1];
rz(-1.1991869) q[2];
sx q[2];
rz(-1.2174165) q[2];
sx q[2];
rz(1.4591188) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4403968) q[1];
sx q[1];
rz(-0.21322589) q[1];
sx q[1];
rz(1.4529983) q[1];
x q[2];
rz(1.791648) q[3];
sx q[3];
rz(-2.0851118) q[3];
sx q[3];
rz(3.0855892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8461385) q[2];
sx q[2];
rz(-1.6590154) q[2];
sx q[2];
rz(-0.42567483) q[2];
rz(-3.1221636) q[3];
sx q[3];
rz(-0.21605505) q[3];
sx q[3];
rz(2.4470636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6762125) q[0];
sx q[0];
rz(-3.1349482) q[0];
sx q[0];
rz(-3.0131969) q[0];
rz(0.68583268) q[1];
sx q[1];
rz(-2.1452955) q[1];
sx q[1];
rz(0.54840666) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6507051) q[0];
sx q[0];
rz(-0.99579949) q[0];
sx q[0];
rz(-0.046396359) q[0];
rz(-pi) q[1];
rz(1.3850645) q[2];
sx q[2];
rz(-0.81132946) q[2];
sx q[2];
rz(-3.082049) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1141495) q[1];
sx q[1];
rz(-0.47361923) q[1];
sx q[1];
rz(1.7042392) q[1];
rz(-pi) q[2];
rz(-3.0382285) q[3];
sx q[3];
rz(-1.169239) q[3];
sx q[3];
rz(2.7523224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2942865) q[2];
sx q[2];
rz(-0.30964482) q[2];
sx q[2];
rz(1.6748641) q[2];
rz(2.0560125) q[3];
sx q[3];
rz(-1.3054566) q[3];
sx q[3];
rz(0.87695688) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10228957) q[0];
sx q[0];
rz(-1.1941432) q[0];
sx q[0];
rz(1.0992682) q[0];
rz(2.8052203) q[1];
sx q[1];
rz(-2.4772494) q[1];
sx q[1];
rz(0.99463314) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69840616) q[0];
sx q[0];
rz(-0.61265677) q[0];
sx q[0];
rz(-1.1820656) q[0];
rz(-pi) q[1];
rz(-1.5671533) q[2];
sx q[2];
rz(-2.074303) q[2];
sx q[2];
rz(2.7267981) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1514046) q[1];
sx q[1];
rz(-1.908704) q[1];
sx q[1];
rz(3.1347472) q[1];
rz(2.6732413) q[3];
sx q[3];
rz(-0.75707179) q[3];
sx q[3];
rz(-0.10145951) q[3];
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
rz(-0.80727243) q[3];
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
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.826236) q[0];
sx q[0];
rz(-2.5671791) q[0];
sx q[0];
rz(0.89609599) q[0];
rz(0.90944666) q[1];
sx q[1];
rz(-0.25032955) q[1];
sx q[1];
rz(0.075008579) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79338193) q[0];
sx q[0];
rz(-2.777034) q[0];
sx q[0];
rz(-1.7478463) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1298112) q[2];
sx q[2];
rz(-0.85306877) q[2];
sx q[2];
rz(0.43874028) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8304886) q[1];
sx q[1];
rz(-1.6357592) q[1];
sx q[1];
rz(0.91598367) q[1];
rz(-0.9644896) q[3];
sx q[3];
rz(-2.4006667) q[3];
sx q[3];
rz(2.9348506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.16671495) q[2];
sx q[2];
rz(-1.1324984) q[2];
sx q[2];
rz(1.9141076) q[2];
rz(-0.04315367) q[3];
sx q[3];
rz(-1.6601325) q[3];
sx q[3];
rz(-2.9149122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
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
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.055187125) q[0];
sx q[0];
rz(-2.1336745) q[0];
sx q[0];
rz(-1.5124209) q[0];
x q[1];
rz(1.5831645) q[2];
sx q[2];
rz(-2.6855199) q[2];
sx q[2];
rz(0.2218483) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.11074556) q[1];
sx q[1];
rz(-1.7304286) q[1];
sx q[1];
rz(2.2975132) q[1];
x q[2];
rz(1.3063752) q[3];
sx q[3];
rz(-1.6248871) q[3];
sx q[3];
rz(-0.63187481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3070613) q[2];
sx q[2];
rz(-2.8246911) q[2];
sx q[2];
rz(-2.0397662) q[2];
rz(-2.7448591) q[3];
sx q[3];
rz(-1.4361897) q[3];
sx q[3];
rz(-0.81469369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48801625) q[0];
sx q[0];
rz(-0.62960136) q[0];
sx q[0];
rz(-0.6189515) q[0];
rz(-2.1221819) q[1];
sx q[1];
rz(-1.5588201) q[1];
sx q[1];
rz(-0.5272665) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0256955) q[0];
sx q[0];
rz(-0.39931116) q[0];
sx q[0];
rz(-2.3621109) q[0];
x q[1];
rz(-2.1617266) q[2];
sx q[2];
rz(-2.1716989) q[2];
sx q[2];
rz(-0.88496937) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.077721715) q[1];
sx q[1];
rz(-1.9943024) q[1];
sx q[1];
rz(1.0788171) q[1];
rz(-pi) q[2];
rz(2.9016568) q[3];
sx q[3];
rz(-0.92471189) q[3];
sx q[3];
rz(-0.59350384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.28425372) q[2];
sx q[2];
rz(-1.3353835) q[2];
sx q[2];
rz(0.93150345) q[2];
rz(-0.81106097) q[3];
sx q[3];
rz(-2.5274726) q[3];
sx q[3];
rz(-1.5739937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60944027) q[0];
sx q[0];
rz(-1.0487707) q[0];
sx q[0];
rz(-0.47927454) q[0];
rz(0.88538623) q[1];
sx q[1];
rz(-0.65941864) q[1];
sx q[1];
rz(-2.9634109) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7893716) q[0];
sx q[0];
rz(-2.2749315) q[0];
sx q[0];
rz(2.2602918) q[0];
x q[1];
rz(2.7257596) q[2];
sx q[2];
rz(-0.73199474) q[2];
sx q[2];
rz(0.40643613) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0308932) q[1];
sx q[1];
rz(-2.5712214) q[1];
sx q[1];
rz(1.13899) q[1];
rz(-pi) q[2];
rz(-0.43042572) q[3];
sx q[3];
rz(-1.6820551) q[3];
sx q[3];
rz(2.984897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6726154) q[2];
sx q[2];
rz(-2.4185833) q[2];
sx q[2];
rz(-2.318312) q[2];
rz(-3.0205884) q[3];
sx q[3];
rz(-0.76366097) q[3];
sx q[3];
rz(0.96414375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
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
rz(-1.490996) q[2];
sx q[2];
rz(-0.36016338) q[2];
sx q[2];
rz(2.2463837) q[2];
rz(-2.0291438) q[3];
sx q[3];
rz(-0.67451285) q[3];
sx q[3];
rz(0.79620517) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
