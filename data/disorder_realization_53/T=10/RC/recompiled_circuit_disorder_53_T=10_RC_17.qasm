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
rz(2.610511) q[1];
sx q[1];
rz(-0.34871066) q[1];
sx q[1];
rz(-1.7928064) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38219163) q[0];
sx q[0];
rz(-0.2022976) q[0];
sx q[0];
rz(2.5281744) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.30795745) q[2];
sx q[2];
rz(-0.27543682) q[2];
sx q[2];
rz(0.76560417) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4781293) q[1];
sx q[1];
rz(-0.98065286) q[1];
sx q[1];
rz(0.66901916) q[1];
rz(-pi) q[2];
x q[2];
rz(0.19282135) q[3];
sx q[3];
rz(-1.3746327) q[3];
sx q[3];
rz(1.7455268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3794043) q[2];
sx q[2];
rz(-2.4193802) q[2];
sx q[2];
rz(-0.34525004) q[2];
rz(2.9521862) q[3];
sx q[3];
rz(-0.49049401) q[3];
sx q[3];
rz(0.96864831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9703366) q[0];
sx q[0];
rz(-0.67983627) q[0];
sx q[0];
rz(2.7217857) q[0];
rz(2.7833815) q[1];
sx q[1];
rz(-2.5102291) q[1];
sx q[1];
rz(1.0158687) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37594721) q[0];
sx q[0];
rz(-2.489466) q[0];
sx q[0];
rz(2.447522) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.83769669) q[2];
sx q[2];
rz(-2.838755) q[2];
sx q[2];
rz(2.8267415) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.42900436) q[1];
sx q[1];
rz(-1.9468369) q[1];
sx q[1];
rz(1.4274548) q[1];
rz(-pi) q[2];
rz(2.2855371) q[3];
sx q[3];
rz(-0.62969172) q[3];
sx q[3];
rz(-0.39470181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.35280716) q[2];
sx q[2];
rz(-0.31987) q[2];
sx q[2];
rz(-2.0324198) q[2];
rz(-0.43168133) q[3];
sx q[3];
rz(-0.3698529) q[3];
sx q[3];
rz(3.0811908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1576841) q[0];
sx q[0];
rz(-1.8439872) q[0];
sx q[0];
rz(0.88009673) q[0];
rz(-1.2616715) q[1];
sx q[1];
rz(-0.59363669) q[1];
sx q[1];
rz(-2.9601011) q[1];
rz(-pi/2) q[2];
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
rz(-2.491465) q[2];
sx q[2];
rz(-1.0368772) q[2];
sx q[2];
rz(2.9549753) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1101802) q[1];
sx q[1];
rz(-1.2978147) q[1];
sx q[1];
rz(3.1168373) q[1];
rz(2.7441032) q[3];
sx q[3];
rz(-2.6263413) q[3];
sx q[3];
rz(1.4504364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1168388) q[2];
sx q[2];
rz(-2.0334058) q[2];
sx q[2];
rz(-1.7987569) q[2];
rz(1.3252307) q[3];
sx q[3];
rz(-0.67458761) q[3];
sx q[3];
rz(2.0480115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4745859) q[0];
sx q[0];
rz(-2.6285567) q[0];
sx q[0];
rz(2.4530607) q[0];
rz(2.9225598) q[1];
sx q[1];
rz(-0.34866798) q[1];
sx q[1];
rz(2.9825488) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.521579) q[0];
sx q[0];
rz(-1.0431246) q[0];
sx q[0];
rz(0.85711993) q[0];
x q[1];
rz(1.9424058) q[2];
sx q[2];
rz(-1.2174165) q[2];
sx q[2];
rz(1.4591188) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.58069431) q[1];
sx q[1];
rz(-1.782522) q[1];
sx q[1];
rz(-3.1161518) q[1];
rz(-0.52491297) q[3];
sx q[3];
rz(-1.7626926) q[3];
sx q[3];
rz(-1.5167985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.29545414) q[2];
sx q[2];
rz(-1.6590154) q[2];
sx q[2];
rz(-2.7159178) q[2];
rz(-3.1221636) q[3];
sx q[3];
rz(-0.21605505) q[3];
sx q[3];
rz(2.4470636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6762125) q[0];
sx q[0];
rz(-3.1349482) q[0];
sx q[0];
rz(0.12839578) q[0];
rz(-2.45576) q[1];
sx q[1];
rz(-2.1452955) q[1];
sx q[1];
rz(0.54840666) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4908876) q[0];
sx q[0];
rz(-0.99579949) q[0];
sx q[0];
rz(3.0951963) q[0];
rz(-pi) q[1];
rz(2.3734599) q[2];
sx q[2];
rz(-1.7051201) q[2];
sx q[2];
rz(-1.3825934) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.17715958) q[1];
sx q[1];
rz(-2.0398643) q[1];
sx q[1];
rz(-3.0735077) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3324276) q[3];
sx q[3];
rz(-0.41394627) q[3];
sx q[3];
rz(-0.64869374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8473062) q[2];
sx q[2];
rz(-0.30964482) q[2];
sx q[2];
rz(-1.6748641) q[2];
rz(-1.0855801) q[3];
sx q[3];
rz(-1.836136) q[3];
sx q[3];
rz(-0.87695688) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0393031) q[0];
sx q[0];
rz(-1.1941432) q[0];
sx q[0];
rz(-2.0423245) q[0];
rz(0.33637235) q[1];
sx q[1];
rz(-0.66434324) q[1];
sx q[1];
rz(-2.1469595) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1625527) q[0];
sx q[0];
rz(-1.0096706) q[0];
sx q[0];
rz(2.8812376) q[0];
rz(-pi) q[1];
rz(3.1349796) q[2];
sx q[2];
rz(-2.6380739) q[2];
sx q[2];
rz(2.7343482) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5632538) q[1];
sx q[1];
rz(-1.5772547) q[1];
sx q[1];
rz(-1.2328813) q[1];
x q[2];
rz(0.70049882) q[3];
sx q[3];
rz(-1.2555712) q[3];
sx q[3];
rz(-1.8217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8424592) q[2];
sx q[2];
rz(-2.2392539) q[2];
sx q[2];
rz(-2.6966406) q[2];
rz(0.80727243) q[3];
sx q[3];
rz(-0.10891309) q[3];
sx q[3];
rz(-0.81594938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31535661) q[0];
sx q[0];
rz(-0.5744136) q[0];
sx q[0];
rz(0.89609599) q[0];
rz(2.232146) q[1];
sx q[1];
rz(-0.25032955) q[1];
sx q[1];
rz(-0.075008579) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98260005) q[0];
sx q[0];
rz(-1.2121965) q[0];
sx q[0];
rz(0.067111777) q[0];
x q[1];
rz(-2.3415065) q[2];
sx q[2];
rz(-1.9817838) q[2];
sx q[2];
rz(-1.5223741) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8320964) q[1];
sx q[1];
rz(-2.2239904) q[1];
sx q[1];
rz(3.0597568) q[1];
x q[2];
rz(2.1771031) q[3];
sx q[3];
rz(-0.74092591) q[3];
sx q[3];
rz(0.20674202) q[3];
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
rz(1.2274851) q[2];
rz(-3.098439) q[3];
sx q[3];
rz(-1.4814601) q[3];
sx q[3];
rz(0.22668049) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1535004) q[0];
sx q[0];
rz(-0.38689125) q[0];
sx q[0];
rz(-2.7283227) q[0];
rz(2.5832672) q[1];
sx q[1];
rz(-2.3007326) q[1];
sx q[1];
rz(-2.4024898) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.055187125) q[0];
sx q[0];
rz(-1.0079181) q[0];
sx q[0];
rz(1.5124209) q[0];
x q[1];
rz(1.5584281) q[2];
sx q[2];
rz(-2.6855199) q[2];
sx q[2];
rz(-0.2218483) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.11074556) q[1];
sx q[1];
rz(-1.4111641) q[1];
sx q[1];
rz(0.84407945) q[1];
x q[2];
rz(-1.3665175) q[3];
sx q[3];
rz(-2.8718227) q[3];
sx q[3];
rz(2.0055874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3070613) q[2];
sx q[2];
rz(-2.8246911) q[2];
sx q[2];
rz(1.1018264) q[2];
rz(2.7448591) q[3];
sx q[3];
rz(-1.705403) q[3];
sx q[3];
rz(-0.81469369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48801625) q[0];
sx q[0];
rz(-0.62960136) q[0];
sx q[0];
rz(2.5226412) q[0];
rz(1.0194107) q[1];
sx q[1];
rz(-1.5588201) q[1];
sx q[1];
rz(-0.5272665) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4370678) q[0];
sx q[0];
rz(-1.2906115) q[0];
sx q[0];
rz(-1.8591451) q[0];
rz(0.68248827) q[2];
sx q[2];
rz(-2.3253369) q[2];
sx q[2];
rz(1.3860821) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.83890115) q[1];
sx q[1];
rz(-0.63759241) q[1];
sx q[1];
rz(-2.332815) q[1];
rz(-pi) q[2];
rz(-2.9016568) q[3];
sx q[3];
rz(-0.92471189) q[3];
sx q[3];
rz(-2.5480888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8573389) q[2];
sx q[2];
rz(-1.8062091) q[2];
sx q[2];
rz(-2.2100892) q[2];
rz(0.81106097) q[3];
sx q[3];
rz(-2.5274726) q[3];
sx q[3];
rz(-1.567599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5321524) q[0];
sx q[0];
rz(-1.0487707) q[0];
sx q[0];
rz(0.47927454) q[0];
rz(0.88538623) q[1];
sx q[1];
rz(-2.482174) q[1];
sx q[1];
rz(2.9634109) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70885926) q[0];
sx q[0];
rz(-2.0769925) q[0];
sx q[0];
rz(2.3082255) q[0];
x q[1];
rz(0.41583305) q[2];
sx q[2];
rz(-0.73199474) q[2];
sx q[2];
rz(-0.40643613) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0308932) q[1];
sx q[1];
rz(-0.57037121) q[1];
sx q[1];
rz(-1.13899) q[1];
rz(-1.448477) q[3];
sx q[3];
rz(-1.9983851) q[3];
sx q[3];
rz(-1.7784255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6726154) q[2];
sx q[2];
rz(-2.4185833) q[2];
sx q[2];
rz(2.318312) q[2];
rz(-3.0205884) q[3];
sx q[3];
rz(-0.76366097) q[3];
sx q[3];
rz(-2.1774489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9085893) q[0];
sx q[0];
rz(-1.0738666) q[0];
sx q[0];
rz(-0.67169541) q[0];
rz(-1.2188777) q[1];
sx q[1];
rz(-1.8119443) q[1];
sx q[1];
rz(-1.3655566) q[1];
rz(1.490996) q[2];
sx q[2];
rz(-2.7814293) q[2];
sx q[2];
rz(-0.89520892) q[2];
rz(0.34006313) q[3];
sx q[3];
rz(-0.9763413) q[3];
sx q[3];
rz(-1.7819596) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];