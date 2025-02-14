OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3545544) q[0];
sx q[0];
rz(-2.3784502) q[0];
sx q[0];
rz(-0.95844498) q[0];
rz(-2.7148442) q[1];
sx q[1];
rz(-2.7172105) q[1];
sx q[1];
rz(0.97270614) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63494191) q[0];
sx q[0];
rz(-1.6641064) q[0];
sx q[0];
rz(2.9397381) q[0];
rz(-1.502566) q[2];
sx q[2];
rz(-1.3989324) q[2];
sx q[2];
rz(1.9910938) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0578459) q[1];
sx q[1];
rz(-1.3120756) q[1];
sx q[1];
rz(-0.73212917) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8593721) q[3];
sx q[3];
rz(-1.422554) q[3];
sx q[3];
rz(2.6765598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0613609) q[2];
sx q[2];
rz(-1.1876567) q[2];
sx q[2];
rz(-0.80002552) q[2];
rz(-0.13036615) q[3];
sx q[3];
rz(-0.76821199) q[3];
sx q[3];
rz(-2.8243294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3404563) q[0];
sx q[0];
rz(-1.2202593) q[0];
sx q[0];
rz(-0.98943797) q[0];
rz(0.89775741) q[1];
sx q[1];
rz(-2.0866626) q[1];
sx q[1];
rz(-2.3602233) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88523524) q[0];
sx q[0];
rz(-1.7634749) q[0];
sx q[0];
rz(-2.0847793) q[0];
rz(-1.6590674) q[2];
sx q[2];
rz(-1.8039743) q[2];
sx q[2];
rz(-2.7172497) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0195535) q[1];
sx q[1];
rz(-1.0940503) q[1];
sx q[1];
rz(2.6470584) q[1];
rz(-2.492401) q[3];
sx q[3];
rz(-0.4684557) q[3];
sx q[3];
rz(2.2476803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.117131) q[2];
sx q[2];
rz(-1.7639561) q[2];
sx q[2];
rz(2.3632939) q[2];
rz(-0.69027573) q[3];
sx q[3];
rz(-2.2199151) q[3];
sx q[3];
rz(-1.5811623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86163259) q[0];
sx q[0];
rz(-2.3065688) q[0];
sx q[0];
rz(2.7296208) q[0];
rz(-2.1906134) q[1];
sx q[1];
rz(-2.488766) q[1];
sx q[1];
rz(-1.2118118) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83203234) q[0];
sx q[0];
rz(-1.6044084) q[0];
sx q[0];
rz(1.5322973) q[0];
x q[1];
rz(-3.0427709) q[2];
sx q[2];
rz(-1.3490236) q[2];
sx q[2];
rz(2.2840471) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.46309323) q[1];
sx q[1];
rz(-0.86113595) q[1];
sx q[1];
rz(1.5343496) q[1];
rz(1.4097628) q[3];
sx q[3];
rz(-1.2098243) q[3];
sx q[3];
rz(-1.7019113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.95436207) q[2];
sx q[2];
rz(-0.53202859) q[2];
sx q[2];
rz(1.6383891) q[2];
rz(3.1014118) q[3];
sx q[3];
rz(-2.2153722) q[3];
sx q[3];
rz(2.9068376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9784341) q[0];
sx q[0];
rz(-1.5143159) q[0];
sx q[0];
rz(0.12983313) q[0];
rz(1.2854598) q[1];
sx q[1];
rz(-1.9655656) q[1];
sx q[1];
rz(2.1066378) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1945788) q[0];
sx q[0];
rz(-1.1348551) q[0];
sx q[0];
rz(2.2313124) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0327243) q[2];
sx q[2];
rz(-2.2316859) q[2];
sx q[2];
rz(0.79273293) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.98945204) q[1];
sx q[1];
rz(-1.3229645) q[1];
sx q[1];
rz(1.6590483) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0134936) q[3];
sx q[3];
rz(-2.3689007) q[3];
sx q[3];
rz(-1.1232291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3639565) q[2];
sx q[2];
rz(-1.7549606) q[2];
sx q[2];
rz(-3.0470972) q[2];
rz(0.42094055) q[3];
sx q[3];
rz(-1.4237483) q[3];
sx q[3];
rz(-1.7054935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41973758) q[0];
sx q[0];
rz(-2.5696745) q[0];
sx q[0];
rz(3.1255334) q[0];
rz(-0.84469604) q[1];
sx q[1];
rz(-0.41929308) q[1];
sx q[1];
rz(0.90763456) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7020589) q[0];
sx q[0];
rz(-2.3204813) q[0];
sx q[0];
rz(1.2653714) q[0];
x q[1];
rz(1.9908484) q[2];
sx q[2];
rz(-1.7909652) q[2];
sx q[2];
rz(3.028462) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6542977) q[1];
sx q[1];
rz(-2.0072486) q[1];
sx q[1];
rz(0.36466332) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5441555) q[3];
sx q[3];
rz(-1.9147885) q[3];
sx q[3];
rz(-2.2709803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5661261) q[2];
sx q[2];
rz(-1.0593654) q[2];
sx q[2];
rz(0.063684138) q[2];
rz(2.5746386) q[3];
sx q[3];
rz(-0.83438116) q[3];
sx q[3];
rz(-0.59757346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8276234) q[0];
sx q[0];
rz(-2.9997928) q[0];
sx q[0];
rz(-1.9837448) q[0];
rz(-1.6572378) q[1];
sx q[1];
rz(-1.6969705) q[1];
sx q[1];
rz(-0.2690014) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6873574) q[0];
sx q[0];
rz(-1.7116065) q[0];
sx q[0];
rz(1.6160377) q[0];
x q[1];
rz(-3.0458955) q[2];
sx q[2];
rz(-1.2996593) q[2];
sx q[2];
rz(2.1094131) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3693302) q[1];
sx q[1];
rz(-0.20552615) q[1];
sx q[1];
rz(-2.0783246) q[1];
rz(-1.7888467) q[3];
sx q[3];
rz(-1.9761128) q[3];
sx q[3];
rz(0.97783711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77702648) q[2];
sx q[2];
rz(-2.624161) q[2];
sx q[2];
rz(-2.2606134) q[2];
rz(-0.12428728) q[3];
sx q[3];
rz(-1.4275987) q[3];
sx q[3];
rz(2.4833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86150259) q[0];
sx q[0];
rz(-2.8741591) q[0];
sx q[0];
rz(-0.68914831) q[0];
rz(-0.46734494) q[1];
sx q[1];
rz(-0.57599774) q[1];
sx q[1];
rz(2.0435832) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3680237) q[0];
sx q[0];
rz(-1.7360285) q[0];
sx q[0];
rz(1.6501897) q[0];
rz(-pi) q[1];
rz(-2.6729692) q[2];
sx q[2];
rz(-0.59983095) q[2];
sx q[2];
rz(0.020102321) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5913691) q[1];
sx q[1];
rz(-2.6003631) q[1];
sx q[1];
rz(-3.079097) q[1];
x q[2];
rz(2.2192248) q[3];
sx q[3];
rz(-0.73966129) q[3];
sx q[3];
rz(1.2303703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.565862) q[2];
sx q[2];
rz(-0.87811676) q[2];
sx q[2];
rz(2.511054) q[2];
rz(-0.60728836) q[3];
sx q[3];
rz(-2.2346965) q[3];
sx q[3];
rz(-1.0926532) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8232515) q[0];
sx q[0];
rz(-0.14597758) q[0];
sx q[0];
rz(-1.3463705) q[0];
rz(-1.3660376) q[1];
sx q[1];
rz(-0.64907688) q[1];
sx q[1];
rz(2.5128561) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7808439) q[0];
sx q[0];
rz(-2.1018902) q[0];
sx q[0];
rz(1.2718334) q[0];
rz(-pi) q[1];
rz(2.628831) q[2];
sx q[2];
rz(-0.96549851) q[2];
sx q[2];
rz(-2.5255741) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6727515) q[1];
sx q[1];
rz(-0.76460014) q[1];
sx q[1];
rz(-2.4048664) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9325065) q[3];
sx q[3];
rz(-2.9291398) q[3];
sx q[3];
rz(0.3750876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8667355) q[2];
sx q[2];
rz(-1.8161512) q[2];
sx q[2];
rz(0.87230116) q[2];
rz(-0.21271475) q[3];
sx q[3];
rz(-1.3958967) q[3];
sx q[3];
rz(-0.58342903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.23671959) q[0];
sx q[0];
rz(-1.6107591) q[0];
sx q[0];
rz(1.2637631) q[0];
rz(2.1243375) q[1];
sx q[1];
rz(-2.2433498) q[1];
sx q[1];
rz(-0.57473007) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1270458) q[0];
sx q[0];
rz(-1.0406063) q[0];
sx q[0];
rz(1.6457902) q[0];
rz(-pi) q[1];
rz(2.9024486) q[2];
sx q[2];
rz(-0.84296339) q[2];
sx q[2];
rz(2.4424038) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.1803429) q[1];
sx q[1];
rz(-1.7579792) q[1];
sx q[1];
rz(-0.3734407) q[1];
rz(-pi) q[2];
rz(2.05288) q[3];
sx q[3];
rz(-1.3441372) q[3];
sx q[3];
rz(0.72666336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7175674) q[2];
sx q[2];
rz(-1.307345) q[2];
sx q[2];
rz(0.032111017) q[2];
rz(-1.1561681) q[3];
sx q[3];
rz(-2.7572032) q[3];
sx q[3];
rz(1.2110075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.1011937) q[0];
sx q[0];
rz(-2.5960584) q[0];
sx q[0];
rz(-1.1241166) q[0];
rz(2.595937) q[1];
sx q[1];
rz(-2.7898495) q[1];
sx q[1];
rz(-0.55145946) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41697793) q[0];
sx q[0];
rz(-1.7584086) q[0];
sx q[0];
rz(-0.12599385) q[0];
x q[1];
rz(-1.3842756) q[2];
sx q[2];
rz(-0.2787481) q[2];
sx q[2];
rz(2.2997625) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.76186) q[1];
sx q[1];
rz(-1.1071536) q[1];
sx q[1];
rz(-0.54370466) q[1];
rz(-pi) q[2];
rz(-1.7343246) q[3];
sx q[3];
rz(-1.8732605) q[3];
sx q[3];
rz(1.0176942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5896899) q[2];
sx q[2];
rz(-0.51108131) q[2];
sx q[2];
rz(1.7299293) q[2];
rz(2.3949413) q[3];
sx q[3];
rz(-1.6801497) q[3];
sx q[3];
rz(0.97486973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.830037) q[0];
sx q[0];
rz(-1.5811601) q[0];
sx q[0];
rz(1.5720221) q[0];
rz(-0.39881067) q[1];
sx q[1];
rz(-1.8067982) q[1];
sx q[1];
rz(-1.6511818) q[1];
rz(-2.0592897) q[2];
sx q[2];
rz(-1.5776967) q[2];
sx q[2];
rz(-2.7803905) q[2];
rz(1.2697826) q[3];
sx q[3];
rz(-2.3475937) q[3];
sx q[3];
rz(0.54318843) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
