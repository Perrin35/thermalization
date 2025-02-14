OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0786809) q[0];
sx q[0];
rz(-0.39830387) q[0];
sx q[0];
rz(-0.31893536) q[0];
rz(-0.4660663) q[1];
sx q[1];
rz(-1.708344) q[1];
sx q[1];
rz(0.27369174) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6382054) q[0];
sx q[0];
rz(-2.2944258) q[0];
sx q[0];
rz(0.47576018) q[0];
rz(-pi) q[1];
rz(-2.3897522) q[2];
sx q[2];
rz(-2.813857) q[2];
sx q[2];
rz(-0.40204266) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4727386) q[1];
sx q[1];
rz(-1.7484807) q[1];
sx q[1];
rz(0.11858003) q[1];
x q[2];
rz(-0.44173794) q[3];
sx q[3];
rz(-2.3328247) q[3];
sx q[3];
rz(1.7925151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8493001) q[2];
sx q[2];
rz(-2.2165074) q[2];
sx q[2];
rz(1.7577457) q[2];
rz(2.2960831) q[3];
sx q[3];
rz(-0.011761646) q[3];
sx q[3];
rz(0.99172926) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94075769) q[0];
sx q[0];
rz(-2.2469914) q[0];
sx q[0];
rz(1.0351329) q[0];
rz(-1.9738522) q[1];
sx q[1];
rz(-2.9393241) q[1];
sx q[1];
rz(0.77441961) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61078888) q[0];
sx q[0];
rz(-2.4770081) q[0];
sx q[0];
rz(-2.2749645) q[0];
rz(-pi) q[1];
rz(2.0032194) q[2];
sx q[2];
rz(-0.70821643) q[2];
sx q[2];
rz(1.1348292) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.32774156) q[1];
sx q[1];
rz(-0.5264414) q[1];
sx q[1];
rz(1.3415643) q[1];
rz(1.6210591) q[3];
sx q[3];
rz(-1.5368965) q[3];
sx q[3];
rz(0.14824319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8059798) q[2];
sx q[2];
rz(-0.077294417) q[2];
sx q[2];
rz(-0.94266164) q[2];
rz(-3.0521657) q[3];
sx q[3];
rz(-0.663203) q[3];
sx q[3];
rz(-2.8360143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5494004) q[0];
sx q[0];
rz(-0.58627272) q[0];
sx q[0];
rz(2.9495682) q[0];
rz(-2.5281455) q[1];
sx q[1];
rz(-0.44812056) q[1];
sx q[1];
rz(0.014785756) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3131994) q[0];
sx q[0];
rz(-1.1737185) q[0];
sx q[0];
rz(-2.3262789) q[0];
rz(-pi) q[1];
rz(-0.68876617) q[2];
sx q[2];
rz(-1.7851794) q[2];
sx q[2];
rz(0.87300473) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5303684) q[1];
sx q[1];
rz(-0.99677982) q[1];
sx q[1];
rz(-2.8635129) q[1];
x q[2];
rz(0.1802894) q[3];
sx q[3];
rz(-3.0766682) q[3];
sx q[3];
rz(-0.29496671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.52745596) q[2];
sx q[2];
rz(-1.5469896) q[2];
sx q[2];
rz(-3.1318437) q[2];
rz(0.57864946) q[3];
sx q[3];
rz(-1.0520244) q[3];
sx q[3];
rz(-0.78154045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11827271) q[0];
sx q[0];
rz(-0.8096205) q[0];
sx q[0];
rz(0.6672346) q[0];
rz(2.1369333) q[1];
sx q[1];
rz(-2.9861082) q[1];
sx q[1];
rz(-1.7612693) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3074493) q[0];
sx q[0];
rz(-0.47795313) q[0];
sx q[0];
rz(-1.5505928) q[0];
rz(-0.99242248) q[2];
sx q[2];
rz(-2.592591) q[2];
sx q[2];
rz(1.3379607) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.39394571) q[1];
sx q[1];
rz(-2.6470199) q[1];
sx q[1];
rz(3.073758) q[1];
x q[2];
rz(2.0949239) q[3];
sx q[3];
rz(-1.7780684) q[3];
sx q[3];
rz(-0.85974271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7669547) q[2];
sx q[2];
rz(-1.1344323) q[2];
sx q[2];
rz(0.053442027) q[2];
rz(0.63353574) q[3];
sx q[3];
rz(-2.7233248) q[3];
sx q[3];
rz(3.1327278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87090129) q[0];
sx q[0];
rz(-0.57796657) q[0];
sx q[0];
rz(-2.0722678) q[0];
rz(-1.8536192) q[1];
sx q[1];
rz(-0.063871495) q[1];
sx q[1];
rz(0.091648253) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.033460334) q[0];
sx q[0];
rz(-3.1164451) q[0];
sx q[0];
rz(1.0914241) q[0];
rz(-pi) q[1];
rz(0.32817082) q[2];
sx q[2];
rz(-2.0053021) q[2];
sx q[2];
rz(-0.19319867) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3507668) q[1];
sx q[1];
rz(-1.9996399) q[1];
sx q[1];
rz(0.40216202) q[1];
rz(-pi) q[2];
rz(0.76988585) q[3];
sx q[3];
rz(-1.5428233) q[3];
sx q[3];
rz(2.6370492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.84534591) q[2];
sx q[2];
rz(-2.3643934) q[2];
sx q[2];
rz(-0.11543154) q[2];
rz(-2.3760065) q[3];
sx q[3];
rz(-1.6439532) q[3];
sx q[3];
rz(2.7287741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5033443) q[0];
sx q[0];
rz(-0.64953506) q[0];
sx q[0];
rz(-0.58393884) q[0];
rz(-3.0931547) q[1];
sx q[1];
rz(-0.2225288) q[1];
sx q[1];
rz(-2.4979874) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9209401) q[0];
sx q[0];
rz(-0.31198129) q[0];
sx q[0];
rz(1.7259773) q[0];
rz(-pi) q[1];
rz(0.90952833) q[2];
sx q[2];
rz(-1.9478746) q[2];
sx q[2];
rz(0.23185767) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.19932718) q[1];
sx q[1];
rz(-2.2761227) q[1];
sx q[1];
rz(-0.17374111) q[1];
rz(-pi) q[2];
rz(2.3619283) q[3];
sx q[3];
rz(-1.7105967) q[3];
sx q[3];
rz(-1.1598905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.863997) q[2];
sx q[2];
rz(-2.3433351) q[2];
sx q[2];
rz(2.4994728) q[2];
rz(-0.13907214) q[3];
sx q[3];
rz(-3.0032447) q[3];
sx q[3];
rz(-1.4596435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36215255) q[0];
sx q[0];
rz(-0.42069778) q[0];
sx q[0];
rz(-0.62533373) q[0];
rz(-0.23896898) q[1];
sx q[1];
rz(-2.9049554) q[1];
sx q[1];
rz(-1.3561358) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55955333) q[0];
sx q[0];
rz(-1.6252717) q[0];
sx q[0];
rz(-3.1097058) q[0];
x q[1];
rz(-0.37650336) q[2];
sx q[2];
rz(-1.2991221) q[2];
sx q[2];
rz(1.6133378) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7069088) q[1];
sx q[1];
rz(-1.7699935) q[1];
sx q[1];
rz(-1.2046709) q[1];
x q[2];
rz(-1.4144355) q[3];
sx q[3];
rz(-1.3605474) q[3];
sx q[3];
rz(0.69867902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.29505342) q[2];
sx q[2];
rz(-1.2650547) q[2];
sx q[2];
rz(0.96578252) q[2];
rz(-0.24671181) q[3];
sx q[3];
rz(-1.8762981) q[3];
sx q[3];
rz(-1.770796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63916373) q[0];
sx q[0];
rz(-2.736709) q[0];
sx q[0];
rz(-3.0054481) q[0];
rz(-2.4038521) q[1];
sx q[1];
rz(-2.9164011) q[1];
sx q[1];
rz(-1.9053316) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.284974) q[0];
sx q[0];
rz(-0.98561937) q[0];
sx q[0];
rz(2.7057458) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8137553) q[2];
sx q[2];
rz(-0.65082538) q[2];
sx q[2];
rz(-1.6215768) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.41676109) q[1];
sx q[1];
rz(-2.6264852) q[1];
sx q[1];
rz(-0.58578844) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2216283) q[3];
sx q[3];
rz(-1.8355838) q[3];
sx q[3];
rz(1.8643606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.69011921) q[2];
sx q[2];
rz(-1.5560919) q[2];
sx q[2];
rz(2.1915009) q[2];
rz(2.7443366) q[3];
sx q[3];
rz(-0.49644956) q[3];
sx q[3];
rz(2.6955786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5600679) q[0];
sx q[0];
rz(-0.6811322) q[0];
sx q[0];
rz(-3.0233622) q[0];
rz(-2.0589578) q[1];
sx q[1];
rz(-1.2018459) q[1];
sx q[1];
rz(1.7892249) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9064241) q[0];
sx q[0];
rz(-1.6733067) q[0];
sx q[0];
rz(-2.7467523) q[0];
rz(-pi) q[1];
rz(0.544205) q[2];
sx q[2];
rz(-1.4883071) q[2];
sx q[2];
rz(1.5832886) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.39555672) q[1];
sx q[1];
rz(-2.6415351) q[1];
sx q[1];
rz(-1.8386503) q[1];
rz(2.1264137) q[3];
sx q[3];
rz(-0.54856442) q[3];
sx q[3];
rz(3.0231904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.81789404) q[2];
sx q[2];
rz(-2.4511621) q[2];
sx q[2];
rz(-0.77198088) q[2];
rz(-0.082993232) q[3];
sx q[3];
rz(-1.8738184) q[3];
sx q[3];
rz(-0.52223372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25987396) q[0];
sx q[0];
rz(-0.75749713) q[0];
sx q[0];
rz(2.3990193) q[0];
rz(1.2965797) q[1];
sx q[1];
rz(-0.80755889) q[1];
sx q[1];
rz(-0.47320941) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4873692) q[0];
sx q[0];
rz(-0.28074902) q[0];
sx q[0];
rz(1.8868179) q[0];
x q[1];
rz(-0.80188216) q[2];
sx q[2];
rz(-2.7919765) q[2];
sx q[2];
rz(-2.5494818) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1291581) q[1];
sx q[1];
rz(-1.9170463) q[1];
sx q[1];
rz(2.064567) q[1];
x q[2];
rz(0.27373154) q[3];
sx q[3];
rz(-2.377284) q[3];
sx q[3];
rz(0.66018644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5903198) q[2];
sx q[2];
rz(-1.4164305) q[2];
sx q[2];
rz(1.9080706) q[2];
rz(2.7963855) q[3];
sx q[3];
rz(-0.39185169) q[3];
sx q[3];
rz(0.83693081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.681916) q[0];
sx q[0];
rz(-1.3688594) q[0];
sx q[0];
rz(-1.3874227) q[0];
rz(3.0826898) q[1];
sx q[1];
rz(-1.420493) q[1];
sx q[1];
rz(-2.4891985) q[1];
rz(-0.8609345) q[2];
sx q[2];
rz(-0.55649011) q[2];
sx q[2];
rz(1.1330806) q[2];
rz(0.98109365) q[3];
sx q[3];
rz(-2.4193939) q[3];
sx q[3];
rz(-1.9388225) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
