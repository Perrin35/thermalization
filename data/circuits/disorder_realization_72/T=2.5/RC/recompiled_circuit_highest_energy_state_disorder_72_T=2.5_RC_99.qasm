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
rz(0.82062757) q[0];
sx q[0];
rz(-0.66250116) q[0];
sx q[0];
rz(2.878046) q[0];
rz(1.1072371) q[1];
sx q[1];
rz(-1.8594445) q[1];
sx q[1];
rz(0.67556226) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1206664) q[0];
sx q[0];
rz(-1.5772737) q[0];
sx q[0];
rz(-1.2351888) q[0];
rz(-pi) q[1];
rz(2.1525429) q[2];
sx q[2];
rz(-2.0768696) q[2];
sx q[2];
rz(-2.842234) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2517667) q[1];
sx q[1];
rz(-1.4205564) q[1];
sx q[1];
rz(1.1778234) q[1];
rz(-0.3278424) q[3];
sx q[3];
rz(-1.0910506) q[3];
sx q[3];
rz(0.2827417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2519309) q[2];
sx q[2];
rz(-1.9193005) q[2];
sx q[2];
rz(2.7737854) q[2];
rz(-2.8273072) q[3];
sx q[3];
rz(-0.29044423) q[3];
sx q[3];
rz(-2.9634326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3322068) q[0];
sx q[0];
rz(-0.089501373) q[0];
sx q[0];
rz(-1.3595164) q[0];
rz(-1.8571732) q[1];
sx q[1];
rz(-1.1651243) q[1];
sx q[1];
rz(3.0899835) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2211014) q[0];
sx q[0];
rz(-1.8360385) q[0];
sx q[0];
rz(-0.18764253) q[0];
rz(2.935595) q[2];
sx q[2];
rz(-1.2075181) q[2];
sx q[2];
rz(-1.3000803) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7588531) q[1];
sx q[1];
rz(-1.5713646) q[1];
sx q[1];
rz(-1.196839) q[1];
x q[2];
rz(-0.47447954) q[3];
sx q[3];
rz(-1.9360376) q[3];
sx q[3];
rz(2.7547835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0525558) q[2];
sx q[2];
rz(-2.8812228) q[2];
sx q[2];
rz(2.5944749) q[2];
rz(-1.268092) q[3];
sx q[3];
rz(-1.781446) q[3];
sx q[3];
rz(-2.8030296) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7608305) q[0];
sx q[0];
rz(-1.0365423) q[0];
sx q[0];
rz(1.8007675) q[0];
rz(-2.2876168) q[1];
sx q[1];
rz(-1.2869336) q[1];
sx q[1];
rz(-0.30240789) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9895565) q[0];
sx q[0];
rz(-1.2947417) q[0];
sx q[0];
rz(-2.2130475) q[0];
rz(-pi) q[1];
rz(2.688681) q[2];
sx q[2];
rz(-0.49395091) q[2];
sx q[2];
rz(-0.043722186) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9871618) q[1];
sx q[1];
rz(-0.87244528) q[1];
sx q[1];
rz(0.72583801) q[1];
rz(-pi) q[2];
rz(-0.27389483) q[3];
sx q[3];
rz(-1.0948101) q[3];
sx q[3];
rz(2.7416503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2053908) q[2];
sx q[2];
rz(-1.8855636) q[2];
sx q[2];
rz(-0.54226792) q[2];
rz(2.1451456) q[3];
sx q[3];
rz(-2.9139329) q[3];
sx q[3];
rz(3.0136133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9484613) q[0];
sx q[0];
rz(-1.6574991) q[0];
sx q[0];
rz(1.2633854) q[0];
rz(-0.43426934) q[1];
sx q[1];
rz(-2.3691005) q[1];
sx q[1];
rz(-1.6901406) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6497191) q[0];
sx q[0];
rz(-1.8611835) q[0];
sx q[0];
rz(2.6650653) q[0];
rz(0.54884882) q[2];
sx q[2];
rz(-1.1804695) q[2];
sx q[2];
rz(3.0496719) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.60148224) q[1];
sx q[1];
rz(-2.9992691) q[1];
sx q[1];
rz(-0.92592923) q[1];
rz(-pi) q[2];
rz(0.94971117) q[3];
sx q[3];
rz(-2.0406282) q[3];
sx q[3];
rz(-2.7787152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0666535) q[2];
sx q[2];
rz(-0.94330072) q[2];
sx q[2];
rz(0.88610506) q[2];
rz(0.0063627176) q[3];
sx q[3];
rz(-1.8013835) q[3];
sx q[3];
rz(-1.1191012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56115443) q[0];
sx q[0];
rz(-0.4564603) q[0];
sx q[0];
rz(1.2503257) q[0];
rz(-2.9087032) q[1];
sx q[1];
rz(-2.3857375) q[1];
sx q[1];
rz(0.09045352) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88921493) q[0];
sx q[0];
rz(-1.6929292) q[0];
sx q[0];
rz(-1.3197164) q[0];
x q[1];
rz(2.6457067) q[2];
sx q[2];
rz(-0.52589881) q[2];
sx q[2];
rz(-2.8484707) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.0098086987) q[1];
sx q[1];
rz(-0.65612853) q[1];
sx q[1];
rz(-3.0387278) q[1];
rz(-pi) q[2];
rz(-3.0556942) q[3];
sx q[3];
rz(-1.2441564) q[3];
sx q[3];
rz(0.24490164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6470486) q[2];
sx q[2];
rz(-0.19379751) q[2];
sx q[2];
rz(-1.2079283) q[2];
rz(-0.9211933) q[3];
sx q[3];
rz(-2.2974206) q[3];
sx q[3];
rz(-0.45879656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10814609) q[0];
sx q[0];
rz(-1.7195846) q[0];
sx q[0];
rz(0.49149996) q[0];
rz(-2.4132701) q[1];
sx q[1];
rz(-2.0305384) q[1];
sx q[1];
rz(0.19457766) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27721632) q[0];
sx q[0];
rz(-0.61232448) q[0];
sx q[0];
rz(1.178647) q[0];
rz(-pi) q[1];
x q[1];
rz(0.79030444) q[2];
sx q[2];
rz(-2.2111597) q[2];
sx q[2];
rz(-1.1932288) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.814333) q[1];
sx q[1];
rz(-2.2047298) q[1];
sx q[1];
rz(-2.4933706) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29816924) q[3];
sx q[3];
rz(-0.67884261) q[3];
sx q[3];
rz(-2.6303675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7877833) q[2];
sx q[2];
rz(-1.0156735) q[2];
sx q[2];
rz(-0.0018399012) q[2];
rz(1.4932102) q[3];
sx q[3];
rz(-0.73072481) q[3];
sx q[3];
rz(-0.28436896) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2581185) q[0];
sx q[0];
rz(-0.26171568) q[0];
sx q[0];
rz(-1.0650241) q[0];
rz(0.26271543) q[1];
sx q[1];
rz(-0.5717259) q[1];
sx q[1];
rz(0.46331847) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.036788851) q[0];
sx q[0];
rz(-2.1533009) q[0];
sx q[0];
rz(0.1881442) q[0];
rz(-pi) q[1];
rz(-2.5538374) q[2];
sx q[2];
rz(-2.3238733) q[2];
sx q[2];
rz(-0.14434926) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9874679) q[1];
sx q[1];
rz(-0.5041669) q[1];
sx q[1];
rz(-2.8567863) q[1];
x q[2];
rz(1.3635395) q[3];
sx q[3];
rz(-3.0240144) q[3];
sx q[3];
rz(-2.1685947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.00039438417) q[2];
sx q[2];
rz(-1.1606777) q[2];
sx q[2];
rz(-0.82994962) q[2];
rz(-2.0937894) q[3];
sx q[3];
rz(-2.6959097) q[3];
sx q[3];
rz(2.9629663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8038427) q[0];
sx q[0];
rz(-2.1512845) q[0];
sx q[0];
rz(2.3867699) q[0];
rz(0.38914514) q[1];
sx q[1];
rz(-1.6837589) q[1];
sx q[1];
rz(2.7438502) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75466418) q[0];
sx q[0];
rz(-2.7542186) q[0];
sx q[0];
rz(-0.8327011) q[0];
rz(-pi) q[1];
rz(-1.3427469) q[2];
sx q[2];
rz(-2.1104276) q[2];
sx q[2];
rz(-1.9613105) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5379679) q[1];
sx q[1];
rz(-1.2984526) q[1];
sx q[1];
rz(2.8773099) q[1];
rz(-pi) q[2];
rz(2.3944309) q[3];
sx q[3];
rz(-2.0748169) q[3];
sx q[3];
rz(0.51591831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8336746) q[2];
sx q[2];
rz(-2.8454744) q[2];
sx q[2];
rz(-2.8274242) q[2];
rz(-1.8999772) q[3];
sx q[3];
rz(-1.5887458) q[3];
sx q[3];
rz(1.9023021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-0.32387963) q[0];
sx q[0];
rz(-0.68809026) q[0];
sx q[0];
rz(2.896198) q[0];
rz(1.9112401) q[1];
sx q[1];
rz(-3.0407258) q[1];
sx q[1];
rz(-2.6255677) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2986261) q[0];
sx q[0];
rz(-0.43382513) q[0];
sx q[0];
rz(2.6515657) q[0];
x q[1];
rz(2.3886959) q[2];
sx q[2];
rz(-2.12924) q[2];
sx q[2];
rz(3.0679997) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.83938767) q[1];
sx q[1];
rz(-0.47014648) q[1];
sx q[1];
rz(-2.4757828) q[1];
x q[2];
rz(-3.0050659) q[3];
sx q[3];
rz(-1.6322005) q[3];
sx q[3];
rz(1.2207954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.55769354) q[2];
sx q[2];
rz(-1.9198753) q[2];
sx q[2];
rz(-1.0830797) q[2];
rz(-0.22100581) q[3];
sx q[3];
rz(-0.14328863) q[3];
sx q[3];
rz(0.7964645) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1204656) q[0];
sx q[0];
rz(-3.0607304) q[0];
sx q[0];
rz(-3.0057111) q[0];
rz(-1.7120301) q[1];
sx q[1];
rz(-1.1212965) q[1];
sx q[1];
rz(2.8543465) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3207239) q[0];
sx q[0];
rz(-0.31021665) q[0];
sx q[0];
rz(-0.91281548) q[0];
rz(-pi) q[1];
rz(-1.2340698) q[2];
sx q[2];
rz(-2.8904466) q[2];
sx q[2];
rz(-0.8539356) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.72037) q[1];
sx q[1];
rz(-2.1541641) q[1];
sx q[1];
rz(1.8349343) q[1];
x q[2];
rz(1.997825) q[3];
sx q[3];
rz(-0.62814071) q[3];
sx q[3];
rz(-2.028156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.91934943) q[2];
sx q[2];
rz(-1.2390077) q[2];
sx q[2];
rz(2.3075721) q[2];
rz(0.30759865) q[3];
sx q[3];
rz(-0.37720507) q[3];
sx q[3];
rz(2.8789177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1766227) q[0];
sx q[0];
rz(-1.3055834) q[0];
sx q[0];
rz(2.1156043) q[0];
rz(2.7550244) q[1];
sx q[1];
rz(-1.3810806) q[1];
sx q[1];
rz(2.5234533) q[1];
rz(2.7170277) q[2];
sx q[2];
rz(-1.4506315) q[2];
sx q[2];
rz(-1.8898024) q[2];
rz(-0.97053075) q[3];
sx q[3];
rz(-0.25401345) q[3];
sx q[3];
rz(0.69380466) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
