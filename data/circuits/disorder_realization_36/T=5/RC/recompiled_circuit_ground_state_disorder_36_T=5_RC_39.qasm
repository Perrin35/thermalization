OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1898243) q[0];
sx q[0];
rz(-0.70280743) q[0];
sx q[0];
rz(1.2220609) q[0];
rz(-0.05161461) q[1];
sx q[1];
rz(-1.5049223) q[1];
sx q[1];
rz(1.4805036) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9400222) q[0];
sx q[0];
rz(-2.4594677) q[0];
sx q[0];
rz(1.8322102) q[0];
rz(-3.0323276) q[2];
sx q[2];
rz(-1.7393629) q[2];
sx q[2];
rz(1.8336465) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.29917078) q[1];
sx q[1];
rz(-1.3963194) q[1];
sx q[1];
rz(-0.4293231) q[1];
x q[2];
rz(-1.0869671) q[3];
sx q[3];
rz(-1.9619521) q[3];
sx q[3];
rz(0.018875518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9606955) q[2];
sx q[2];
rz(-0.81520671) q[2];
sx q[2];
rz(2.1107471) q[2];
rz(2.893462) q[3];
sx q[3];
rz(-1.3458359) q[3];
sx q[3];
rz(2.8743751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2613075) q[0];
sx q[0];
rz(-1.945865) q[0];
sx q[0];
rz(-1.1241359) q[0];
rz(-1.0719489) q[1];
sx q[1];
rz(-1.5644466) q[1];
sx q[1];
rz(0.42676485) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4845571) q[0];
sx q[0];
rz(-0.74321514) q[0];
sx q[0];
rz(-1.5667746) q[0];
x q[1];
rz(-0.34554225) q[2];
sx q[2];
rz(-1.3239408) q[2];
sx q[2];
rz(0.051356476) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.93226953) q[1];
sx q[1];
rz(-2.435782) q[1];
sx q[1];
rz(2.9747444) q[1];
rz(-pi) q[2];
rz(-2.9396966) q[3];
sx q[3];
rz(-1.1888973) q[3];
sx q[3];
rz(0.61743067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4830604) q[2];
sx q[2];
rz(-2.6250562) q[2];
sx q[2];
rz(0.10022441) q[2];
rz(-1.077486) q[3];
sx q[3];
rz(-1.5764377) q[3];
sx q[3];
rz(1.235435) q[3];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60233068) q[0];
sx q[0];
rz(-0.97719181) q[0];
sx q[0];
rz(-1.4343028) q[0];
rz(-2.3151248) q[1];
sx q[1];
rz(-1.6441328) q[1];
sx q[1];
rz(-0.43063393) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9481107) q[0];
sx q[0];
rz(-1.8405582) q[0];
sx q[0];
rz(1.7408235) q[0];
rz(1.8437064) q[2];
sx q[2];
rz(-1.4452955) q[2];
sx q[2];
rz(0.93965209) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.34738008) q[1];
sx q[1];
rz(-2.2092704) q[1];
sx q[1];
rz(0.49974167) q[1];
x q[2];
rz(-2.1770559) q[3];
sx q[3];
rz(-2.5068847) q[3];
sx q[3];
rz(-0.026594435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.87381252) q[2];
sx q[2];
rz(-1.9766821) q[2];
sx q[2];
rz(-2.5679892) q[2];
rz(2.7576647) q[3];
sx q[3];
rz(-1.3150747) q[3];
sx q[3];
rz(-2.4887776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4189932) q[0];
sx q[0];
rz(-2.3498131) q[0];
sx q[0];
rz(-2.0844039) q[0];
rz(-2.3253697) q[1];
sx q[1];
rz(-2.1482601) q[1];
sx q[1];
rz(2.9073471) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.809991) q[0];
sx q[0];
rz(-1.1561285) q[0];
sx q[0];
rz(0.11964397) q[0];
x q[1];
rz(-0.42452888) q[2];
sx q[2];
rz(-1.0840975) q[2];
sx q[2];
rz(2.8667712) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2709759) q[1];
sx q[1];
rz(-2.2249195) q[1];
sx q[1];
rz(-2.4821217) q[1];
rz(-pi) q[2];
rz(1.0660421) q[3];
sx q[3];
rz(-1.1467993) q[3];
sx q[3];
rz(0.056886176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8546042) q[2];
sx q[2];
rz(-1.333326) q[2];
sx q[2];
rz(-0.90281478) q[2];
rz(1.3826987) q[3];
sx q[3];
rz(-2.8185676) q[3];
sx q[3];
rz(-0.84669101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.46120241) q[0];
sx q[0];
rz(-1.2795376) q[0];
sx q[0];
rz(-0.70988208) q[0];
rz(-2.6704125) q[1];
sx q[1];
rz(-1.2268365) q[1];
sx q[1];
rz(0.064362854) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1854537) q[0];
sx q[0];
rz(-0.50489932) q[0];
sx q[0];
rz(1.8527777) q[0];
rz(-2.0662599) q[2];
sx q[2];
rz(-1.7771942) q[2];
sx q[2];
rz(-1.0546233) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2379265) q[1];
sx q[1];
rz(-1.9896547) q[1];
sx q[1];
rz(0.55540696) q[1];
rz(-pi) q[2];
rz(3.0547906) q[3];
sx q[3];
rz(-1.5613114) q[3];
sx q[3];
rz(2.529151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.430696) q[2];
sx q[2];
rz(-0.5527834) q[2];
sx q[2];
rz(2.1882679) q[2];
rz(-0.55771762) q[3];
sx q[3];
rz(-1.6744813) q[3];
sx q[3];
rz(-2.6127156) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3144658) q[0];
sx q[0];
rz(-1.0449266) q[0];
sx q[0];
rz(2.1449828) q[0];
rz(3.0820471) q[1];
sx q[1];
rz(-2.153502) q[1];
sx q[1];
rz(-1.7893808) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9080862) q[0];
sx q[0];
rz(-2.3829989) q[0];
sx q[0];
rz(2.3367466) q[0];
rz(-pi) q[1];
rz(2.3699721) q[2];
sx q[2];
rz(-1.811027) q[2];
sx q[2];
rz(-0.26773237) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2595166) q[1];
sx q[1];
rz(-1.8456371) q[1];
sx q[1];
rz(-0.3104574) q[1];
x q[2];
rz(-0.34525235) q[3];
sx q[3];
rz(-2.1998029) q[3];
sx q[3];
rz(-2.0382413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.754564) q[2];
sx q[2];
rz(-1.2754385) q[2];
sx q[2];
rz(-0.55956364) q[2];
rz(2.3061421) q[3];
sx q[3];
rz(-2.7159034) q[3];
sx q[3];
rz(1.3028418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.68592042) q[0];
sx q[0];
rz(-2.4113825) q[0];
sx q[0];
rz(-0.36668229) q[0];
rz(0.8815676) q[1];
sx q[1];
rz(-2.5036948) q[1];
sx q[1];
rz(-1.7104023) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87302215) q[0];
sx q[0];
rz(-2.2013469) q[0];
sx q[0];
rz(-2.0591169) q[0];
rz(-pi) q[1];
rz(-2.4553599) q[2];
sx q[2];
rz(-1.8424818) q[2];
sx q[2];
rz(-1.2154798) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8579432) q[1];
sx q[1];
rz(-1.6873296) q[1];
sx q[1];
rz(1.5999419) q[1];
x q[2];
rz(-2.8411631) q[3];
sx q[3];
rz(-0.68044956) q[3];
sx q[3];
rz(-0.013338683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.35393474) q[2];
sx q[2];
rz(-0.5160318) q[2];
sx q[2];
rz(-0.77989522) q[2];
rz(-0.53906131) q[3];
sx q[3];
rz(-1.5122248) q[3];
sx q[3];
rz(-2.5299634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1270545) q[0];
sx q[0];
rz(-0.67033613) q[0];
sx q[0];
rz(2.3624453) q[0];
rz(2.6226793) q[1];
sx q[1];
rz(-1.2823558) q[1];
sx q[1];
rz(0.40294495) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76501827) q[0];
sx q[0];
rz(-0.85581644) q[0];
sx q[0];
rz(0.9299703) q[0];
x q[1];
rz(2.8267118) q[2];
sx q[2];
rz(-1.1078896) q[2];
sx q[2];
rz(2.7969282) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0190474) q[1];
sx q[1];
rz(-1.1752988) q[1];
sx q[1];
rz(2.4707153) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.99767606) q[3];
sx q[3];
rz(-1.3460288) q[3];
sx q[3];
rz(-2.5483957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.11401033) q[2];
sx q[2];
rz(-1.9738013) q[2];
sx q[2];
rz(1.1172509) q[2];
rz(2.4452325) q[3];
sx q[3];
rz(-1.1098692) q[3];
sx q[3];
rz(-1.8541568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6663412) q[0];
sx q[0];
rz(-0.23463686) q[0];
sx q[0];
rz(-2.7269205) q[0];
rz(1.0617725) q[1];
sx q[1];
rz(-2.2551408) q[1];
sx q[1];
rz(2.3019703) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9794724) q[0];
sx q[0];
rz(-1.9764672) q[0];
sx q[0];
rz(1.0825054) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4213441) q[2];
sx q[2];
rz(-2.0145085) q[2];
sx q[2];
rz(-2.4202079) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.95123467) q[1];
sx q[1];
rz(-2.2273185) q[1];
sx q[1];
rz(-1.3912398) q[1];
rz(0.59323816) q[3];
sx q[3];
rz(-2.7194033) q[3];
sx q[3];
rz(1.7144904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7825369) q[2];
sx q[2];
rz(-1.7724719) q[2];
sx q[2];
rz(-3.0699406) q[2];
rz(-0.045529384) q[3];
sx q[3];
rz(-1.9436911) q[3];
sx q[3];
rz(2.6825452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5422106) q[0];
sx q[0];
rz(-0.57294232) q[0];
sx q[0];
rz(-1.0546767) q[0];
rz(3.1090464) q[1];
sx q[1];
rz(-2.1701505) q[1];
sx q[1];
rz(0.5221101) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73682937) q[0];
sx q[0];
rz(-1.0008345) q[0];
sx q[0];
rz(1.9362157) q[0];
x q[1];
rz(-2.1565151) q[2];
sx q[2];
rz(-0.92364468) q[2];
sx q[2];
rz(0.070731846) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5473816) q[1];
sx q[1];
rz(-1.949233) q[1];
sx q[1];
rz(-2.7345279) q[1];
rz(-pi) q[2];
rz(-2.9154112) q[3];
sx q[3];
rz(-2.1070679) q[3];
sx q[3];
rz(2.8482343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1471499) q[2];
sx q[2];
rz(-2.9394737) q[2];
sx q[2];
rz(-1.7333376) q[2];
rz(-1.5857006) q[3];
sx q[3];
rz(-2.9097911) q[3];
sx q[3];
rz(2.2304992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0178575) q[0];
sx q[0];
rz(-1.3029079) q[0];
sx q[0];
rz(-2.2444176) q[0];
rz(-0.97902117) q[1];
sx q[1];
rz(-1.9985825) q[1];
sx q[1];
rz(-0.40252007) q[1];
rz(1.2382837) q[2];
sx q[2];
rz(-1.3855743) q[2];
sx q[2];
rz(-0.34045548) q[2];
rz(1.7579982) q[3];
sx q[3];
rz(-0.91619195) q[3];
sx q[3];
rz(-0.34294101) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
