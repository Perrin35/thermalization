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
rz(-2.6920707) q[0];
sx q[0];
rz(-1.7188526) q[0];
sx q[0];
rz(2.0916405) q[0];
rz(2.6226251) q[1];
sx q[1];
rz(-1.608404) q[1];
sx q[1];
rz(-0.050921507) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.19456) q[0];
sx q[0];
rz(-1.2480535) q[0];
sx q[0];
rz(-2.6850924) q[0];
rz(-2.9257183) q[2];
sx q[2];
rz(-1.2435438) q[2];
sx q[2];
rz(2.808771) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9999034) q[1];
sx q[1];
rz(-1.8898655) q[1];
sx q[1];
rz(1.2525108) q[1];
x q[2];
rz(2.9031119) q[3];
sx q[3];
rz(-2.2422504) q[3];
sx q[3];
rz(-2.8136498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.27429399) q[2];
sx q[2];
rz(-2.2653502) q[2];
sx q[2];
rz(1.743861) q[2];
rz(2.6869669) q[3];
sx q[3];
rz(-2.2635098) q[3];
sx q[3];
rz(-1.4286058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4947263) q[0];
sx q[0];
rz(-2.2541101) q[0];
sx q[0];
rz(2.5383762) q[0];
rz(-1.2546722) q[1];
sx q[1];
rz(-2.5277977) q[1];
sx q[1];
rz(-0.57993728) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86444002) q[0];
sx q[0];
rz(-1.7371558) q[0];
sx q[0];
rz(-2.600738) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0236778) q[2];
sx q[2];
rz(-1.0476026) q[2];
sx q[2];
rz(-2.4459237) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9442288) q[1];
sx q[1];
rz(-1.6104638) q[1];
sx q[1];
rz(-2.1751389) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.5671366) q[3];
sx q[3];
rz(-1.8923762) q[3];
sx q[3];
rz(1.086832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6569528) q[2];
sx q[2];
rz(-2.3886949) q[2];
sx q[2];
rz(0.37772712) q[2];
rz(1.7083302) q[3];
sx q[3];
rz(-1.5092311) q[3];
sx q[3];
rz(1.1908971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3483873) q[0];
sx q[0];
rz(-0.054940104) q[0];
sx q[0];
rz(1.6435664) q[0];
rz(1.161423) q[1];
sx q[1];
rz(-1.252251) q[1];
sx q[1];
rz(-0.84322554) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7347618) q[0];
sx q[0];
rz(-2.1929682) q[0];
sx q[0];
rz(-0.0733331) q[0];
rz(-pi) q[1];
rz(-1.1528422) q[2];
sx q[2];
rz(-1.3713297) q[2];
sx q[2];
rz(2.0355088) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.49065152) q[1];
sx q[1];
rz(-2.2434714) q[1];
sx q[1];
rz(2.9789799) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2349818) q[3];
sx q[3];
rz(-0.89794176) q[3];
sx q[3];
rz(-1.8115753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2735542) q[2];
sx q[2];
rz(-0.66323438) q[2];
sx q[2];
rz(2.3466477) q[2];
rz(0.76977175) q[3];
sx q[3];
rz(-0.92307463) q[3];
sx q[3];
rz(-0.15708378) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5539826) q[0];
sx q[0];
rz(-1.3251323) q[0];
sx q[0];
rz(0.28833589) q[0];
rz(-0.68710697) q[1];
sx q[1];
rz(-1.4930875) q[1];
sx q[1];
rz(1.7126602) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16609678) q[0];
sx q[0];
rz(-1.5642484) q[0];
sx q[0];
rz(-0.01207821) q[0];
rz(0.32055118) q[2];
sx q[2];
rz(-0.42141576) q[2];
sx q[2];
rz(-1.5326064) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1667119) q[1];
sx q[1];
rz(-0.26435164) q[1];
sx q[1];
rz(0.2522142) q[1];
x q[2];
rz(-2.3219548) q[3];
sx q[3];
rz(-2.313531) q[3];
sx q[3];
rz(-2.6520906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5700506) q[2];
sx q[2];
rz(-2.1752581) q[2];
sx q[2];
rz(-0.16560444) q[2];
rz(-0.064519493) q[3];
sx q[3];
rz(-0.12972984) q[3];
sx q[3];
rz(1.8701514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22873838) q[0];
sx q[0];
rz(-1.1930635) q[0];
sx q[0];
rz(-2.2304529) q[0];
rz(-0.97081026) q[1];
sx q[1];
rz(-1.3263005) q[1];
sx q[1];
rz(2.2784746) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1939094) q[0];
sx q[0];
rz(-1.860431) q[0];
sx q[0];
rz(-2.0783483) q[0];
rz(-pi) q[1];
rz(-1.4865033) q[2];
sx q[2];
rz(-0.26940036) q[2];
sx q[2];
rz(2.0790554) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6533901) q[1];
sx q[1];
rz(-2.0987256) q[1];
sx q[1];
rz(1.3009682) q[1];
rz(-pi) q[2];
rz(3.0961995) q[3];
sx q[3];
rz(-1.8828585) q[3];
sx q[3];
rz(1.454221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.054964885) q[2];
sx q[2];
rz(-1.2530155) q[2];
sx q[2];
rz(2.6524554) q[2];
rz(-2.1431811) q[3];
sx q[3];
rz(-0.17397927) q[3];
sx q[3];
rz(0.099099549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84422207) q[0];
sx q[0];
rz(-0.64592823) q[0];
sx q[0];
rz(-2.5883664) q[0];
rz(-0.79752254) q[1];
sx q[1];
rz(-2.1452466) q[1];
sx q[1];
rz(1.1513938) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0259375) q[0];
sx q[0];
rz(-2.135072) q[0];
sx q[0];
rz(1.7274117) q[0];
rz(-pi) q[1];
rz(0.55091605) q[2];
sx q[2];
rz(-1.4247494) q[2];
sx q[2];
rz(-0.82376007) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4321255) q[1];
sx q[1];
rz(-1.8399939) q[1];
sx q[1];
rz(-1.3666735) q[1];
x q[2];
rz(-0.90535935) q[3];
sx q[3];
rz(-1.2631772) q[3];
sx q[3];
rz(0.24505982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.79661757) q[2];
sx q[2];
rz(-1.7285708) q[2];
sx q[2];
rz(-2.3100992) q[2];
rz(1.3384532) q[3];
sx q[3];
rz(-2.0667388) q[3];
sx q[3];
rz(-0.89053806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.7688585) q[0];
sx q[0];
rz(-1.9356198) q[0];
sx q[0];
rz(-0.15705577) q[0];
rz(1.318469) q[1];
sx q[1];
rz(-0.942197) q[1];
sx q[1];
rz(-0.35776055) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6899932) q[0];
sx q[0];
rz(-2.4781961) q[0];
sx q[0];
rz(-2.3729352) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0249457) q[2];
sx q[2];
rz(-2.3966925) q[2];
sx q[2];
rz(-0.91889436) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6727698) q[1];
sx q[1];
rz(-0.92255083) q[1];
sx q[1];
rz(-2.286377) q[1];
x q[2];
rz(1.7292496) q[3];
sx q[3];
rz(-2.1395349) q[3];
sx q[3];
rz(0.67711745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1498108) q[2];
sx q[2];
rz(-1.2573743) q[2];
sx q[2];
rz(3.120976) q[2];
rz(1.2425544) q[3];
sx q[3];
rz(-0.81762448) q[3];
sx q[3];
rz(2.8500565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6108516) q[0];
sx q[0];
rz(-1.956097) q[0];
sx q[0];
rz(-0.32671842) q[0];
rz(-0.84699455) q[1];
sx q[1];
rz(-2.0920483) q[1];
sx q[1];
rz(-2.0832031) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3922472) q[0];
sx q[0];
rz(-1.6405182) q[0];
sx q[0];
rz(2.0850943) q[0];
rz(1.7472668) q[2];
sx q[2];
rz(-2.2024849) q[2];
sx q[2];
rz(1.8827918) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2742334) q[1];
sx q[1];
rz(-1.6489442) q[1];
sx q[1];
rz(1.1876501) q[1];
rz(-pi) q[2];
rz(-2.136738) q[3];
sx q[3];
rz(-1.8396688) q[3];
sx q[3];
rz(0.0026800935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1191001) q[2];
sx q[2];
rz(-2.7969226) q[2];
sx q[2];
rz(1.9692839) q[2];
rz(-3.1398224) q[3];
sx q[3];
rz(-0.5032731) q[3];
sx q[3];
rz(0.9790023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86972791) q[0];
sx q[0];
rz(-0.80369049) q[0];
sx q[0];
rz(-1.8029689) q[0];
rz(1.9610693) q[1];
sx q[1];
rz(-2.0484643) q[1];
sx q[1];
rz(-0.48042935) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6201037) q[0];
sx q[0];
rz(-2.0248981) q[0];
sx q[0];
rz(0.2257077) q[0];
rz(-2.9715638) q[2];
sx q[2];
rz(-1.6126093) q[2];
sx q[2];
rz(-0.60838503) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1611966) q[1];
sx q[1];
rz(-0.56884407) q[1];
sx q[1];
rz(-1.7561164) q[1];
x q[2];
rz(0.35828405) q[3];
sx q[3];
rz(-1.7833774) q[3];
sx q[3];
rz(-2.6019118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.43634811) q[2];
sx q[2];
rz(-2.1465325) q[2];
sx q[2];
rz(-2.056541) q[2];
rz(1.1121701) q[3];
sx q[3];
rz(-2.2472436) q[3];
sx q[3];
rz(2.3280242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28806624) q[0];
sx q[0];
rz(-0.94678322) q[0];
sx q[0];
rz(1.0119525) q[0];
rz(-2.2400253) q[1];
sx q[1];
rz(-0.26600599) q[1];
sx q[1];
rz(2.576135) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6265427) q[0];
sx q[0];
rz(-1.3626271) q[0];
sx q[0];
rz(2.2295537) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1419292) q[2];
sx q[2];
rz(-1.326401) q[2];
sx q[2];
rz(1.7606869) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9931075) q[1];
sx q[1];
rz(-0.7229745) q[1];
sx q[1];
rz(1.430368) q[1];
rz(-pi) q[2];
rz(-1.8826671) q[3];
sx q[3];
rz(-0.33167517) q[3];
sx q[3];
rz(-1.498095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3507583) q[2];
sx q[2];
rz(-1.5459205) q[2];
sx q[2];
rz(-1.2606384) q[2];
rz(0.44220051) q[3];
sx q[3];
rz(-1.0298157) q[3];
sx q[3];
rz(-1.376576) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5653391) q[0];
sx q[0];
rz(-2.2591142) q[0];
sx q[0];
rz(2.2478065) q[0];
rz(1.5743938) q[1];
sx q[1];
rz(-1.4566474) q[1];
sx q[1];
rz(1.7172071) q[1];
rz(1.8509751) q[2];
sx q[2];
rz(-1.1600672) q[2];
sx q[2];
rz(1.5008012) q[2];
rz(-2.9775053) q[3];
sx q[3];
rz(-1.9174977) q[3];
sx q[3];
rz(1.4227572) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
