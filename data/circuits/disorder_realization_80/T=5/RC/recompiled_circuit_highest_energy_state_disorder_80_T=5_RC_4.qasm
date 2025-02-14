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
rz(3.5340478) q[0];
sx q[0];
rz(3.4253149) q[0];
sx q[0];
rz(14.520159) q[0];
rz(-0.44161931) q[1];
sx q[1];
rz(-1.8752357) q[1];
sx q[1];
rz(-1.9904354) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8609223) q[0];
sx q[0];
rz(-0.09083561) q[0];
sx q[0];
rz(-1.3263741) q[0];
rz(-pi) q[1];
rz(-0.63090905) q[2];
sx q[2];
rz(-2.563347) q[2];
sx q[2];
rz(1.1743682) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.85271231) q[1];
sx q[1];
rz(-1.4433644) q[1];
sx q[1];
rz(2.7136449) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.092716) q[3];
sx q[3];
rz(-1.8303695) q[3];
sx q[3];
rz(-2.5775016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.89897951) q[2];
sx q[2];
rz(-2.2448764) q[2];
sx q[2];
rz(2.8420281) q[2];
rz(0.35563955) q[3];
sx q[3];
rz(-2.1483597) q[3];
sx q[3];
rz(-2.7049844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3886609) q[0];
sx q[0];
rz(-0.54838538) q[0];
sx q[0];
rz(0.27632982) q[0];
rz(-0.07946864) q[1];
sx q[1];
rz(-2.6092968) q[1];
sx q[1];
rz(-0.4963378) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4014028) q[0];
sx q[0];
rz(-1.3200134) q[0];
sx q[0];
rz(2.0987505) q[0];
rz(1.6480108) q[2];
sx q[2];
rz(-1.9081915) q[2];
sx q[2];
rz(1.5959306) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.73442921) q[1];
sx q[1];
rz(-1.8968655) q[1];
sx q[1];
rz(0.77496518) q[1];
x q[2];
rz(-1.7326991) q[3];
sx q[3];
rz(-1.5739958) q[3];
sx q[3];
rz(-0.42409431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3695099) q[2];
sx q[2];
rz(-1.3630867) q[2];
sx q[2];
rz(2.962964) q[2];
rz(-1.6102839) q[3];
sx q[3];
rz(-2.2154) q[3];
sx q[3];
rz(2.3791651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71132201) q[0];
sx q[0];
rz(-2.0549759) q[0];
sx q[0];
rz(1.4104728) q[0];
rz(-1.7544282) q[1];
sx q[1];
rz(-1.6491363) q[1];
sx q[1];
rz(-14/(3*pi)) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1992599) q[0];
sx q[0];
rz(-0.6845277) q[0];
sx q[0];
rz(-3.0897762) q[0];
rz(-pi) q[1];
x q[1];
rz(0.40992592) q[2];
sx q[2];
rz(-0.84884825) q[2];
sx q[2];
rz(-0.97463911) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8285813) q[1];
sx q[1];
rz(-1.269956) q[1];
sx q[1];
rz(-1.1238696) q[1];
rz(-pi) q[2];
rz(-0.36538045) q[3];
sx q[3];
rz(-0.93910223) q[3];
sx q[3];
rz(1.2022579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2014655) q[2];
sx q[2];
rz(-1.0041693) q[2];
sx q[2];
rz(-2.6489769) q[2];
rz(-0.41871873) q[3];
sx q[3];
rz(-2.1188354) q[3];
sx q[3];
rz(-0.59188265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.017460499) q[0];
sx q[0];
rz(-5/(13*pi)) q[0];
sx q[0];
rz(-0.19805743) q[0];
rz(-0.65912229) q[1];
sx q[1];
rz(-2.5810869) q[1];
sx q[1];
rz(-2.9445599) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2360508) q[0];
sx q[0];
rz(-1.6044378) q[0];
sx q[0];
rz(-0.5307401) q[0];
rz(1.8577099) q[2];
sx q[2];
rz(-1.3409145) q[2];
sx q[2];
rz(2.2429466) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9777364) q[1];
sx q[1];
rz(-1.2173554) q[1];
sx q[1];
rz(-2.5720235) q[1];
x q[2];
rz(2.1157589) q[3];
sx q[3];
rz(-1.539388) q[3];
sx q[3];
rz(1.9592913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.78843242) q[2];
sx q[2];
rz(-1.3015231) q[2];
sx q[2];
rz(-1.2737459) q[2];
rz(1.5984009) q[3];
sx q[3];
rz(-1.1929932) q[3];
sx q[3];
rz(2.8852692) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0850328) q[0];
sx q[0];
rz(-1.7340478) q[0];
sx q[0];
rz(-0.10277596) q[0];
rz(-0.13725266) q[1];
sx q[1];
rz(-2.6922701) q[1];
sx q[1];
rz(-1.4385673) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5202839) q[0];
sx q[0];
rz(-0.94693333) q[0];
sx q[0];
rz(-2.413723) q[0];
x q[1];
rz(-2.6784727) q[2];
sx q[2];
rz(-1.6119405) q[2];
sx q[2];
rz(-2.9477313) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9341497) q[1];
sx q[1];
rz(-0.29950073) q[1];
sx q[1];
rz(-2.092176) q[1];
x q[2];
rz(2.4047671) q[3];
sx q[3];
rz(-2.9758436) q[3];
sx q[3];
rz(-0.98893569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.30456257) q[2];
sx q[2];
rz(-1.4745159) q[2];
sx q[2];
rz(-2.4414818) q[2];
rz(2.2434798) q[3];
sx q[3];
rz(-2.5105295) q[3];
sx q[3];
rz(-1.496544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9333618) q[0];
sx q[0];
rz(-0.58657402) q[0];
sx q[0];
rz(1.5775648) q[0];
rz(0.6598407) q[1];
sx q[1];
rz(-2.3520062) q[1];
sx q[1];
rz(-2.1671364) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75622065) q[0];
sx q[0];
rz(-1.5819183) q[0];
sx q[0];
rz(-1.5626989) q[0];
rz(-pi) q[1];
rz(1.7436557) q[2];
sx q[2];
rz(-2.0522016) q[2];
sx q[2];
rz(1.0499894) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2657915) q[1];
sx q[1];
rz(-0.50170021) q[1];
sx q[1];
rz(-0.022345101) q[1];
x q[2];
rz(1.5100689) q[3];
sx q[3];
rz(-1.0332881) q[3];
sx q[3];
rz(0.93858657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.66880354) q[2];
sx q[2];
rz(-0.69837022) q[2];
sx q[2];
rz(-0.93639708) q[2];
rz(-2.2514553) q[3];
sx q[3];
rz(-1.7809159) q[3];
sx q[3];
rz(0.20588017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6137921) q[0];
sx q[0];
rz(-2.9449936) q[0];
sx q[0];
rz(-0.11372456) q[0];
rz(-1.3971036) q[1];
sx q[1];
rz(-1.0662268) q[1];
sx q[1];
rz(3.1166335) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1393046) q[0];
sx q[0];
rz(-1.6253259) q[0];
sx q[0];
rz(1.789458) q[0];
rz(-pi) q[1];
rz(-0.64348797) q[2];
sx q[2];
rz(-1.7805791) q[2];
sx q[2];
rz(-1.6580559) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.61103454) q[1];
sx q[1];
rz(-1.5054346) q[1];
sx q[1];
rz(-2.2863131) q[1];
rz(-pi) q[2];
rz(-0.19784854) q[3];
sx q[3];
rz(-1.5157082) q[3];
sx q[3];
rz(-0.69869374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6156442) q[2];
sx q[2];
rz(-1.7082461) q[2];
sx q[2];
rz(-1.0053267) q[2];
rz(-2.2510236) q[3];
sx q[3];
rz(-1.4054479) q[3];
sx q[3];
rz(-2.4864206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(1.5263851) q[0];
sx q[0];
rz(-0.40281519) q[0];
sx q[0];
rz(2.030754) q[0];
rz(3.0530744) q[1];
sx q[1];
rz(-2.7169777) q[1];
sx q[1];
rz(-1.7875338) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7294069) q[0];
sx q[0];
rz(-0.78347396) q[0];
sx q[0];
rz(-1.7307348) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0646787) q[2];
sx q[2];
rz(-2.8324443) q[2];
sx q[2];
rz(-2.791601) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4072666) q[1];
sx q[1];
rz(-1.3248492) q[1];
sx q[1];
rz(1.6489563) q[1];
rz(0.99127524) q[3];
sx q[3];
rz(-2.2439828) q[3];
sx q[3];
rz(1.1330825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.84305772) q[2];
sx q[2];
rz(-2.3044846) q[2];
sx q[2];
rz(-0.38541547) q[2];
rz(-2.391583) q[3];
sx q[3];
rz(-0.95732006) q[3];
sx q[3];
rz(3.0756522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1239922) q[0];
sx q[0];
rz(-1.0753068) q[0];
sx q[0];
rz(-2.3093834) q[0];
rz(-1.0819134) q[1];
sx q[1];
rz(-2.3257906) q[1];
sx q[1];
rz(1.5218511) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6523217) q[0];
sx q[0];
rz(-1.1025363) q[0];
sx q[0];
rz(-0.44083198) q[0];
rz(-pi) q[1];
rz(2.887313) q[2];
sx q[2];
rz(-1.8449718) q[2];
sx q[2];
rz(2.5667896) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9389006) q[1];
sx q[1];
rz(-2.1317516) q[1];
sx q[1];
rz(-2.2748018) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.91401871) q[3];
sx q[3];
rz(-1.7631712) q[3];
sx q[3];
rz(-3.0466444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.99721432) q[2];
sx q[2];
rz(-1.245456) q[2];
sx q[2];
rz(1.9027556) q[2];
rz(-1.8501806) q[3];
sx q[3];
rz(-1.8700799) q[3];
sx q[3];
rz(-3.0102357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0111897) q[0];
sx q[0];
rz(-2.1881396) q[0];
sx q[0];
rz(2.3977872) q[0];
rz(-0.047608308) q[1];
sx q[1];
rz(-1.1046537) q[1];
sx q[1];
rz(-1.1415175) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6358546) q[0];
sx q[0];
rz(-0.72162745) q[0];
sx q[0];
rz(2.3139335) q[0];
rz(2.52454) q[2];
sx q[2];
rz(-1.2484372) q[2];
sx q[2];
rz(0.53215357) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1117259) q[1];
sx q[1];
rz(-2.2891785) q[1];
sx q[1];
rz(2.5938382) q[1];
rz(0.64472736) q[3];
sx q[3];
rz(-1.7162185) q[3];
sx q[3];
rz(2.4032088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6514674) q[2];
sx q[2];
rz(-0.23568025) q[2];
sx q[2];
rz(0.30533314) q[2];
rz(-1.4281645) q[3];
sx q[3];
rz(-1.3729264) q[3];
sx q[3];
rz(1.8491687) q[3];
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
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3931428) q[0];
sx q[0];
rz(-0.23101692) q[0];
sx q[0];
rz(1.1363181) q[0];
rz(-2.582386) q[1];
sx q[1];
rz(-1.8103841) q[1];
sx q[1];
rz(0.24787535) q[1];
rz(-2.7087173) q[2];
sx q[2];
rz(-0.71472157) q[2];
sx q[2];
rz(1.0403487) q[2];
rz(-2.2918626) q[3];
sx q[3];
rz(-0.17912514) q[3];
sx q[3];
rz(2.9368243) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
