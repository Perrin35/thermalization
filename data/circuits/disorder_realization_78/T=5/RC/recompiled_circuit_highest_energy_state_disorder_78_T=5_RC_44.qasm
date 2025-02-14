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
rz(1.5155091) q[0];
sx q[0];
rz(-2.8647487) q[0];
sx q[0];
rz(-0.11685637) q[0];
rz(-0.46861831) q[1];
sx q[1];
rz(-2.7873971) q[1];
sx q[1];
rz(2.6503704) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.011176414) q[0];
sx q[0];
rz(-2.6493086) q[0];
sx q[0];
rz(1.4468132) q[0];
rz(-pi) q[1];
rz(-0.38552706) q[2];
sx q[2];
rz(-2.6221199) q[2];
sx q[2];
rz(2.9105718) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3957805) q[1];
sx q[1];
rz(-2.2015101) q[1];
sx q[1];
rz(-1.6935901) q[1];
rz(0.025229021) q[3];
sx q[3];
rz(-1.2275891) q[3];
sx q[3];
rz(-1.1017088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3789856) q[2];
sx q[2];
rz(-0.28756046) q[2];
sx q[2];
rz(0.22988764) q[2];
rz(3.0843132) q[3];
sx q[3];
rz(-0.66904896) q[3];
sx q[3];
rz(0.91276401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3463992) q[0];
sx q[0];
rz(-2.754358) q[0];
sx q[0];
rz(-2.9246395) q[0];
rz(0.03288658) q[1];
sx q[1];
rz(-2.2078728) q[1];
sx q[1];
rz(-2.4864973) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74467662) q[0];
sx q[0];
rz(-1.7707877) q[0];
sx q[0];
rz(2.295664) q[0];
rz(-pi) q[1];
rz(2.6397471) q[2];
sx q[2];
rz(-1.7463533) q[2];
sx q[2];
rz(-1.6278536) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.82258525) q[1];
sx q[1];
rz(-2.2862541) q[1];
sx q[1];
rz(-0.69147488) q[1];
x q[2];
rz(1.7294077) q[3];
sx q[3];
rz(-0.74071032) q[3];
sx q[3];
rz(-0.22601688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.32404262) q[2];
sx q[2];
rz(-0.62411672) q[2];
sx q[2];
rz(-1.8435271) q[2];
rz(1.921418) q[3];
sx q[3];
rz(-1.4379359) q[3];
sx q[3];
rz(0.65742457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87797457) q[0];
sx q[0];
rz(-2.2442696) q[0];
sx q[0];
rz(2.5275912) q[0];
rz(-1.3677431) q[1];
sx q[1];
rz(-2.4847993) q[1];
sx q[1];
rz(-2.6460389) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14392631) q[0];
sx q[0];
rz(-0.94519061) q[0];
sx q[0];
rz(-0.99220522) q[0];
rz(-pi) q[1];
x q[1];
rz(0.96769963) q[2];
sx q[2];
rz(-0.7843547) q[2];
sx q[2];
rz(-1.7627782) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.55807251) q[1];
sx q[1];
rz(-2.8360543) q[1];
sx q[1];
rz(2.2347441) q[1];
rz(-pi) q[2];
rz(2.5088247) q[3];
sx q[3];
rz(-0.55857165) q[3];
sx q[3];
rz(0.66890034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.76501781) q[2];
sx q[2];
rz(-2.0707264) q[2];
sx q[2];
rz(-2.5490254) q[2];
rz(-0.23558922) q[3];
sx q[3];
rz(-2.3841136) q[3];
sx q[3];
rz(-0.45430115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6465103) q[0];
sx q[0];
rz(-0.80146924) q[0];
sx q[0];
rz(0.0058280514) q[0];
rz(-0.5799154) q[1];
sx q[1];
rz(-0.35570759) q[1];
sx q[1];
rz(-2.6204806) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2297008) q[0];
sx q[0];
rz(-1.5661514) q[0];
sx q[0];
rz(1.5779431) q[0];
rz(-2.1717511) q[2];
sx q[2];
rz(-1.5950987) q[2];
sx q[2];
rz(-1.5497335) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.24758928) q[1];
sx q[1];
rz(-0.62452468) q[1];
sx q[1];
rz(-0.7288868) q[1];
rz(-pi) q[2];
rz(-1.6507876) q[3];
sx q[3];
rz(-1.8820401) q[3];
sx q[3];
rz(-2.1224006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2516605) q[2];
sx q[2];
rz(-0.53091383) q[2];
sx q[2];
rz(-0.31822515) q[2];
rz(2.4288154) q[3];
sx q[3];
rz(-2.8409499) q[3];
sx q[3];
rz(-1.6596863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7707959) q[0];
sx q[0];
rz(-0.21368055) q[0];
sx q[0];
rz(0.036291432) q[0];
rz(2.6151784) q[1];
sx q[1];
rz(-1.4585835) q[1];
sx q[1];
rz(-2.859419) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9045893) q[0];
sx q[0];
rz(-0.95247277) q[0];
sx q[0];
rz(0.34070054) q[0];
rz(-2.2660115) q[2];
sx q[2];
rz(-1.1291847) q[2];
sx q[2];
rz(-1.832806) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5314922) q[1];
sx q[1];
rz(-0.86750114) q[1];
sx q[1];
rz(1.478757) q[1];
rz(-2.3822278) q[3];
sx q[3];
rz(-0.96624422) q[3];
sx q[3];
rz(-0.60455632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.37461773) q[2];
sx q[2];
rz(-0.66156113) q[2];
sx q[2];
rz(-3.0401163) q[2];
rz(-0.96595079) q[3];
sx q[3];
rz(-0.92065293) q[3];
sx q[3];
rz(-0.11400338) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5517956) q[0];
sx q[0];
rz(-2.9647201) q[0];
sx q[0];
rz(1.9675323) q[0];
rz(-2.7575098) q[1];
sx q[1];
rz(-1.1840772) q[1];
sx q[1];
rz(-0.039965872) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6718681) q[0];
sx q[0];
rz(-1.3626119) q[0];
sx q[0];
rz(1.2428268) q[0];
rz(0.33127516) q[2];
sx q[2];
rz(-0.9563947) q[2];
sx q[2];
rz(1.5671052) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0607073) q[1];
sx q[1];
rz(-0.97700143) q[1];
sx q[1];
rz(-1.1546561) q[1];
x q[2];
rz(0.21585606) q[3];
sx q[3];
rz(-2.3036968) q[3];
sx q[3];
rz(-0.49666109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.93969718) q[2];
sx q[2];
rz(-2.5494718) q[2];
sx q[2];
rz(-0.8197909) q[2];
rz(0.78153265) q[3];
sx q[3];
rz(-2.235409) q[3];
sx q[3];
rz(1.3801105) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49067295) q[0];
sx q[0];
rz(-1.0395721) q[0];
sx q[0];
rz(1.8248722) q[0];
rz(-1.3680178) q[1];
sx q[1];
rz(-1.9098234) q[1];
sx q[1];
rz(-2.7046611) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5812213) q[0];
sx q[0];
rz(-0.19393714) q[0];
sx q[0];
rz(-1.3790129) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11796911) q[2];
sx q[2];
rz(-1.6802537) q[2];
sx q[2];
rz(0.80889672) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1679042) q[1];
sx q[1];
rz(-1.2120795) q[1];
sx q[1];
rz(-2.3547947) q[1];
x q[2];
rz(-2.223001) q[3];
sx q[3];
rz(-0.64213412) q[3];
sx q[3];
rz(2.0854502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.74066073) q[2];
sx q[2];
rz(-2.0621982) q[2];
sx q[2];
rz(3.0069922) q[2];
rz(1.2234737) q[3];
sx q[3];
rz(-1.7569907) q[3];
sx q[3];
rz(1.1659762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80158919) q[0];
sx q[0];
rz(-2.8497301) q[0];
sx q[0];
rz(-2.959751) q[0];
rz(-2.8002411) q[1];
sx q[1];
rz(-1.1266484) q[1];
sx q[1];
rz(-0.15886074) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0358064) q[0];
sx q[0];
rz(-1.6393568) q[0];
sx q[0];
rz(0.44802702) q[0];
x q[1];
rz(2.7188825) q[2];
sx q[2];
rz(-1.3101995) q[2];
sx q[2];
rz(2.8963331) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3594209) q[1];
sx q[1];
rz(-1.8295994) q[1];
sx q[1];
rz(1.3464298) q[1];
x q[2];
rz(-1.4182866) q[3];
sx q[3];
rz(-1.0220064) q[3];
sx q[3];
rz(-1.2330518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.86026704) q[2];
sx q[2];
rz(-2.211326) q[2];
sx q[2];
rz(-0.88480985) q[2];
rz(1.4191215) q[3];
sx q[3];
rz(-1.0249187) q[3];
sx q[3];
rz(-2.6917246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.849843) q[0];
sx q[0];
rz(-1.2717286) q[0];
sx q[0];
rz(-0.9935317) q[0];
rz(-1.7811071) q[1];
sx q[1];
rz(-0.88881701) q[1];
sx q[1];
rz(-0.55307585) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38163227) q[0];
sx q[0];
rz(-2.6928718) q[0];
sx q[0];
rz(-0.54872175) q[0];
rz(-2.4031248) q[2];
sx q[2];
rz(-1.5212125) q[2];
sx q[2];
rz(2.2253583) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.40402624) q[1];
sx q[1];
rz(-1.8282923) q[1];
sx q[1];
rz(-2.5016258) q[1];
rz(1.9286128) q[3];
sx q[3];
rz(-0.67104895) q[3];
sx q[3];
rz(-0.34085654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9370344) q[2];
sx q[2];
rz(-2.5356346) q[2];
sx q[2];
rz(2.9045612) q[2];
rz(-1.6009181) q[3];
sx q[3];
rz(-0.76434869) q[3];
sx q[3];
rz(-1.9717533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.071526214) q[0];
sx q[0];
rz(-1.5109753) q[0];
sx q[0];
rz(-0.1189098) q[0];
rz(-2.7023244) q[1];
sx q[1];
rz(-1.2989137) q[1];
sx q[1];
rz(-3.0781504) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6574888) q[0];
sx q[0];
rz(-2.2229337) q[0];
sx q[0];
rz(-0.87638292) q[0];
rz(-pi) q[1];
x q[1];
rz(0.62427743) q[2];
sx q[2];
rz(-1.7221525) q[2];
sx q[2];
rz(2.4454012) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0902151) q[1];
sx q[1];
rz(-1.8974638) q[1];
sx q[1];
rz(3.1326381) q[1];
rz(-pi) q[2];
rz(2.1367705) q[3];
sx q[3];
rz(-1.1558954) q[3];
sx q[3];
rz(1.3674629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0110317) q[2];
sx q[2];
rz(-2.06879) q[2];
sx q[2];
rz(1.6434742) q[2];
rz(-0.9324075) q[3];
sx q[3];
rz(-1.1199718) q[3];
sx q[3];
rz(1.2137265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3907923) q[0];
sx q[0];
rz(-1.9600497) q[0];
sx q[0];
rz(2.0586769) q[0];
rz(2.2687268) q[1];
sx q[1];
rz(-2.0582336) q[1];
sx q[1];
rz(2.2179926) q[1];
rz(2.4571669) q[2];
sx q[2];
rz(-1.050907) q[2];
sx q[2];
rz(2.4466865) q[2];
rz(-0.57742248) q[3];
sx q[3];
rz(-1.631177) q[3];
sx q[3];
rz(-2.9384818) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
