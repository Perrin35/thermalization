OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.6289829) q[0];
sx q[0];
rz(4.1689685) q[0];
sx q[0];
rz(12.203759) q[0];
rz(0.91712046) q[1];
sx q[1];
rz(-0.4904823) q[1];
sx q[1];
rz(-0.34163707) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5255614) q[0];
sx q[0];
rz(-1.6934868) q[0];
sx q[0];
rz(0.37567715) q[0];
x q[1];
rz(-0.59092893) q[2];
sx q[2];
rz(-1.8701631) q[2];
sx q[2];
rz(2.9654944) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.29400723) q[1];
sx q[1];
rz(-2.5348186) q[1];
sx q[1];
rz(-1.0311544) q[1];
rz(-pi) q[2];
rz(1.0757252) q[3];
sx q[3];
rz(-0.44701156) q[3];
sx q[3];
rz(1.3017553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5477156) q[2];
sx q[2];
rz(-1.1527858) q[2];
sx q[2];
rz(-2.1851052) q[2];
rz(2.9603413) q[3];
sx q[3];
rz(-0.58877188) q[3];
sx q[3];
rz(-1.6199002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0797121) q[0];
sx q[0];
rz(-0.065608874) q[0];
sx q[0];
rz(1.1454426) q[0];
rz(-1.068813) q[1];
sx q[1];
rz(-0.83199465) q[1];
sx q[1];
rz(2.4157445) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9024591) q[0];
sx q[0];
rz(-2.4665138) q[0];
sx q[0];
rz(0.59160561) q[0];
rz(2.0306573) q[2];
sx q[2];
rz(-2.1150555) q[2];
sx q[2];
rz(2.0366675) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5686381) q[1];
sx q[1];
rz(-1.8796088) q[1];
sx q[1];
rz(-3.1152578) q[1];
rz(-0.29048357) q[3];
sx q[3];
rz(-2.4986914) q[3];
sx q[3];
rz(0.089240616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.48250616) q[2];
sx q[2];
rz(-1.7616452) q[2];
sx q[2];
rz(2.144311) q[2];
rz(-1.7287792) q[3];
sx q[3];
rz(-1.0935067) q[3];
sx q[3];
rz(0.93878448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(2.7230351) q[0];
sx q[0];
rz(-2.1756873) q[0];
sx q[0];
rz(2.0130656) q[0];
rz(-1.8380802) q[1];
sx q[1];
rz(-0.729527) q[1];
sx q[1];
rz(2.2361141) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3826686) q[0];
sx q[0];
rz(-1.7610468) q[0];
sx q[0];
rz(-2.8959136) q[0];
rz(-pi) q[1];
rz(-2.682914) q[2];
sx q[2];
rz(-2.3346666) q[2];
sx q[2];
rz(0.57070953) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9943774) q[1];
sx q[1];
rz(-0.71502393) q[1];
sx q[1];
rz(1.6836402) q[1];
rz(-2.2842992) q[3];
sx q[3];
rz(-1.5382574) q[3];
sx q[3];
rz(-1.1710222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3811615) q[2];
sx q[2];
rz(-2.9734549) q[2];
sx q[2];
rz(-0.9872438) q[2];
rz(-0.045771249) q[3];
sx q[3];
rz(-1.2922623) q[3];
sx q[3];
rz(-2.9614017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10953294) q[0];
sx q[0];
rz(-0.68234545) q[0];
sx q[0];
rz(3.0155244) q[0];
rz(2.9571422) q[1];
sx q[1];
rz(-1.0686864) q[1];
sx q[1];
rz(-1.9741845) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6895034) q[0];
sx q[0];
rz(-2.3514682) q[0];
sx q[0];
rz(0.28496123) q[0];
rz(-pi) q[1];
rz(2.8235769) q[2];
sx q[2];
rz(-0.52581767) q[2];
sx q[2];
rz(1.9981245) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3525225) q[1];
sx q[1];
rz(-1.8348502) q[1];
sx q[1];
rz(-2.4697484) q[1];
x q[2];
rz(-0.91058369) q[3];
sx q[3];
rz(-2.108421) q[3];
sx q[3];
rz(1.1066574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2771153) q[2];
sx q[2];
rz(-0.41431674) q[2];
sx q[2];
rz(0.63151044) q[2];
rz(2.946092) q[3];
sx q[3];
rz(-1.2771527) q[3];
sx q[3];
rz(2.7235532) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9773848) q[0];
sx q[0];
rz(-0.13810869) q[0];
sx q[0];
rz(0.50267977) q[0];
rz(-2.0282822) q[1];
sx q[1];
rz(-2.138442) q[1];
sx q[1];
rz(-0.46708333) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57468092) q[0];
sx q[0];
rz(-1.3710638) q[0];
sx q[0];
rz(-1.9499705) q[0];
rz(-pi) q[1];
rz(2.5555243) q[2];
sx q[2];
rz(-1.35891) q[2];
sx q[2];
rz(0.52473247) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0568697) q[1];
sx q[1];
rz(-1.9922171) q[1];
sx q[1];
rz(0.23878581) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.763837) q[3];
sx q[3];
rz(-0.69460624) q[3];
sx q[3];
rz(-2.4823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6749394) q[2];
sx q[2];
rz(-1.249524) q[2];
sx q[2];
rz(-2.6494027) q[2];
rz(-0.28218937) q[3];
sx q[3];
rz(-1.7084833) q[3];
sx q[3];
rz(1.5757489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2714587) q[0];
sx q[0];
rz(-0.54015714) q[0];
sx q[0];
rz(3.0555994) q[0];
rz(1.7135235) q[1];
sx q[1];
rz(-1.0079039) q[1];
sx q[1];
rz(1.6451947) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7487885) q[0];
sx q[0];
rz(-1.7668056) q[0];
sx q[0];
rz(-1.0402354) q[0];
x q[1];
rz(1.1804579) q[2];
sx q[2];
rz(-0.57363735) q[2];
sx q[2];
rz(-0.2804642) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3277668) q[1];
sx q[1];
rz(-1.5944949) q[1];
sx q[1];
rz(1.2162672) q[1];
rz(-pi) q[2];
rz(-3.1268678) q[3];
sx q[3];
rz(-0.71299362) q[3];
sx q[3];
rz(-1.7409489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.66203403) q[2];
sx q[2];
rz(-2.5223314) q[2];
sx q[2];
rz(2.7453444) q[2];
rz(0.59404343) q[3];
sx q[3];
rz(-2.0500573) q[3];
sx q[3];
rz(1.5007639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2632161) q[0];
sx q[0];
rz(-0.57290572) q[0];
sx q[0];
rz(-1.0923882) q[0];
rz(1.7842402) q[1];
sx q[1];
rz(-1.4211632) q[1];
sx q[1];
rz(-3.1299652) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8893338) q[0];
sx q[0];
rz(-2.6703435) q[0];
sx q[0];
rz(2.528119) q[0];
x q[1];
rz(1.5989223) q[2];
sx q[2];
rz(-1.5453494) q[2];
sx q[2];
rz(2.5487119) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7203622) q[1];
sx q[1];
rz(-1.0248529) q[1];
sx q[1];
rz(-0.027688428) q[1];
rz(-pi) q[2];
x q[2];
rz(2.743268) q[3];
sx q[3];
rz(-2.4279804) q[3];
sx q[3];
rz(-2.6834965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6860883) q[2];
sx q[2];
rz(-0.65767613) q[2];
sx q[2];
rz(2.8536076) q[2];
rz(1.1503495) q[3];
sx q[3];
rz(-1.8201613) q[3];
sx q[3];
rz(2.5434125) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8312254) q[0];
sx q[0];
rz(-2.6103525) q[0];
sx q[0];
rz(1.9294552) q[0];
rz(0.37777004) q[1];
sx q[1];
rz(-1.2021844) q[1];
sx q[1];
rz(-1.4935965) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0069107) q[0];
sx q[0];
rz(-1.784076) q[0];
sx q[0];
rz(2.8365305) q[0];
x q[1];
rz(1.1475032) q[2];
sx q[2];
rz(-0.14483843) q[2];
sx q[2];
rz(2.7763979) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.32840604) q[1];
sx q[1];
rz(-1.9771736) q[1];
sx q[1];
rz(-0.74192585) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7842403) q[3];
sx q[3];
rz(-2.6297914) q[3];
sx q[3];
rz(-1.7970999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3528072) q[2];
sx q[2];
rz(-2.0264503) q[2];
sx q[2];
rz(1.2517694) q[2];
rz(-2.2552323) q[3];
sx q[3];
rz(-0.76615196) q[3];
sx q[3];
rz(3.0322976) q[3];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4733474) q[0];
sx q[0];
rz(-2.1003523) q[0];
sx q[0];
rz(2.9633203) q[0];
rz(3.0687304) q[1];
sx q[1];
rz(-2.5477414) q[1];
sx q[1];
rz(1.9624306) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.060886325) q[0];
sx q[0];
rz(-1.0651677) q[0];
sx q[0];
rz(-1.9067184) q[0];
rz(-pi) q[1];
rz(-2.6369008) q[2];
sx q[2];
rz(-1.0070649) q[2];
sx q[2];
rz(1.6187514) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.8930248) q[1];
sx q[1];
rz(-2.5556459) q[1];
sx q[1];
rz(-2.663661) q[1];
rz(-1.0696899) q[3];
sx q[3];
rz(-1.2341098) q[3];
sx q[3];
rz(1.6285553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.69958413) q[2];
sx q[2];
rz(-2.1506491) q[2];
sx q[2];
rz(-2.1098095) q[2];
rz(-1.7992841) q[3];
sx q[3];
rz(-1.7248025) q[3];
sx q[3];
rz(-0.22527307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.853302) q[0];
sx q[0];
rz(-1.0422215) q[0];
sx q[0];
rz(-0.061696079) q[0];
rz(1.8348947) q[1];
sx q[1];
rz(-1.6625762) q[1];
sx q[1];
rz(-2.0526989) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93931224) q[0];
sx q[0];
rz(-1.5081524) q[0];
sx q[0];
rz(-1.5059727) q[0];
rz(-0.18386545) q[2];
sx q[2];
rz(-0.92391787) q[2];
sx q[2];
rz(-1.136214) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.1354462) q[1];
sx q[1];
rz(-1.0137614) q[1];
sx q[1];
rz(2.1250493) q[1];
x q[2];
rz(1.0763361) q[3];
sx q[3];
rz(-0.38443243) q[3];
sx q[3];
rz(1.2151375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6178199) q[2];
sx q[2];
rz(-2.4536295) q[2];
sx q[2];
rz(-3.0715959) q[2];
rz(-2.2648515) q[3];
sx q[3];
rz(-1.8165959) q[3];
sx q[3];
rz(2.783412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4927647) q[0];
sx q[0];
rz(-1.6236826) q[0];
sx q[0];
rz(0.59326011) q[0];
rz(-1.6246673) q[1];
sx q[1];
rz(-2.0943874) q[1];
sx q[1];
rz(2.692683) q[1];
rz(2.2286227) q[2];
sx q[2];
rz(-1.2428478) q[2];
sx q[2];
rz(-1.7075677) q[2];
rz(-1.9061448) q[3];
sx q[3];
rz(-1.5318346) q[3];
sx q[3];
rz(-0.58605015) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
