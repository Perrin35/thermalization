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
rz(0.76437104) q[0];
sx q[0];
rz(1.8005014) q[0];
sx q[0];
rz(10.330248) q[0];
rz(1.5068997) q[1];
sx q[1];
rz(-2.1807179) q[1];
sx q[1];
rz(1.5009872) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.577271) q[0];
sx q[0];
rz(-0.73350302) q[0];
sx q[0];
rz(-1.2726239) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1314279) q[2];
sx q[2];
rz(-1.4110392) q[2];
sx q[2];
rz(-1.8607832) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5792306) q[1];
sx q[1];
rz(-0.048431245) q[1];
sx q[1];
rz(0.50244759) q[1];
x q[2];
rz(-1.0085868) q[3];
sx q[3];
rz(-1.2744181) q[3];
sx q[3];
rz(2.5135771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.91244873) q[2];
sx q[2];
rz(-1.7488166) q[2];
sx q[2];
rz(2.8196715) q[2];
rz(-0.95494444) q[3];
sx q[3];
rz(-2.3117282) q[3];
sx q[3];
rz(1.9979075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8973273) q[0];
sx q[0];
rz(-1.0538415) q[0];
sx q[0];
rz(-2.6296997) q[0];
rz(-1.5533252) q[1];
sx q[1];
rz(-2.6542108) q[1];
sx q[1];
rz(2.0194676) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13817912) q[0];
sx q[0];
rz(-1.2597359) q[0];
sx q[0];
rz(1.4640693) q[0];
x q[1];
rz(-1.6585412) q[2];
sx q[2];
rz(-0.43102396) q[2];
sx q[2];
rz(-1.6602248) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1074688) q[1];
sx q[1];
rz(-1.144088) q[1];
sx q[1];
rz(-0.36860768) q[1];
rz(2.4347794) q[3];
sx q[3];
rz(-1.9403096) q[3];
sx q[3];
rz(2.564424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.7684218) q[2];
sx q[2];
rz(-1.3139407) q[2];
sx q[2];
rz(-0.46305099) q[2];
rz(-2.566346) q[3];
sx q[3];
rz(-1.3653711) q[3];
sx q[3];
rz(-0.099253207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4082044) q[0];
sx q[0];
rz(-0.84771228) q[0];
sx q[0];
rz(-2.7114765) q[0];
rz(0.46562132) q[1];
sx q[1];
rz(-2.4211113) q[1];
sx q[1];
rz(-0.63708416) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67733902) q[0];
sx q[0];
rz(-1.6116983) q[0];
sx q[0];
rz(-0.014057191) q[0];
rz(-pi) q[1];
rz(0.026224296) q[2];
sx q[2];
rz(-2.3560239) q[2];
sx q[2];
rz(2.3303967) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8608537) q[1];
sx q[1];
rz(-2.2855084) q[1];
sx q[1];
rz(1.2398941) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6227642) q[3];
sx q[3];
rz(-1.227855) q[3];
sx q[3];
rz(-0.9615306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.78274) q[2];
sx q[2];
rz(-1.0461067) q[2];
sx q[2];
rz(-2.4165238) q[2];
rz(-1.8105761) q[3];
sx q[3];
rz(-1.015181) q[3];
sx q[3];
rz(2.5031808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26938874) q[0];
sx q[0];
rz(-1.6870455) q[0];
sx q[0];
rz(1.697502) q[0];
rz(3.0929502) q[1];
sx q[1];
rz(-1.8154058) q[1];
sx q[1];
rz(-2.8526502) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5015857) q[0];
sx q[0];
rz(-1.3292392) q[0];
sx q[0];
rz(1.9511186) q[0];
rz(-pi) q[1];
x q[1];
rz(0.028325105) q[2];
sx q[2];
rz(-1.7212756) q[2];
sx q[2];
rz(0.17720824) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1565981) q[1];
sx q[1];
rz(-0.20594507) q[1];
sx q[1];
rz(-1.0099645) q[1];
rz(-1.6252609) q[3];
sx q[3];
rz(-0.90622444) q[3];
sx q[3];
rz(1.2913454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.284953) q[2];
sx q[2];
rz(-2.608947) q[2];
sx q[2];
rz(-2.7613769) q[2];
rz(-0.63816655) q[3];
sx q[3];
rz(-1.2297945) q[3];
sx q[3];
rz(1.1285454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83288348) q[0];
sx q[0];
rz(-0.55279151) q[0];
sx q[0];
rz(1.0943476) q[0];
rz(0.87567466) q[1];
sx q[1];
rz(-2.5119669) q[1];
sx q[1];
rz(-2.0424776) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6051145) q[0];
sx q[0];
rz(-2.7301359) q[0];
sx q[0];
rz(-0.34939975) q[0];
rz(-pi) q[1];
rz(0.34277041) q[2];
sx q[2];
rz(-1.8963328) q[2];
sx q[2];
rz(2.8936762) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2420038) q[1];
sx q[1];
rz(-2.5676641) q[1];
sx q[1];
rz(-1.276488) q[1];
rz(-pi) q[2];
x q[2];
rz(0.91644561) q[3];
sx q[3];
rz(-1.6977842) q[3];
sx q[3];
rz(0.45230745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2431474) q[2];
sx q[2];
rz(-1.7869608) q[2];
sx q[2];
rz(0.97664991) q[2];
rz(0.53168932) q[3];
sx q[3];
rz(-0.44638005) q[3];
sx q[3];
rz(-0.71438742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0831182) q[0];
sx q[0];
rz(-2.0396905) q[0];
sx q[0];
rz(-1.2716768) q[0];
rz(2.7187128) q[1];
sx q[1];
rz(-1.6115178) q[1];
sx q[1];
rz(-1.6082825) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5137702) q[0];
sx q[0];
rz(-1.5585527) q[0];
sx q[0];
rz(-1.6374541) q[0];
rz(-pi) q[1];
rz(-1.4727696) q[2];
sx q[2];
rz(-1.5614126) q[2];
sx q[2];
rz(0.90017747) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.530513) q[1];
sx q[1];
rz(-2.8332498) q[1];
sx q[1];
rz(-0.18402305) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9098784) q[3];
sx q[3];
rz(-2.3149256) q[3];
sx q[3];
rz(-2.3208914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5292458) q[2];
sx q[2];
rz(-1.8491448) q[2];
sx q[2];
rz(-2.7563654) q[2];
rz(-3.0569844) q[3];
sx q[3];
rz(-2.707983) q[3];
sx q[3];
rz(2.2658074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92152921) q[0];
sx q[0];
rz(-2.0863057) q[0];
sx q[0];
rz(0.39749843) q[0];
rz(2.6037604) q[1];
sx q[1];
rz(-2.718524) q[1];
sx q[1];
rz(0.0040815512) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5631112) q[0];
sx q[0];
rz(-0.77881506) q[0];
sx q[0];
rz(-1.7527197) q[0];
rz(-0.88433684) q[2];
sx q[2];
rz(-1.8914273) q[2];
sx q[2];
rz(-2.8291631) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.34985182) q[1];
sx q[1];
rz(-2.428973) q[1];
sx q[1];
rz(2.5373759) q[1];
rz(0.8936196) q[3];
sx q[3];
rz(-2.1359518) q[3];
sx q[3];
rz(-0.86097417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5092545) q[2];
sx q[2];
rz(-1.961144) q[2];
sx q[2];
rz(0.21305591) q[2];
rz(-0.80900711) q[3];
sx q[3];
rz(-0.041497858) q[3];
sx q[3];
rz(1.8607148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1236561) q[0];
sx q[0];
rz(-1.9391215) q[0];
sx q[0];
rz(1.1630455) q[0];
rz(0.2991547) q[1];
sx q[1];
rz(-1.2048293) q[1];
sx q[1];
rz(-0.97602731) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9508787) q[0];
sx q[0];
rz(-1.4602721) q[0];
sx q[0];
rz(0.24675639) q[0];
rz(-1.7485745) q[2];
sx q[2];
rz(-1.4206829) q[2];
sx q[2];
rz(-2.3689388) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.73607549) q[1];
sx q[1];
rz(-1.3323405) q[1];
sx q[1];
rz(1.9333436) q[1];
rz(-1.8194514) q[3];
sx q[3];
rz(-1.9832522) q[3];
sx q[3];
rz(1.6909042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.30190793) q[2];
sx q[2];
rz(-0.32543108) q[2];
sx q[2];
rz(-0.1304661) q[2];
rz(-1.3236375) q[3];
sx q[3];
rz(-1.8925083) q[3];
sx q[3];
rz(0.29449335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94006222) q[0];
sx q[0];
rz(-0.37684965) q[0];
sx q[0];
rz(1.4991722) q[0];
rz(2.6014853) q[1];
sx q[1];
rz(-0.79743782) q[1];
sx q[1];
rz(-1.7971719) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8136052) q[0];
sx q[0];
rz(-1.179197) q[0];
sx q[0];
rz(-2.8285976) q[0];
rz(-2.1189519) q[2];
sx q[2];
rz(-1.022199) q[2];
sx q[2];
rz(-2.5238069) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0522103) q[1];
sx q[1];
rz(-1.3213514) q[1];
sx q[1];
rz(0.98813842) q[1];
rz(-pi) q[2];
x q[2];
rz(0.88560652) q[3];
sx q[3];
rz(-2.2937498) q[3];
sx q[3];
rz(0.9566488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6092047) q[2];
sx q[2];
rz(-2.8049073) q[2];
sx q[2];
rz(-1.4538291) q[2];
rz(-0.42803556) q[3];
sx q[3];
rz(-1.5513709) q[3];
sx q[3];
rz(-1.9868896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55353272) q[0];
sx q[0];
rz(-0.45192161) q[0];
sx q[0];
rz(1.6499299) q[0];
rz(-0.55117575) q[1];
sx q[1];
rz(-2.1291514) q[1];
sx q[1];
rz(2.5490882) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4649749) q[0];
sx q[0];
rz(-1.3335557) q[0];
sx q[0];
rz(-0.28740164) q[0];
rz(-pi) q[1];
rz(-1.0414174) q[2];
sx q[2];
rz(-1.0905078) q[2];
sx q[2];
rz(2.4872371) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4372651) q[1];
sx q[1];
rz(-0.47050899) q[1];
sx q[1];
rz(-0.60948845) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5207401) q[3];
sx q[3];
rz(-1.0721237) q[3];
sx q[3];
rz(0.071144516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.42980117) q[2];
sx q[2];
rz(-2.1000803) q[2];
sx q[2];
rz(3.0832624) q[2];
rz(0.37060261) q[3];
sx q[3];
rz(-2.8648418) q[3];
sx q[3];
rz(3.1378194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8534828) q[0];
sx q[0];
rz(-1.3567038) q[0];
sx q[0];
rz(0.10467341) q[0];
rz(-1.3660322) q[1];
sx q[1];
rz(-0.80148253) q[1];
sx q[1];
rz(-2.186224) q[1];
rz(2.7544207) q[2];
sx q[2];
rz(-0.6091112) q[2];
sx q[2];
rz(2.1412639) q[2];
rz(-0.52363734) q[3];
sx q[3];
rz(-2.0515531) q[3];
sx q[3];
rz(2.6691347) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
