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
rz(4.6484923) q[1];
sx q[1];
rz(5.3223106) q[1];
sx q[1];
rz(7.9237908) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1849821) q[0];
sx q[0];
rz(-0.8763823) q[0];
sx q[0];
rz(0.25882369) q[0];
rz(-pi) q[1];
rz(-1.1314279) q[2];
sx q[2];
rz(-1.4110392) q[2];
sx q[2];
rz(1.2808094) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5792306) q[1];
sx q[1];
rz(-3.0931614) q[1];
sx q[1];
rz(-0.50244759) q[1];
rz(-pi) q[2];
rz(-2.1330058) q[3];
sx q[3];
rz(-1.2744181) q[3];
sx q[3];
rz(0.62801559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.91244873) q[2];
sx q[2];
rz(-1.392776) q[2];
sx q[2];
rz(-2.8196715) q[2];
rz(-2.1866482) q[3];
sx q[3];
rz(-0.82986444) q[3];
sx q[3];
rz(-1.1436852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2442653) q[0];
sx q[0];
rz(-1.0538415) q[0];
sx q[0];
rz(0.51189297) q[0];
rz(-1.5533252) q[1];
sx q[1];
rz(-2.6542108) q[1];
sx q[1];
rz(-1.1221251) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47488362) q[0];
sx q[0];
rz(-0.32829744) q[0];
sx q[0];
rz(0.31995456) q[0];
rz(1.1412337) q[2];
sx q[2];
rz(-1.5341752) q[2];
sx q[2];
rz(3.1319194) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0341238) q[1];
sx q[1];
rz(-1.144088) q[1];
sx q[1];
rz(0.36860768) q[1];
rz(2.6038325) q[3];
sx q[3];
rz(-2.3590292) q[3];
sx q[3];
rz(1.3937221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3731709) q[2];
sx q[2];
rz(-1.3139407) q[2];
sx q[2];
rz(-0.46305099) q[2];
rz(0.57524663) q[3];
sx q[3];
rz(-1.7762215) q[3];
sx q[3];
rz(-3.0423394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7333882) q[0];
sx q[0];
rz(-2.2938804) q[0];
sx q[0];
rz(-0.43011618) q[0];
rz(0.46562132) q[1];
sx q[1];
rz(-0.72048134) q[1];
sx q[1];
rz(0.63708416) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2475605) q[0];
sx q[0];
rz(-1.5848418) q[0];
sx q[0];
rz(-1.5298903) q[0];
rz(-pi) q[1];
rz(0.78539679) q[2];
sx q[2];
rz(-1.5893418) q[2];
sx q[2];
rz(0.77814276) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8608537) q[1];
sx q[1];
rz(-2.2855084) q[1];
sx q[1];
rz(-1.2398941) q[1];
rz(0.62406059) q[3];
sx q[3];
rz(-2.5284323) q[3];
sx q[3];
rz(2.000119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.78274) q[2];
sx q[2];
rz(-2.095486) q[2];
sx q[2];
rz(2.4165238) q[2];
rz(-1.3310165) q[3];
sx q[3];
rz(-1.015181) q[3];
sx q[3];
rz(0.63841188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26938874) q[0];
sx q[0];
rz(-1.6870455) q[0];
sx q[0];
rz(-1.4440906) q[0];
rz(3.0929502) q[1];
sx q[1];
rz(-1.3261869) q[1];
sx q[1];
rz(2.8526502) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3061482) q[0];
sx q[0];
rz(-1.939538) q[0];
sx q[0];
rz(-0.25935092) q[0];
rz(-0.028325105) q[2];
sx q[2];
rz(-1.420317) q[2];
sx q[2];
rz(0.17720824) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.41449857) q[1];
sx q[1];
rz(-1.3967522) q[1];
sx q[1];
rz(-0.11066173) q[1];
rz(-1.5163317) q[3];
sx q[3];
rz(-2.2353682) q[3];
sx q[3];
rz(-1.8502473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.85663969) q[2];
sx q[2];
rz(-2.608947) q[2];
sx q[2];
rz(0.3802158) q[2];
rz(-0.63816655) q[3];
sx q[3];
rz(-1.9117982) q[3];
sx q[3];
rz(-1.1285454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3087092) q[0];
sx q[0];
rz(-2.5888011) q[0];
sx q[0];
rz(-2.047245) q[0];
rz(0.87567466) q[1];
sx q[1];
rz(-0.62962571) q[1];
sx q[1];
rz(2.0424776) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35661429) q[0];
sx q[0];
rz(-1.7081424) q[0];
sx q[0];
rz(-2.7524968) q[0];
x q[1];
rz(2.7988222) q[2];
sx q[2];
rz(-1.8963328) q[2];
sx q[2];
rz(-2.8936762) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2420038) q[1];
sx q[1];
rz(-2.5676641) q[1];
sx q[1];
rz(-1.8651047) q[1];
rz(-0.91644561) q[3];
sx q[3];
rz(-1.4438085) q[3];
sx q[3];
rz(0.45230745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2431474) q[2];
sx q[2];
rz(-1.3546319) q[2];
sx q[2];
rz(-0.97664991) q[2];
rz(-2.6099033) q[3];
sx q[3];
rz(-2.6952126) q[3];
sx q[3];
rz(-2.4272052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0831182) q[0];
sx q[0];
rz(-2.0396905) q[0];
sx q[0];
rz(-1.8699159) q[0];
rz(2.7187128) q[1];
sx q[1];
rz(-1.6115178) q[1];
sx q[1];
rz(-1.6082825) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23841183) q[0];
sx q[0];
rz(-3.0738214) q[0];
sx q[0];
rz(1.3890024) q[0];
x q[1];
rz(1.6688231) q[2];
sx q[2];
rz(-1.5614126) q[2];
sx q[2];
rz(0.90017747) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6110797) q[1];
sx q[1];
rz(-0.30834282) q[1];
sx q[1];
rz(2.9575696) q[1];
rz(1.815237) q[3];
sx q[3];
rz(-2.3688753) q[3];
sx q[3];
rz(0.48549197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5292458) q[2];
sx q[2];
rz(-1.2924478) q[2];
sx q[2];
rz(0.38522729) q[2];
rz(-3.0569844) q[3];
sx q[3];
rz(-2.707983) q[3];
sx q[3];
rz(2.2658074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92152921) q[0];
sx q[0];
rz(-2.0863057) q[0];
sx q[0];
rz(0.39749843) q[0];
rz(-2.6037604) q[1];
sx q[1];
rz(-2.718524) q[1];
sx q[1];
rz(3.1375111) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83140552) q[0];
sx q[0];
rz(-2.3334529) q[0];
sx q[0];
rz(-2.9648997) q[0];
rz(2.7360544) q[2];
sx q[2];
rz(-2.2161762) q[2];
sx q[2];
rz(2.1359512) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.38975484) q[1];
sx q[1];
rz(-2.1389277) q[1];
sx q[1];
rz(1.1144494) q[1];
rz(-2.3621906) q[3];
sx q[3];
rz(-2.2891697) q[3];
sx q[3];
rz(-1.8441594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6323382) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0179366) q[0];
sx q[0];
rz(-1.2024711) q[0];
sx q[0];
rz(-1.1630455) q[0];
rz(2.842438) q[1];
sx q[1];
rz(-1.2048293) q[1];
sx q[1];
rz(0.97602731) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7892864) q[0];
sx q[0];
rz(-1.8160161) q[0];
sx q[0];
rz(1.6847436) q[0];
rz(-pi) q[1];
rz(2.2783504) q[2];
sx q[2];
rz(-0.23216557) q[2];
sx q[2];
rz(-1.6492998) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8638226) q[1];
sx q[1];
rz(-0.43102095) q[1];
sx q[1];
rz(2.1716539) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6292442) q[3];
sx q[3];
rz(-2.6636925) q[3];
sx q[3];
rz(2.0153139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8396847) q[2];
sx q[2];
rz(-0.32543108) q[2];
sx q[2];
rz(0.1304661) q[2];
rz(1.3236375) q[3];
sx q[3];
rz(-1.8925083) q[3];
sx q[3];
rz(-0.29449335) q[3];
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
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94006222) q[0];
sx q[0];
rz(-0.37684965) q[0];
sx q[0];
rz(-1.4991722) q[0];
rz(0.54010737) q[1];
sx q[1];
rz(-2.3441548) q[1];
sx q[1];
rz(1.3444208) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5169056) q[0];
sx q[0];
rz(-2.6453827) q[0];
sx q[0];
rz(0.93006706) q[0];
rz(-pi) q[1];
rz(-2.4355679) q[2];
sx q[2];
rz(-0.75504061) q[2];
sx q[2];
rz(-1.6598827) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0522103) q[1];
sx q[1];
rz(-1.3213514) q[1];
sx q[1];
rz(0.98813842) q[1];
x q[2];
rz(2.2559861) q[3];
sx q[3];
rz(-0.84784283) q[3];
sx q[3];
rz(0.9566488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6092047) q[2];
sx q[2];
rz(-2.8049073) q[2];
sx q[2];
rz(-1.4538291) q[2];
rz(2.7135571) q[3];
sx q[3];
rz(-1.5513709) q[3];
sx q[3];
rz(-1.9868896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55353272) q[0];
sx q[0];
rz(-0.45192161) q[0];
sx q[0];
rz(1.6499299) q[0];
rz(0.55117575) q[1];
sx q[1];
rz(-1.0124413) q[1];
sx q[1];
rz(-0.59250441) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77757072) q[0];
sx q[0];
rz(-0.3705655) q[0];
sx q[0];
rz(2.4353566) q[0];
rz(-2.3717688) q[2];
sx q[2];
rz(-2.4425659) q[2];
sx q[2];
rz(1.5848643) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4514102) q[1];
sx q[1];
rz(-1.3082778) q[1];
sx q[1];
rz(2.7464944) q[1];
x q[2];
rz(-3.0499712) q[3];
sx q[3];
rz(-2.6406248) q[3];
sx q[3];
rz(2.9660781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7117915) q[2];
sx q[2];
rz(-2.1000803) q[2];
sx q[2];
rz(0.05833021) q[2];
rz(-0.37060261) q[3];
sx q[3];
rz(-0.27675089) q[3];
sx q[3];
rz(-0.0037732865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2881099) q[0];
sx q[0];
rz(-1.7848889) q[0];
sx q[0];
rz(-3.0369192) q[0];
rz(-1.7755605) q[1];
sx q[1];
rz(-2.3401101) q[1];
sx q[1];
rz(0.95536864) q[1];
rz(1.3132533) q[2];
sx q[2];
rz(-2.1291485) q[2];
sx q[2];
rz(2.6026405) q[2];
rz(2.6179553) q[3];
sx q[3];
rz(-2.0515531) q[3];
sx q[3];
rz(2.6691347) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
