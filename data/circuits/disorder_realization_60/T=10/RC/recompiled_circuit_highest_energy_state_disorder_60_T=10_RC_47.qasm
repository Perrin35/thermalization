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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95661058) q[0];
sx q[0];
rz(-0.8763823) q[0];
sx q[0];
rz(-0.25882369) q[0];
rz(-pi) q[1];
rz(2.9653984) q[2];
sx q[2];
rz(-2.004188) q[2];
sx q[2];
rz(-2.7769763) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.49351827) q[1];
sx q[1];
rz(-1.5474802) q[1];
sx q[1];
rz(-3.0991395) q[1];
x q[2];
rz(-2.7952173) q[3];
sx q[3];
rz(-1.0358255) q[3];
sx q[3];
rz(-0.76081027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2291439) q[2];
sx q[2];
rz(-1.392776) q[2];
sx q[2];
rz(0.32192117) q[2];
rz(0.95494444) q[3];
sx q[3];
rz(-2.3117282) q[3];
sx q[3];
rz(-1.9979075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2442653) q[0];
sx q[0];
rz(-2.0877512) q[0];
sx q[0];
rz(-0.51189297) q[0];
rz(-1.5533252) q[1];
sx q[1];
rz(-2.6542108) q[1];
sx q[1];
rz(2.0194676) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7417541) q[0];
sx q[0];
rz(-1.4692093) q[0];
sx q[0];
rz(-0.31272696) q[0];
rz(3.1013158) q[2];
sx q[2];
rz(-1.1415408) q[2];
sx q[2];
rz(-1.5778936) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0341238) q[1];
sx q[1];
rz(-1.144088) q[1];
sx q[1];
rz(2.772985) q[1];
rz(-2.6038325) q[3];
sx q[3];
rz(-0.7825635) q[3];
sx q[3];
rz(-1.7478706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.7684218) q[2];
sx q[2];
rz(-1.827652) q[2];
sx q[2];
rz(-0.46305099) q[2];
rz(0.57524663) q[3];
sx q[3];
rz(-1.3653711) q[3];
sx q[3];
rz(3.0423394) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4082044) q[0];
sx q[0];
rz(-2.2938804) q[0];
sx q[0];
rz(-2.7114765) q[0];
rz(-0.46562132) q[1];
sx q[1];
rz(-2.4211113) q[1];
sx q[1];
rz(-2.5045085) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34619949) q[0];
sx q[0];
rz(-0.043248873) q[0];
sx q[0];
rz(1.9016483) q[0];
rz(1.5970206) q[2];
sx q[2];
rz(-0.78557149) q[2];
sx q[2];
rz(-0.77411133) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3780669) q[1];
sx q[1];
rz(-2.3664306) q[1];
sx q[1];
rz(-0.35825348) q[1];
x q[2];
rz(-2.6227642) q[3];
sx q[3];
rz(-1.227855) q[3];
sx q[3];
rz(2.1800621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.78274) q[2];
sx q[2];
rz(-1.0461067) q[2];
sx q[2];
rz(0.72506881) q[2];
rz(-1.8105761) q[3];
sx q[3];
rz(-2.1264117) q[3];
sx q[3];
rz(-2.5031808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8722039) q[0];
sx q[0];
rz(-1.4545472) q[0];
sx q[0];
rz(1.697502) q[0];
rz(0.048642453) q[1];
sx q[1];
rz(-1.8154058) q[1];
sx q[1];
rz(-0.28894249) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83544448) q[0];
sx q[0];
rz(-1.2020547) q[0];
sx q[0];
rz(0.25935092) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.028325105) q[2];
sx q[2];
rz(-1.420317) q[2];
sx q[2];
rz(0.17720824) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.98499456) q[1];
sx q[1];
rz(-2.9356476) q[1];
sx q[1];
rz(-1.0099645) q[1];
rz(-pi) q[2];
rz(-2.4763002) q[3];
sx q[3];
rz(-1.5279309) q[3];
sx q[3];
rz(-0.24584082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.85663969) q[2];
sx q[2];
rz(-2.608947) q[2];
sx q[2];
rz(-2.7613769) q[2];
rz(2.5034261) q[3];
sx q[3];
rz(-1.9117982) q[3];
sx q[3];
rz(2.0130472) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83288348) q[0];
sx q[0];
rz(-2.5888011) q[0];
sx q[0];
rz(-1.0943476) q[0];
rz(-2.265918) q[1];
sx q[1];
rz(-2.5119669) q[1];
sx q[1];
rz(1.099115) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35661429) q[0];
sx q[0];
rz(-1.7081424) q[0];
sx q[0];
rz(-2.7524968) q[0];
rz(2.7988222) q[2];
sx q[2];
rz(-1.2452599) q[2];
sx q[2];
rz(-0.24791644) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.079541072) q[1];
sx q[1];
rz(-1.4126443) q[1];
sx q[1];
rz(-2.1248716) q[1];
rz(-pi) q[2];
rz(1.3640251) q[3];
sx q[3];
rz(-2.4768157) q[3];
sx q[3];
rz(2.1867276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8984453) q[2];
sx q[2];
rz(-1.7869608) q[2];
sx q[2];
rz(-0.97664991) q[2];
rz(-2.6099033) q[3];
sx q[3];
rz(-2.6952126) q[3];
sx q[3];
rz(0.71438742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0584745) q[0];
sx q[0];
rz(-1.1019022) q[0];
sx q[0];
rz(-1.8699159) q[0];
rz(0.42287982) q[1];
sx q[1];
rz(-1.6115178) q[1];
sx q[1];
rz(-1.5333102) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5137702) q[0];
sx q[0];
rz(-1.5830399) q[0];
sx q[0];
rz(1.6374541) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6663867) q[2];
sx q[2];
rz(-0.098473452) q[2];
sx q[2];
rz(2.5661039) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.530513) q[1];
sx q[1];
rz(-0.30834282) q[1];
sx q[1];
rz(-2.9575696) q[1];
x q[2];
rz(2.3284376) q[3];
sx q[3];
rz(-1.7405563) q[3];
sx q[3];
rz(-0.90857279) q[3];
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
rz(0.084608229) q[3];
sx q[3];
rz(-2.707983) q[3];
sx q[3];
rz(-0.87578526) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
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
rz(-0.42306867) q[1];
sx q[1];
rz(-3.1375111) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3101871) q[0];
sx q[0];
rz(-2.3334529) q[0];
sx q[0];
rz(2.9648997) q[0];
rz(-2.7360544) q[2];
sx q[2];
rz(-0.92541646) q[2];
sx q[2];
rz(2.1359512) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7023125) q[1];
sx q[1];
rz(-1.951362) q[1];
sx q[1];
rz(2.5234532) q[1];
rz(0.68304045) q[3];
sx q[3];
rz(-1.0131825) q[3];
sx q[3];
rz(-0.30323365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5092545) q[2];
sx q[2];
rz(-1.1804487) q[2];
sx q[2];
rz(-2.9285367) q[2];
rz(0.80900711) q[3];
sx q[3];
rz(-0.041497858) q[3];
sx q[3];
rz(-1.8607148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1236561) q[0];
sx q[0];
rz(-1.2024711) q[0];
sx q[0];
rz(-1.9785471) q[0];
rz(-0.2991547) q[1];
sx q[1];
rz(-1.2048293) q[1];
sx q[1];
rz(0.97602731) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3487745) q[0];
sx q[0];
rz(-2.8716757) q[0];
sx q[0];
rz(-2.7151373) q[0];
x q[1];
rz(-2.2783504) q[2];
sx q[2];
rz(-2.9094271) q[2];
sx q[2];
rz(-1.6492998) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2175156) q[1];
sx q[1];
rz(-1.922632) q[1];
sx q[1];
rz(2.887243) q[1];
rz(1.3221413) q[3];
sx q[3];
rz(-1.1583405) q[3];
sx q[3];
rz(1.4506884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8396847) q[2];
sx q[2];
rz(-0.32543108) q[2];
sx q[2];
rz(3.0111266) q[2];
rz(1.3236375) q[3];
sx q[3];
rz(-1.2490844) q[3];
sx q[3];
rz(0.29449335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(2.2015304) q[0];
sx q[0];
rz(-2.764743) q[0];
sx q[0];
rz(-1.4991722) q[0];
rz(-2.6014853) q[1];
sx q[1];
rz(-0.79743782) q[1];
sx q[1];
rz(1.7971719) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5169056) q[0];
sx q[0];
rz(-2.6453827) q[0];
sx q[0];
rz(2.2115256) q[0];
rz(0.62144582) q[2];
sx q[2];
rz(-2.0315731) q[2];
sx q[2];
rz(-2.4968392) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1229748) q[1];
sx q[1];
rz(-2.5135165) q[1];
sx q[1];
rz(2.0043892) q[1];
rz(-0.88560652) q[3];
sx q[3];
rz(-0.84784283) q[3];
sx q[3];
rz(-2.1849439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6092047) q[2];
sx q[2];
rz(-0.33668533) q[2];
sx q[2];
rz(-1.6877635) q[2];
rz(2.7135571) q[3];
sx q[3];
rz(-1.5513709) q[3];
sx q[3];
rz(-1.9868896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5880599) q[0];
sx q[0];
rz(-2.689671) q[0];
sx q[0];
rz(1.6499299) q[0];
rz(0.55117575) q[1];
sx q[1];
rz(-2.1291514) q[1];
sx q[1];
rz(0.59250441) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6766178) q[0];
sx q[0];
rz(-1.8080369) q[0];
sx q[0];
rz(2.854191) q[0];
x q[1];
rz(2.1001753) q[2];
sx q[2];
rz(-1.0905078) q[2];
sx q[2];
rz(-0.65435556) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.77280918) q[1];
sx q[1];
rz(-1.9516489) q[1];
sx q[1];
rz(1.2874777) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0499712) q[3];
sx q[3];
rz(-2.6406248) q[3];
sx q[3];
rz(2.9660781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7117915) q[2];
sx q[2];
rz(-1.0415123) q[2];
sx q[2];
rz(-0.05833021) q[2];
rz(-0.37060261) q[3];
sx q[3];
rz(-2.8648418) q[3];
sx q[3];
rz(0.0037732865) q[3];
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
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8534828) q[0];
sx q[0];
rz(-1.3567038) q[0];
sx q[0];
rz(0.10467341) q[0];
rz(1.3660322) q[1];
sx q[1];
rz(-2.3401101) q[1];
sx q[1];
rz(0.95536864) q[1];
rz(-1.8283394) q[2];
sx q[2];
rz(-2.1291485) q[2];
sx q[2];
rz(2.6026405) q[2];
rz(-2.1128863) q[3];
sx q[3];
rz(-1.1114612) q[3];
sx q[3];
rz(-2.3041861) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
