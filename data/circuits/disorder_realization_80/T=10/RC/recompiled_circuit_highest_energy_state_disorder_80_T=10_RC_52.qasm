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
rz(-1.9266204) q[0];
sx q[0];
rz(-1.5705234) q[0];
sx q[0];
rz(-2.7745752) q[0];
rz(3.1103599) q[1];
sx q[1];
rz(-2.2366958) q[1];
sx q[1];
rz(-2.4651405) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.518153) q[0];
sx q[0];
rz(-1.4777967) q[0];
sx q[0];
rz(2.8661845) q[0];
rz(-pi) q[1];
rz(-2.7713369) q[2];
sx q[2];
rz(-0.91846648) q[2];
sx q[2];
rz(2.6877162) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1156173) q[1];
sx q[1];
rz(-1.3968803) q[1];
sx q[1];
rz(2.0439201) q[1];
rz(-pi) q[2];
rz(1.4047926) q[3];
sx q[3];
rz(-0.91445427) q[3];
sx q[3];
rz(2.9350077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5260432) q[2];
sx q[2];
rz(-0.56576663) q[2];
sx q[2];
rz(2.2541798) q[2];
rz(3.060107) q[3];
sx q[3];
rz(-1.8946276) q[3];
sx q[3];
rz(2.2639349) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74251974) q[0];
sx q[0];
rz(-0.44624534) q[0];
sx q[0];
rz(-1.2267858) q[0];
rz(2.1555105) q[1];
sx q[1];
rz(-1.4619275) q[1];
sx q[1];
rz(-0.90423924) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3021224) q[0];
sx q[0];
rz(-1.4581437) q[0];
sx q[0];
rz(2.95999) q[0];
rz(2.9071525) q[2];
sx q[2];
rz(-0.16142217) q[2];
sx q[2];
rz(-1.8983253) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0457669) q[1];
sx q[1];
rz(-0.44039044) q[1];
sx q[1];
rz(-0.38118036) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1911325) q[3];
sx q[3];
rz(-2.2588552) q[3];
sx q[3];
rz(-2.9241274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6282689) q[2];
sx q[2];
rz(-0.29418918) q[2];
sx q[2];
rz(2.7163556) q[2];
rz(0.50906316) q[3];
sx q[3];
rz(-1.1278917) q[3];
sx q[3];
rz(1.4396671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2929448) q[0];
sx q[0];
rz(-3.0631298) q[0];
sx q[0];
rz(-1.5899832) q[0];
rz(-1.6665005) q[1];
sx q[1];
rz(-1.0868797) q[1];
sx q[1];
rz(-2.1045254) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7860216) q[0];
sx q[0];
rz(-1.6986587) q[0];
sx q[0];
rz(0.95604817) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.19107341) q[2];
sx q[2];
rz(-0.5454692) q[2];
sx q[2];
rz(2.7857162) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0817226) q[1];
sx q[1];
rz(-1.3716501) q[1];
sx q[1];
rz(-0.14330602) q[1];
rz(2.7855139) q[3];
sx q[3];
rz(-0.48569187) q[3];
sx q[3];
rz(-0.49957032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.40917778) q[2];
sx q[2];
rz(-2.7517509) q[2];
sx q[2];
rz(-1.2624435) q[2];
rz(0.08517313) q[3];
sx q[3];
rz(-0.44819918) q[3];
sx q[3];
rz(-2.0195154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20268102) q[0];
sx q[0];
rz(-2.1124463) q[0];
sx q[0];
rz(0.19126782) q[0];
rz(0.8910886) q[1];
sx q[1];
rz(-2.8137408) q[1];
sx q[1];
rz(-1.0506312) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3252661) q[0];
sx q[0];
rz(-0.83503316) q[0];
sx q[0];
rz(-0.28006552) q[0];
rz(2.2810535) q[2];
sx q[2];
rz(-1.3293174) q[2];
sx q[2];
rz(-2.2975486) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5393319) q[1];
sx q[1];
rz(-2.6289399) q[1];
sx q[1];
rz(0.67474483) q[1];
rz(0.13728085) q[3];
sx q[3];
rz(-1.0044668) q[3];
sx q[3];
rz(-0.3504741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8656859) q[2];
sx q[2];
rz(-0.17912093) q[2];
sx q[2];
rz(-0.56037819) q[2];
rz(0.066702453) q[3];
sx q[3];
rz(-1.3542465) q[3];
sx q[3];
rz(-2.1393447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10750833) q[0];
sx q[0];
rz(-1.1488687) q[0];
sx q[0];
rz(-2.4565571) q[0];
rz(-0.3512474) q[1];
sx q[1];
rz(-2.5076187) q[1];
sx q[1];
rz(-1.3397071) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9046341) q[0];
sx q[0];
rz(-1.5642421) q[0];
sx q[0];
rz(0.48906869) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1451911) q[2];
sx q[2];
rz(-2.1972924) q[2];
sx q[2];
rz(0.76001924) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9205528) q[1];
sx q[1];
rz(-1.5410081) q[1];
sx q[1];
rz(1.4090621) q[1];
rz(-pi) q[2];
x q[2];
rz(0.48974006) q[3];
sx q[3];
rz(-1.1818172) q[3];
sx q[3];
rz(1.31969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.54905218) q[2];
sx q[2];
rz(-1.668674) q[2];
sx q[2];
rz(2.9259017) q[2];
rz(2.9933764) q[3];
sx q[3];
rz(-0.18311466) q[3];
sx q[3];
rz(-0.7557925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.53629476) q[0];
sx q[0];
rz(-1.7114102) q[0];
sx q[0];
rz(-0.17912616) q[0];
rz(2.1419549) q[1];
sx q[1];
rz(-2.3536286) q[1];
sx q[1];
rz(2.9771908) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9793376) q[0];
sx q[0];
rz(-1.7712403) q[0];
sx q[0];
rz(-2.3810099) q[0];
rz(2.908488) q[2];
sx q[2];
rz(-2.6118661) q[2];
sx q[2];
rz(0.41981439) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.33949125) q[1];
sx q[1];
rz(-0.65971148) q[1];
sx q[1];
rz(-0.011899634) q[1];
x q[2];
rz(-1.3429759) q[3];
sx q[3];
rz(-1.9285255) q[3];
sx q[3];
rz(1.2219971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1228235) q[2];
sx q[2];
rz(-1.5888701) q[2];
sx q[2];
rz(1.7079879) q[2];
rz(-2.150599) q[3];
sx q[3];
rz(-1.3905448) q[3];
sx q[3];
rz(1.2436793) q[3];
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
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7551512) q[0];
sx q[0];
rz(-3.0720818) q[0];
sx q[0];
rz(2.884927) q[0];
rz(3.0811884) q[1];
sx q[1];
rz(-0.81788617) q[1];
sx q[1];
rz(-2.7108257) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11837932) q[0];
sx q[0];
rz(-0.97985615) q[0];
sx q[0];
rz(0.33765461) q[0];
rz(-pi) q[1];
x q[1];
rz(0.68773261) q[2];
sx q[2];
rz(-1.0842676) q[2];
sx q[2];
rz(-0.4358049) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3045945) q[1];
sx q[1];
rz(-0.82436383) q[1];
sx q[1];
rz(-2.5489105) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7458785) q[3];
sx q[3];
rz(-2.4444067) q[3];
sx q[3];
rz(0.28466636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7595235) q[2];
sx q[2];
rz(-2.6117595) q[2];
sx q[2];
rz(1.0158553) q[2];
rz(1.2782512) q[3];
sx q[3];
rz(-1.9081554) q[3];
sx q[3];
rz(2.5281233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.627219) q[0];
sx q[0];
rz(-2.6173499) q[0];
sx q[0];
rz(-2.4216477) q[0];
rz(2.4607957) q[1];
sx q[1];
rz(-0.38377181) q[1];
sx q[1];
rz(-2.0489571) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2075296) q[0];
sx q[0];
rz(-2.060256) q[0];
sx q[0];
rz(-1.6013157) q[0];
rz(-pi) q[1];
rz(-0.60694867) q[2];
sx q[2];
rz(-1.8766093) q[2];
sx q[2];
rz(2.1987777) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5515308) q[1];
sx q[1];
rz(-2.2023337) q[1];
sx q[1];
rz(0.88330357) q[1];
rz(-pi) q[2];
rz(0.44050782) q[3];
sx q[3];
rz(-2.8603209) q[3];
sx q[3];
rz(-0.59556304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0138268) q[2];
sx q[2];
rz(-2.3769145) q[2];
sx q[2];
rz(-2.428425) q[2];
rz(-1.3753043) q[3];
sx q[3];
rz(-1.2262552) q[3];
sx q[3];
rz(-1.3885952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16667287) q[0];
sx q[0];
rz(-2.9627934) q[0];
sx q[0];
rz(1.4270225) q[0];
rz(0.40766454) q[1];
sx q[1];
rz(-1.0382321) q[1];
sx q[1];
rz(2.8020249) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0660865) q[0];
sx q[0];
rz(-1.3097449) q[0];
sx q[0];
rz(-2.2737696) q[0];
rz(-pi) q[1];
rz(2.7251516) q[2];
sx q[2];
rz(-1.526262) q[2];
sx q[2];
rz(0.85516667) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0915874) q[1];
sx q[1];
rz(-1.2500804) q[1];
sx q[1];
rz(1.0343814) q[1];
rz(-pi) q[2];
rz(0.154687) q[3];
sx q[3];
rz(-2.5167234) q[3];
sx q[3];
rz(-2.0578602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0582383) q[2];
sx q[2];
rz(-2.3149172) q[2];
sx q[2];
rz(-0.44462407) q[2];
rz(-2.5104163) q[3];
sx q[3];
rz(-1.3381713) q[3];
sx q[3];
rz(1.9806503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4525962) q[0];
sx q[0];
rz(-1.8190374) q[0];
sx q[0];
rz(3.1153862) q[0];
rz(-1.7474489) q[1];
sx q[1];
rz(-1.988764) q[1];
sx q[1];
rz(-2.0004415) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25309207) q[0];
sx q[0];
rz(-1.9733175) q[0];
sx q[0];
rz(-2.5293468) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0442002) q[2];
sx q[2];
rz(-0.44008128) q[2];
sx q[2];
rz(0.68005622) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2403912) q[1];
sx q[1];
rz(-1.9878042) q[1];
sx q[1];
rz(-2.1621125) q[1];
x q[2];
rz(2.5716173) q[3];
sx q[3];
rz(-1.4432657) q[3];
sx q[3];
rz(0.018477378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.074284) q[2];
sx q[2];
rz(-1.1003541) q[2];
sx q[2];
rz(-1.9791774) q[2];
rz(0.31636604) q[3];
sx q[3];
rz(-1.1917944) q[3];
sx q[3];
rz(1.4828064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3027073) q[0];
sx q[0];
rz(-0.067262983) q[0];
sx q[0];
rz(2.5433232) q[0];
rz(0.0089946714) q[1];
sx q[1];
rz(-0.79575494) q[1];
sx q[1];
rz(1.5473821) q[1];
rz(-1.4228504) q[2];
sx q[2];
rz(-2.2276217) q[2];
sx q[2];
rz(-1.2256636) q[2];
rz(0.056445382) q[3];
sx q[3];
rz(-0.79462449) q[3];
sx q[3];
rz(0.36804646) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
