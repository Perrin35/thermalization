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
rz(0.79701841) q[0];
sx q[0];
rz(-2.2198644) q[0];
sx q[0];
rz(-1.8812688) q[0];
rz(-2.4951275) q[1];
sx q[1];
rz(-0.87352455) q[1];
sx q[1];
rz(0.33198196) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8409451) q[0];
sx q[0];
rz(-1.993538) q[0];
sx q[0];
rz(3.1147248) q[0];
rz(-pi) q[1];
rz(-2.7012965) q[2];
sx q[2];
rz(-0.74661359) q[2];
sx q[2];
rz(-0.25737112) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.60377893) q[1];
sx q[1];
rz(-0.75944501) q[1];
sx q[1];
rz(3.095597) q[1];
rz(-pi) q[2];
rz(-2.6993887) q[3];
sx q[3];
rz(-1.5118268) q[3];
sx q[3];
rz(-1.9249431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.476764) q[2];
sx q[2];
rz(-0.9330743) q[2];
sx q[2];
rz(-0.19922166) q[2];
rz(-2.9557989) q[3];
sx q[3];
rz(-1.4150861) q[3];
sx q[3];
rz(-1.7004405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81714565) q[0];
sx q[0];
rz(-2.6528093) q[0];
sx q[0];
rz(2.4031438) q[0];
rz(-1.6959408) q[1];
sx q[1];
rz(-1.7414469) q[1];
sx q[1];
rz(-2.6587528) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33979499) q[0];
sx q[0];
rz(-1.4301824) q[0];
sx q[0];
rz(1.4184345) q[0];
rz(-pi) q[1];
rz(2.0245527) q[2];
sx q[2];
rz(-2.7575709) q[2];
sx q[2];
rz(-2.6817317) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9571217) q[1];
sx q[1];
rz(-2.1163242) q[1];
sx q[1];
rz(-1.0218744) q[1];
rz(2.9779997) q[3];
sx q[3];
rz(-1.0181659) q[3];
sx q[3];
rz(1.9528509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8740251) q[2];
sx q[2];
rz(-1.4635307) q[2];
sx q[2];
rz(1.2919424) q[2];
rz(1.1235631) q[3];
sx q[3];
rz(-0.70190391) q[3];
sx q[3];
rz(-1.4152214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8181151) q[0];
sx q[0];
rz(-0.98387843) q[0];
sx q[0];
rz(-1.6997319) q[0];
rz(-1.2280751) q[1];
sx q[1];
rz(-1.2308729) q[1];
sx q[1];
rz(-1.0736116) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.613205) q[0];
sx q[0];
rz(-0.79314418) q[0];
sx q[0];
rz(2.6499477) q[0];
rz(0.91543897) q[2];
sx q[2];
rz(-1.4129593) q[2];
sx q[2];
rz(-2.137568) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1086928) q[1];
sx q[1];
rz(-1.6039023) q[1];
sx q[1];
rz(0.86825235) q[1];
rz(-pi) q[2];
rz(-2.8183636) q[3];
sx q[3];
rz(-2.6628837) q[3];
sx q[3];
rz(-1.1532825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.29232612) q[2];
sx q[2];
rz(-2.8385415) q[2];
sx q[2];
rz(-2.1185875) q[2];
rz(-0.060221378) q[3];
sx q[3];
rz(-1.9950208) q[3];
sx q[3];
rz(-2.4165418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66134727) q[0];
sx q[0];
rz(-2.9954973) q[0];
sx q[0];
rz(0.51937854) q[0];
rz(-0.6005148) q[1];
sx q[1];
rz(-1.0888211) q[1];
sx q[1];
rz(2.3868938) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31904991) q[0];
sx q[0];
rz(-0.62915914) q[0];
sx q[0];
rz(-0.99790093) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2738931) q[2];
sx q[2];
rz(-1.461173) q[2];
sx q[2];
rz(-0.90849691) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0097451) q[1];
sx q[1];
rz(-2.5452217) q[1];
sx q[1];
rz(-0.98747298) q[1];
x q[2];
rz(0.56086991) q[3];
sx q[3];
rz(-1.9114693) q[3];
sx q[3];
rz(2.7531227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.34919136) q[2];
sx q[2];
rz(-2.3442522) q[2];
sx q[2];
rz(0.1758197) q[2];
rz(1.1261806) q[3];
sx q[3];
rz(-1.3577941) q[3];
sx q[3];
rz(-0.064182909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.675451) q[0];
sx q[0];
rz(-1.6459246) q[0];
sx q[0];
rz(-2.5597036) q[0];
rz(1.9582845) q[1];
sx q[1];
rz(-2.4300523) q[1];
sx q[1];
rz(1.2540832) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7039827) q[0];
sx q[0];
rz(-2.2334198) q[0];
sx q[0];
rz(-0.41351701) q[0];
rz(-pi) q[1];
rz(0.78127258) q[2];
sx q[2];
rz(-0.63129497) q[2];
sx q[2];
rz(-2.6897813) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5595793) q[1];
sx q[1];
rz(-1.2767643) q[1];
sx q[1];
rz(-0.21802417) q[1];
rz(-pi) q[2];
x q[2];
rz(0.33361907) q[3];
sx q[3];
rz(-0.52894178) q[3];
sx q[3];
rz(1.2490619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1648272) q[2];
sx q[2];
rz(-0.81563121) q[2];
sx q[2];
rz(-2.7867479) q[2];
rz(-1.9918848) q[3];
sx q[3];
rz(-2.0668273) q[3];
sx q[3];
rz(-0.84158516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-0.069816291) q[0];
sx q[0];
rz(-2.8700097) q[0];
sx q[0];
rz(0.040104453) q[0];
rz(-1.6768203) q[1];
sx q[1];
rz(-2.1720839) q[1];
sx q[1];
rz(-0.47223314) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67481638) q[0];
sx q[0];
rz(-0.32186723) q[0];
sx q[0];
rz(0.9573413) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5750999) q[2];
sx q[2];
rz(-0.65281463) q[2];
sx q[2];
rz(-3.0497361) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.66989723) q[1];
sx q[1];
rz(-1.9006839) q[1];
sx q[1];
rz(-2.8389588) q[1];
rz(1.8566441) q[3];
sx q[3];
rz(-2.3116391) q[3];
sx q[3];
rz(-2.2835116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8670696) q[2];
sx q[2];
rz(-0.7663061) q[2];
sx q[2];
rz(-1.4806032) q[2];
rz(-3.1001422) q[3];
sx q[3];
rz(-2.4292414) q[3];
sx q[3];
rz(-0.97431549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90254766) q[0];
sx q[0];
rz(-2.3426265) q[0];
sx q[0];
rz(0.7582742) q[0];
rz(0.43765086) q[1];
sx q[1];
rz(-1.9270908) q[1];
sx q[1];
rz(-0.95380107) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8361369) q[0];
sx q[0];
rz(-1.9528332) q[0];
sx q[0];
rz(-3.0925095) q[0];
rz(-pi) q[1];
rz(2.4720947) q[2];
sx q[2];
rz(-0.84247103) q[2];
sx q[2];
rz(-2.5491722) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.055102417) q[1];
sx q[1];
rz(-2.2728572) q[1];
sx q[1];
rz(0.081900077) q[1];
x q[2];
rz(2.8607868) q[3];
sx q[3];
rz(-1.1698517) q[3];
sx q[3];
rz(1.795648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7268251) q[2];
sx q[2];
rz(-2.813952) q[2];
sx q[2];
rz(2.7437575) q[2];
rz(-1.2658524) q[3];
sx q[3];
rz(-1.419302) q[3];
sx q[3];
rz(0.46441594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-1.9809113) q[0];
sx q[0];
rz(-0.92781624) q[0];
sx q[0];
rz(2.8371147) q[0];
rz(1.132384) q[1];
sx q[1];
rz(-0.67190036) q[1];
sx q[1];
rz(2.0679881) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.96655) q[0];
sx q[0];
rz(-1.5538408) q[0];
sx q[0];
rz(1.5477563) q[0];
rz(-1.9307641) q[2];
sx q[2];
rz(-1.9641293) q[2];
sx q[2];
rz(2.8434812) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4661387) q[1];
sx q[1];
rz(-2.4917648) q[1];
sx q[1];
rz(-1.1625579) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6210446) q[3];
sx q[3];
rz(-2.2138688) q[3];
sx q[3];
rz(-2.1724043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.73211804) q[2];
sx q[2];
rz(-1.858859) q[2];
sx q[2];
rz(3.0916302) q[2];
rz(0.91160715) q[3];
sx q[3];
rz(-2.215569) q[3];
sx q[3];
rz(0.91134206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6066345) q[0];
sx q[0];
rz(-1.3953403) q[0];
sx q[0];
rz(0.362679) q[0];
rz(0.72049385) q[1];
sx q[1];
rz(-0.81169218) q[1];
sx q[1];
rz(1.9627176) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6828109) q[0];
sx q[0];
rz(-2.1401554) q[0];
sx q[0];
rz(-2.1853133) q[0];
rz(1.6495398) q[2];
sx q[2];
rz(-1.8905565) q[2];
sx q[2];
rz(-1.3516689) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9423811) q[1];
sx q[1];
rz(-1.4742715) q[1];
sx q[1];
rz(1.7418899) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5963327) q[3];
sx q[3];
rz(-1.2125748) q[3];
sx q[3];
rz(-2.6521366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.0044340213) q[2];
sx q[2];
rz(-2.4538071) q[2];
sx q[2];
rz(-1.9384109) q[2];
rz(-0.016599003) q[3];
sx q[3];
rz(-1.6417475) q[3];
sx q[3];
rz(-1.871292) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5430629) q[0];
sx q[0];
rz(-0.5235343) q[0];
sx q[0];
rz(-1.6981) q[0];
rz(-0.47830018) q[1];
sx q[1];
rz(-2.4595478) q[1];
sx q[1];
rz(-1.9313448) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2624274) q[0];
sx q[0];
rz(-2.0602003) q[0];
sx q[0];
rz(2.4995245) q[0];
x q[1];
rz(2.4912676) q[2];
sx q[2];
rz(-1.4286094) q[2];
sx q[2];
rz(-1.9498169) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0134971) q[1];
sx q[1];
rz(-1.6229318) q[1];
sx q[1];
rz(-0.32491046) q[1];
rz(-pi) q[2];
x q[2];
rz(0.5860255) q[3];
sx q[3];
rz(-0.3908433) q[3];
sx q[3];
rz(1.217921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.18834194) q[2];
sx q[2];
rz(-2.5276999) q[2];
sx q[2];
rz(-0.71315145) q[2];
rz(0.6330511) q[3];
sx q[3];
rz(-0.96708599) q[3];
sx q[3];
rz(2.2658394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0093832) q[0];
sx q[0];
rz(-1.8108981) q[0];
sx q[0];
rz(-0.43999408) q[0];
rz(-2.412759) q[1];
sx q[1];
rz(-0.76580096) q[1];
sx q[1];
rz(-2.0892807) q[1];
rz(-1.960878) q[2];
sx q[2];
rz(-1.4643535) q[2];
sx q[2];
rz(2.6753078) q[2];
rz(1.9464372) q[3];
sx q[3];
rz(-0.79938625) q[3];
sx q[3];
rz(1.1342794) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
