OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8785414) q[0];
sx q[0];
rz(-0.6113373) q[0];
sx q[0];
rz(1.4933458) q[0];
rz(-2.2892294) q[1];
sx q[1];
rz(-1.747793) q[1];
sx q[1];
rz(1.2300904) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6540007) q[0];
sx q[0];
rz(-1.0941668) q[0];
sx q[0];
rz(-2.0311902) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9622494) q[2];
sx q[2];
rz(-1.9916996) q[2];
sx q[2];
rz(-2.1564623) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7609357) q[1];
sx q[1];
rz(-2.0208997) q[1];
sx q[1];
rz(-0.37827079) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3268747) q[3];
sx q[3];
rz(-2.8496242) q[3];
sx q[3];
rz(2.1413745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7386231) q[2];
sx q[2];
rz(-1.3961184) q[2];
sx q[2];
rz(1.9782861) q[2];
rz(2.2033384) q[3];
sx q[3];
rz(-0.80621755) q[3];
sx q[3];
rz(1.7136542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52629483) q[0];
sx q[0];
rz(-1.9087003) q[0];
sx q[0];
rz(-1.9799318) q[0];
rz(1.42234) q[1];
sx q[1];
rz(-1.6484345) q[1];
sx q[1];
rz(-0.086440451) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0924226) q[0];
sx q[0];
rz(-2.0960652) q[0];
sx q[0];
rz(-2.9591333) q[0];
rz(1.3474355) q[2];
sx q[2];
rz(-1.0901951) q[2];
sx q[2];
rz(-1.640412) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0965496) q[1];
sx q[1];
rz(-3.109906) q[1];
sx q[1];
rz(-2.6331054) q[1];
x q[2];
rz(1.2041041) q[3];
sx q[3];
rz(-2.9584998) q[3];
sx q[3];
rz(2.2568373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.57261673) q[2];
sx q[2];
rz(-0.47875753) q[2];
sx q[2];
rz(3.0968481) q[2];
rz(-0.71988002) q[3];
sx q[3];
rz(-1.6241112) q[3];
sx q[3];
rz(-1.450479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8986664) q[0];
sx q[0];
rz(-1.9654322) q[0];
sx q[0];
rz(-1.2258919) q[0];
rz(2.2236845) q[1];
sx q[1];
rz(-1.8501015) q[1];
sx q[1];
rz(-1.2845405) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51581406) q[0];
sx q[0];
rz(-0.5042133) q[0];
sx q[0];
rz(1.4674241) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.011083885) q[2];
sx q[2];
rz(-2.5199118) q[2];
sx q[2];
rz(-1.7116837) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0513902) q[1];
sx q[1];
rz(-1.0388684) q[1];
sx q[1];
rz(2.4064785) q[1];
rz(-pi) q[2];
rz(2.7906832) q[3];
sx q[3];
rz(-0.65334645) q[3];
sx q[3];
rz(0.86337435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7104177) q[2];
sx q[2];
rz(-1.9937036) q[2];
sx q[2];
rz(0.61719027) q[2];
rz(-0.91444531) q[3];
sx q[3];
rz(-0.72943288) q[3];
sx q[3];
rz(0.35703737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0964088) q[0];
sx q[0];
rz(-3.0953188) q[0];
sx q[0];
rz(-2.1531877) q[0];
rz(-1.660396) q[1];
sx q[1];
rz(-2.587187) q[1];
sx q[1];
rz(-2.5344417) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5051712) q[0];
sx q[0];
rz(-1.8431516) q[0];
sx q[0];
rz(0.26557458) q[0];
rz(-pi) q[1];
rz(0.64035236) q[2];
sx q[2];
rz(-1.8812064) q[2];
sx q[2];
rz(-2.7057757) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1478575) q[1];
sx q[1];
rz(-0.58489913) q[1];
sx q[1];
rz(-0.92928504) q[1];
x q[2];
rz(-1.3903343) q[3];
sx q[3];
rz(-1.4752098) q[3];
sx q[3];
rz(-1.4845899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0561169) q[2];
sx q[2];
rz(-1.7530707) q[2];
sx q[2];
rz(0.2307387) q[2];
rz(-0.76656109) q[3];
sx q[3];
rz(-0.14403382) q[3];
sx q[3];
rz(-1.8421596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6799927) q[0];
sx q[0];
rz(-1.5138641) q[0];
sx q[0];
rz(-3.0738714) q[0];
rz(-1.8374775) q[1];
sx q[1];
rz(-1.291052) q[1];
sx q[1];
rz(-1.7052604) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50698603) q[0];
sx q[0];
rz(-1.9546967) q[0];
sx q[0];
rz(-2.3070154) q[0];
x q[1];
rz(-2.0642989) q[2];
sx q[2];
rz(-1.4736325) q[2];
sx q[2];
rz(-0.0070126931) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.98187145) q[1];
sx q[1];
rz(-1.2471732) q[1];
sx q[1];
rz(1.0713599) q[1];
rz(-1.0338289) q[3];
sx q[3];
rz(-2.7269533) q[3];
sx q[3];
rz(-1.4495759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0962254) q[2];
sx q[2];
rz(-2.2787091) q[2];
sx q[2];
rz(-2.6048342) q[2];
rz(-1.2196352) q[3];
sx q[3];
rz(-2.2370971) q[3];
sx q[3];
rz(1.9513244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9959975) q[0];
sx q[0];
rz(-0.97766101) q[0];
sx q[0];
rz(-1.4033432) q[0];
rz(-2.6373236) q[1];
sx q[1];
rz(-1.206617) q[1];
sx q[1];
rz(-0.23955841) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58695108) q[0];
sx q[0];
rz(-1.3901781) q[0];
sx q[0];
rz(-0.15968542) q[0];
rz(-pi) q[1];
rz(-0.94953612) q[2];
sx q[2];
rz(-0.64709787) q[2];
sx q[2];
rz(-2.4280649) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.18232432) q[1];
sx q[1];
rz(-2.2367918) q[1];
sx q[1];
rz(-1.4592501) q[1];
rz(-pi) q[2];
x q[2];
rz(0.40294874) q[3];
sx q[3];
rz(-2.0411754) q[3];
sx q[3];
rz(1.345548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2704894) q[2];
sx q[2];
rz(-1.6051382) q[2];
sx q[2];
rz(-0.25320539) q[2];
rz(-2.1829055) q[3];
sx q[3];
rz(-2.2264693) q[3];
sx q[3];
rz(-1.2747214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36509982) q[0];
sx q[0];
rz(-0.72117844) q[0];
sx q[0];
rz(1.7505296) q[0];
rz(-2.5513388) q[1];
sx q[1];
rz(-1.2140467) q[1];
sx q[1];
rz(-2.6388157) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25994008) q[0];
sx q[0];
rz(-1.5852196) q[0];
sx q[0];
rz(1.4653066) q[0];
rz(-2.4847772) q[2];
sx q[2];
rz(-1.8699099) q[2];
sx q[2];
rz(0.5123891) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.16044438) q[1];
sx q[1];
rz(-2.2878621) q[1];
sx q[1];
rz(-2.8129548) q[1];
rz(-pi) q[2];
rz(-0.54347968) q[3];
sx q[3];
rz(-1.8414652) q[3];
sx q[3];
rz(1.9668996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7949924) q[2];
sx q[2];
rz(-1.9575926) q[2];
sx q[2];
rz(2.6896175) q[2];
rz(2.82708) q[3];
sx q[3];
rz(-1.9032685) q[3];
sx q[3];
rz(-0.048129169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5101584) q[0];
sx q[0];
rz(-2.2619814) q[0];
sx q[0];
rz(-1.6360224) q[0];
rz(0.10558852) q[1];
sx q[1];
rz(-2.0629203) q[1];
sx q[1];
rz(0.67965913) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1435232) q[0];
sx q[0];
rz(-2.1532602) q[0];
sx q[0];
rz(-0.67654718) q[0];
rz(-1.0318884) q[2];
sx q[2];
rz(-2.4499345) q[2];
sx q[2];
rz(-2.0641278) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.41808203) q[1];
sx q[1];
rz(-0.46996221) q[1];
sx q[1];
rz(1.5638007) q[1];
rz(-pi) q[2];
x q[2];
rz(1.345349) q[3];
sx q[3];
rz(-1.7525008) q[3];
sx q[3];
rz(-0.34442201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.52123657) q[2];
sx q[2];
rz(-2.607589) q[2];
sx q[2];
rz(-2.3109069) q[2];
rz(-2.1031117) q[3];
sx q[3];
rz(-2.4887648) q[3];
sx q[3];
rz(-1.7415107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0251544) q[0];
sx q[0];
rz(-1.9898299) q[0];
sx q[0];
rz(1.0333767) q[0];
rz(-1.0721463) q[1];
sx q[1];
rz(-1.1786345) q[1];
sx q[1];
rz(-2.5796366) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28918761) q[0];
sx q[0];
rz(-2.4187208) q[0];
sx q[0];
rz(1.7815352) q[0];
x q[1];
rz(-1.1216036) q[2];
sx q[2];
rz(-0.70708924) q[2];
sx q[2];
rz(2.0769151) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.6889557) q[1];
sx q[1];
rz(-1.9105043) q[1];
sx q[1];
rz(0.79381408) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4212332) q[3];
sx q[3];
rz(-3.0834167) q[3];
sx q[3];
rz(-0.10567372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0434887) q[2];
sx q[2];
rz(-1.6489112) q[2];
sx q[2];
rz(-1.4081504) q[2];
rz(3.0812541) q[3];
sx q[3];
rz(-1.8294168) q[3];
sx q[3];
rz(-1.8797125) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14335808) q[0];
sx q[0];
rz(-2.3682605) q[0];
sx q[0];
rz(1.7711357) q[0];
rz(-2.1342318) q[1];
sx q[1];
rz(-1.7528563) q[1];
sx q[1];
rz(3.0000684) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25989448) q[0];
sx q[0];
rz(-2.0260677) q[0];
sx q[0];
rz(-1.0819525) q[0];
rz(0.74574377) q[2];
sx q[2];
rz(-0.088938449) q[2];
sx q[2];
rz(0.76751417) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.255573) q[1];
sx q[1];
rz(-1.951978) q[1];
sx q[1];
rz(-2.0440434) q[1];
x q[2];
rz(0.022079682) q[3];
sx q[3];
rz(-2.0661269) q[3];
sx q[3];
rz(-1.9324404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0274028) q[2];
sx q[2];
rz(-1.9582615) q[2];
sx q[2];
rz(-1.2664504) q[2];
rz(3.0053084) q[3];
sx q[3];
rz(-2.2253021) q[3];
sx q[3];
rz(0.19792476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4895353) q[0];
sx q[0];
rz(-1.0132402) q[0];
sx q[0];
rz(-1.5860438) q[0];
rz(0.96792211) q[1];
sx q[1];
rz(-2.617901) q[1];
sx q[1];
rz(-1.0760457) q[1];
rz(3.0332312) q[2];
sx q[2];
rz(-0.81365267) q[2];
sx q[2];
rz(-1.1968422) q[2];
rz(-3.089681) q[3];
sx q[3];
rz(-1.3811227) q[3];
sx q[3];
rz(1.5863252) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
