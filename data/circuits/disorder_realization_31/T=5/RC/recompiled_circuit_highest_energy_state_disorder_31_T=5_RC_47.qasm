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
rz(1.2827058) q[0];
sx q[0];
rz(-2.1252706) q[0];
sx q[0];
rz(1.8860201) q[0];
rz(1.7805055) q[1];
sx q[1];
rz(-2.1828916) q[1];
sx q[1];
rz(-0.60671848) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9656096) q[0];
sx q[0];
rz(-2.2874766) q[0];
sx q[0];
rz(0.59881439) q[0];
x q[1];
rz(1.6902906) q[2];
sx q[2];
rz(-0.26675376) q[2];
sx q[2];
rz(0.96386749) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7156034) q[1];
sx q[1];
rz(-1.2066325) q[1];
sx q[1];
rz(3.0028572) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9806271) q[3];
sx q[3];
rz(-1.4518629) q[3];
sx q[3];
rz(-1.87784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2029734) q[2];
sx q[2];
rz(-1.9846658) q[2];
sx q[2];
rz(0.17091664) q[2];
rz(-2.0140698) q[3];
sx q[3];
rz(-2.5850962) q[3];
sx q[3];
rz(-3.0881622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4776066) q[0];
sx q[0];
rz(-2.8150616) q[0];
sx q[0];
rz(-2.4531181) q[0];
rz(-1.7636048) q[1];
sx q[1];
rz(-2.1207899) q[1];
sx q[1];
rz(0.14150208) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3162075) q[0];
sx q[0];
rz(-2.0658464) q[0];
sx q[0];
rz(-1.113167) q[0];
rz(-2.9269518) q[2];
sx q[2];
rz(-2.7032489) q[2];
sx q[2];
rz(0.12210309) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5192109) q[1];
sx q[1];
rz(-1.8528717) q[1];
sx q[1];
rz(0.36954986) q[1];
x q[2];
rz(-3.0612437) q[3];
sx q[3];
rz(-1.9544722) q[3];
sx q[3];
rz(-1.2017045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7576393) q[2];
sx q[2];
rz(-2.6944104) q[2];
sx q[2];
rz(0.028707061) q[2];
rz(0.38255295) q[3];
sx q[3];
rz(-1.5111204) q[3];
sx q[3];
rz(-2.9401275) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89556995) q[0];
sx q[0];
rz(-1.0423132) q[0];
sx q[0];
rz(2.7699455) q[0];
rz(2.9314575) q[1];
sx q[1];
rz(-0.34966436) q[1];
sx q[1];
rz(-1.9336611) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71825224) q[0];
sx q[0];
rz(-2.1499942) q[0];
sx q[0];
rz(0.96286003) q[0];
x q[1];
rz(3.0591427) q[2];
sx q[2];
rz(-1.5918613) q[2];
sx q[2];
rz(-2.9631071) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.722586) q[1];
sx q[1];
rz(-1.0204781) q[1];
sx q[1];
rz(2.2645878) q[1];
x q[2];
rz(0.11318993) q[3];
sx q[3];
rz(-0.34750156) q[3];
sx q[3];
rz(2.6828535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.62260425) q[2];
sx q[2];
rz(-3.0486139) q[2];
sx q[2];
rz(0.63817111) q[2];
rz(2.7291164) q[3];
sx q[3];
rz(-1.6715489) q[3];
sx q[3];
rz(-2.0392058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.07688044) q[0];
sx q[0];
rz(-0.91016722) q[0];
sx q[0];
rz(-1.3077211) q[0];
rz(0.67963302) q[1];
sx q[1];
rz(-2.4145587) q[1];
sx q[1];
rz(-1.9576498) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7789202) q[0];
sx q[0];
rz(-2.2244172) q[0];
sx q[0];
rz(-1.9833167) q[0];
x q[1];
rz(-0.45071256) q[2];
sx q[2];
rz(-0.26077429) q[2];
sx q[2];
rz(-0.63577543) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7412926) q[1];
sx q[1];
rz(-1.502599) q[1];
sx q[1];
rz(-0.3505716) q[1];
rz(-pi) q[2];
rz(0.77149646) q[3];
sx q[3];
rz(-2.960254) q[3];
sx q[3];
rz(-2.4490631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1188941) q[2];
sx q[2];
rz(-0.1476295) q[2];
sx q[2];
rz(-1.0559319) q[2];
rz(-2.8583156) q[3];
sx q[3];
rz(-1.8746459) q[3];
sx q[3];
rz(1.7578846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0530171) q[0];
sx q[0];
rz(-0.96240369) q[0];
sx q[0];
rz(-1.7330633) q[0];
rz(-0.22964302) q[1];
sx q[1];
rz(-0.85669986) q[1];
sx q[1];
rz(-1.1913258) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4588683) q[0];
sx q[0];
rz(-2.7153569) q[0];
sx q[0];
rz(1.2221673) q[0];
x q[1];
rz(-1.1896594) q[2];
sx q[2];
rz(-0.32763619) q[2];
sx q[2];
rz(2.0608847) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3670259) q[1];
sx q[1];
rz(-1.8205678) q[1];
sx q[1];
rz(-0.43065177) q[1];
rz(2.2681885) q[3];
sx q[3];
rz(-2.3583989) q[3];
sx q[3];
rz(2.7627373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4952937) q[2];
sx q[2];
rz(-0.36094347) q[2];
sx q[2];
rz(1.6634644) q[2];
rz(-0.16119257) q[3];
sx q[3];
rz(-1.5321923) q[3];
sx q[3];
rz(0.78981367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-1.5021055) q[0];
sx q[0];
rz(-0.4244856) q[0];
sx q[0];
rz(-2.2232527) q[0];
rz(-0.25686747) q[1];
sx q[1];
rz(-2.145642) q[1];
sx q[1];
rz(1.1161944) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6300121) q[0];
sx q[0];
rz(-2.0697547) q[0];
sx q[0];
rz(-0.60006882) q[0];
rz(-pi) q[1];
rz(2.844238) q[2];
sx q[2];
rz(-2.8562244) q[2];
sx q[2];
rz(-3.0653846) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5815365) q[1];
sx q[1];
rz(-1.7651084) q[1];
sx q[1];
rz(-2.8981853) q[1];
rz(-pi) q[2];
rz(-1.7446872) q[3];
sx q[3];
rz(-1.8106204) q[3];
sx q[3];
rz(0.88527977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9548698) q[2];
sx q[2];
rz(-0.8064417) q[2];
sx q[2];
rz(-2.7346129) q[2];
rz(2.5278029) q[3];
sx q[3];
rz(-1.5301306) q[3];
sx q[3];
rz(-0.49155244) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75519049) q[0];
sx q[0];
rz(-2.3004005) q[0];
sx q[0];
rz(2.2156773) q[0];
rz(2.6689802) q[1];
sx q[1];
rz(-2.0874529) q[1];
sx q[1];
rz(2.6714163) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42379728) q[0];
sx q[0];
rz(-1.4873624) q[0];
sx q[0];
rz(2.3071737) q[0];
x q[1];
rz(2.7514303) q[2];
sx q[2];
rz(-0.96218006) q[2];
sx q[2];
rz(0.89649342) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2963841) q[1];
sx q[1];
rz(-1.5344347) q[1];
sx q[1];
rz(2.3567289) q[1];
x q[2];
rz(-2.3005465) q[3];
sx q[3];
rz(-1.5435092) q[3];
sx q[3];
rz(-2.5430665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.65972313) q[2];
sx q[2];
rz(-2.0707776) q[2];
sx q[2];
rz(0.83079633) q[2];
rz(3.0610436) q[3];
sx q[3];
rz(-1.0815257) q[3];
sx q[3];
rz(-0.95675937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.047121) q[0];
sx q[0];
rz(-1.5249277) q[0];
sx q[0];
rz(-0.15810529) q[0];
rz(-2.415601) q[1];
sx q[1];
rz(-2.7114365) q[1];
sx q[1];
rz(-0.37011883) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67664941) q[0];
sx q[0];
rz(-1.9029362) q[0];
sx q[0];
rz(1.667883) q[0];
rz(-pi) q[1];
rz(0.49741715) q[2];
sx q[2];
rz(-0.82225613) q[2];
sx q[2];
rz(2.9679839) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.41851317) q[1];
sx q[1];
rz(-1.9620158) q[1];
sx q[1];
rz(1.9926975) q[1];
rz(-2.1158989) q[3];
sx q[3];
rz(-1.926356) q[3];
sx q[3];
rz(0.73425435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1607754) q[2];
sx q[2];
rz(-1.2951916) q[2];
sx q[2];
rz(-0.13713169) q[2];
rz(1.49617) q[3];
sx q[3];
rz(-0.6936332) q[3];
sx q[3];
rz(0.48367286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5982323) q[0];
sx q[0];
rz(-1.0022663) q[0];
sx q[0];
rz(0.12635669) q[0];
rz(2.5670746) q[1];
sx q[1];
rz(-0.9205598) q[1];
sx q[1];
rz(-2.394302) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80602801) q[0];
sx q[0];
rz(-1.3845516) q[0];
sx q[0];
rz(-1.1769017) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4392082) q[2];
sx q[2];
rz(-2.0364072) q[2];
sx q[2];
rz(1.3872176) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0382749) q[1];
sx q[1];
rz(-2.5468198) q[1];
sx q[1];
rz(-0.86869356) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7036311) q[3];
sx q[3];
rz(-1.3107302) q[3];
sx q[3];
rz(-0.65937926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1025461) q[2];
sx q[2];
rz(-1.1209844) q[2];
sx q[2];
rz(-2.5940564) q[2];
rz(1.3012137) q[3];
sx q[3];
rz(-1.1373212) q[3];
sx q[3];
rz(-3.014452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.309677) q[0];
sx q[0];
rz(-1.7597821) q[0];
sx q[0];
rz(2.5610913) q[0];
rz(-0.25513395) q[1];
sx q[1];
rz(-0.83108035) q[1];
sx q[1];
rz(1.9560122) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0017173926) q[0];
sx q[0];
rz(-2.1809362) q[0];
sx q[0];
rz(2.3185179) q[0];
rz(-pi) q[1];
rz(-0.0015179356) q[2];
sx q[2];
rz(-1.9337144) q[2];
sx q[2];
rz(1.6281106) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6839468) q[1];
sx q[1];
rz(-0.47858979) q[1];
sx q[1];
rz(-2.9475888) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3316292) q[3];
sx q[3];
rz(-0.73430919) q[3];
sx q[3];
rz(-1.7074761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.4127976) q[2];
sx q[2];
rz(-2.2788861) q[2];
sx q[2];
rz(2.1346788) q[2];
rz(1.5236731) q[3];
sx q[3];
rz(-0.98098743) q[3];
sx q[3];
rz(1.9699008) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11378743) q[0];
sx q[0];
rz(-1.2600949) q[0];
sx q[0];
rz(-2.8366198) q[0];
rz(1.8995359) q[1];
sx q[1];
rz(-1.8546974) q[1];
sx q[1];
rz(0.7484662) q[1];
rz(1.1390252) q[2];
sx q[2];
rz(-1.4130269) q[2];
sx q[2];
rz(-2.4561981) q[2];
rz(-0.83132311) q[3];
sx q[3];
rz(-1.1047805) q[3];
sx q[3];
rz(-1.4428231) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
