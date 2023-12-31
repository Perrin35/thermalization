OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5390227) q[0];
sx q[0];
rz(-2.5780926) q[0];
sx q[0];
rz(2.6846057) q[0];
rz(-0.62178388) q[1];
sx q[1];
rz(-0.68067247) q[1];
sx q[1];
rz(-1.8656123) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2109283) q[0];
sx q[0];
rz(-1.352248) q[0];
sx q[0];
rz(-1.4320089) q[0];
rz(-3.0511191) q[2];
sx q[2];
rz(-1.4960559) q[2];
sx q[2];
rz(0.59404101) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9717325) q[1];
sx q[1];
rz(-1.5105376) q[1];
sx q[1];
rz(1.7737284) q[1];
rz(2.071322) q[3];
sx q[3];
rz(-0.71159092) q[3];
sx q[3];
rz(0.69858944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1093381) q[2];
sx q[2];
rz(-2.2685752) q[2];
sx q[2];
rz(-2.9764552) q[2];
rz(1.5103546) q[3];
sx q[3];
rz(-2.813208) q[3];
sx q[3];
rz(0.74222773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6363268) q[0];
sx q[0];
rz(-2.9716182) q[0];
sx q[0];
rz(-3.0965565) q[0];
rz(-2.8024407) q[1];
sx q[1];
rz(-1.3749326) q[1];
sx q[1];
rz(0.80274686) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0720285) q[0];
sx q[0];
rz(-1.9575798) q[0];
sx q[0];
rz(-2.8974429) q[0];
x q[1];
rz(-3.0375508) q[2];
sx q[2];
rz(-1.0696971) q[2];
sx q[2];
rz(0.83545557) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.75623679) q[1];
sx q[1];
rz(-1.5192351) q[1];
sx q[1];
rz(2.1897584) q[1];
rz(1.734415) q[3];
sx q[3];
rz(-2.46582) q[3];
sx q[3];
rz(2.8676652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7654483) q[2];
sx q[2];
rz(-2.5120698) q[2];
sx q[2];
rz(1.7155898) q[2];
rz(2.5358893) q[3];
sx q[3];
rz(-2.1798445) q[3];
sx q[3];
rz(-0.081993016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95214343) q[0];
sx q[0];
rz(-1.9316297) q[0];
sx q[0];
rz(2.8564575) q[0];
rz(-0.51672283) q[1];
sx q[1];
rz(-1.6500094) q[1];
sx q[1];
rz(1.3396938) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.08683603) q[0];
sx q[0];
rz(-1.1491547) q[0];
sx q[0];
rz(-0.7936759) q[0];
rz(-pi) q[1];
rz(2.8324098) q[2];
sx q[2];
rz(-0.90484607) q[2];
sx q[2];
rz(0.14112976) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9078319) q[1];
sx q[1];
rz(-1.3843752) q[1];
sx q[1];
rz(-0.16252653) q[1];
rz(-2.6684746) q[3];
sx q[3];
rz(-0.74672532) q[3];
sx q[3];
rz(2.8603922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.67119917) q[2];
sx q[2];
rz(-0.931804) q[2];
sx q[2];
rz(1.7335256) q[2];
rz(-2.9584598) q[3];
sx q[3];
rz(-2.1268842) q[3];
sx q[3];
rz(-0.071921913) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6614439) q[0];
sx q[0];
rz(-2.792795) q[0];
sx q[0];
rz(2.9273422) q[0];
rz(-2.7125773) q[1];
sx q[1];
rz(-2.0206101) q[1];
sx q[1];
rz(2.4750211) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4144856) q[0];
sx q[0];
rz(-1.3654728) q[0];
sx q[0];
rz(1.7617102) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33118403) q[2];
sx q[2];
rz(-2.3996155) q[2];
sx q[2];
rz(-0.77510288) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2506634) q[1];
sx q[1];
rz(-0.86859413) q[1];
sx q[1];
rz(2.3198747) q[1];
x q[2];
rz(-1.7188844) q[3];
sx q[3];
rz(-1.7184966) q[3];
sx q[3];
rz(1.5642017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.2864705) q[2];
sx q[2];
rz(-0.21229599) q[2];
sx q[2];
rz(-2.4812223) q[2];
rz(1.1874229) q[3];
sx q[3];
rz(-1.238845) q[3];
sx q[3];
rz(0.56046265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51263556) q[0];
sx q[0];
rz(-1.9218788) q[0];
sx q[0];
rz(-0.83706013) q[0];
rz(2.191026) q[1];
sx q[1];
rz(-2.4730813) q[1];
sx q[1];
rz(1.1594835) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54407185) q[0];
sx q[0];
rz(-0.87771767) q[0];
sx q[0];
rz(-1.1580628) q[0];
x q[1];
rz(-1.7972838) q[2];
sx q[2];
rz(-0.61033568) q[2];
sx q[2];
rz(-3.0938873) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0598462) q[1];
sx q[1];
rz(-1.0330327) q[1];
sx q[1];
rz(0.50077011) q[1];
rz(-pi) q[2];
x q[2];
rz(2.106835) q[3];
sx q[3];
rz(-1.0039819) q[3];
sx q[3];
rz(-2.00622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.435047) q[2];
sx q[2];
rz(-1.1112301) q[2];
sx q[2];
rz(-1.8699899) q[2];
rz(-2.0909677) q[3];
sx q[3];
rz(-1.2579066) q[3];
sx q[3];
rz(0.064855382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5746675) q[0];
sx q[0];
rz(-2.8820679) q[0];
sx q[0];
rz(2.0264453) q[0];
rz(0.043958157) q[1];
sx q[1];
rz(-1.3479439) q[1];
sx q[1];
rz(3.040722) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2642919) q[0];
sx q[0];
rz(-1.722324) q[0];
sx q[0];
rz(-0.63440462) q[0];
x q[1];
rz(0.51668824) q[2];
sx q[2];
rz(-2.5942205) q[2];
sx q[2];
rz(-1.9043497) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.069177376) q[1];
sx q[1];
rz(-2.8848856) q[1];
sx q[1];
rz(-2.2290345) q[1];
rz(-pi) q[2];
rz(0.50562596) q[3];
sx q[3];
rz(-2.6856832) q[3];
sx q[3];
rz(-1.0491919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4042525) q[2];
sx q[2];
rz(-1.6254144) q[2];
sx q[2];
rz(-0.92528382) q[2];
rz(2.960079) q[3];
sx q[3];
rz(-2.4003024) q[3];
sx q[3];
rz(1.293175) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7589384) q[0];
sx q[0];
rz(-2.721334) q[0];
sx q[0];
rz(-0.96310258) q[0];
rz(1.5902279) q[1];
sx q[1];
rz(-2.1236877) q[1];
sx q[1];
rz(2.4538453) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15935414) q[0];
sx q[0];
rz(-1.9040477) q[0];
sx q[0];
rz(-2.9100145) q[0];
x q[1];
rz(-1.2528) q[2];
sx q[2];
rz(-1.7858292) q[2];
sx q[2];
rz(-3.0766578) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.63989418) q[1];
sx q[1];
rz(-1.5940897) q[1];
sx q[1];
rz(0.1293907) q[1];
x q[2];
rz(-1.5958105) q[3];
sx q[3];
rz(-2.5731312) q[3];
sx q[3];
rz(-0.75077885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8035651) q[2];
sx q[2];
rz(-1.4146718) q[2];
sx q[2];
rz(-0.93758279) q[2];
rz(0.98085105) q[3];
sx q[3];
rz(-0.3347446) q[3];
sx q[3];
rz(-2.3576095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5279609) q[0];
sx q[0];
rz(-0.95454803) q[0];
sx q[0];
rz(3.0623073) q[0];
rz(2.0203363) q[1];
sx q[1];
rz(-0.93833485) q[1];
sx q[1];
rz(1.1987196) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0872902) q[0];
sx q[0];
rz(-1.1371005) q[0];
sx q[0];
rz(1.4951697) q[0];
x q[1];
rz(3.0555658) q[2];
sx q[2];
rz(-1.2089529) q[2];
sx q[2];
rz(-1.8334243) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.63003507) q[1];
sx q[1];
rz(-0.50543284) q[1];
sx q[1];
rz(0.10314421) q[1];
rz(-pi) q[2];
rz(-1.8317322) q[3];
sx q[3];
rz(-2.0203777) q[3];
sx q[3];
rz(-0.98228067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.85598677) q[2];
sx q[2];
rz(-0.14582835) q[2];
sx q[2];
rz(0.78906995) q[2];
rz(-0.35813913) q[3];
sx q[3];
rz(-1.5538244) q[3];
sx q[3];
rz(-1.8607128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.566074) q[0];
sx q[0];
rz(-1.386336) q[0];
sx q[0];
rz(0.58832204) q[0];
rz(-0.026780216) q[1];
sx q[1];
rz(-0.90590042) q[1];
sx q[1];
rz(2.2081597) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3073472) q[0];
sx q[0];
rz(-1.630162) q[0];
sx q[0];
rz(3.0349915) q[0];
rz(-pi) q[1];
x q[1];
rz(2.80745) q[2];
sx q[2];
rz(-1.9311937) q[2];
sx q[2];
rz(2.6956218) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.48298207) q[1];
sx q[1];
rz(-2.254199) q[1];
sx q[1];
rz(-1.5233836) q[1];
rz(-pi) q[2];
rz(0.45456072) q[3];
sx q[3];
rz(-0.70485669) q[3];
sx q[3];
rz(2.9398033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0788706) q[2];
sx q[2];
rz(-2.3214919) q[2];
sx q[2];
rz(2.771647) q[2];
rz(2.2423819) q[3];
sx q[3];
rz(-1.4827385) q[3];
sx q[3];
rz(-0.95329154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9545492) q[0];
sx q[0];
rz(-1.8148913) q[0];
sx q[0];
rz(-0.6533587) q[0];
rz(1.5006784) q[1];
sx q[1];
rz(-1.4705642) q[1];
sx q[1];
rz(-0.18383372) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.375658) q[0];
sx q[0];
rz(-1.4544393) q[0];
sx q[0];
rz(0.12136202) q[0];
rz(-pi) q[1];
x q[1];
rz(0.83689883) q[2];
sx q[2];
rz(-0.95015929) q[2];
sx q[2];
rz(-1.9377973) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.373133) q[1];
sx q[1];
rz(-0.8994973) q[1];
sx q[1];
rz(0.21233227) q[1];
x q[2];
rz(3.1179908) q[3];
sx q[3];
rz(-2.5581103) q[3];
sx q[3];
rz(-0.49707801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.9188149) q[2];
sx q[2];
rz(-0.99514014) q[2];
sx q[2];
rz(-0.20239057) q[2];
rz(-1.0732132) q[3];
sx q[3];
rz(-1.9066633) q[3];
sx q[3];
rz(-1.7741268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68207537) q[0];
sx q[0];
rz(-0.85953241) q[0];
sx q[0];
rz(-0.55534242) q[0];
rz(-2.7783685) q[1];
sx q[1];
rz(-1.3365311) q[1];
sx q[1];
rz(2.8844759) q[1];
rz(2.6828962) q[2];
sx q[2];
rz(-1.9445322) q[2];
sx q[2];
rz(-1.4945488) q[2];
rz(-1.8300874) q[3];
sx q[3];
rz(-2.5355831) q[3];
sx q[3];
rz(-2.2021962) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
