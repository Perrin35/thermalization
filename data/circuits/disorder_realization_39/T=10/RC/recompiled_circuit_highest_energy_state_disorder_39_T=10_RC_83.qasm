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
rz(-2.5290639) q[0];
sx q[0];
rz(-0.6331442) q[0];
sx q[0];
rz(-0.36344114) q[0];
rz(-0.58703047) q[1];
sx q[1];
rz(-2.6122004) q[1];
sx q[1];
rz(1.5076948) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76985659) q[0];
sx q[0];
rz(-1.5205859) q[0];
sx q[0];
rz(-1.4771853) q[0];
x q[1];
rz(-2.426067) q[2];
sx q[2];
rz(-1.160485) q[2];
sx q[2];
rz(-0.35718756) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4775042) q[1];
sx q[1];
rz(-1.3574151) q[1];
sx q[1];
rz(2.1413436) q[1];
rz(-0.81256688) q[3];
sx q[3];
rz(-1.1538394) q[3];
sx q[3];
rz(-1.5244694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2734566) q[2];
sx q[2];
rz(-3.0507705) q[2];
sx q[2];
rz(0.4503797) q[2];
rz(0.069933794) q[3];
sx q[3];
rz(-2.0262597) q[3];
sx q[3];
rz(2.6993338) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.062334199) q[0];
sx q[0];
rz(-2.313995) q[0];
sx q[0];
rz(-0.29912478) q[0];
rz(-2.4807855) q[1];
sx q[1];
rz(-2.0313163) q[1];
sx q[1];
rz(0.92346907) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8689821) q[0];
sx q[0];
rz(-1.8911757) q[0];
sx q[0];
rz(2.9584209) q[0];
rz(-0.23552509) q[2];
sx q[2];
rz(-1.7851637) q[2];
sx q[2];
rz(1.216824) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0638173) q[1];
sx q[1];
rz(-0.37250054) q[1];
sx q[1];
rz(0.17613303) q[1];
x q[2];
rz(-0.98387444) q[3];
sx q[3];
rz(-1.8070584) q[3];
sx q[3];
rz(2.7191702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.48351651) q[2];
sx q[2];
rz(-2.7407756) q[2];
sx q[2];
rz(3.0576903) q[2];
rz(0.99533254) q[3];
sx q[3];
rz(-1.1570802) q[3];
sx q[3];
rz(-1.5108494) q[3];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8453269) q[0];
sx q[0];
rz(-1.9060598) q[0];
sx q[0];
rz(-0.96257019) q[0];
rz(-0.30749097) q[1];
sx q[1];
rz(-0.88491416) q[1];
sx q[1];
rz(0.21893758) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5318659) q[0];
sx q[0];
rz(-1.5012739) q[0];
sx q[0];
rz(-0.34644925) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.89246558) q[2];
sx q[2];
rz(-1.7545059) q[2];
sx q[2];
rz(3.122729) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1477222) q[1];
sx q[1];
rz(-1.9703373) q[1];
sx q[1];
rz(3.1385204) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2943324) q[3];
sx q[3];
rz(-0.81786441) q[3];
sx q[3];
rz(0.62608214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4585877) q[2];
sx q[2];
rz(-1.3938437) q[2];
sx q[2];
rz(-1.5031507) q[2];
rz(-0.25034869) q[3];
sx q[3];
rz(-0.72385794) q[3];
sx q[3];
rz(2.943128) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4918764) q[0];
sx q[0];
rz(-2.8174077) q[0];
sx q[0];
rz(-2.8247996) q[0];
rz(0.016238796) q[1];
sx q[1];
rz(-2.3490678) q[1];
sx q[1];
rz(-0.24831717) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1386928) q[0];
sx q[0];
rz(-1.0637245) q[0];
sx q[0];
rz(0.72082989) q[0];
rz(-pi) q[1];
rz(0.14713471) q[2];
sx q[2];
rz(-2.6700041) q[2];
sx q[2];
rz(1.6560645) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.59245306) q[1];
sx q[1];
rz(-0.35129476) q[1];
sx q[1];
rz(0.37254803) q[1];
x q[2];
rz(2.2894182) q[3];
sx q[3];
rz(-0.2874473) q[3];
sx q[3];
rz(1.3369657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.877964) q[2];
sx q[2];
rz(-0.10136494) q[2];
sx q[2];
rz(1.9975837) q[2];
rz(2.3794203) q[3];
sx q[3];
rz(-0.52505571) q[3];
sx q[3];
rz(-2.150056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85631973) q[0];
sx q[0];
rz(-3.0506436) q[0];
sx q[0];
rz(-1.0524622) q[0];
rz(-2.6513124) q[1];
sx q[1];
rz(-1.0888381) q[1];
sx q[1];
rz(3.0163684) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1220684) q[0];
sx q[0];
rz(-0.53036371) q[0];
sx q[0];
rz(0.85774647) q[0];
rz(-pi) q[1];
rz(0.43951054) q[2];
sx q[2];
rz(-0.59328967) q[2];
sx q[2];
rz(0.35542929) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.050712) q[1];
sx q[1];
rz(-1.4987603) q[1];
sx q[1];
rz(1.8211831) q[1];
x q[2];
rz(2.9216705) q[3];
sx q[3];
rz(-0.87080216) q[3];
sx q[3];
rz(0.31693566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.48309717) q[2];
sx q[2];
rz(-0.53762895) q[2];
sx q[2];
rz(-1.5527234) q[2];
rz(-0.056564864) q[3];
sx q[3];
rz(-1.5141234) q[3];
sx q[3];
rz(0.27730832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5796984) q[0];
sx q[0];
rz(-0.46439207) q[0];
sx q[0];
rz(0.42771801) q[0];
rz(-0.5411717) q[1];
sx q[1];
rz(-2.7310557) q[1];
sx q[1];
rz(-1.7778832) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83374524) q[0];
sx q[0];
rz(-2.5302093) q[0];
sx q[0];
rz(-1.3788542) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.12597106) q[2];
sx q[2];
rz(-0.4301542) q[2];
sx q[2];
rz(-0.90829721) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7040055) q[1];
sx q[1];
rz(-2.2477337) q[1];
sx q[1];
rz(2.747405) q[1];
x q[2];
rz(1.1579944) q[3];
sx q[3];
rz(-1.5607474) q[3];
sx q[3];
rz(-1.5749941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0987739) q[2];
sx q[2];
rz(-0.51977485) q[2];
sx q[2];
rz(2.7322312) q[2];
rz(2.9218819) q[3];
sx q[3];
rz(-1.5922056) q[3];
sx q[3];
rz(-1.6519914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4396502) q[0];
sx q[0];
rz(-0.42124978) q[0];
sx q[0];
rz(-0.14081328) q[0];
rz(-0.5367865) q[1];
sx q[1];
rz(-1.7861007) q[1];
sx q[1];
rz(2.7649194) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46464065) q[0];
sx q[0];
rz(-2.8429705) q[0];
sx q[0];
rz(-0.99891164) q[0];
x q[1];
rz(-2.4666653) q[2];
sx q[2];
rz(-1.2528193) q[2];
sx q[2];
rz(-0.392158) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.46493775) q[1];
sx q[1];
rz(-2.3466517) q[1];
sx q[1];
rz(1.388768) q[1];
rz(-2.7984058) q[3];
sx q[3];
rz(-0.43794855) q[3];
sx q[3];
rz(-1.7476976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.10304) q[2];
sx q[2];
rz(-1.1662177) q[2];
sx q[2];
rz(2.340359) q[2];
rz(1.0771105) q[3];
sx q[3];
rz(-0.25682768) q[3];
sx q[3];
rz(-0.1424772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(0.24256663) q[0];
sx q[0];
rz(-0.16066708) q[0];
sx q[0];
rz(2.2005079) q[0];
rz(-0.89648992) q[1];
sx q[1];
rz(-1.4184378) q[1];
sx q[1];
rz(-2.3268708) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.488637) q[0];
sx q[0];
rz(-1.9359447) q[0];
sx q[0];
rz(0.40819755) q[0];
rz(-pi) q[1];
x q[1];
rz(0.28622921) q[2];
sx q[2];
rz(-2.9636418) q[2];
sx q[2];
rz(-0.19831698) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.43489698) q[1];
sx q[1];
rz(-2.416947) q[1];
sx q[1];
rz(-1.0452634) q[1];
x q[2];
rz(-0.37121935) q[3];
sx q[3];
rz(-1.0346078) q[3];
sx q[3];
rz(0.49193383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.61646378) q[2];
sx q[2];
rz(-0.63462555) q[2];
sx q[2];
rz(-2.5647707) q[2];
rz(0.47354627) q[3];
sx q[3];
rz(-1.1487995) q[3];
sx q[3];
rz(0.043341652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-0.43805495) q[0];
sx q[0];
rz(-2.3860768) q[0];
sx q[0];
rz(-3.0536861) q[0];
rz(-1.5880623) q[1];
sx q[1];
rz(-0.5843662) q[1];
sx q[1];
rz(2.8839674) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15730298) q[0];
sx q[0];
rz(-2.7098795) q[0];
sx q[0];
rz(1.2832416) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.070775971) q[2];
sx q[2];
rz(-2.1969866) q[2];
sx q[2];
rz(-2.6232121) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0637197) q[1];
sx q[1];
rz(-1.8053836) q[1];
sx q[1];
rz(-0.44198702) q[1];
rz(-2.1472503) q[3];
sx q[3];
rz(-1.0691158) q[3];
sx q[3];
rz(0.27015206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0660925) q[2];
sx q[2];
rz(-1.6801715) q[2];
sx q[2];
rz(3.063805) q[2];
rz(1.0545688) q[3];
sx q[3];
rz(-2.4013897) q[3];
sx q[3];
rz(-2.4906702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2442378) q[0];
sx q[0];
rz(-0.033441823) q[0];
sx q[0];
rz(0.43575409) q[0];
rz(3.1089973) q[1];
sx q[1];
rz(-2.4041912) q[1];
sx q[1];
rz(-1.3479412) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.01719) q[0];
sx q[0];
rz(-1.5637056) q[0];
sx q[0];
rz(-0.75049863) q[0];
rz(-0.52705898) q[2];
sx q[2];
rz(-1.3599476) q[2];
sx q[2];
rz(-1.859889) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0465728) q[1];
sx q[1];
rz(-3.1029178) q[1];
sx q[1];
rz(-1.2016618) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3147821) q[3];
sx q[3];
rz(-1.7932995) q[3];
sx q[3];
rz(3.0526628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2160448) q[2];
sx q[2];
rz(-2.1074769) q[2];
sx q[2];
rz(-0.75630581) q[2];
rz(-0.082152724) q[3];
sx q[3];
rz(-0.3923471) q[3];
sx q[3];
rz(0.69123417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3638851) q[0];
sx q[0];
rz(-2.4019699) q[0];
sx q[0];
rz(-0.85516047) q[0];
rz(2.0289519) q[1];
sx q[1];
rz(-2.009195) q[1];
sx q[1];
rz(2.7809273) q[1];
rz(-1.821221) q[2];
sx q[2];
rz(-1.741521) q[2];
sx q[2];
rz(2.3157816) q[2];
rz(-2.6769911) q[3];
sx q[3];
rz(-2.0809485) q[3];
sx q[3];
rz(2.8551024) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
