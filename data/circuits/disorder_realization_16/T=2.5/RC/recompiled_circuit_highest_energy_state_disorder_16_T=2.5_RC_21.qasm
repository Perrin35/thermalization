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
rz(0.056869153) q[0];
sx q[0];
rz(-0.19357227) q[0];
sx q[0];
rz(-0.40785664) q[0];
rz(-0.017539311) q[1];
sx q[1];
rz(-1.8899625) q[1];
sx q[1];
rz(-1.6012021) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6787036) q[0];
sx q[0];
rz(-0.4400095) q[0];
sx q[0];
rz(1.6848887) q[0];
rz(-pi) q[1];
x q[1];
rz(0.27215927) q[2];
sx q[2];
rz(-2.6725884) q[2];
sx q[2];
rz(2.9763165) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7471024) q[1];
sx q[1];
rz(-1.1848406) q[1];
sx q[1];
rz(-1.3420022) q[1];
rz(-pi) q[2];
rz(-0.82117053) q[3];
sx q[3];
rz(-1.0945012) q[3];
sx q[3];
rz(-2.0408292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5876329) q[2];
sx q[2];
rz(-2.0822058) q[2];
sx q[2];
rz(-2.9239192) q[2];
rz(2.830128) q[3];
sx q[3];
rz(-2.5296827) q[3];
sx q[3];
rz(1.9757087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4199453) q[0];
sx q[0];
rz(-2.8657275) q[0];
sx q[0];
rz(-0.69197792) q[0];
rz(2.3852589) q[1];
sx q[1];
rz(-0.5223918) q[1];
sx q[1];
rz(-1.5501032) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78849906) q[0];
sx q[0];
rz(-2.278557) q[0];
sx q[0];
rz(1.4097296) q[0];
x q[1];
rz(1.814957) q[2];
sx q[2];
rz(-0.40882822) q[2];
sx q[2];
rz(-0.93016184) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2427502) q[1];
sx q[1];
rz(-1.4489343) q[1];
sx q[1];
rz(0.16074462) q[1];
rz(-pi) q[2];
rz(-0.85568537) q[3];
sx q[3];
rz(-2.4484903) q[3];
sx q[3];
rz(0.86831304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.99593607) q[2];
sx q[2];
rz(-2.9684976) q[2];
sx q[2];
rz(-2.9023329) q[2];
rz(-1.4907106) q[3];
sx q[3];
rz(-1.2942634) q[3];
sx q[3];
rz(-1.9132805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2759129) q[0];
sx q[0];
rz(-0.90621197) q[0];
sx q[0];
rz(-0.014634125) q[0];
rz(2.2485661) q[1];
sx q[1];
rz(-2.1664186) q[1];
sx q[1];
rz(-2.9258974) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53447002) q[0];
sx q[0];
rz(-1.5661582) q[0];
sx q[0];
rz(1.5195373) q[0];
rz(-pi) q[1];
rz(-0.65962445) q[2];
sx q[2];
rz(-1.2953892) q[2];
sx q[2];
rz(3.1153452) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0552935) q[1];
sx q[1];
rz(-0.26565105) q[1];
sx q[1];
rz(2.6363027) q[1];
x q[2];
rz(-0.87743296) q[3];
sx q[3];
rz(-2.139353) q[3];
sx q[3];
rz(1.7753778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0594242) q[2];
sx q[2];
rz(-1.1801722) q[2];
sx q[2];
rz(0.99349418) q[2];
rz(0.70837402) q[3];
sx q[3];
rz(-1.1185442) q[3];
sx q[3];
rz(0.12349252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74743903) q[0];
sx q[0];
rz(-0.49462947) q[0];
sx q[0];
rz(2.4476449) q[0];
rz(-0.28678647) q[1];
sx q[1];
rz(-1.6096121) q[1];
sx q[1];
rz(2.0538816) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74831731) q[0];
sx q[0];
rz(-1.4469742) q[0];
sx q[0];
rz(2.078915) q[0];
x q[1];
rz(1.0469646) q[2];
sx q[2];
rz(-1.0817384) q[2];
sx q[2];
rz(-0.2917052) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.775271) q[1];
sx q[1];
rz(-1.873981) q[1];
sx q[1];
rz(-0.94667411) q[1];
rz(-pi) q[2];
rz(-1.5252211) q[3];
sx q[3];
rz(-2.2846095) q[3];
sx q[3];
rz(2.8136106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9168758) q[2];
sx q[2];
rz(-1.1563053) q[2];
sx q[2];
rz(1.7737596) q[2];
rz(-1.2380838) q[3];
sx q[3];
rz(-1.5690469) q[3];
sx q[3];
rz(-2.9253173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5991768) q[0];
sx q[0];
rz(-1.2681862) q[0];
sx q[0];
rz(-0.61781484) q[0];
rz(1.1082331) q[1];
sx q[1];
rz(-0.30312678) q[1];
sx q[1];
rz(2.2584426) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9882433) q[0];
sx q[0];
rz(-0.81479615) q[0];
sx q[0];
rz(-1.5606461) q[0];
rz(-pi) q[1];
rz(0.3381392) q[2];
sx q[2];
rz(-0.45881841) q[2];
sx q[2];
rz(-0.79566075) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0269491) q[1];
sx q[1];
rz(-1.0806688) q[1];
sx q[1];
rz(3.0672706) q[1];
x q[2];
rz(0.95728504) q[3];
sx q[3];
rz(-0.64910474) q[3];
sx q[3];
rz(1.4278936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9250179) q[2];
sx q[2];
rz(-2.8897132) q[2];
sx q[2];
rz(1.3612755) q[2];
rz(1.2733634) q[3];
sx q[3];
rz(-1.7073771) q[3];
sx q[3];
rz(-0.98289615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.044416044) q[0];
sx q[0];
rz(-1.8955078) q[0];
sx q[0];
rz(0.62028766) q[0];
rz(2.7445131) q[1];
sx q[1];
rz(-2.670791) q[1];
sx q[1];
rz(1.1995859) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.967554) q[0];
sx q[0];
rz(-0.10835526) q[0];
sx q[0];
rz(-1.3757785) q[0];
rz(-2.994147) q[2];
sx q[2];
rz(-1.7803528) q[2];
sx q[2];
rz(3.0990019) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0850272) q[1];
sx q[1];
rz(-1.051184) q[1];
sx q[1];
rz(1.7934133) q[1];
rz(2.7971917) q[3];
sx q[3];
rz(-1.6823497) q[3];
sx q[3];
rz(-1.9007746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.24641307) q[2];
sx q[2];
rz(-1.2440888) q[2];
sx q[2];
rz(-0.67071521) q[2];
rz(-2.8504168) q[3];
sx q[3];
rz(-2.5305735) q[3];
sx q[3];
rz(-2.5197855) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7156242) q[0];
sx q[0];
rz(-1.427303) q[0];
sx q[0];
rz(-0.17844644) q[0];
rz(2.2715691) q[1];
sx q[1];
rz(-0.88302892) q[1];
sx q[1];
rz(-0.80284405) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69235301) q[0];
sx q[0];
rz(-2.521551) q[0];
sx q[0];
rz(0.043278261) q[0];
rz(-pi) q[1];
rz(2.5576791) q[2];
sx q[2];
rz(-1.7758992) q[2];
sx q[2];
rz(-2.5896304) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9446148) q[1];
sx q[1];
rz(-2.3715545) q[1];
sx q[1];
rz(-0.47702392) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25814806) q[3];
sx q[3];
rz(-2.1827563) q[3];
sx q[3];
rz(-1.985509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.37658438) q[2];
sx q[2];
rz(-1.8439801) q[2];
sx q[2];
rz(0.89361781) q[2];
rz(-1.5997959) q[3];
sx q[3];
rz(-1.4897646) q[3];
sx q[3];
rz(2.5347575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1641418) q[0];
sx q[0];
rz(-2.7315388) q[0];
sx q[0];
rz(2.9264911) q[0];
rz(0.84456259) q[1];
sx q[1];
rz(-1.3203878) q[1];
sx q[1];
rz(0.79427687) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1298908) q[0];
sx q[0];
rz(-3.112535) q[0];
sx q[0];
rz(-0.41883166) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6849943) q[2];
sx q[2];
rz(-1.2933044) q[2];
sx q[2];
rz(3.1268529) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8819184) q[1];
sx q[1];
rz(-1.4800023) q[1];
sx q[1];
rz(-2.1114248) q[1];
rz(2.8135962) q[3];
sx q[3];
rz(-0.80910002) q[3];
sx q[3];
rz(-2.3878333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8636785) q[2];
sx q[2];
rz(-1.2028376) q[2];
sx q[2];
rz(-1.4097144) q[2];
rz(-2.388741) q[3];
sx q[3];
rz(-1.9949621) q[3];
sx q[3];
rz(2.1394155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3442605) q[0];
sx q[0];
rz(-2.7334038) q[0];
sx q[0];
rz(-2.1687188) q[0];
rz(-1.5219888) q[1];
sx q[1];
rz(-1.3643967) q[1];
sx q[1];
rz(2.5111228) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3469543) q[0];
sx q[0];
rz(-1.3890059) q[0];
sx q[0];
rz(-1.7111045) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5233598) q[2];
sx q[2];
rz(-0.15711297) q[2];
sx q[2];
rz(2.4542798) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3296632) q[1];
sx q[1];
rz(-2.3721937) q[1];
sx q[1];
rz(-3.0698464) q[1];
x q[2];
rz(1.6338324) q[3];
sx q[3];
rz(-1.3975289) q[3];
sx q[3];
rz(1.2931371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1260599) q[2];
sx q[2];
rz(-1.1522747) q[2];
sx q[2];
rz(-1.3014334) q[2];
rz(1.556501) q[3];
sx q[3];
rz(-0.82587487) q[3];
sx q[3];
rz(-0.41527709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5179317) q[0];
sx q[0];
rz(-0.15075891) q[0];
sx q[0];
rz(0.49322042) q[0];
rz(-1.8863691) q[1];
sx q[1];
rz(-1.8938277) q[1];
sx q[1];
rz(-0.36466041) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4987558) q[0];
sx q[0];
rz(-2.7793573) q[0];
sx q[0];
rz(1.9693768) q[0];
rz(-2.1388059) q[2];
sx q[2];
rz(-2.4052561) q[2];
sx q[2];
rz(1.9970837) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.030487) q[1];
sx q[1];
rz(-2.3685622) q[1];
sx q[1];
rz(1.4902601) q[1];
x q[2];
rz(-1.8871898) q[3];
sx q[3];
rz(-1.5324161) q[3];
sx q[3];
rz(-0.61172133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6898592) q[2];
sx q[2];
rz(-2.1351337) q[2];
sx q[2];
rz(-0.9922007) q[2];
rz(0.36710468) q[3];
sx q[3];
rz(-1.4398984) q[3];
sx q[3];
rz(2.4632857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51919666) q[0];
sx q[0];
rz(-1.665103) q[0];
sx q[0];
rz(-1.7892224) q[0];
rz(2.2182111) q[1];
sx q[1];
rz(-2.2877749) q[1];
sx q[1];
rz(1.9705082) q[1];
rz(2.1869356) q[2];
sx q[2];
rz(-2.8980394) q[2];
sx q[2];
rz(0.63799636) q[2];
rz(0.97040776) q[3];
sx q[3];
rz(-2.1526489) q[3];
sx q[3];
rz(-2.7027705) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
