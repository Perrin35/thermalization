OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5698009) q[0];
sx q[0];
rz(-2.2324012) q[0];
sx q[0];
rz(-1.5067014) q[0];
rz(1.2070967) q[1];
sx q[1];
rz(6.7758898) q[1];
sx q[1];
rz(13.125782) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.089559) q[0];
sx q[0];
rz(-1.8843643) q[0];
sx q[0];
rz(1.3960209) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.71601358) q[2];
sx q[2];
rz(-0.56173827) q[2];
sx q[2];
rz(1.5519976) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.31025793) q[1];
sx q[1];
rz(-1.2911011) q[1];
sx q[1];
rz(3.1344942) q[1];
x q[2];
rz(0.97752882) q[3];
sx q[3];
rz(-1.1479605) q[3];
sx q[3];
rz(-2.5893167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.59014615) q[2];
sx q[2];
rz(-0.85942736) q[2];
sx q[2];
rz(-2.0476511) q[2];
rz(-1.236773) q[3];
sx q[3];
rz(-1.1636795) q[3];
sx q[3];
rz(2.8804603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7673489) q[0];
sx q[0];
rz(-1.2232895) q[0];
sx q[0];
rz(-2.4236524) q[0];
rz(1.9473437) q[1];
sx q[1];
rz(-2.7337044) q[1];
sx q[1];
rz(1.6494707) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8555174) q[0];
sx q[0];
rz(-1.6950894) q[0];
sx q[0];
rz(3.0086424) q[0];
rz(-1.1928333) q[2];
sx q[2];
rz(-1.1587669) q[2];
sx q[2];
rz(-2.2343324) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1528984) q[1];
sx q[1];
rz(-1.1225268) q[1];
sx q[1];
rz(-0.50181915) q[1];
rz(-pi) q[2];
rz(0.82573311) q[3];
sx q[3];
rz(-1.4108218) q[3];
sx q[3];
rz(2.3478594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0933928) q[2];
sx q[2];
rz(-0.77646774) q[2];
sx q[2];
rz(0.36273599) q[2];
rz(-2.2705966) q[3];
sx q[3];
rz(-1.3716776) q[3];
sx q[3];
rz(-1.8179651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
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
rz(0.094630346) q[0];
sx q[0];
rz(-1.3032664) q[0];
sx q[0];
rz(2.5808425) q[0];
rz(2.1669855) q[1];
sx q[1];
rz(-2.2798996) q[1];
sx q[1];
rz(2.9240756) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91470766) q[0];
sx q[0];
rz(-0.98641005) q[0];
sx q[0];
rz(1.547772) q[0];
rz(1.752408) q[2];
sx q[2];
rz(-2.6208502) q[2];
sx q[2];
rz(0.40042675) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8629093) q[1];
sx q[1];
rz(-2.8562198) q[1];
sx q[1];
rz(-1.1527658) q[1];
rz(1.6337745) q[3];
sx q[3];
rz(-2.1788886) q[3];
sx q[3];
rz(0.62880317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4983623) q[2];
sx q[2];
rz(-0.78780323) q[2];
sx q[2];
rz(0.95334774) q[2];
rz(-2.0638454) q[3];
sx q[3];
rz(-1.6809623) q[3];
sx q[3];
rz(-0.47422153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.137602) q[0];
sx q[0];
rz(-2.9676262) q[0];
sx q[0];
rz(2.2250788) q[0];
rz(0.42463955) q[1];
sx q[1];
rz(-1.3820796) q[1];
sx q[1];
rz(2.0133846) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5461499) q[0];
sx q[0];
rz(-2.2227123) q[0];
sx q[0];
rz(-2.2007887) q[0];
rz(-pi) q[1];
rz(2.9229972) q[2];
sx q[2];
rz(-2.3724883) q[2];
sx q[2];
rz(-1.1111271) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.90626806) q[1];
sx q[1];
rz(-0.61329326) q[1];
sx q[1];
rz(-2.6032531) q[1];
x q[2];
rz(2.1682285) q[3];
sx q[3];
rz(-1.7586244) q[3];
sx q[3];
rz(1.1115766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.26587036) q[2];
sx q[2];
rz(-1.6758726) q[2];
sx q[2];
rz(-1.9423368) q[2];
rz(-1.9843598) q[3];
sx q[3];
rz(-0.80409378) q[3];
sx q[3];
rz(2.6257302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-2.1542094) q[0];
sx q[0];
rz(-0.043954285) q[0];
sx q[0];
rz(-2.2162345) q[0];
rz(2.5788653) q[1];
sx q[1];
rz(-1.0147164) q[1];
sx q[1];
rz(-2.7664807) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5432388) q[0];
sx q[0];
rz(-1.028201) q[0];
sx q[0];
rz(1.6981237) q[0];
x q[1];
rz(-3.0882224) q[2];
sx q[2];
rz(-2.0056021) q[2];
sx q[2];
rz(-0.31440266) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7081333) q[1];
sx q[1];
rz(-2.9707751) q[1];
sx q[1];
rz(2.9950735) q[1];
x q[2];
rz(0.45467037) q[3];
sx q[3];
rz(-1.934762) q[3];
sx q[3];
rz(2.2157247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.90317) q[2];
sx q[2];
rz(-2.2245202) q[2];
sx q[2];
rz(2.6440716) q[2];
rz(-0.43429747) q[3];
sx q[3];
rz(-1.4562166) q[3];
sx q[3];
rz(0.98880497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76181805) q[0];
sx q[0];
rz(-0.86123818) q[0];
sx q[0];
rz(-0.27989835) q[0];
rz(-2.1474536) q[1];
sx q[1];
rz(-0.86052624) q[1];
sx q[1];
rz(2.5631189) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3247447) q[0];
sx q[0];
rz(-1.9473796) q[0];
sx q[0];
rz(-2.5927932) q[0];
rz(-pi) q[1];
x q[1];
rz(0.70273383) q[2];
sx q[2];
rz(-1.2652072) q[2];
sx q[2];
rz(-1.0440799) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.76982564) q[1];
sx q[1];
rz(-2.3845362) q[1];
sx q[1];
rz(0.28626059) q[1];
x q[2];
rz(1.469645) q[3];
sx q[3];
rz(-1.3037852) q[3];
sx q[3];
rz(1.284541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0824739) q[2];
sx q[2];
rz(-0.59933496) q[2];
sx q[2];
rz(0.88258755) q[2];
rz(0.62697083) q[3];
sx q[3];
rz(-0.65492237) q[3];
sx q[3];
rz(-2.2466808) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51666981) q[0];
sx q[0];
rz(-0.94180095) q[0];
sx q[0];
rz(-3.140977) q[0];
rz(0.90283886) q[1];
sx q[1];
rz(-0.87264624) q[1];
sx q[1];
rz(3.0965064) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0233399) q[0];
sx q[0];
rz(-1.0824507) q[0];
sx q[0];
rz(-1.4834542) q[0];
rz(1.3602518) q[2];
sx q[2];
rz(-2.3017052) q[2];
sx q[2];
rz(2.108253) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.0087456044) q[1];
sx q[1];
rz(-2.3373621) q[1];
sx q[1];
rz(-2.7102317) q[1];
x q[2];
rz(-0.99211971) q[3];
sx q[3];
rz(-2.1014155) q[3];
sx q[3];
rz(-1.19095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8467329) q[2];
sx q[2];
rz(-2.5856057) q[2];
sx q[2];
rz(-2.2247458) q[2];
rz(0.89546853) q[3];
sx q[3];
rz(-1.6296891) q[3];
sx q[3];
rz(1.2112613) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.004892) q[0];
sx q[0];
rz(-3.0558375) q[0];
sx q[0];
rz(-0.46491796) q[0];
rz(-1.1888986) q[1];
sx q[1];
rz(-1.3465954) q[1];
sx q[1];
rz(2.7704923) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2087323) q[0];
sx q[0];
rz(-1.1602872) q[0];
sx q[0];
rz(-0.44637827) q[0];
rz(3.1041652) q[2];
sx q[2];
rz(-2.3212395) q[2];
sx q[2];
rz(-1.6714753) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5587204) q[1];
sx q[1];
rz(-2.0077188) q[1];
sx q[1];
rz(-0.82804273) q[1];
rz(-pi) q[2];
x q[2];
rz(0.38406541) q[3];
sx q[3];
rz(-1.3017941) q[3];
sx q[3];
rz(-2.254738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4798639) q[2];
sx q[2];
rz(-0.78833818) q[2];
sx q[2];
rz(-2.7313477) q[2];
rz(1.5069626) q[3];
sx q[3];
rz(-0.40532902) q[3];
sx q[3];
rz(0.043225616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0931382) q[0];
sx q[0];
rz(-2.6148836) q[0];
sx q[0];
rz(2.9658537) q[0];
rz(-2.3161092) q[1];
sx q[1];
rz(-1.8800507) q[1];
sx q[1];
rz(0.82297355) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9776445) q[0];
sx q[0];
rz(-2.2744482) q[0];
sx q[0];
rz(0.73935853) q[0];
x q[1];
rz(-2.1954567) q[2];
sx q[2];
rz(-0.51057928) q[2];
sx q[2];
rz(0.389314) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8205679) q[1];
sx q[1];
rz(-1.9181644) q[1];
sx q[1];
rz(-0.72183164) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4863344) q[3];
sx q[3];
rz(-0.41197398) q[3];
sx q[3];
rz(-2.1503445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2814111) q[2];
sx q[2];
rz(-1.9965636) q[2];
sx q[2];
rz(-2.6004876) q[2];
rz(-2.7152854) q[3];
sx q[3];
rz(-0.46190327) q[3];
sx q[3];
rz(-2.2947252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68778872) q[0];
sx q[0];
rz(-2.332088) q[0];
sx q[0];
rz(2.8059106) q[0];
rz(-3.1188534) q[1];
sx q[1];
rz(-0.72240654) q[1];
sx q[1];
rz(0.67363277) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27972066) q[0];
sx q[0];
rz(-2.4422944) q[0];
sx q[0];
rz(-2.5847816) q[0];
rz(2.6725476) q[2];
sx q[2];
rz(-1.6613591) q[2];
sx q[2];
rz(-0.064571206) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2362263) q[1];
sx q[1];
rz(-2.1520808) q[1];
sx q[1];
rz(-2.8984215) q[1];
rz(-pi) q[2];
rz(-1.4309747) q[3];
sx q[3];
rz(-1.2396024) q[3];
sx q[3];
rz(0.14023031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.444904) q[2];
sx q[2];
rz(-0.99385571) q[2];
sx q[2];
rz(2.135684) q[2];
rz(2.3575947) q[3];
sx q[3];
rz(-2.8204212) q[3];
sx q[3];
rz(-1.706749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1226596) q[0];
sx q[0];
rz(-1.4228595) q[0];
sx q[0];
rz(2.0056437) q[0];
rz(-2.7574273) q[1];
sx q[1];
rz(-1.1687678) q[1];
sx q[1];
rz(1.6484177) q[1];
rz(-2.3758985) q[2];
sx q[2];
rz(-2.8509344) q[2];
sx q[2];
rz(0.40727587) q[2];
rz(0.85687153) q[3];
sx q[3];
rz(-1.8663434) q[3];
sx q[3];
rz(0.28874884) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
