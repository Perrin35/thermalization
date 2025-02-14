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
rz(0.87585706) q[0];
sx q[0];
rz(-0.96867222) q[0];
sx q[0];
rz(-0.92585603) q[0];
rz(-2.5246188) q[1];
sx q[1];
rz(-2.4763835) q[1];
sx q[1];
rz(1.8170504) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2331794) q[0];
sx q[0];
rz(-2.1696013) q[0];
sx q[0];
rz(1.6480867) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1466897) q[2];
sx q[2];
rz(-2.5518199) q[2];
sx q[2];
rz(-1.1663933) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.032925978) q[1];
sx q[1];
rz(-2.1582728) q[1];
sx q[1];
rz(-1.0125748) q[1];
x q[2];
rz(-2.376287) q[3];
sx q[3];
rz(-2.0767143) q[3];
sx q[3];
rz(-1.7047061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2609743) q[2];
sx q[2];
rz(-1.3857434) q[2];
sx q[2];
rz(-0.79208881) q[2];
rz(-2.2657307) q[3];
sx q[3];
rz(-0.10871092) q[3];
sx q[3];
rz(-0.66332269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51139128) q[0];
sx q[0];
rz(-1.317861) q[0];
sx q[0];
rz(-2.200101) q[0];
rz(2.1425653) q[1];
sx q[1];
rz(-0.92728725) q[1];
sx q[1];
rz(-3.0587382) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6951286) q[0];
sx q[0];
rz(-2.3020491) q[0];
sx q[0];
rz(-1.7056998) q[0];
rz(-pi) q[1];
rz(0.98495667) q[2];
sx q[2];
rz(-2.7212866) q[2];
sx q[2];
rz(2.676385) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1260316) q[1];
sx q[1];
rz(-2.0210316) q[1];
sx q[1];
rz(-2.5712988) q[1];
rz(2.4965246) q[3];
sx q[3];
rz(-2.0760025) q[3];
sx q[3];
rz(-0.92407214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2443709) q[2];
sx q[2];
rz(-1.3542078) q[2];
sx q[2];
rz(-1.2574035) q[2];
rz(-2.2952378) q[3];
sx q[3];
rz(-2.9823163) q[3];
sx q[3];
rz(0.67811051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1988679) q[0];
sx q[0];
rz(-1.6833479) q[0];
sx q[0];
rz(2.9123836) q[0];
rz(2.9275059) q[1];
sx q[1];
rz(-0.87132088) q[1];
sx q[1];
rz(0.3814989) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70951916) q[0];
sx q[0];
rz(-0.95147485) q[0];
sx q[0];
rz(-1.1362057) q[0];
rz(-0.83865954) q[2];
sx q[2];
rz(-1.7146972) q[2];
sx q[2];
rz(0.14140192) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4134564) q[1];
sx q[1];
rz(-1.6237139) q[1];
sx q[1];
rz(3.1127597) q[1];
x q[2];
rz(0.014147357) q[3];
sx q[3];
rz(-2.7944744) q[3];
sx q[3];
rz(-1.0747758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4182338) q[2];
sx q[2];
rz(-0.25622076) q[2];
sx q[2];
rz(-0.56824938) q[2];
rz(1.4189643) q[3];
sx q[3];
rz(-1.4293554) q[3];
sx q[3];
rz(1.0770146) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.028246183) q[0];
sx q[0];
rz(-1.1320817) q[0];
sx q[0];
rz(2.8908253) q[0];
rz(-0.084550683) q[1];
sx q[1];
rz(-1.0944347) q[1];
sx q[1];
rz(1.9409404) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4290743) q[0];
sx q[0];
rz(-1.9299283) q[0];
sx q[0];
rz(-2.284689) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8941706) q[2];
sx q[2];
rz(-2.0004002) q[2];
sx q[2];
rz(1.0216433) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8507217) q[1];
sx q[1];
rz(-2.4230036) q[1];
sx q[1];
rz(2.9903837) q[1];
rz(-pi) q[2];
rz(2.0968998) q[3];
sx q[3];
rz(-1.3402437) q[3];
sx q[3];
rz(-1.5073204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4667929) q[2];
sx q[2];
rz(-1.1226706) q[2];
sx q[2];
rz(1.2910845) q[2];
rz(2.71991) q[3];
sx q[3];
rz(-0.45160523) q[3];
sx q[3];
rz(1.5765367) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52294937) q[0];
sx q[0];
rz(-2.4186501) q[0];
sx q[0];
rz(-1.500754) q[0];
rz(-0.067642637) q[1];
sx q[1];
rz(-1.5539955) q[1];
sx q[1];
rz(1.1281475) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2826506) q[0];
sx q[0];
rz(-2.9500486) q[0];
sx q[0];
rz(2.9454977) q[0];
x q[1];
rz(-2.0968151) q[2];
sx q[2];
rz(-1.8785718) q[2];
sx q[2];
rz(0.36164944) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.060202816) q[1];
sx q[1];
rz(-3.0160025) q[1];
sx q[1];
rz(-2.0141898) q[1];
rz(0.23388548) q[3];
sx q[3];
rz(-1.2986168) q[3];
sx q[3];
rz(-0.50258499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.03380123) q[2];
sx q[2];
rz(-1.1148323) q[2];
sx q[2];
rz(-0.6768674) q[2];
rz(2.233861) q[3];
sx q[3];
rz(-1.9151442) q[3];
sx q[3];
rz(2.7270253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21861976) q[0];
sx q[0];
rz(-0.75160471) q[0];
sx q[0];
rz(0.76939097) q[0];
rz(1.2409302) q[1];
sx q[1];
rz(-0.85378328) q[1];
sx q[1];
rz(-2.9262537) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8148635) q[0];
sx q[0];
rz(-1.0007326) q[0];
sx q[0];
rz(0.22320052) q[0];
x q[1];
rz(-3.1245232) q[2];
sx q[2];
rz(-1.9015549) q[2];
sx q[2];
rz(-1.1137373) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3591738) q[1];
sx q[1];
rz(-0.69643387) q[1];
sx q[1];
rz(-0.38621186) q[1];
rz(-pi) q[2];
rz(-2.4438079) q[3];
sx q[3];
rz(-1.8337733) q[3];
sx q[3];
rz(-1.7587723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.26663366) q[2];
sx q[2];
rz(-1.5077488) q[2];
sx q[2];
rz(-0.61666644) q[2];
rz(-0.2002317) q[3];
sx q[3];
rz(-0.73035208) q[3];
sx q[3];
rz(2.0088137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.013833372) q[0];
sx q[0];
rz(-0.56202373) q[0];
sx q[0];
rz(-3.1264937) q[0];
rz(0.12241441) q[1];
sx q[1];
rz(-1.2950803) q[1];
sx q[1];
rz(2.1133568) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0434179) q[0];
sx q[0];
rz(-2.2255996) q[0];
sx q[0];
rz(1.5006939) q[0];
x q[1];
rz(0.21034849) q[2];
sx q[2];
rz(-2.1783683) q[2];
sx q[2];
rz(1.4841532) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9837357) q[1];
sx q[1];
rz(-2.7089556) q[1];
sx q[1];
rz(3.0745201) q[1];
x q[2];
rz(2.8260303) q[3];
sx q[3];
rz(-1.8286108) q[3];
sx q[3];
rz(2.1699048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.60468173) q[2];
sx q[2];
rz(-1.2102419) q[2];
sx q[2];
rz(0.49986419) q[2];
rz(2.8258421) q[3];
sx q[3];
rz(-2.4707268) q[3];
sx q[3];
rz(-1.0189112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4278118) q[0];
sx q[0];
rz(-1.918387) q[0];
sx q[0];
rz(-0.94386238) q[0];
rz(-1.7367412) q[1];
sx q[1];
rz(-1.2662788) q[1];
sx q[1];
rz(2.2199383) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93225828) q[0];
sx q[0];
rz(-0.32296041) q[0];
sx q[0];
rz(2.9684116) q[0];
rz(-pi) q[1];
rz(1.2723075) q[2];
sx q[2];
rz(-1.3278074) q[2];
sx q[2];
rz(-1.9858433) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.72295318) q[1];
sx q[1];
rz(-2.2721842) q[1];
sx q[1];
rz(0.69461125) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0579487) q[3];
sx q[3];
rz(-1.5745224) q[3];
sx q[3];
rz(-1.3386352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9965808) q[2];
sx q[2];
rz(-1.8700446) q[2];
sx q[2];
rz(-1.5865631) q[2];
rz(-1.0541213) q[3];
sx q[3];
rz(-1.328822) q[3];
sx q[3];
rz(1.740295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15040511) q[0];
sx q[0];
rz(-1.0973955) q[0];
sx q[0];
rz(0.64252585) q[0];
rz(-1.7036899) q[1];
sx q[1];
rz(-2.3846886) q[1];
sx q[1];
rz(0.14370758) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8929307) q[0];
sx q[0];
rz(-0.65048671) q[0];
sx q[0];
rz(2.6269795) q[0];
rz(0.92387284) q[2];
sx q[2];
rz(-1.1387741) q[2];
sx q[2];
rz(-2.1763755) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4776996) q[1];
sx q[1];
rz(-1.8346922) q[1];
sx q[1];
rz(-0.66910918) q[1];
rz(-pi) q[2];
rz(-0.71227422) q[3];
sx q[3];
rz(-0.55432075) q[3];
sx q[3];
rz(-1.8795151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6045427) q[2];
sx q[2];
rz(-1.9436676) q[2];
sx q[2];
rz(0.26270467) q[2];
rz(0.25775868) q[3];
sx q[3];
rz(-2.7596605) q[3];
sx q[3];
rz(0.8530544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2951374) q[0];
sx q[0];
rz(-1.4895804) q[0];
sx q[0];
rz(1.2204131) q[0];
rz(-0.48183164) q[1];
sx q[1];
rz(-2.316663) q[1];
sx q[1];
rz(-0.89734546) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7677444) q[0];
sx q[0];
rz(-1.4486827) q[0];
sx q[0];
rz(2.6227345) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4608896) q[2];
sx q[2];
rz(-1.8898003) q[2];
sx q[2];
rz(1.6750592) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7792977) q[1];
sx q[1];
rz(-0.26302281) q[1];
sx q[1];
rz(1.3785326) q[1];
x q[2];
rz(-0.96199022) q[3];
sx q[3];
rz(-1.655984) q[3];
sx q[3];
rz(-0.60768581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.64409488) q[2];
sx q[2];
rz(-2.2150025) q[2];
sx q[2];
rz(0.11478718) q[2];
rz(2.3980906) q[3];
sx q[3];
rz(-1.8758352) q[3];
sx q[3];
rz(0.44873294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7262309) q[0];
sx q[0];
rz(-2.3804433) q[0];
sx q[0];
rz(-0.95809715) q[0];
rz(1.6356946) q[1];
sx q[1];
rz(-1.1264569) q[1];
sx q[1];
rz(1.1503848) q[1];
rz(0.74581292) q[2];
sx q[2];
rz(-1.0619166) q[2];
sx q[2];
rz(-0.3084736) q[2];
rz(-0.43396797) q[3];
sx q[3];
rz(-1.4186191) q[3];
sx q[3];
rz(-0.71604244) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
