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
rz(0.61697382) q[1];
sx q[1];
rz(-0.66520912) q[1];
sx q[1];
rz(1.3245423) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0966372) q[0];
sx q[0];
rz(-2.5384266) q[0];
sx q[0];
rz(-3.0289193) q[0];
x q[1];
rz(-0.1466897) q[2];
sx q[2];
rz(-2.5518199) q[2];
sx q[2];
rz(1.1663933) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1086667) q[1];
sx q[1];
rz(-0.98331988) q[1];
sx q[1];
rz(1.0125748) q[1];
rz(-pi) q[2];
rz(0.67456986) q[3];
sx q[3];
rz(-2.2534182) q[3];
sx q[3];
rz(2.5404524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2609743) q[2];
sx q[2];
rz(-1.3857434) q[2];
sx q[2];
rz(2.3495038) q[2];
rz(0.875862) q[3];
sx q[3];
rz(-3.0328817) q[3];
sx q[3];
rz(0.66332269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6302014) q[0];
sx q[0];
rz(-1.317861) q[0];
sx q[0];
rz(2.200101) q[0];
rz(-0.9990274) q[1];
sx q[1];
rz(-0.92728725) q[1];
sx q[1];
rz(0.082854465) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1076528) q[0];
sx q[0];
rz(-1.4705188) q[0];
sx q[0];
rz(-2.4058008) q[0];
x q[1];
rz(1.9272958) q[2];
sx q[2];
rz(-1.7983603) q[2];
sx q[2];
rz(2.5806732) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.99020834) q[1];
sx q[1];
rz(-0.71077222) q[1];
sx q[1];
rz(-0.73020331) q[1];
x q[2];
rz(0.96534713) q[3];
sx q[3];
rz(-2.1248528) q[3];
sx q[3];
rz(-0.29747552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2443709) q[2];
sx q[2];
rz(-1.3542078) q[2];
sx q[2];
rz(1.2574035) q[2];
rz(-0.8463549) q[3];
sx q[3];
rz(-0.1592764) q[3];
sx q[3];
rz(-2.4634821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9427247) q[0];
sx q[0];
rz(-1.6833479) q[0];
sx q[0];
rz(2.9123836) q[0];
rz(-0.21408679) q[1];
sx q[1];
rz(-0.87132088) q[1];
sx q[1];
rz(0.3814989) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4320735) q[0];
sx q[0];
rz(-2.1901178) q[0];
sx q[0];
rz(1.1362057) q[0];
x q[1];
rz(0.83865954) q[2];
sx q[2];
rz(-1.7146972) q[2];
sx q[2];
rz(3.0001907) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2974071) q[1];
sx q[1];
rz(-1.5995889) q[1];
sx q[1];
rz(1.6237358) q[1];
x q[2];
rz(-1.5656785) q[3];
sx q[3];
rz(-1.9178784) q[3];
sx q[3];
rz(-2.0818613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4182338) q[2];
sx q[2];
rz(-0.25622076) q[2];
sx q[2];
rz(-0.56824938) q[2];
rz(-1.4189643) q[3];
sx q[3];
rz(-1.4293554) q[3];
sx q[3];
rz(-1.0770146) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1133465) q[0];
sx q[0];
rz(-1.1320817) q[0];
sx q[0];
rz(0.25076732) q[0];
rz(0.084550683) q[1];
sx q[1];
rz(-2.0471579) q[1];
sx q[1];
rz(1.9409404) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47287175) q[0];
sx q[0];
rz(-0.78470147) q[0];
sx q[0];
rz(-2.0913823) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0798825) q[2];
sx q[2];
rz(-2.6497095) q[2];
sx q[2];
rz(1.5668004) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.975864) q[1];
sx q[1];
rz(-1.670125) q[1];
sx q[1];
rz(-0.71290675) q[1];
rz(-2.8765466) q[3];
sx q[3];
rz(-1.0599905) q[3];
sx q[3];
rz(-3.0731415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4667929) q[2];
sx q[2];
rz(-1.1226706) q[2];
sx q[2];
rz(1.8505081) q[2];
rz(0.42168266) q[3];
sx q[3];
rz(-0.45160523) q[3];
sx q[3];
rz(-1.5765367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6186433) q[0];
sx q[0];
rz(-0.72294253) q[0];
sx q[0];
rz(-1.6408386) q[0];
rz(-3.07395) q[1];
sx q[1];
rz(-1.5875971) q[1];
sx q[1];
rz(1.1281475) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65928946) q[0];
sx q[0];
rz(-1.7586252) q[0];
sx q[0];
rz(1.5330305) q[0];
rz(-1.0063926) q[2];
sx q[2];
rz(-0.60205209) q[2];
sx q[2];
rz(2.4133701) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.50666821) q[1];
sx q[1];
rz(-1.6841869) q[1];
sx q[1];
rz(3.0874814) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2637783) q[3];
sx q[3];
rz(-0.35696176) q[3];
sx q[3];
rz(-1.9138543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.03380123) q[2];
sx q[2];
rz(-1.1148323) q[2];
sx q[2];
rz(2.4647253) q[2];
rz(-2.233861) q[3];
sx q[3];
rz(-1.9151442) q[3];
sx q[3];
rz(0.41456732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9229729) q[0];
sx q[0];
rz(-0.75160471) q[0];
sx q[0];
rz(2.3722017) q[0];
rz(-1.9006624) q[1];
sx q[1];
rz(-0.85378328) q[1];
sx q[1];
rz(0.21533899) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0701987) q[0];
sx q[0];
rz(-0.60766534) q[0];
sx q[0];
rz(-1.2383226) q[0];
rz(-pi) q[1];
rz(-1.5211283) q[2];
sx q[2];
rz(-2.81041) q[2];
sx q[2];
rz(-2.0803723) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0907545) q[1];
sx q[1];
rz(-1.8148481) q[1];
sx q[1];
rz(0.65907101) q[1];
x q[2];
rz(1.2329383) q[3];
sx q[3];
rz(-0.90150276) q[3];
sx q[3];
rz(3.11495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.26663366) q[2];
sx q[2];
rz(-1.6338438) q[2];
sx q[2];
rz(2.5249262) q[2];
rz(0.2002317) q[3];
sx q[3];
rz(-0.73035208) q[3];
sx q[3];
rz(1.1327789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.013833372) q[0];
sx q[0];
rz(-2.5795689) q[0];
sx q[0];
rz(3.1264937) q[0];
rz(0.12241441) q[1];
sx q[1];
rz(-1.8465123) q[1];
sx q[1];
rz(-2.1133568) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.016617386) q[0];
sx q[0];
rz(-2.4835971) q[0];
sx q[0];
rz(-3.0506177) q[0];
rz(-pi) q[1];
x q[1];
rz(0.21034849) q[2];
sx q[2];
rz(-0.96322434) q[2];
sx q[2];
rz(1.6574394) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.084000951) q[1];
sx q[1];
rz(-1.1391974) q[1];
sx q[1];
rz(-1.6017385) q[1];
x q[2];
rz(0.31556231) q[3];
sx q[3];
rz(-1.3129819) q[3];
sx q[3];
rz(2.1699048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5369109) q[2];
sx q[2];
rz(-1.2102419) q[2];
sx q[2];
rz(-0.49986419) q[2];
rz(-0.31575051) q[3];
sx q[3];
rz(-2.4707268) q[3];
sx q[3];
rz(2.1226814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71378088) q[0];
sx q[0];
rz(-1.918387) q[0];
sx q[0];
rz(-2.1977303) q[0];
rz(-1.4048514) q[1];
sx q[1];
rz(-1.8753139) q[1];
sx q[1];
rz(2.2199383) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2093344) q[0];
sx q[0];
rz(-0.32296041) q[0];
sx q[0];
rz(-2.9684116) q[0];
rz(2.8878288) q[2];
sx q[2];
rz(-1.2813338) q[2];
sx q[2];
rz(0.34115215) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4186395) q[1];
sx q[1];
rz(-2.2721842) q[1];
sx q[1];
rz(-0.69461125) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1373164) q[3];
sx q[3];
rz(-2.0836401) q[3];
sx q[3];
rz(-2.9073334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9965808) q[2];
sx q[2];
rz(-1.8700446) q[2];
sx q[2];
rz(-1.5550295) q[2];
rz(-1.0541213) q[3];
sx q[3];
rz(-1.328822) q[3];
sx q[3];
rz(1.740295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15040511) q[0];
sx q[0];
rz(-1.0973955) q[0];
sx q[0];
rz(-0.64252585) q[0];
rz(-1.7036899) q[1];
sx q[1];
rz(-0.75690401) q[1];
sx q[1];
rz(-0.14370758) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7449581) q[0];
sx q[0];
rz(-1.2681343) q[0];
sx q[0];
rz(-0.58505262) q[0];
rz(-pi) q[1];
rz(-2.6176378) q[2];
sx q[2];
rz(-2.1499976) q[2];
sx q[2];
rz(-2.8423345) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4776996) q[1];
sx q[1];
rz(-1.8346922) q[1];
sx q[1];
rz(2.4724835) q[1];
rz(-0.71227422) q[3];
sx q[3];
rz(-0.55432075) q[3];
sx q[3];
rz(-1.8795151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5370499) q[2];
sx q[2];
rz(-1.9436676) q[2];
sx q[2];
rz(2.878888) q[2];
rz(-2.883834) q[3];
sx q[3];
rz(-0.38193211) q[3];
sx q[3];
rz(-0.8530544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8464552) q[0];
sx q[0];
rz(-1.6520123) q[0];
sx q[0];
rz(1.2204131) q[0];
rz(-0.48183164) q[1];
sx q[1];
rz(-2.316663) q[1];
sx q[1];
rz(-0.89734546) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37384826) q[0];
sx q[0];
rz(-1.4486827) q[0];
sx q[0];
rz(2.6227345) q[0];
rz(-0.32081066) q[2];
sx q[2];
rz(-1.4664553) q[2];
sx q[2];
rz(-3.0719245) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5612302) q[1];
sx q[1];
rz(-1.8288611) q[1];
sx q[1];
rz(0.051405426) q[1];
rz(-pi) q[2];
rz(0.10372377) q[3];
sx q[3];
rz(-0.96451603) q[3];
sx q[3];
rz(1.0223573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4974978) q[2];
sx q[2];
rz(-0.92659014) q[2];
sx q[2];
rz(0.11478718) q[2];
rz(-0.74350205) q[3];
sx q[3];
rz(-1.2657575) q[3];
sx q[3];
rz(2.6928597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.92123564) q[2];
sx q[2];
rz(-0.9365281) q[2];
sx q[2];
rz(-2.3021883) q[2];
rz(1.7382449) q[3];
sx q[3];
rz(-1.1421775) q[3];
sx q[3];
rz(0.92489064) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
