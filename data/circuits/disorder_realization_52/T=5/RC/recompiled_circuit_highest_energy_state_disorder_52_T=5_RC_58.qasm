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
rz(-1.4327383) q[0];
sx q[0];
rz(-0.32810768) q[0];
sx q[0];
rz(-2.3018667) q[0];
rz(-1.7307164) q[1];
sx q[1];
rz(-1.6003992) q[1];
sx q[1];
rz(2.3720001) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1043848) q[0];
sx q[0];
rz(-1.6055371) q[0];
sx q[0];
rz(-1.4912692) q[0];
rz(-pi) q[1];
rz(0.02564557) q[2];
sx q[2];
rz(-2.316595) q[2];
sx q[2];
rz(-2.1767669) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.283764) q[1];
sx q[1];
rz(-0.35175465) q[1];
sx q[1];
rz(-1.5726388) q[1];
rz(-pi) q[2];
rz(2.493164) q[3];
sx q[3];
rz(-0.75639137) q[3];
sx q[3];
rz(-1.3677011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8771693) q[2];
sx q[2];
rz(-1.2494272) q[2];
sx q[2];
rz(-1.1733615) q[2];
rz(-0.42826432) q[3];
sx q[3];
rz(-2.5731125) q[3];
sx q[3];
rz(2.4885524) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1332755) q[0];
sx q[0];
rz(-2.7011217) q[0];
sx q[0];
rz(1.0622729) q[0];
rz(1.5509037) q[1];
sx q[1];
rz(-2.5804602) q[1];
sx q[1];
rz(-1.0947469) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0355201) q[0];
sx q[0];
rz(-2.8815334) q[0];
sx q[0];
rz(-1.7709269) q[0];
rz(-pi) q[1];
rz(-2.1065478) q[2];
sx q[2];
rz(-2.3887815) q[2];
sx q[2];
rz(-2.8967146) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.96910406) q[1];
sx q[1];
rz(-1.2840233) q[1];
sx q[1];
rz(-2.7412834) q[1];
rz(0.21778222) q[3];
sx q[3];
rz(-2.3541303) q[3];
sx q[3];
rz(-1.2125208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.67794472) q[2];
sx q[2];
rz(-2.8539168) q[2];
sx q[2];
rz(0.90052432) q[2];
rz(2.1053704) q[3];
sx q[3];
rz(-0.13768727) q[3];
sx q[3];
rz(0.010995939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4765332) q[0];
sx q[0];
rz(-2.0158975) q[0];
sx q[0];
rz(1.4680468) q[0];
rz(0.92848575) q[1];
sx q[1];
rz(-1.0037582) q[1];
sx q[1];
rz(-1.4220062) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6674468) q[0];
sx q[0];
rz(-0.26413879) q[0];
sx q[0];
rz(2.788782) q[0];
x q[1];
rz(0.85365752) q[2];
sx q[2];
rz(-2.2480534) q[2];
sx q[2];
rz(-0.028606107) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.068068935) q[1];
sx q[1];
rz(-1.0367107) q[1];
sx q[1];
rz(1.1492386) q[1];
rz(-pi) q[2];
rz(2.7169796) q[3];
sx q[3];
rz(-1.8311005) q[3];
sx q[3];
rz(-1.9003968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.15098393) q[2];
sx q[2];
rz(-1.7093806) q[2];
sx q[2];
rz(-0.81986156) q[2];
rz(0.12623434) q[3];
sx q[3];
rz(-1.1773959) q[3];
sx q[3];
rz(-1.646515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4594629) q[0];
sx q[0];
rz(-1.7403025) q[0];
sx q[0];
rz(-1.6023741) q[0];
rz(1.0914717) q[1];
sx q[1];
rz(-0.83408728) q[1];
sx q[1];
rz(1.9713255) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4919969) q[0];
sx q[0];
rz(-1.2536712) q[0];
sx q[0];
rz(1.1742422) q[0];
rz(-pi) q[1];
rz(2.290904) q[2];
sx q[2];
rz(-1.3733092) q[2];
sx q[2];
rz(-2.6732004) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5398769) q[1];
sx q[1];
rz(-1.6605494) q[1];
sx q[1];
rz(1.490277) q[1];
x q[2];
rz(0.50962944) q[3];
sx q[3];
rz(-2.7197803) q[3];
sx q[3];
rz(3.075656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9200661) q[2];
sx q[2];
rz(-2.2130794) q[2];
sx q[2];
rz(-0.083219223) q[2];
rz(-0.65256882) q[3];
sx q[3];
rz(-3.1196399) q[3];
sx q[3];
rz(2.2799802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49587747) q[0];
sx q[0];
rz(-0.28384122) q[0];
sx q[0];
rz(-0.19677095) q[0];
rz(1.1605877) q[1];
sx q[1];
rz(-2.6996758) q[1];
sx q[1];
rz(1.7534076) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3236448) q[0];
sx q[0];
rz(-2.9716688) q[0];
sx q[0];
rz(3.1325045) q[0];
x q[1];
rz(0.20241995) q[2];
sx q[2];
rz(-1.4223863) q[2];
sx q[2];
rz(-2.9925516) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.95294556) q[1];
sx q[1];
rz(-1.7428696) q[1];
sx q[1];
rz(0.91305542) q[1];
rz(-pi) q[2];
rz(1.1673683) q[3];
sx q[3];
rz(-0.77047548) q[3];
sx q[3];
rz(1.6516754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.56472003) q[2];
sx q[2];
rz(-1.9331965) q[2];
sx q[2];
rz(-1.8604856) q[2];
rz(-2.7992904) q[3];
sx q[3];
rz(-3.0908995) q[3];
sx q[3];
rz(-2.2127693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.36011919) q[0];
sx q[0];
rz(-2.0442648) q[0];
sx q[0];
rz(2.1088364) q[0];
rz(-0.67689854) q[1];
sx q[1];
rz(-2.758226) q[1];
sx q[1];
rz(2.3597609) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1843518) q[0];
sx q[0];
rz(-2.4515984) q[0];
sx q[0];
rz(2.5646567) q[0];
rz(3.1091613) q[2];
sx q[2];
rz(-1.6219181) q[2];
sx q[2];
rz(-2.0699208) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1345152) q[1];
sx q[1];
rz(-2.482153) q[1];
sx q[1];
rz(-0.51523988) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6138462) q[3];
sx q[3];
rz(-1.2331748) q[3];
sx q[3];
rz(-1.7204325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3013762) q[2];
sx q[2];
rz(-0.8679114) q[2];
sx q[2];
rz(-0.8737348) q[2];
rz(-0.78911632) q[3];
sx q[3];
rz(-1.5566166) q[3];
sx q[3];
rz(0.48621392) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9588722) q[0];
sx q[0];
rz(-0.046165753) q[0];
sx q[0];
rz(-0.090855457) q[0];
rz(0.60984045) q[1];
sx q[1];
rz(-1.5515386) q[1];
sx q[1];
rz(-0.59744936) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0598568) q[0];
sx q[0];
rz(-1.6993521) q[0];
sx q[0];
rz(-2.9946694) q[0];
x q[1];
rz(1.9436854) q[2];
sx q[2];
rz(-1.882236) q[2];
sx q[2];
rz(-2.5562037) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6521921) q[1];
sx q[1];
rz(-2.0058549) q[1];
sx q[1];
rz(-1.0580741) q[1];
rz(-pi) q[2];
rz(-3.1037056) q[3];
sx q[3];
rz(-1.1514613) q[3];
sx q[3];
rz(1.9333206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.37683836) q[2];
sx q[2];
rz(-2.59616) q[2];
sx q[2];
rz(-3.0010014) q[2];
rz(1.8600672) q[3];
sx q[3];
rz(-1.8257273) q[3];
sx q[3];
rz(0.52711058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22290467) q[0];
sx q[0];
rz(-1.0162901) q[0];
sx q[0];
rz(-1.6131529) q[0];
rz(-0.81360045) q[1];
sx q[1];
rz(-3.0366812) q[1];
sx q[1];
rz(2.334107) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59283108) q[0];
sx q[0];
rz(-0.27540311) q[0];
sx q[0];
rz(2.3345678) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0139461) q[2];
sx q[2];
rz(-1.6880182) q[2];
sx q[2];
rz(1.6991311) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.85881104) q[1];
sx q[1];
rz(-0.75560299) q[1];
sx q[1];
rz(2.1367461) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.42515884) q[3];
sx q[3];
rz(-2.0704554) q[3];
sx q[3];
rz(2.5806346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3537102) q[2];
sx q[2];
rz(-1.3545802) q[2];
sx q[2];
rz(-2.9969969) q[2];
rz(2.0374129) q[3];
sx q[3];
rz(-2.927533) q[3];
sx q[3];
rz(-1.0415227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60177326) q[0];
sx q[0];
rz(-2.0105392) q[0];
sx q[0];
rz(2.5850776) q[0];
rz(-0.76983184) q[1];
sx q[1];
rz(-0.12431215) q[1];
sx q[1];
rz(0.46129033) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88474482) q[0];
sx q[0];
rz(-2.1408653) q[0];
sx q[0];
rz(0.056572551) q[0];
x q[1];
rz(-0.83314216) q[2];
sx q[2];
rz(-2.1066446) q[2];
sx q[2];
rz(-1.7876724) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0374679) q[1];
sx q[1];
rz(-2.3437683) q[1];
sx q[1];
rz(-0.41607694) q[1];
rz(-0.56311816) q[3];
sx q[3];
rz(-0.65059911) q[3];
sx q[3];
rz(-1.4103149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7028659) q[2];
sx q[2];
rz(-0.30539572) q[2];
sx q[2];
rz(-2.3766282) q[2];
rz(1.9276098) q[3];
sx q[3];
rz(-0.58911222) q[3];
sx q[3];
rz(-2.7629619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42534378) q[0];
sx q[0];
rz(-3.0916164) q[0];
sx q[0];
rz(-2.7783527) q[0];
rz(-1.4092457) q[1];
sx q[1];
rz(-2.1556518) q[1];
sx q[1];
rz(0.22663103) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0186414) q[0];
sx q[0];
rz(-1.567807) q[0];
sx q[0];
rz(0.13521533) q[0];
x q[1];
rz(-0.66618528) q[2];
sx q[2];
rz(-2.2720784) q[2];
sx q[2];
rz(1.2291069) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5635101) q[1];
sx q[1];
rz(-2.3704766) q[1];
sx q[1];
rz(-1.095849) q[1];
rz(0.33785401) q[3];
sx q[3];
rz(-1.5607335) q[3];
sx q[3];
rz(0.73051605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5994485) q[2];
sx q[2];
rz(-2.6207974) q[2];
sx q[2];
rz(-0.64860541) q[2];
rz(-0.058622807) q[3];
sx q[3];
rz(-0.62871814) q[3];
sx q[3];
rz(0.64529836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36515737) q[0];
sx q[0];
rz(-1.2189652) q[0];
sx q[0];
rz(-0.49268876) q[0];
rz(1.0538712) q[1];
sx q[1];
rz(-2.5851879) q[1];
sx q[1];
rz(-2.7256706) q[1];
rz(0.024439288) q[2];
sx q[2];
rz(-1.8069488) q[2];
sx q[2];
rz(-1.9843742) q[2];
rz(-2.2661187) q[3];
sx q[3];
rz(-1.8322104) q[3];
sx q[3];
rz(-1.9598243) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
