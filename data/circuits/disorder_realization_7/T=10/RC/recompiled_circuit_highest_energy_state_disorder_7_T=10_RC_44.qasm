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
rz(-1.7718908) q[0];
sx q[0];
rz(-2.7112024) q[0];
sx q[0];
rz(-2.9769467) q[0];
rz(1.4511664) q[1];
sx q[1];
rz(-2.7632406) q[1];
sx q[1];
rz(-0.70986706) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9116316) q[0];
sx q[0];
rz(-2.8459186) q[0];
sx q[0];
rz(1.1082746) q[0];
x q[1];
rz(0.69605445) q[2];
sx q[2];
rz(-2.5230319) q[2];
sx q[2];
rz(1.5355009) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8566411) q[1];
sx q[1];
rz(-1.2177693) q[1];
sx q[1];
rz(0.052327144) q[1];
rz(-pi) q[2];
rz(-0.48810256) q[3];
sx q[3];
rz(-2.8781366) q[3];
sx q[3];
rz(-1.4266222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5300488) q[2];
sx q[2];
rz(-0.350746) q[2];
sx q[2];
rz(-1.2190399) q[2];
rz(-0.48398316) q[3];
sx q[3];
rz(-0.84570208) q[3];
sx q[3];
rz(-1.774196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71500635) q[0];
sx q[0];
rz(-2.8026717) q[0];
sx q[0];
rz(-1.4434848) q[0];
rz(1.273607) q[1];
sx q[1];
rz(-1.0304291) q[1];
sx q[1];
rz(-0.66676203) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.043676959) q[0];
sx q[0];
rz(-2.7975492) q[0];
sx q[0];
rz(-2.6187569) q[0];
rz(-pi) q[1];
rz(-0.09632907) q[2];
sx q[2];
rz(-1.6821292) q[2];
sx q[2];
rz(-1.1679315) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.34109182) q[1];
sx q[1];
rz(-1.1665163) q[1];
sx q[1];
rz(1.0914108) q[1];
rz(-0.81979378) q[3];
sx q[3];
rz(-1.3043376) q[3];
sx q[3];
rz(0.72270061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4296809) q[2];
sx q[2];
rz(-1.1744171) q[2];
sx q[2];
rz(-2.5453117) q[2];
rz(-1.0239673) q[3];
sx q[3];
rz(-1.3919316) q[3];
sx q[3];
rz(-1.1202687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73151082) q[0];
sx q[0];
rz(-2.5937268) q[0];
sx q[0];
rz(1.4672853) q[0];
rz(0.038854988) q[1];
sx q[1];
rz(-1.424574) q[1];
sx q[1];
rz(0.88599667) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4873427) q[0];
sx q[0];
rz(-2.4375101) q[0];
sx q[0];
rz(0.6939597) q[0];
x q[1];
rz(-1.6396693) q[2];
sx q[2];
rz(-0.67126545) q[2];
sx q[2];
rz(-1.4466937) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8140602) q[1];
sx q[1];
rz(-2.7814354) q[1];
sx q[1];
rz(3.0292757) q[1];
rz(-2.233611) q[3];
sx q[3];
rz(-3.040979) q[3];
sx q[3];
rz(0.20357547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2140865) q[2];
sx q[2];
rz(-1.6897886) q[2];
sx q[2];
rz(0.56373325) q[2];
rz(2.2659414) q[3];
sx q[3];
rz(-2.0893658) q[3];
sx q[3];
rz(2.0503069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5617274) q[0];
sx q[0];
rz(-2.2932597) q[0];
sx q[0];
rz(-2.0735829) q[0];
rz(2.7470398) q[1];
sx q[1];
rz(-2.6119472) q[1];
sx q[1];
rz(-1.669917) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1354692) q[0];
sx q[0];
rz(-2.5803225) q[0];
sx q[0];
rz(-2.4563588) q[0];
x q[1];
rz(-1.285423) q[2];
sx q[2];
rz(-1.4625408) q[2];
sx q[2];
rz(1.6433604) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2100239) q[1];
sx q[1];
rz(-1.9940901) q[1];
sx q[1];
rz(-2.2970437) q[1];
rz(-2.7427892) q[3];
sx q[3];
rz(-1.670518) q[3];
sx q[3];
rz(-2.698363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7431405) q[2];
sx q[2];
rz(-1.3116216) q[2];
sx q[2];
rz(0.12942806) q[2];
rz(-2.5982924) q[3];
sx q[3];
rz(-2.048251) q[3];
sx q[3];
rz(-1.3885385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.038079) q[0];
sx q[0];
rz(-1.0032126) q[0];
sx q[0];
rz(-2.6858618) q[0];
rz(-2.575846) q[1];
sx q[1];
rz(-0.63876286) q[1];
sx q[1];
rz(2.9260213) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0545313) q[0];
sx q[0];
rz(-2.2534459) q[0];
sx q[0];
rz(2.9222617) q[0];
rz(-pi) q[1];
rz(0.59026123) q[2];
sx q[2];
rz(-2.1356842) q[2];
sx q[2];
rz(-1.60162) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.15102303) q[1];
sx q[1];
rz(-2.8672979) q[1];
sx q[1];
rz(1.3132233) q[1];
rz(-pi) q[2];
rz(-1.2108675) q[3];
sx q[3];
rz(-2.8382819) q[3];
sx q[3];
rz(-0.018050628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4244708) q[2];
sx q[2];
rz(-2.8346546) q[2];
sx q[2];
rz(-1.4324987) q[2];
rz(-2.0004382) q[3];
sx q[3];
rz(-1.2204095) q[3];
sx q[3];
rz(0.46877638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6345217) q[0];
sx q[0];
rz(-0.85352007) q[0];
sx q[0];
rz(-0.61990196) q[0];
rz(1.9339405) q[1];
sx q[1];
rz(-2.4687605) q[1];
sx q[1];
rz(2.8501453) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6648367) q[0];
sx q[0];
rz(-2.5546722) q[0];
sx q[0];
rz(-1.8700325) q[0];
rz(2.7619128) q[2];
sx q[2];
rz(-0.34647339) q[2];
sx q[2];
rz(-2.068067) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4174704) q[1];
sx q[1];
rz(-0.63878585) q[1];
sx q[1];
rz(-2.7330509) q[1];
x q[2];
rz(-1.5602371) q[3];
sx q[3];
rz(-2.227475) q[3];
sx q[3];
rz(-2.9845723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.42346272) q[2];
sx q[2];
rz(-2.4726157) q[2];
sx q[2];
rz(0.97406975) q[2];
rz(1.918321) q[3];
sx q[3];
rz(-1.7134824) q[3];
sx q[3];
rz(-2.1035002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.85422) q[0];
sx q[0];
rz(-1.9725476) q[0];
sx q[0];
rz(2.8337692) q[0];
rz(-0.27990118) q[1];
sx q[1];
rz(-0.24903909) q[1];
sx q[1];
rz(-1.8797967) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7578106) q[0];
sx q[0];
rz(-1.5368665) q[0];
sx q[0];
rz(1.4365804) q[0];
x q[1];
rz(-1.1818188) q[2];
sx q[2];
rz(-2.2643467) q[2];
sx q[2];
rz(2.3626568) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9692291) q[1];
sx q[1];
rz(-1.5595655) q[1];
sx q[1];
rz(1.2757236) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.070075017) q[3];
sx q[3];
rz(-1.106602) q[3];
sx q[3];
rz(-1.7572559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5032924) q[2];
sx q[2];
rz(-1.952848) q[2];
sx q[2];
rz(-0.80257455) q[2];
rz(-1.5447626) q[3];
sx q[3];
rz(-0.39512008) q[3];
sx q[3];
rz(1.4748658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0751727) q[0];
sx q[0];
rz(-1.5936699) q[0];
sx q[0];
rz(-2.5286034) q[0];
rz(2.0892443) q[1];
sx q[1];
rz(-0.78951013) q[1];
sx q[1];
rz(1.2552415) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4583115) q[0];
sx q[0];
rz(-2.0052687) q[0];
sx q[0];
rz(0.47554728) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0381211) q[2];
sx q[2];
rz(-1.3824859) q[2];
sx q[2];
rz(-1.7323158) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.44574983) q[1];
sx q[1];
rz(-0.90409213) q[1];
sx q[1];
rz(-2.1882344) q[1];
rz(-2.8792198) q[3];
sx q[3];
rz(-1.7011569) q[3];
sx q[3];
rz(-2.8487132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5455948) q[2];
sx q[2];
rz(-2.2277446) q[2];
sx q[2];
rz(2.9929898) q[2];
rz(-2.4871067) q[3];
sx q[3];
rz(-1.196967) q[3];
sx q[3];
rz(1.436208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8504836) q[0];
sx q[0];
rz(-1.5330667) q[0];
sx q[0];
rz(-0.0012375687) q[0];
rz(1.2507863) q[1];
sx q[1];
rz(-0.93235278) q[1];
sx q[1];
rz(2.1900182) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5244999) q[0];
sx q[0];
rz(-1.4513591) q[0];
sx q[0];
rz(2.9725084) q[0];
rz(-pi) q[1];
x q[1];
rz(0.96828332) q[2];
sx q[2];
rz(-0.70455019) q[2];
sx q[2];
rz(2.8271443) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.34500319) q[1];
sx q[1];
rz(-3.0751347) q[1];
sx q[1];
rz(2.4065325) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6660085) q[3];
sx q[3];
rz(-1.9692752) q[3];
sx q[3];
rz(1.0131032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1452267) q[2];
sx q[2];
rz(-1.9731015) q[2];
sx q[2];
rz(-2.1053947) q[2];
rz(0.869831) q[3];
sx q[3];
rz(-1.4709604) q[3];
sx q[3];
rz(-0.74365348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0319801) q[0];
sx q[0];
rz(-0.59015048) q[0];
sx q[0];
rz(-1.1085283) q[0];
rz(0.96500665) q[1];
sx q[1];
rz(-1.7644278) q[1];
sx q[1];
rz(2.5722497) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78380221) q[0];
sx q[0];
rz(-2.0817502) q[0];
sx q[0];
rz(-1.820867) q[0];
rz(-2.0351719) q[2];
sx q[2];
rz(-1.9074252) q[2];
sx q[2];
rz(-2.6576633) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3778622) q[1];
sx q[1];
rz(-1.5625349) q[1];
sx q[1];
rz(-2.0935835) q[1];
x q[2];
rz(0.17579349) q[3];
sx q[3];
rz(-0.7657649) q[3];
sx q[3];
rz(-0.24336963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1278594) q[2];
sx q[2];
rz(-0.49488417) q[2];
sx q[2];
rz(-3.0950756) q[2];
rz(0.30139309) q[3];
sx q[3];
rz(-1.9207585) q[3];
sx q[3];
rz(-2.9197599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90060577) q[0];
sx q[0];
rz(-1.2099246) q[0];
sx q[0];
rz(-1.8151617) q[0];
rz(-1.2251414) q[1];
sx q[1];
rz(-2.6381208) q[1];
sx q[1];
rz(-2.0874964) q[1];
rz(-0.27041783) q[2];
sx q[2];
rz(-2.5106988) q[2];
sx q[2];
rz(-0.0028263447) q[2];
rz(2.2598738) q[3];
sx q[3];
rz(-0.93986012) q[3];
sx q[3];
rz(-2.8840251) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
