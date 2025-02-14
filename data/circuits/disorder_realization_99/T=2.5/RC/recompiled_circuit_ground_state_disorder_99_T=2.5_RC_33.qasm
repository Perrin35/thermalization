OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.34243256) q[0];
sx q[0];
rz(-0.39781308) q[0];
sx q[0];
rz(1.8737268) q[0];
rz(-0.59029382) q[1];
sx q[1];
rz(2.3550912) q[1];
sx q[1];
rz(13.785706) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3800664) q[0];
sx q[0];
rz(-0.60193578) q[0];
sx q[0];
rz(-1.8695037) q[0];
rz(-2.8706495) q[2];
sx q[2];
rz(-2.5763955) q[2];
sx q[2];
rz(2.924356) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.87675512) q[1];
sx q[1];
rz(-1.7746468) q[1];
sx q[1];
rz(-1.8556103) q[1];
rz(0.23363913) q[3];
sx q[3];
rz(-2.2147182) q[3];
sx q[3];
rz(-1.131191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5775602) q[2];
sx q[2];
rz(-1.5215678) q[2];
sx q[2];
rz(1.7731898) q[2];
rz(0.91156256) q[3];
sx q[3];
rz(-1.3699968) q[3];
sx q[3];
rz(-3.1172359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1034705) q[0];
sx q[0];
rz(-0.39762527) q[0];
sx q[0];
rz(0.11904112) q[0];
rz(1.1301522) q[1];
sx q[1];
rz(-2.1739013) q[1];
sx q[1];
rz(1.4281323) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6551681) q[0];
sx q[0];
rz(-1.055724) q[0];
sx q[0];
rz(-1.0977655) q[0];
x q[1];
rz(-0.4006673) q[2];
sx q[2];
rz(-1.9611437) q[2];
sx q[2];
rz(2.0290749) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2547639) q[1];
sx q[1];
rz(-2.6568251) q[1];
sx q[1];
rz(-1.280275) q[1];
rz(-pi) q[2];
rz(0.30562206) q[3];
sx q[3];
rz(-1.8942617) q[3];
sx q[3];
rz(0.68438578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.52593645) q[2];
sx q[2];
rz(-1.4639414) q[2];
sx q[2];
rz(-0.76457912) q[2];
rz(2.6584451) q[3];
sx q[3];
rz(-1.8254447) q[3];
sx q[3];
rz(-0.84883261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0371542) q[0];
sx q[0];
rz(-0.26532441) q[0];
sx q[0];
rz(1.4720488) q[0];
rz(-2.5813591) q[1];
sx q[1];
rz(-1.3872223) q[1];
sx q[1];
rz(1.4405506) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.154523) q[0];
sx q[0];
rz(-1.3366342) q[0];
sx q[0];
rz(-0.20407853) q[0];
rz(-pi) q[1];
rz(-3.1344536) q[2];
sx q[2];
rz(-2.4666365) q[2];
sx q[2];
rz(-1.0085886) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8941108) q[1];
sx q[1];
rz(-1.9097345) q[1];
sx q[1];
rz(-2.1268658) q[1];
rz(-0.76466842) q[3];
sx q[3];
rz(-0.32161412) q[3];
sx q[3];
rz(2.4177616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4855839) q[2];
sx q[2];
rz(-1.6501082) q[2];
sx q[2];
rz(0.34240016) q[2];
rz(-1.7229236) q[3];
sx q[3];
rz(-2.4125621) q[3];
sx q[3];
rz(-2.5820144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90958428) q[0];
sx q[0];
rz(-0.37264687) q[0];
sx q[0];
rz(-1.5268071) q[0];
rz(-2.3341663) q[1];
sx q[1];
rz(-1.5734943) q[1];
sx q[1];
rz(-1.2015013) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3897409) q[0];
sx q[0];
rz(-0.83503113) q[0];
sx q[0];
rz(2.7324008) q[0];
rz(-pi) q[1];
rz(-0.64862675) q[2];
sx q[2];
rz(-1.0264215) q[2];
sx q[2];
rz(-0.2444765) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0302194) q[1];
sx q[1];
rz(-0.89240701) q[1];
sx q[1];
rz(-2.6256843) q[1];
x q[2];
rz(-2.9462847) q[3];
sx q[3];
rz(-2.8805519) q[3];
sx q[3];
rz(-1.9633499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.80369192) q[2];
sx q[2];
rz(-2.1521229) q[2];
sx q[2];
rz(1.9039512) q[2];
rz(-1.8303653) q[3];
sx q[3];
rz(-1.6236191) q[3];
sx q[3];
rz(-1.9633912) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.118498) q[0];
sx q[0];
rz(-1.3730405) q[0];
sx q[0];
rz(2.232724) q[0];
rz(-0.88242775) q[1];
sx q[1];
rz(-0.71417037) q[1];
sx q[1];
rz(-0.73208255) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86504146) q[0];
sx q[0];
rz(-2.2892729) q[0];
sx q[0];
rz(0.9263692) q[0];
x q[1];
rz(-1.2901487) q[2];
sx q[2];
rz(-0.82260231) q[2];
sx q[2];
rz(1.4952687) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8632224) q[1];
sx q[1];
rz(-1.9626106) q[1];
sx q[1];
rz(2.6696221) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74589269) q[3];
sx q[3];
rz(-2.159366) q[3];
sx q[3];
rz(2.4323127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.070179209) q[2];
sx q[2];
rz(-0.81255239) q[2];
sx q[2];
rz(1.2893691) q[2];
rz(1.6728801) q[3];
sx q[3];
rz(-1.9332935) q[3];
sx q[3];
rz(2.7220791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5656723) q[0];
sx q[0];
rz(-1.9603632) q[0];
sx q[0];
rz(-1.1015724) q[0];
rz(2.9167602) q[1];
sx q[1];
rz(-1.79554) q[1];
sx q[1];
rz(2.1563931) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6533642) q[0];
sx q[0];
rz(-2.4765827) q[0];
sx q[0];
rz(-1.8652161) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2217789) q[2];
sx q[2];
rz(-0.98528457) q[2];
sx q[2];
rz(1.7988009) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0538865) q[1];
sx q[1];
rz(-0.83982491) q[1];
sx q[1];
rz(2.5165416) q[1];
x q[2];
rz(0.52805488) q[3];
sx q[3];
rz(-1.5925358) q[3];
sx q[3];
rz(-2.4491212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.78099403) q[2];
sx q[2];
rz(-0.66086078) q[2];
sx q[2];
rz(-1.9471656) q[2];
rz(3.0418975) q[3];
sx q[3];
rz(-1.3705148) q[3];
sx q[3];
rz(-2.1626507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7128971) q[0];
sx q[0];
rz(-0.45658699) q[0];
sx q[0];
rz(-2.0275443) q[0];
rz(-1.3975337) q[1];
sx q[1];
rz(-1.7618529) q[1];
sx q[1];
rz(0.31375113) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6614051) q[0];
sx q[0];
rz(-1.620694) q[0];
sx q[0];
rz(-1.2779092) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4391104) q[2];
sx q[2];
rz(-1.734231) q[2];
sx q[2];
rz(-1.861546) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.75668979) q[1];
sx q[1];
rz(-0.80763615) q[1];
sx q[1];
rz(-2.7212423) q[1];
rz(-2.8242495) q[3];
sx q[3];
rz(-2.4306524) q[3];
sx q[3];
rz(-1.8985572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.067001192) q[2];
sx q[2];
rz(-2.7374697) q[2];
sx q[2];
rz(-0.60603777) q[2];
rz(-1.0780942) q[3];
sx q[3];
rz(-0.86229101) q[3];
sx q[3];
rz(-1.3585453) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9692877) q[0];
sx q[0];
rz(-2.2242039) q[0];
sx q[0];
rz(1.1267927) q[0];
rz(2.2276095) q[1];
sx q[1];
rz(-0.4387478) q[1];
sx q[1];
rz(-0.21805683) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1153961) q[0];
sx q[0];
rz(-0.65803981) q[0];
sx q[0];
rz(-0.81320073) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.70773196) q[2];
sx q[2];
rz(-0.92365217) q[2];
sx q[2];
rz(0.92491481) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.14971241) q[1];
sx q[1];
rz(-2.212095) q[1];
sx q[1];
rz(2.9845357) q[1];
rz(-pi) q[2];
rz(-3.0239437) q[3];
sx q[3];
rz(-1.434064) q[3];
sx q[3];
rz(-2.4916388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5128936) q[2];
sx q[2];
rz(-0.67535526) q[2];
sx q[2];
rz(2.8524354) q[2];
rz(-2.1469927) q[3];
sx q[3];
rz(-1.9528439) q[3];
sx q[3];
rz(2.6836256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7116123) q[0];
sx q[0];
rz(-2.0531605) q[0];
sx q[0];
rz(-2.3751538) q[0];
rz(-1.537716) q[1];
sx q[1];
rz(-1.8285373) q[1];
sx q[1];
rz(0.91805735) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9452451) q[0];
sx q[0];
rz(-2.1736988) q[0];
sx q[0];
rz(2.5640998) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7154022) q[2];
sx q[2];
rz(-1.7376126) q[2];
sx q[2];
rz(1.6289935) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.80528211) q[1];
sx q[1];
rz(-0.86617058) q[1];
sx q[1];
rz(1.8971838) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3206743) q[3];
sx q[3];
rz(-1.6733352) q[3];
sx q[3];
rz(-0.96270442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.90091577) q[2];
sx q[2];
rz(-0.46611163) q[2];
sx q[2];
rz(0.051699836) q[2];
rz(-0.85203552) q[3];
sx q[3];
rz(-1.4308735) q[3];
sx q[3];
rz(2.8813072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.32605115) q[0];
sx q[0];
rz(-1.4651848) q[0];
sx q[0];
rz(0.091212243) q[0];
rz(-2.3444029) q[1];
sx q[1];
rz(-1.7557095) q[1];
sx q[1];
rz(-2.7104654) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6043458) q[0];
sx q[0];
rz(-1.6291999) q[0];
sx q[0];
rz(1.6184774) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2195039) q[2];
sx q[2];
rz(-1.6493634) q[2];
sx q[2];
rz(2.9804413) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8517837) q[1];
sx q[1];
rz(-2.0153592) q[1];
sx q[1];
rz(-1.3291994) q[1];
rz(-0.43041269) q[3];
sx q[3];
rz(-1.0135883) q[3];
sx q[3];
rz(0.9924953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5690696) q[2];
sx q[2];
rz(-1.5670245) q[2];
sx q[2];
rz(1.265906) q[2];
rz(-1.2931394) q[3];
sx q[3];
rz(-2.0796516) q[3];
sx q[3];
rz(-2.7354447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.010392808) q[0];
sx q[0];
rz(-0.69080234) q[0];
sx q[0];
rz(0.77558415) q[0];
rz(2.5776183) q[1];
sx q[1];
rz(-1.0697983) q[1];
sx q[1];
rz(-1.0615798) q[1];
rz(-1.1376913) q[2];
sx q[2];
rz(-1.644293) q[2];
sx q[2];
rz(-2.1183011) q[2];
rz(0.33071409) q[3];
sx q[3];
rz(-1.1138108) q[3];
sx q[3];
rz(2.6295233) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
