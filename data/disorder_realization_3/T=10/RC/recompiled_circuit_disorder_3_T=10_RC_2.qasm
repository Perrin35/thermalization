OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.91575032) q[0];
sx q[0];
rz(-3.1103818) q[0];
sx q[0];
rz(-2.6565235) q[0];
rz(-2.354061) q[1];
sx q[1];
rz(-2.1252316) q[1];
sx q[1];
rz(0.41419849) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.480455) q[0];
sx q[0];
rz(-1.9248795) q[0];
sx q[0];
rz(0.35868355) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0010927) q[2];
sx q[2];
rz(-1.4909407) q[2];
sx q[2];
rz(0.18730883) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2990883) q[1];
sx q[1];
rz(-0.43049225) q[1];
sx q[1];
rz(1.4166142) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8768086) q[3];
sx q[3];
rz(-1.4720819) q[3];
sx q[3];
rz(-2.7717154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9238613) q[2];
sx q[2];
rz(-1.2552746) q[2];
sx q[2];
rz(0.031575354) q[2];
rz(1.8850373) q[3];
sx q[3];
rz(-0.50285181) q[3];
sx q[3];
rz(0.52662915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8388222) q[0];
sx q[0];
rz(-1.6844203) q[0];
sx q[0];
rz(-0.1698499) q[0];
rz(0.70392144) q[1];
sx q[1];
rz(-2.070065) q[1];
sx q[1];
rz(-0.53952113) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9233421) q[0];
sx q[0];
rz(-1.5241511) q[0];
sx q[0];
rz(2.0204861) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.72210724) q[2];
sx q[2];
rz(-1.2566084) q[2];
sx q[2];
rz(-2.9177641) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.19903781) q[1];
sx q[1];
rz(-2.6032762) q[1];
sx q[1];
rz(-1.946279) q[1];
x q[2];
rz(0.94445618) q[3];
sx q[3];
rz(-2.3271022) q[3];
sx q[3];
rz(0.32523793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.66723055) q[2];
sx q[2];
rz(-1.2381866) q[2];
sx q[2];
rz(-0.24307069) q[2];
rz(2.4754751) q[3];
sx q[3];
rz(-2.5770498) q[3];
sx q[3];
rz(1.8977785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3617525) q[0];
sx q[0];
rz(-3.0267974) q[0];
sx q[0];
rz(2.6932122) q[0];
rz(-1.386863) q[1];
sx q[1];
rz(-1.9877537) q[1];
sx q[1];
rz(2.8853436) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.019779531) q[0];
sx q[0];
rz(-2.461713) q[0];
sx q[0];
rz(-2.7133184) q[0];
rz(-pi) q[1];
rz(-0.82661144) q[2];
sx q[2];
rz(-2.5275143) q[2];
sx q[2];
rz(2.5342864) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.2910034) q[1];
sx q[1];
rz(-0.91445078) q[1];
sx q[1];
rz(-2.7235051) q[1];
x q[2];
rz(-1.5330629) q[3];
sx q[3];
rz(-1.3940666) q[3];
sx q[3];
rz(2.7023774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3391352) q[2];
sx q[2];
rz(-2.4352303) q[2];
sx q[2];
rz(2.5668872) q[2];
rz(1.7859219) q[3];
sx q[3];
rz(-1.1698497) q[3];
sx q[3];
rz(-1.2333966) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.213585) q[0];
sx q[0];
rz(-1.428823) q[0];
sx q[0];
rz(0.25948778) q[0];
rz(1.9909987) q[1];
sx q[1];
rz(-1.3535627) q[1];
sx q[1];
rz(2.4096699) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4703092) q[0];
sx q[0];
rz(-2.0844315) q[0];
sx q[0];
rz(-1.0259823) q[0];
x q[1];
rz(-1.6875793) q[2];
sx q[2];
rz(-2.2846966) q[2];
sx q[2];
rz(1.2920979) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1720097) q[1];
sx q[1];
rz(-1.828555) q[1];
sx q[1];
rz(1.7715363) q[1];
rz(-pi) q[2];
rz(2.6287574) q[3];
sx q[3];
rz(-2.8898015) q[3];
sx q[3];
rz(0.73392111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6966454) q[2];
sx q[2];
rz(-1.2735294) q[2];
sx q[2];
rz(-0.15110061) q[2];
rz(-0.54667306) q[3];
sx q[3];
rz(-1.0497382) q[3];
sx q[3];
rz(-0.37724075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54995173) q[0];
sx q[0];
rz(-2.0786091) q[0];
sx q[0];
rz(-1.2623825) q[0];
rz(-1.4683912) q[1];
sx q[1];
rz(-0.60931283) q[1];
sx q[1];
rz(2.343822) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7933465) q[0];
sx q[0];
rz(-0.12013809) q[0];
sx q[0];
rz(-1.5500463) q[0];
rz(-2.5227929) q[2];
sx q[2];
rz(-0.36703645) q[2];
sx q[2];
rz(0.4180846) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3623912) q[1];
sx q[1];
rz(-2.5070842) q[1];
sx q[1];
rz(-1.7290551) q[1];
rz(-pi) q[2];
rz(-0.4425211) q[3];
sx q[3];
rz(-1.2708775) q[3];
sx q[3];
rz(0.72344852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.795934) q[2];
sx q[2];
rz(-0.63085932) q[2];
sx q[2];
rz(-0.30203715) q[2];
rz(-1.1473514) q[3];
sx q[3];
rz(-1.4619504) q[3];
sx q[3];
rz(2.5938477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40925947) q[0];
sx q[0];
rz(-0.091826037) q[0];
sx q[0];
rz(-1.9858032) q[0];
rz(1.0844768) q[1];
sx q[1];
rz(-0.98025727) q[1];
sx q[1];
rz(3.0715122) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.018774059) q[0];
sx q[0];
rz(-2.9111324) q[0];
sx q[0];
rz(-2.3019058) q[0];
rz(1.0191392) q[2];
sx q[2];
rz(-2.2746804) q[2];
sx q[2];
rz(-0.64955074) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1706714) q[1];
sx q[1];
rz(-1.0075924) q[1];
sx q[1];
rz(1.3794273) q[1];
rz(-0.035404215) q[3];
sx q[3];
rz(-1.6453214) q[3];
sx q[3];
rz(0.41985928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.30248102) q[2];
sx q[2];
rz(-1.182686) q[2];
sx q[2];
rz(-2.5202259) q[2];
rz(1.7012043) q[3];
sx q[3];
rz(-2.6337603) q[3];
sx q[3];
rz(-0.5733718) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086534111) q[0];
sx q[0];
rz(-1.4165514) q[0];
sx q[0];
rz(-0.85987464) q[0];
rz(-1.2043918) q[1];
sx q[1];
rz(-0.87221968) q[1];
sx q[1];
rz(3.133657) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82848362) q[0];
sx q[0];
rz(-0.61921739) q[0];
sx q[0];
rz(2.6909268) q[0];
rz(0.87848778) q[2];
sx q[2];
rz(-1.2307067) q[2];
sx q[2];
rz(1.0483339) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8994645) q[1];
sx q[1];
rz(-0.90485307) q[1];
sx q[1];
rz(-2.4813586) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8273724) q[3];
sx q[3];
rz(-1.946297) q[3];
sx q[3];
rz(0.34208959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.69592151) q[2];
sx q[2];
rz(-1.9182703) q[2];
sx q[2];
rz(2.725214) q[2];
rz(1.773206) q[3];
sx q[3];
rz(-1.2973283) q[3];
sx q[3];
rz(-2.2369475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96173441) q[0];
sx q[0];
rz(-2.8680153) q[0];
sx q[0];
rz(2.7767048) q[0];
rz(-2.2015613) q[1];
sx q[1];
rz(-0.53661984) q[1];
sx q[1];
rz(-1.6392802) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4125) q[0];
sx q[0];
rz(-1.3206498) q[0];
sx q[0];
rz(1.581122) q[0];
x q[1];
rz(2.084923) q[2];
sx q[2];
rz(-2.3377315) q[2];
sx q[2];
rz(-0.67509292) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2316206) q[1];
sx q[1];
rz(-0.93712229) q[1];
sx q[1];
rz(0.85068591) q[1];
x q[2];
rz(-1.3650465) q[3];
sx q[3];
rz(-1.4806517) q[3];
sx q[3];
rz(0.84116018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.49729785) q[2];
sx q[2];
rz(-2.6361894) q[2];
sx q[2];
rz(-2.0765182) q[2];
rz(0.30125695) q[3];
sx q[3];
rz(-1.732429) q[3];
sx q[3];
rz(1.5208972) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97312462) q[0];
sx q[0];
rz(-3.0597866) q[0];
sx q[0];
rz(-0.43564963) q[0];
rz(1.7565953) q[1];
sx q[1];
rz(-0.47043097) q[1];
sx q[1];
rz(-0.41697821) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96796658) q[0];
sx q[0];
rz(-2.2244503) q[0];
sx q[0];
rz(-0.047441479) q[0];
x q[1];
rz(-2.4531104) q[2];
sx q[2];
rz(-0.79198972) q[2];
sx q[2];
rz(2.288523) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0028487) q[1];
sx q[1];
rz(-1.4321623) q[1];
sx q[1];
rz(-2.0744051) q[1];
rz(-pi) q[2];
rz(0.74420332) q[3];
sx q[3];
rz(-2.9508698) q[3];
sx q[3];
rz(-0.71803367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.10432648) q[2];
sx q[2];
rz(-1.4929079) q[2];
sx q[2];
rz(-2.9157675) q[2];
rz(0.2078235) q[3];
sx q[3];
rz(-2.4184629) q[3];
sx q[3];
rz(-0.58661714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41480961) q[0];
sx q[0];
rz(-2.8746958) q[0];
sx q[0];
rz(-1.6171932) q[0];
rz(2.1879451) q[1];
sx q[1];
rz(-1.2748268) q[1];
sx q[1];
rz(-1.8189925) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1117301) q[0];
sx q[0];
rz(-0.39806453) q[0];
sx q[0];
rz(-1.2533305) q[0];
x q[1];
rz(-2.0199213) q[2];
sx q[2];
rz(-1.7977062) q[2];
sx q[2];
rz(0.34044468) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9433371) q[1];
sx q[1];
rz(-0.88611929) q[1];
sx q[1];
rz(-2.4556922) q[1];
x q[2];
rz(-1.3721458) q[3];
sx q[3];
rz(-2.47654) q[3];
sx q[3];
rz(-2.8952451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3502729) q[2];
sx q[2];
rz(-2.0952756) q[2];
sx q[2];
rz(2.1255169) q[2];
rz(1.919205) q[3];
sx q[3];
rz(-0.17761579) q[3];
sx q[3];
rz(0.55541742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14324698) q[0];
sx q[0];
rz(-0.9220985) q[0];
sx q[0];
rz(-2.0621598) q[0];
rz(1.7779508) q[1];
sx q[1];
rz(-1.2294055) q[1];
sx q[1];
rz(1.3235863) q[1];
rz(-0.50921847) q[2];
sx q[2];
rz(-1.5621395) q[2];
sx q[2];
rz(-0.12315673) q[2];
rz(0.090311269) q[3];
sx q[3];
rz(-1.3406546) q[3];
sx q[3];
rz(-1.8484074) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
